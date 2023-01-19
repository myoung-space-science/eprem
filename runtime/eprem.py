import argparse
import collections.abc
import contextlib
import datetime
import json
import os
import pathlib
import shutil
import subprocess
import types
import typing

try:
    import yaml
    _HAVE_YAML = True
except ModuleNotFoundError:
    _HAVE_YAML = False


PathLike = typing.Union[str, os.PathLike]


def fullpath(p: PathLike):
    """Expand and resolve the given path."""
    return pathlib.Path(p).expanduser().resolve()


class LogKeyError(KeyError):
    """The log has no entry for a given key."""


class RunKeyError(Exception):
    """The log entry does not include a given key."""


class ReadTypeError(Exception):
    """There is no support for reading a given file type."""


RunLogType = typing.TypeVar('RunLogType', bound='RunLog')


class RunLog(collections.abc.Mapping):
    """Mapping-based interface to an EPREM project log."""

    def __init__(
        self,
        path: PathLike,
        branches: typing.Iterable[str]=None,
        **common
    ) -> None:
        """Create a new project log.
        
        Parameters
        ----------
        path : path-like
            The path at which to create the log file. May be relative to the
            current directory.

        branches : iterable of strings, optional
            The branch names, if any, in the project.

        common
            Key-value pairs of attributes that are common to all runs.
        """
        self._path = fullpath(path)
        self._branches = tuple(branches or [])
        self.dump(common)

    def __len__(self) -> int:
        """Called for len(self)."""
        return len(self._asdict)

    def __iter__(self) -> typing.Iterator[str]:
        """Called for iter(self)."""
        return iter(self._asdict)

    def __getitem__(self, __k: str):
        """Get metadata for a run."""
        if __k in self._asdict:
            return self._asdict[__k]
        raise LogKeyError(f"Unknown run {__k!r}")

    @typing.overload
    def create(
        self: RunLogType,
        key: str,
        mapping: typing.Mapping
    ) -> RunLogType:
        """Create a new entry from a mapping."""

    @typing.overload
    def create(
        self: RunLogType,
        key: str,
        **items
    ) -> RunLogType:
        """Create a new entry from key-value pairs."""

    def create(self, key: PathLike, *mapping: typing.Mapping, **items):
        contents = self._asdict.copy()
        if mapping and items:
            raise TypeError(
                f"{self.__class__.__qualname__}.create accepts"
                " a single mapping or multiple key-value pairs"
                ", but not both"
            ) from None
        new = dict(*mapping) or items
        contents[str(key)] = new
        self.dump(contents)
        return self

    def load(self, key: PathLike, source: PathLike, filetype: str=None):
        """Create a new entry in this log file from a file."""
        contents = self._asdict.copy()
        contents[str(key)] = self._read_from_file(source, filetype)
        self.dump(contents)
        return self

    def _read_from_file(self, source, filetype):
        """Update the current contents from the `source` file."""
        contents = self._asdict.copy()
        loader = self._source_loader(filetype or pathlib.Path(source).suffix)
        with pathlib.Path(source).open('r') as fp:
            if loaded := loader(fp):
                contents.update(loaded)
        return contents

    def _source_loader(self, filetype: str):
        """Get a format-specific file-reader."""
        if filetype.lower().lstrip('.') == 'yaml':
            if _HAVE_YAML:
                return yaml.safe_load
            raise ReadTypeError("No support for reading YAML") from None
        if filetype.lower().lstrip('.') == 'json':
            return json.load
        raise ValueError(
            f"Unknown file type: {filetype!r}"
        ) from None

    def append(self, target: str, key: str, metadata):
        """Append metadata to `target`."""
        contents = self._asdict.copy()
        try:
            record = contents[target]
        except KeyError as err:
            raise LogKeyError(
                f"Cannot append to unknown run {target!r}"
            ) from err
        record[key] = metadata
        self.dump(contents)
        return self

    def mv(self, source: PathLike, target: PathLike):
        """Rename `source` to `target` in this log file."""
        old = str(source)
        new = str(target)
        if old not in self:
            raise LogKeyError(
                f"Cannot rename unknown run {old!r}"
            ) from None
        updated = {k: v for k, v in self._asdict.items() if k != old}
        updated[new] = self._asdict[old]
        self.dump(updated)
        return self

    def rm(self, target: PathLike):
        """Remove the target run from this log file."""
        updated = {k: v for k, v in self._asdict.items() if k != str(target)}
        self.dump(updated)
        return self

    @property
    def _asdict(self) -> typing.Dict[str, typing.Dict[str, typing.Any]]:
        """Internal dictionary representing the current contents."""
        with self.path.open('r') as fp:
            return dict(json.load(fp))

    def dump(self, contents):
        """Write `contents` to this log file."""
        with self.path.open('w') as fp:
            json.dump(contents, fp, indent=4, sort_keys=True)

    @property
    def name(self):
        """The name of this log file.
        
        Same as `RunLog.path.name`.
        """
        return self.path.name

    @property
    def path(self):
        """The path to this log file."""
        return self._path

    def __str__(self) -> str:
        """A simplified representation of this object."""
        return str(self._asdict)

    def __repr__(self) -> str:
        """An unambiguous representation of this object."""
        return f"{self.__class__.__qualname__}({self.path})"


class PathOperationError(Exception):
    """This operation is not allowed on the given path(s)."""


class ProjectExistsError(Exception):
    """A project with this name already exists."""


_P = typing.TypeVar('_P', bound=pathlib.Path)


class _ProjectInit(typing.Mapping):
    """A mapping of `~Project` initialization attributes."""

    _kwargs = {
        'branches': {'type': tuple, 'default': ()},
        'config': {'type': str, 'default': 'eprem.cfg'},
        'output': {'type': str, 'default': 'eprem.log'},
        'rundir': {'type': str, 'default': 'runs'},
        'logstem': {'type': str, 'default': 'runs'},
    }

    def __init__(self, root, **kwargs) -> None:
        """Create a new instance."""
        self._attrs = {
            key: this['type'](kwargs.get(key) or this['default'])
            for key, this in self._kwargs.items()
        }
        self._attrs['root'] = str(root)
        self._path = None
        self._root = None
        self._branches = None
        self._config = None
        self._output = None
        self._rundir = None
        self._logname = None
        self._logstem = None

    @property
    def path(self):
        """A fully qualified path to the project root directory."""
        if self._path is None:
            self._path = fullpath(self.root)
        return self._path

    @property
    def root(self):
        """The project root directory."""
        if self._root is None:
            self._root = str(self._attrs['root'])
        return self._root

    @property
    def branches(self) -> typing.Tuple[str]:
        """The distinct project branches, if any."""
        if self._branches is None:
            self._branches = tuple(self._attrs['branches'])
        return self._branches

    @property
    def config(self):
        """The name of the standard project configuration file."""
        if self._config is None:
            self._config = str(self._attrs['config'])
        return self._config

    @property
    def output(self):
        """The name of the standard project output log."""
        if self._output is None:
            # Ensure a string file name, even if the initialization argument was
            # a path. We don't want to write to an arbitrary location on disk!
            self._output = pathlib.Path(self._attrs['output']).name
        return self._output

    @property
    def rundir(self):
        """The name of the standard project run directory."""
        if self._rundir is None:
            self._rundir = str(self._attrs['rundir'])
        return self._rundir

    @property
    def logname(self):
        """The name of the project-wide log file."""
        if self._logname is None:
            self._logname = pathlib.Path(
                self.logstem
            ).with_suffix('.json').name
        return self._logname

    @property
    def logstem(self):
        """The name (without suffix) of the project-wide log file."""
        if self._logstem is None:
            self._logstem = str(self._attrs['logstem'])
        return self._logstem

    def __len__(self) -> int:
        """Called for len(self)."""
        return len(self._attrs)

    def __iter__(self) -> typing.Iterable[str]:
        """Called for iter(self)."""
        return iter(self._attrs)

    def __getitem__(self, __k: str):
        """Key-based access to attribute values."""
        if __k in self._attrs:
            return self._attrs[__k]
        raise KeyError(f"Unknown attribute {__k!r}")

    def __repr__(self) -> str:
        """An unambiguous representation of this object."""
        return '\n'.join(f"{k}: {v}" for k, v in self.items())


ProjectType = typing.TypeVar('ProjectType', bound='Project')


class Project:
    """Interface to an EPREM runtime project."""

    @typing.overload
    def __init__(
        self,
        root: typing.Union[str, pathlib.Path],
    ) -> None: ...

    @typing.overload
    def __init__(
        self,
        root: typing.Union[str, pathlib.Path],
        branches: typing.Iterable[str]=None,
        config: str=None,
        output: str=None,
        rundir: str=None,
        logname: str=None,
    ) -> None: ...

    database = pathlib.Path('.eprem-runtime.json')

    def __init__(self, root, **kwargs):
        """Initialize a new project."""
        self._isvalid = False
        attrs = self._init_attrs(root, kwargs)
        self._log = None
        self._name = None
        directories = [
            attrs.path / branch / attrs.rundir
            for branch in attrs.branches or ['']
        ]
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
        self._log = RunLog(
            attrs.path / attrs.logname,
            branches=attrs.branches,
            config=attrs.config,
            output=attrs.output,
        )
        self._attrs = attrs
        self._directories = directories
        self._isvalid = True

    def _init_attrs(self, root: pathlib.Path, kwargs: dict):
        """Initialize arguments from input or the database."""
        path = fullpath(root)
        if path.exists() and kwargs:
            existing = (
                f"{self.__class__.__qualname__}"
                f"({os.path.relpath(path)!r})"
            )
            raise ProjectExistsError(
                f"The project {path.name!r} already exists in {path.parent}. "
                f"You can access the existing project via {existing}"
            )
        key = str(path)
        if path.exists():
            with self.database.open('r') as fp:
                existing = dict(json.load(fp))
            return _ProjectInit(**existing[key])
        path.mkdir(parents=True)
        init = _ProjectInit(root=path, **kwargs)
        with self.database.open('w') as fp:
            json.dump({key: dict(init)}, fp, indent=4, sort_keys=True)
        return init

    def show(self: ProjectType, *runs: str):
        """Display information about this project or the named run(s)."""
        if not runs:
            self._show_project()
        if len(runs) == 1:
            self._show_run(run)
        for run in runs:
            underline(run)
            self._show_run(run)

    def _show_project(self):
        """Display information about this project."""
        underline("Project")
        print(self._attrs)
        if not self.branches:
            underline("runs")
            print('\n'.join(self.runs))
            return
        for branch, runs in self.branches.items():
            underline(f"Branch {branch}")
            print('\n'.join(runs))

    def _show_run(self, run: str):
        """Display information about the named run."""
        try:
            branches = self.runs[run]
        except KeyError:
            print(f"No run named {run!r}")
            return
        if not branches:
            print(self.root / self._attrs.rundir / run)
        for branch in self.runs[run]:
            print(self.root / branch / self._attrs.rundir / run)

    def reset(self, force: bool=False, silent: bool=False):
        """Reset this project to its initial state."""
        return self.rm('*', errors=(not force), silent=silent)

    def remove(self, force: bool=False, silent: bool=False):
        """Delete this project."""
        shutil.rmtree(self.root, ignore_errors=force)
        with self.database.open('r') as fp:
            current = dict(json.load(fp))
        updated = {
            k: v for k, v in current.items()
            if k != str(self.root)
        }
        with self.database.open('w') as fp:
            json.dump(updated, fp, indent=4, sort_keys=True)
        self._isvalid = False
        if not silent:
            print(f"Removed project at {self.root}")

    def run(
        self: ProjectType,
        config: PathLike,
        name: str=None,
        branches: typing.Union[str, typing.Iterable[str]]=None,
        nproc: int=None,
        environment: typing.Dict[str, str]=None,
        errors: bool=False,
        silent: bool=False,
    ) -> ProjectType:
        """Set up and execute a new EPREM run within this project."""
        for path in self._build_run_paths(name, branches):
            self._create_run(
                config,
                path,
                nproc=nproc,
                environment=environment,
                errors=errors,
                silent=silent,
            )
        return self

    def _build_run_paths(self, name: str, branches):
        """Build a list of paths in which to run the simulation."""
        rundirs = self._get_rundirs(branches)
        run = name or datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S.%f')
        return [rundir / run for rundir in rundirs]

    def _create_run(
        self: ProjectType,
        config: PathLike,
        path: pathlib.Path,
        nproc: int=None,
        environment: typing.Dict[str, str]=None,
        errors: bool=False,
        silent: bool=False,
    ) -> None:
        """Create a single run."""
        if error := self._try_to_make(path):
            if errors:
                raise PathOperationError(error)
            if not silent:
                print(error)
            return
        shutil.copy(config, path / self._attrs.config)
        branch = path.parent.parent
        mpirun = _locate('mpirun', branch, environment or {})
        eprem = _locate('eprem', branch, environment or {})
        command = (
            "nice -n 10 ionice -c 2 -n 3 "
            f"{mpirun} --mca btl_base_warn_component_unused 0 "
            f"-n {nproc or 1} {eprem} eprem.cfg"
        )
        output = path / self._attrs.output
        now = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        with output.open('w') as stdout:
            process = subprocess.Popen(
                command,
                shell=True,
                cwd=path,
                stdout=stdout,
                stderr=subprocess.STDOUT,
            )
            if not silent:
                print(f"\n[{process.pid}]")
                print(f"Started at {now}")
            process.wait()
            if not silent:
                print(f"Created {path.name} in branch {branch.name!r}")
        if process.returncode == 0:
            logentry = {
                'mpirun': str(mpirun),
                'eprem': str(eprem),
                'time': now,
            }
            self.log.create(str(path), logentry)
        elif not silent:
            print(
                f"WARNING: Process exited with {process.returncode}",
                end='\n\n',
            )

    def _try_to_make(self, this: pathlib.Path):
        """Create `this` only if safe to do so."""
        if this.exists():
            return f"Cannot create {this}: already exists"
        this.mkdir(parents=True)

    def mv(
        self: ProjectType,
        source: str,
        target: str,
        branches: typing.Union[str, typing.Iterable[str]]=None,
        errors: bool=False,
        silent: bool=False,
    ) -> ProjectType:
        """Rename an existing EPREM run within this project."""
        pairs = self._build_mv_pairs(source, target, branches)
        if not pairs:
            if not silent:
                print(f"Nothing to rename for {source!r}")
            return
        for (run, new) in pairs:
            self._rename_run(run, new, errors=errors, silent=silent)
        return self

    def _build_mv_pairs(self, source: str, target: str, branches):
        """Build a list of path pairs for renaming."""
        rundirs = self._get_rundirs(branches)
        return [(rundir / source, rundir / target) for rundir in rundirs]

    def _rename_run(
        self: ProjectType,
        run: pathlib.Path,
        new: pathlib.Path,
        errors: bool=False,
        silent: bool=False,
    ) -> None:
        """Rename a single run."""
        if error := self._try_to_rename(run, new):
            if errors:
                raise PathOperationError(error)
            if not silent:
                print(error)
            return
        self.log.mv(run, new)
        if not silent:
            branch = self._get_branch_name(run)
            base = f"Renamed {run.name!r} to {new.name!r}"
            print(f"{base} in branch {branch!r}" if branch else base)

    def _try_to_rename(self, this: pathlib.Path, that: pathlib.Path):
        """Rename `this` to `that` only if safe to do so."""
        if not this.exists():
            return f"Cannot rename {this}: does not exist"
        if not this.is_dir():
            return f"Cannot rename {this}: not a directory"
        if that.exists():
            return (
                f"Renaming {this.name!r} to {that.name!r} would "
                f"overwrite {that}."
            )
        this.rename(that)

    def rm(
        self: ProjectType,
        run: str,
        branches: typing.Union[str, typing.Iterable[str]]=None,
        errors: bool=False,
        silent: bool=False,
    ) -> ProjectType:
        """Remove an existing EPREM run from this project."""
        paths = self._build_rm_paths(run, branches)
        if not paths:
            if not silent:
                print(f"Nothing to remove for {run!r}")
            return
        for path in paths:
            self._remove_run(path, errors=errors, silent=silent)
        return self

    def _build_rm_paths(self, target: str, branches):
        """Build a list of paths to remove."""
        rundirs = self._get_rundirs(branches)
        return [
            path for rundir in rundirs
            for path in rundir.glob(target)
        ]

    def _remove_run(
        self: ProjectType,
        run: pathlib.Path,
        errors: bool=False,
        silent: bool=False,
    ) -> None:
        """Remove a single run."""
        if error := self._try_to_remove(run):
            if errors:
                raise PathOperationError(error)
            if not silent:
                print(error)
            return
        self.log.rm(run)
        if not silent:
            base = f"Removed {run.name!r}"
            branch = self._get_branch_name(run)
            print(f"{base} from branch {branch!r}" if branch else base)

    def _try_to_remove(self, this: pathlib.Path):
        """Remove `this` only if safe to do so."""
        if not this.exists():
            return f"Cannot remove {this}: does not exist"
        if not this.is_dir():
            return f"Cannot remove {this}: not a directory"
        shutil.rmtree(this)

    def _get_branch_name(self, path: pathlib.Path):
        """Get the project branch name, if any, of `path`."""
        parents = (str(p) for p in path.relative_to(self.root).parents)
        with contextlib.suppress(StopIteration):
            if this := next(p for p in parents if p in self.branches):
                return this

    def _get_rundirs(
        self,
        branches: typing.Union[str, typing.Iterable[str]]=None,
    ) -> typing.List[pathlib.Path]:
        """Build an appropriate collection of runtime directories."""
        if not branches:
            return self.directories
        subset = (
            {branches} if isinstance(branches, str)
            else set(branches)
        )
        return [d for d in self.directories if d.parent.name in subset]

    @property
    def log(self):
        """The log of runs in this project."""
        return self._log

    @property
    def name(self):
        """The name of this project.
        
        Same as `Project.root.name` or `Project.base.name`.
        """
        if self._name is None:
            self._name = self.root.name
        return self._name

    @property
    def runs(self):
        """The available runs, and owning branches, if any."""
        available = {
            run
            for runs in self.branches.values()
            for run in runs
        }
        this = {run: set() for run in available}
        for branch, runs in self.branches.items():
            for run in runs:
                this[run] |= {branch}
        return this

    @property
    def branches(self):
        """The project branches, if any, and their available runs."""
        if not self._attrs.branches:
            return {}
        return {
            branch: {
                path.name for d in self._get_rundirs([branch])
                for path in d.glob('*')
            }
            for branch in self._attrs.branches
        }

    @property
    def directories(self):
        """The full path to each run directory."""
        return self._directories

    @property
    def base(self):
        """Alias for `Project.root`."""
        return self.root

    @property
    def root(self):
        """The top-level directory of this project."""
        return self._attrs.path

    def __bool__(self) -> bool:
        """True if this is a valid project."""
        return self._isvalid

    def __eq__(self, other) -> bool:
        """True if two projects have the same initializing attributes."""
        if isinstance(other, Project):
            return self._attrs == other._attrs
        return NotImplemented

    def __str__(self) -> str:
        """A simplified representation of this object."""
        return self.name

    def __repr__(self) -> str:
        """An unambiguous representation of this object."""
        return f"{self.__class__.__qualname__}({self.root})"


def _locate(
    name: str,
    path: pathlib.Path,
    environment: typing.Dict[str, str]
) -> pathlib.Path:
    """Compute an appropriate path to the named element.

    Notes
    -----
    * Intended for use by `~eprem.Project`.
    * This function will attempt to create a full path (resolving links as
      necessary) based on `environment` or from `path / name`. If neither exist,
      it will return `name` as-is, thereby allowing calling code to default to
      the searching the system path.
    """
    location = environment.get(name) or path / name
    it = fullpath(os.path.realpath(location))
    return it if it.exists() else pathlib.Path(shutil.which(name))


def underline(text: str):
    """Print underlined text."""
    dashes = '-' * len(text)
    print(f"\n{text}")
    print(dashes)


def spawn():
    """Set up and execute a new EPREM run."""


def rename():
    """Rename an existing EPREM run."""


def remove():
    """Remove an existing EPREM run."""


def doc2help(func: types.FunctionType):
    """Convert a function docstring to CLI help text."""
    doclines = func.__doc__.split('\n')
    summary = doclines[0]
    return summary[0].lower() + summary[1:]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Support for operations on EPREM runs.",
    )
    parser.add_argument(
        '-l',
        '--logfile',
        help=(
            "path to the relevant log file"
            "\n(default: runs.json)"
        ),
    )
    parser.add_argument(
        '-v',
        '--verbose',
        help="print runtime messages",
        action='store_true',
    )
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        '--spawn',
        help=doc2help(spawn),
        nargs=2,
        metavar=('CONFIG', 'TARGET'),
    )
    mode_group.add_argument(
        '--rename',
        help=doc2help(rename),
        nargs=2,
        metavar=('SOURCE', 'TARGET'),
    )
    mode_group.add_argument(
        '--remove',
        help=doc2help(remove),
        metavar='TARGET',
    )
    cli = vars(parser.parse_args())
