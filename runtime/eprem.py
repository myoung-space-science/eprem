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


class RunLog(collections.abc.Mapping):
    """Mapping-based interface to an EPREM project log."""

    def __init__(self, path: PathLike, **common) -> None:
        """Create a new project log.
        
        Parameters
        ----------
        path : path-like
            The path at which to create the log file. May be relative to the
            current directory.

        common
            Key-value pairs of attributes that are common to all runs.
        """
        self._path = fullpath(path)
        self._common = common

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

    def create(self, key: str, source=None, filetype: str=None):
        """Create a new entry in this log file."""
        contents = self._asdict.copy()
        if source is None:
            contents[key] = {}
        contents.update({key: self._normalize_source(source, filetype)})
        self.dump(contents)
        return self

    def _normalize_source(self, source, filetype):
        """Convert `source` into an appropriate dictionary."""
        if source is None:
            return {}
        if isinstance(source, typing.Mapping):
            return dict(source)
        if isinstance(source, str):
            return self._read_from_file(source, filetype)
        raise TypeError(f"Unrecognized source type: {type(source)}")

    def load(self, source: str, filetype: str):
        """Create a new entry in this log file from a file."""
        self.dump(self._read_from_file(source, filetype))
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

    def rename(self, source: str, target: str, subdir: str=None):
        """Rename `source` to `target` in this log file."""
        current = self._asdict.copy()
        try:
            record = current[source]
        except KeyError as err:
            raise LogKeyError(
                f"Cannot rename unknown run {source!r}"
            ) from err
        if subdir:
            original = record[subdir]
            new = pathlib.Path(original['directory']).parent / target
            tmp = {
                k: v if k != 'directory' else str(new)
                for k, v in original.items()
            }
            renamed = {subdir: tmp}
        else:
            new = pathlib.Path(record['directory']).parent / target
            renamed = {
                k: v if k != 'directory' else str(new)
                for k, v in record.items()
            }
        updated = {k: v for k, v in current.items() if k != source}
        updated[target] = renamed
        self.dump(updated)
        return self

    def remove(self, *targets: str, subdir: str=None):
        """Remove the target run(s) from this log file."""
        current = self._asdict.copy()
        if target := next((t for t in targets if t not in current), None):
            raise LogKeyError(f"Cannot remove unknown run {target!r}")
        updated = {k: v for k, v in current.items() if k not in targets}
        if subdir:
            for run, info in current.items():
                if run in targets:
                    updated[run] = {
                        k: v for k, v in info.items() if k != subdir
                    }
        self.dump(updated)
        return self

    @property
    def _asdict(self) -> typing.Dict[str, typing.Dict[str, typing.Any]]:
        """Internal dictionary representing the current contents."""
        try:
            with self.path.open('r') as fp:
                return dict(json.load(fp))
        except FileNotFoundError:
            self.dump(self._common)
            return self._asdict

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
    def branches(self):
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
        rundir: str=None,
        logname: str=None,
    ) -> None: ...

    database = pathlib.Path('.eprem-runtime.json')

    def __init__(self, root, **kwargs):
        """Initialize a new project."""
        attrs = self._init_attrs(root, kwargs)
        self._log = None
        self._name = None
        directories = [
            attrs.path / branch / attrs.rundir
            for branch in attrs.branches or ['']
        ]
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
        self._attrs = attrs
        self._directories = directories

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

    def spawn(
        self: ProjectType,
        config: str,
        name: str=None,
        subset: typing.Union[str, typing.Iterable[str]]=None,
        nprocs: int=None,
        runlog: str='eprem.log',
        environment: typing.Dict[str, str]=None,
        silent: bool=False,
    ) -> ProjectType:
        """Set up and execute a new EPREM run within this project."""
        directories = (
            {subset} if isinstance(subset, str)
            else set(subset or ())
        )
        for path in self._make_paths(name, directories):
            branch = path.parent.parent
            mpirun = self._locate('mpirun', branch, environment or {})
            eprem = self._locate('eprem', branch, environment or {})
            shutil.copy(config, path / self._attrs.config)
            command = (
                "nice -n 10 ionice -c 2 -n 3 "
                f"{mpirun} --mca btl_base_warn_component_unused 0 "
                f"-n {nprocs or 1} {eprem} eprem.cfg"
            )
            output = path / runlog
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
                    print(
                        f"[PID: {process.pid} @ {now}] "
                        f"{name} --> {branch.name}/"
                    )
                process.wait()
            if not silent and process.returncode:
                print(f"WARNING: Process exited with {process.returncode}\n")
            logentry = {
                'mpirun': str(mpirun),
                'eprem': str(eprem),
                'directory': str(path),
                'time': datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
            }
            self.log.create(branch.name, logentry)

    def _locate(
        self,
        name: str,
        path: pathlib.Path,
        environment: typing.Dict[str, str]
    ) -> typing.Union[str, pathlib.Path]:
        """Compute an appropriate path to the named element.
        
        Intended for internal use by `~eprem.Project`.

        This method will attempt to create a full path (resolving links as
        necessary) from `environment` or from `path / name`. If neither
        exist, it will return `name` as-is, thereby allowing calling code to
        default to the searching the system path.
        """
        location = environment.get(name) or path / name
        it = fullpath(os.path.realpath(location))
        return it if it.exists() else name

    def rename(
        self: ProjectType,
        source: str,
        target: str,
        subset: typing.Union[str, typing.Iterable[str]]=None,
        silent: bool=False,
    ) -> ProjectType:
        """Rename an existing EPREM run within this project."""
        directories = (
            {subset} if isinstance(subset, str)
            else set(subset or ())
        )
        paths = self._rename_paths(source, target, directories)
        if not paths:
            if not silent:
                print(f"Nothing to rename for {source!r}")
            return
        branches = [path.parent.parent.name for path in paths]
        for branch in branches:
            self.log.rename(source, target, branch)
        if not silent:
            print(f"Updated {self.log.path}")
            underline(f"{self.root}")
            for branch in branches:
                print(f"[{source} -> {target}] in {branch}")

    def remove(
        self: ProjectType,
        run: str,
        subset: typing.Union[str, typing.Iterable[str]]=None,
        silent: bool=False,
    ) -> ProjectType:
        """Remove an existing EPREM run from this project."""
        directories = (
            {subset} if isinstance(subset, str)
            else set(subset or ())
        )
        paths = self._remove_paths(run, directories)
        if not paths:
            if not silent:
                print(f"Nothing to remove for {run!r}")
        branches = [path.parent.parent.name for path in paths]
        targets = [path.name for path in paths]
        for branch in branches:
            self.log.remove(*targets, branch)
        if not silent:
            print(f"Updated {self.log.path}")
            underline(f"{self.root}")
            for branch in branches:
                for target in targets:
                    print(f"removed {target} from {branch}")

    def show(self: ProjectType, *runs: str):
        """Display information about this project or the named run(s)."""

    def _make_paths(self, name: str, subset: typing.Set[str]):
        """Create the target subdirectory in each branch."""
        rundirs = self._get_rundirs(subset)
        paths = [rundir / name for rundir in rundirs]
        action = self._make_paths_check(*paths)
        for path in paths:
            action(path)
        return paths

    def _make_paths_check(self, *paths: pathlib.Path):
        """Return a function to create paths only if safe to do so."""
        def action(path: pathlib.Path):
            path.mkdir(parents=True)
        for path in paths:
            if path.exists():
                raise PathOperationError(
                    f"Cannot create {path}: already exists"
                ) from None
        return action

    def _rename_paths(self, src: str, dst: str, subset: typing.Set[str]):
        """Rename `src` to `dst` in all subdirectories."""
        rundirs = self._get_rundirs(subset)
        pairs = [(rundir / src, rundir / dst) for rundir in rundirs]
        action = self._rename_paths_check(*pairs)
        for pair in pairs:
            action(*pair)
        return rundirs

    def _rename_paths_check(self, *pairs: typing.Tuple[_P, _P]):
        """Return a function to rename paths only if safe to do so."""
        def action(old: _P, new: _P):
            old.rename(new)
        for old, new in pairs:
            if not old.exists():
                raise PathOperationError(
                    f"Cannot rename {old}: does not exist"
                ) from None
            if not old.is_dir():
                raise PathOperationError(
                    f"Cannot rename {old}: not a directory"
                ) from None
            if new.exists():
                raise PathOperationError(
                    f"Renaming {old.name!r} to {new.name!r} would "
                    f"overwrite {new}."
                ) from None
        return action

    def _remove_paths(self, name: str, subset: typing.Set[str]):
        """Remove `target` from all subdirectories."""
        rundirs = self._get_rundirs(subset)
        paths = [
            path for rundir in rundirs
            for path in rundir.glob(name)
        ]
        action = self._remove_paths_check(*paths)
        for path in paths:
            action(path)
        return rundirs

    def _remove_paths_check(self, *paths: pathlib.Path):
        """Return a function to remove paths only if safe to do so."""
        def action(path: pathlib.Path):
            shutil.rmtree(path)
        for path in paths:
            if not path.exists():
                raise PathOperationError(
                    f"Cannot remove {path}: does not exist"
                ) from None
            if not path.is_dir():
                raise PathOperationError(
                    f"Cannot remove {path}: not a directory"
                ) from None
        return action

    def _get_rundirs(self, subset: typing.Iterable[str]):
        """Build an appropriate collection of runtime directories."""
        if not subset:
            return self.directories
        return [d for d in self.directories if d.name in subset]

    @property
    def log(self):
        """The log of runs in this project."""
        if self._log is None:
            self._log = RunLog(
                self.root / self._attrs.logname,
                config=self._attrs.config,
            )
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
