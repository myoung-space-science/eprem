"""
Support for creating and modifying EPREM runs with a defined project.
"""

import argparse
import collections.abc
import contextlib
import datetime
import functools
import json
import os
import pathlib
import shutil
import subprocess
import textwrap
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


class PathTypeError(Exception):
    """The path is the wrong type."""


class ReadTypeError(Exception):
    """There is no support for reading a given file type."""


RunLogType = typing.TypeVar('RunLogType', bound='RunLog')


class RunLog(collections.abc.Mapping):
    """Mapping-based interface to an EPREM project log."""

    def __init__(self, path: PathLike, **common) -> None:
        """Create a new project log.
        
        Parameters
        ----------
        path : path-like
            The path to the log file, including the file name. If the file does
            not exist, this class will create a new one; otherwise, this class
            will use the existing file. The path may be relative to the current
            directory.

        common
            Key-value pairs of attributes that are common to all runs. Not
            allowed when retrieving an existing file.

        Raises
        ------
        TypeError
            Caller attempted to set common attributes on an existing log.
        """
        full = fullpath(path)
        if full.exists():
            self._path = full
            if common:
                raise TypeError(
                    "Cannot change common attributes "
                    "of an existing log"
                ) from None
        else:
            self._path = full
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

    def update_directory(self, target: PathLike):
        """Update the directory path to this file.

        Parameters
        ----------
        target : path-like
            The directory to which to move this file. May be relative.

        Raises
        ------
        PathTypeError
            The target path does not exist or is not a directory.

        Notes
        -----
        * This method assumes that the target path exists.
        * This method does not allow changes to the file name.
        """
        new = fullpath(target)
        if not new.exists():
            raise PathTypeError("The new path must exist") from None
        if not new.is_dir():
            raise PathTypeError("Cannot rename log file") from None
        self._path = new / self.name
        source = str(self.path.parent)
        for key in self:
            self.mv(key, key.replace(source, str(new)))
        return self

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

    def append(self, target: PathLike, key: str, value):
        """Append metadata to `target`."""
        contents = self._asdict.copy()
        run = str(target)
        try:
            record = contents[run]
        except KeyError as err:
            raise LogKeyError(
                f"Cannot append to unknown run {run!r}"
            ) from err
        record[key] = value
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
        display = {k: v for k, v in self.items() if k != 'branches'}
        display['branches'] = self['branches'] or None
        return '\n'.join(f"{k}: {v}" for k, v in display.items())


_Paths = typing.TypeVar(
    '_Paths',
    pathlib.Path,
    typing.Tuple[pathlib.Path, pathlib.Path],
)


class _PathOperation(collections.abc.Iterable):
    """Helper for path-based operations."""

    @typing.overload
    def __init__(
        self,
        operator: typing.Callable[
            [pathlib.Path],
            typing.Union[pathlib.Path, PathOperationError]
        ],
        paths: typing.Iterable[pathlib.Path]=None,
    ) -> None: ...

    @typing.overload
    def __init__(
        self,
        operator: typing.Callable[
            [pathlib.Path, pathlib.Path],
            typing.Union[pathlib.Path, PathOperationError]
        ],
        paths: typing.Iterable[typing.Tuple[pathlib.Path, pathlib.Path]]=None,
    ) -> None: ...

    def __init__(self, operator, paths: _Paths=None) -> None:
        """Create a new instance."""
        self.operator = operator
        self._paths = paths
        self._force = None
        self._silent = None

    def __bool__(self) -> bool:
        """Called for bool(self)."""
        return bool(self._paths)

    def __iter__(self):
        """Iteratively apply the operation, if possible."""
        error = None
        for path in self._paths:
            result = self.apply(path)
            if not isinstance(result, PathOperationError):
                yield result
                error = False
            else:
                error = True
            if error:
                if not self._force:
                    raise result
                if not self._silent:
                    print(result)

    def errors(self, force: bool=False, silent: bool=False):
        """Update how this instance handles errors and warnings."""
        self._force = force
        self._silent = silent
        return self

    def apply(self, this) -> typing.Union[pathlib.Path, typing.Tuple[pathlib.Path, pathlib.Path]]:
        """Apply the operator to a`this`, if possible."""
        if isinstance(this, pathlib.Path):
            return self.operator(this)
        try:
            allpaths = all(isinstance(i, pathlib.Path) for i in this)
        except TypeError:
            allpaths = False
        if allpaths:
            return self.operator(*this)


RunPathsType = typing.TypeVar('RunPathsType', bound='RunPaths')


class RunPaths(collections.abc.Collection):
    """Manager for path operations on runtime directories."""

    def __init__(
        self,
        root: PathLike,
        branches: typing.Iterable[str]=None,
        base: str=None,
    ) -> None:
        self._root = fullpath(root)
        self._base = base or 'runs'
        self._branches = branches
        self._listing = None
        for directory in self.listing:
            directory.mkdir(parents=True, exist_ok=True)

    def update(
        self: RunPathsType,
        root: PathLike=None,
        branches: typing.Iterable[str]=None,
        base: str=None,
    ) -> RunPathsType:
        """Update path components."""
        self._listing = None
        if root:
            path = fullpath(root)
            if path.exists():
                raise PathTypeError(
                    f"Renaming {self.root.name!r} to {path.name!r} would "
                    f"overwrite {path}."
                )
            self.root.rename(path)
            self._root = path
        if branches:
            self._branches = branches
        if base:
            self._base = base
        return self

    def __contains__(self, __x: PathLike) -> bool:
        """Called for __x in self."""
        return __x in self.listing

    def __iter__(self) -> typing.Iterator[pathlib.Path]:
        """Called for iter(self)."""
        return iter(self.listing)

    def __len__(self) -> int:
        """Called for len(self)."""
        return len(self.listing)

    def mkdir(
        self,
        target: PathLike,
        branches: typing.Union[str, typing.Iterable[str]]=None,
    ) -> _PathOperation:
        """Create target run(s) if possible."""
        paths = self.define(target, branches=branches)
        return _PathOperation(self._mkdir, paths=paths)

    def _mkdir(self, target: PathLike):
        """Create the target run if safe to do so."""
        this = fullpath(target)
        if this.exists():
            return PathOperationError(
                f"Cannot create {str(this)!r}: already exists"
            )
        this.mkdir(parents=True)
        return this

    def mv(
        self,
        source: PathLike,
        target: PathLike,
        branches: typing.Union[str, typing.Iterable[str]]=None,
    ) -> _PathOperation:
        """Rename `source` to `target` where possible."""
        pairs = self.define(source, target, branches=branches)
        return _PathOperation(self._mv, paths=pairs)

    def _mv(self, source: PathLike, target: PathLike):
        """Rename `source` to `target` only if safe to do so."""
        this = fullpath(source)
        that = fullpath(target)
        if not this.exists():
            return PathOperationError(
                f"Cannot rename {this}: does not exist"
            )
        if not this.is_dir():
            return PathOperationError(
                f"Cannot rename {this}: not a directory"
            )
        if that.exists():
            return PathOperationError(
                f"Renaming {this.name!r} to {that.name!r} would "
                f"overwrite {that}."
            )
        this.rename(that)
        return (this, that)

    def rm(
        self,
        target: PathLike,
        branches: typing.Union[str, typing.Iterable[str]]=None,
    ) -> _PathOperation:
        """Remove target run(s), if possible."""
        paths = self.define(pattern=target, branches=branches)
        return _PathOperation(self._rm, paths=paths)

    def _rm(self, target: PathLike):
        """Remove the target path only if safe to do so."""
        this = fullpath(target)
        if not this.exists():
            return PathOperationError(
                f"Cannot remove {this}: does not exist"
            )
        if not this.is_dir():
            return PathOperationError(
                f"Cannot remove {this}: not a directory"
            )
        shutil.rmtree(this)
        return this

    # TODO: Python >= 3.9
    # - list[typing.Union[pathlib.Path, tuple[pathlib.Path, ...]]]

    # TODO: Python >= 3.10
    # - list[pathlib.Path | tuple[pathlib.Path, ...]]

    ListOfPathsOrPathTuples = typing.List[
        typing.Union[pathlib.Path, typing.Tuple[pathlib.Path, ...]]
    ]

    def define(
        self,
        *names: str,
        pattern: str=None,
        branches: typing.Union[str, typing.Iterable[str]]=None,
    ) -> ListOfPathsOrPathTuples:
        """Build a list of paths or tuples of paths.
        
        Parameters
        ----------
        *names : string
            Zero or more names of runtime directories to create.
        """
        rundirs = self._get_rundirs(branches)
        if pattern:
            return [
                path for rundir in rundirs
                for path in rundir.glob(pattern)
            ]
        if len(names) == 1:
            name = names[0]
            return [rundir / name for rundir in rundirs]
        if not names:
            return self.define(pattern='*', branches=branches)
        return [tuple(rundir / name for name in names) for rundir in rundirs]

    def resolve(self, target: str, branches: typing.Iterable[str]=None):
        """Compute the full path(s) to the target run."""
        if not branches:
            return [self.root / self.base / target]
        return [
            self.root / branch / self.base / target
            for branch in branches
        ]

    @property
    def runs(self) -> typing.Dict[str, typing.Set[str]]:
        """The available runs, and owning branches, if any."""
        if not self._branches:
            return {
                run.name: set()
                for run in (self.root / self.base).glob('*')
            }
        base = {
            run.name: set()
            for branch in self._branches
            for run in (self.root / branch / self.base).glob('*')
        }
        for run, branches in base.items():
            branches.update(
                branch for branch in self._branches
                if (self.root / branch / self.base / run).is_dir()
            )
        return base

    @property
    def branches(self):
        """The project branches, if any, and their available runs."""
        if not self._branches:
            return {}
        base = {branch: set() for branch in self._branches}
        for run, branches in self.runs.items():
            for branch in branches:
                base[branch].update({run})
        return base

    def _get_rundirs(
        self,
        branches: typing.Union[str, typing.Iterable[str]]=None,
    ) -> typing.List[pathlib.Path]:
        """Build an appropriate collection of runtime directories."""
        if not branches:
            return self.listing
        subset = (
            {branches} if isinstance(branches, str)
            else set(branches)
        )
        return [d for d in self.listing if d.parent.name in subset]

    @property
    def listing(self):
        """Full paths to all runtime directories.
        
        Notes
        -----
        This property does not use the `RunPaths.branches` property to iterate
        over available branches because that would create an infinite recursion
        via `RunPaths._get_rundirs`.
        """
        if self._listing is None:
            self._listing = [
                self.root / branch / self.base
                for branch in self._branches
            ] if self._branches else [self.root / self.base]
        return self._listing

    @property
    def root(self):
        """The root directory."""
        return self._root

    @property
    def base(self):
        """The name of the common base runtime directory."""
        return self._base


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
        self._root = None
        self._log = self._get_log(attrs)
        self._attrs = attrs
        self._directories = RunPaths(attrs.path, attrs.branches, attrs.rundir)
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
        if not self.database.is_file():
            with self.database.open('w') as fp:
                json.dump({}, fp)
        with self.database.open('r') as fp:
            existing = dict(json.load(fp))
        if path.exists():
            return _ProjectInit(**existing[key])
        path.mkdir(parents=True)
        init = _ProjectInit(root=path, **kwargs)
        updated = {**existing, key: dict(init)}
        with self.database.open('w') as fp:
            json.dump(updated, fp, indent=4, sort_keys=True)
        return init

    def _get_log(self, attrs: _ProjectInit):
        """Create or retrieve the log of runs."""
        path = attrs.path / attrs.logname
        if path.exists():
            return RunLog(path)
        return RunLog(
            path,
            config=attrs.config,
            output=attrs.output,
        )

    def show(self: ProjectType, *runs: str):
        """Display information about this project or the named run(s).
        
        Parameters
        ----------
        *runs : string
            The named run(s) to display. This method also accepts `'*'`, which
            will cause it to display information about all available runs.
        """
        if not runs:
            self._show_project()
        requested = (
            tuple(self.runs) if len(runs) == 1 and runs[0] == '*'
            else runs
        )
        for run in requested:
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
        for path in self.directories.resolve(run, branches):
            print(path)

    def reset(self, force: bool=False, silent: bool=False):
        """Reset this project to its initial state."""
        return self.rm('*', errors=(not force), silent=silent)

    def rename(self, target: str, force: bool=False, silent: bool=False):
        """Rename this project to `target`."""
        # Save the current root directory since it will change.
        old = self.root
        # Convert the target to a full path.
        new = fullpath(target)
        # Handle an existing target path.
        if new.exists():
            if not force:
                raise PathOperationError(
                    f"Renaming this project to {target!r} would "
                    f"overwrite {new}."
                )
            if not silent:
                print(f"Overwriting {new}")
        # Update the runtime paths on disk.
        self._directories.update(root=new)
        # Update the log of runtime paths.
        self.log.update_directory(new)
        # Update the project database.
        self._update_database(str(old), str(new))
        # Echo success, if applicable.
        if not silent:
            if self.root.parent == old.parent == pathlib.Path.cwd():
                result = target
            else:
                result = str(self.root)
            print(f"Renamed project to {result!r}")

    def _update_database(self, source: str, target: str):
        """Rename this project's path in the project database."""
        with self.database.open('r') as fp:
            current = dict(json.load(fp))
        updated = {k: v for k, v in current.items() if k != source}
        this = current[source].copy()
        this['root'] = target
        updated[target] = this
        with self.database.open('w') as fp:
            json.dump(updated, fp, indent=4, sort_keys=True)

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
        run = name or datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S.%f')
        paths = self.directories.mkdir(run, branches=branches)
        for path in paths.errors(force=(not errors), silent=silent):
            self._create_run(
                config=config,
                path=path,
                nproc=nproc,
                environment=environment,
                silent=silent,
            )

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
                print(f"[{process.pid}]")
                print(f"Started {path.name!r} at {now}")
            process.wait()
            if not silent:
                base = f"Finished {path.name!r}"
                message = (
                    f"{base} in branch {branch.name!r}" if self.branches
                    else base
                )
                print(message)
        if process.returncode == 0:
            logentry = {
                'command': command,
                'time': now,
            }
            self.log.create(str(path), logentry)
        elif not silent:
            print(
                f"WARNING: Process exited with {process.returncode}",
                end='\n\n',
            )

    def mv(
        self: ProjectType,
        source: str,
        target: str,
        branches: typing.Union[str, typing.Iterable[str]]=None,
        errors: bool=False,
        silent: bool=False,
    ) -> ProjectType:
        """Rename an existing EPREM run within this project."""
        pairs = self.directories.mv(source, target, branches=branches)
        if not pairs:
            if not silent:
                print(f"Nothing to rename for {source!r}")
            return
        for (run, new) in pairs.errors(force=(not errors), silent=silent):
            self.log.mv(run, new)
            if not silent:
                branch = self._get_branch_name(run)
                base = f"Renamed {source!r} to {target!r}"
                print(f"{base} in branch {branch!r}" if branch else base)

    def rm(
        self: ProjectType,
        run: str,
        branches: typing.Union[str, typing.Iterable[str]]=None,
        errors: bool=False,
        silent: bool=False,
    ) -> ProjectType:
        """Remove an existing EPREM run from this project."""
        paths = self.directories.rm(run, branches=branches)
        if not paths:
            if not silent:
                print(f"Nothing to remove for {run!r}")
            return
        for path in paths.errors(force=(not errors), silent=silent):
            self.log.rm(path)
            if not silent:
                base = f"Removed {path.name!r}"
                branch = self._get_branch_name(path)
                print(f"{base} from branch {branch!r}" if branch else base)

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
        
        This property is an alias for `~Project.root.name`.
        """
        if self._name is None:
            self._name = self.root.name
        return self._name

    @property
    def runs(self):
        """The available runs, and owning branches, if any."""
        return self.directories.runs

    @property
    def branches(self):
        """The project branches, if any, and their available runs."""
        return self.directories.branches

    @property
    def root(self):
        """The top-level directory of this project."""
        return self.directories.root

    @property
    def directories(self):
        """The full path to each run directory."""
        return self._directories

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


def doc2help(__x) -> None:
    """Convert a function docstring to CLI help text."""
    try:
        target = __x if isinstance(__x, str) else str(__x.__doc__)
    except AttributeError:
        raise TypeError(f"Cannot create help text from {__x!r}") from None
    doclines = target.lstrip('\n').split('\n')
    summary = doclines[0].rstrip('.')
    return summary[0].lower() + summary[1:]


class CLI(typing.Mapping):
    """A custom command-line parser for EPREM simulations."""

    def __init__(self, *args, **kwargs):
        """Create a new instance."""
        self._args = args
        self._kwargs = kwargs
        self._parser = argparse.ArgumentParser(*self._args, **self._kwargs)
        self._subparsers = None
        self._subcommands = None

    def __len__(self) -> int:
        """Called for len(self)."""
        return len(self.subcommands)

    def __iter__(self) -> typing.Iterator[str]:
        """Called for iter(self)."""
        return iter(self.subcommands)

    def __getitem__(self, __k: str) -> argparse.ArgumentParser:
        """Access subcommands by key."""
        if __k in self.subcommands:
            return self.subcommands[__k]
        raise KeyError(f"No subcommand for {__k!r}")

    def include(self, _func=None, **meta):
        """Register a subcommand."""
        def cli_action(func: types.FunctionType):
            """Decorate `func` as a command-line action."""
            key = func.__name__
            subparser = self.subparsers.add_parser(
                key,
                help=doc2help(func),
                description=doc2help(func),
                formatter_class=argparse.RawTextHelpFormatter,
                **meta
            )
            subparser.set_defaults(func=func)
            self.subcommands[key] = subparser
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                return func(*args, **kwargs)
            return wrapper
        if _func is None:
            return cli_action
        return cli_action(_func)

    def run(self):
        """Execute operations based on command-like arguments."""
        parsed = vars(self.parser.parse_args())
        func = parsed.pop('func')
        func(**parsed)

    @property
    def parser(self):
        """The main argument parser."""
        return self._parser

    @property
    def subparsers(self):
        """The parser for each subcommand."""
        if self._subparsers is None:
            self._subparsers = self.parser.add_subparsers(
                title="supported sub-commands",
            )
        return self._subparsers

    @property
    def subcommands(self) -> typing.Dict[str, argparse.ArgumentParser]:
        """The registered subcommands."""
        if self._subcommands is None:
            self._subcommands = {}
        return self._subcommands


cli = CLI(
    formatter_class=argparse.RawTextHelpFormatter,
    description=doc2help(__doc__),
)


@cli.include
def create(
    path: PathLike,
    branches: typing.Union[str, typing.Iterable[str]]=None,
    config: str=None,
    output: str=None,
    rundir: str=None,
    logname: str=None,
    verbose: bool=False,
) -> None:
    """Create a new project."""
    project = Project(
        path,
        branches=branches,
        config=config,
        output=output,
        rundir=rundir,
        logname=logname,
    )
    if verbose:
        parts = [
            f"Created project {project.name!r} in",
            f"{project.root}",
        ]
        single = ' '.join(parts)
        message = single if len(single) < 70 else '\n'.join(parts)
        print(message)
cli.subcommands['create'].add_argument(
    '-b',
    '--branches',
    help="names of project branches\n(default: none)",
    nargs='*',
    metavar=('A', 'B'),
)
cli.subcommands['create'].add_argument(
    '-c',
    '--config',
    help=(
        "the name to give runtime config files"
        "\n(default: 'eprem.cfg')"
    ),
)
cli.subcommands['create'].add_argument(
    '-o',
    '--output',
    help=(
        "the name to give runtime output logs"
        "\n(default: 'eprem.log')"
    ),
)
cli.subcommands['create'].add_argument(
    '-d',
    '--rundir',
    help=(
        "the name of the directory that will contain simulation runs"
        "\n(default: 'runs')"
    ),
)
cli.subcommands['create'].add_argument(
    '-l',
    '--logname',
    help=(
        "the name of the file to which to log run metadata"
        "\n(default: 'runs.json')"
    ),
)


@cli.include
def reset(
    path: PathLike,
    force: bool=False,
    verbose: bool=False,
) -> None:
    """Reset an existing project."""
    project = Project(path)
    project.reset(force=force, silent=(not verbose))
cli.subcommands['reset'].add_argument(
    '-f',
    '--force',
    help="do not raise exceptions when resetting a project",
    action='store_true',
)


@cli.include
def rename(
    path: PathLike,
    target: PathLike,
    force: bool=False,
    verbose: bool=False,
) -> None:
    """Rename an existing EPREM project."""
    project = Project(path)
    project.rename(target, force=force, silent=(not verbose))
cli.subcommands['rename'].add_argument(
    'target',
    help="the new name of the project",
)
cli.subcommands['rename'].add_argument(
    '-f',
    '--force',
    help="do not raise exceptions when renaming a project",
    action='store_true',
)


@cli.include
def remove(
    path: PathLike,
    force: bool=False,
    verbose: bool=False,
) -> None:
    """Remove an existing EPREM project."""
    project = Project(path)
    project.remove(force=force, silent=(not verbose))
cli.subcommands['remove'].add_argument(
    '-f',
    '--force',
    help="do not raise exceptions when a removing project",
    action='store_true',
)


@cli.include
def run(
    path: PathLike,
    config: PathLike,
    target: str=None,
    branches: typing.Union[str, typing.Iterable[str]]=None,
    nproc: int=None,
    force: bool=False,
    verbose: bool=False,
    **environment
) -> None:
    """Set up and execute a new EPREM run in all applicable branches."""
    project = Project(path)
    project.run(
        config=config,
        name=target,
        branches=branches,
        nproc=nproc,
        environment=environment,
        errors=(not force),
        silent=(not verbose),
    )
cli.subcommands['run'].add_argument(
    'config',
    help="path to the EPREM config file to use",
)
cli.subcommands['run'].add_argument(
    '-t',
    '--target',
    help="the name to the new run\n(default: created from date and time)",
)
cli.subcommands['run'].add_argument(
    '-b',
    '--branches',
    help="the affected branches\n(default: all)",
    nargs='*',
    metavar=('A', 'B'),
)
cli.subcommands['run'].add_argument(
    '-n',
    '--nproc',
    help="the number of parallel processes to use\n(default: 1)",
    type=int,
    default=1,
    metavar='N',
)
cli.subcommands['run'].add_argument(
    '-f',
    '--force',
    help="do not raise exceptions when creating a new run",
    action='store_true',
)
cli.subcommands['run'].add_argument(
    '-m',
    '--mpirun',
    help="path to the MPI binary to use\n(default: $PATH)",
)
cli.subcommands['run'].add_argument(
    '-e',
    '--eprem',
    help="path to the EPREM executable to run\n(default: $PATH)",
)


@cli.include
def mv(
    path: PathLike,
    source: str,
    target: str,
    branches: typing.Union[str, typing.Iterable[str]]=None,
    force: bool=False,
    verbose: bool=False,
) -> None:
    """Rename an existing EPREM run in all applicable branches."""
    project = Project(path)
    project.mv(
        source,
        target,
        branches=branches,
        errors=(not force),
        silent=(not verbose),
    )
cli.subcommands['mv'].add_argument(
    'source',
    help="the current name of the run(s)",
)
cli.subcommands['mv'].add_argument(
    'target',
    help="the new name of the run(s)",
)
cli.subcommands['mv'].add_argument(
    '-b',
    '--branches',
    help="the affected branches\n(default: all)",
    nargs='*',
    metavar=('A', 'B'),
)
cli.subcommands['mv'].add_argument(
    '-f',
    '--force',
    help="do not raise exceptions when renaming a file",
    action='store_true',
)


@cli.include
def rm(
    path: PathLike,
    target: str,
    branches: typing.Union[str, typing.Iterable[str]]=None,
    force: bool=False,
    verbose: bool=False,
) -> None:
    """Remove an existing EPREM run in all applicable branches."""
    project = Project(path)
    project.rm(
        target,
        branches=branches,
        errors=(not force),
        silent=(not verbose),
    )
cli.subcommands['rm'].add_argument(
    'target',
    help="the name of the run(s) to remove",
)
cli.subcommands['rm'].add_argument(
    '-b',
    '--branches',
    help="the affected branches\n(default: all)",
    nargs='*',
    metavar=('A', 'B'),
)
cli.subcommands['rm'].add_argument(
    '-f',
    '--force',
    help="do not raise exceptions when removing a file",
    action='store_true',
)


@cli.include
def show(
    path: PathLike,
    *runs: str,
    **extra # HACK: Ignore erroneous keywords.
) -> None:
    """Display information about an existing project."""
    project = Project(path)
    project.show(*runs)
showruns = cli.subcommands['show'].add_mutually_exclusive_group()
showruns.add_argument(
    '-r',
    '--run',
    dest='runs',
    help=(
        "display information about the named run"
        ";\nspecify multiple times to include multiple runs"
    ),
    action='extend',
)
showruns.add_argument(
    '-a',
    '--all',
    dest='runs',
    help="display information for all runs",
    action='store_const',
    const='*',
)


if __name__ == "__main__":
    cli.parser.add_argument(
        '-p',
        '--project',
        dest='path',
        required=True,
        help="the path to the target project; may be relative",
    )
    cli.parser.add_argument(
        '-v',
        '--verbose',
        help="print runtime messages",
        action='store_true',
    )
    cli.run()


