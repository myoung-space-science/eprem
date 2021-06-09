import argparse
from pathlib import Path
import subprocess
from typing import *


class Logger:
    """Log events and print them, if requested."""
    def __init__(self, filename: str=None) -> None:
        try:
            self.filename = Path(filename).expanduser().resolve()
        except TypeError:
            self.filename = None
        self._fp = None
        self._logging = False

    @property
    def start(self) -> None:
        """Start logging."""
        self._fp = self.filename.open('w') if self.filename else None
        self._logging = True

    def log(self, message: str, **kwargs):
        """Log a single user-defined event."""
        if self._logging:
            print(message, file=self._fp, **kwargs)

    @property
    def stop(self) -> None:
        """Stop logging."""
        if self._fp:
            self._fp.close()
        self._logging = False

    def __bool__(self) -> bool:
        """True if the logger is currently logging."""
        return self._logging

    def __enter__(self) -> 'Logger':
        """Enter a logging context."""
        self.start
        return self

    def __exit__(self, _type, _value, _traceback) -> bool:
        """Exit the current logging context."""
        self.stop
        return False


class BuildError(Exception):
    """There was an error while building EPREM."""
    def __init__(self, stage: str) -> None:
        self.stage = stage

    def __str__(self) -> str:
        return f"Build failed during '{self.stage}'"


class BuildStep:
    """A class to represent a step in the build process."""
    def __init__(
        self,
        name: str,
        description: str,
        method: Callable[..., bool],
        *args,
        **kwargs,
    ) -> None:
        self.name = name
        self.description = description
        self._method = method
        self._args = args,
        self._kwargs = kwargs

    @property
    def succeeded(self):
        """"""
        return self._method()


class EPREMBuilder:
    """A class to manage EPREM building commands."""
    def __init__(
        self,
        target: str=None,
        ext_deps: str=None,
        mode: str=None,
    ) -> None:
        self._target = make_path(target, '.')
        self._build = self._target / 'source' / 'eprem'
        self._ext_deps = make_path(ext_deps, '~') / 'deps'
        self._mode = mode
        self._logger = Logger('build.log')
        self._iteration = None
        self.steps = ['prepare', 'configure', 'make']
        self._recipes = {
            'prepare': {
                'method': self._prepare,
                'description': "Preparing dependencies",
            },
            'configure': {
                'method': self._configure,
                'description': "Configuring EPREM",
                # 'deps': ['libconfig', 'hdf4', 'netcdf'],
            },
            'make': {
                'method': self._make,
                'description': "Building EPREM",
            },
        }

    @property
    def _deps(self) -> List[str]:
        """The dependencies for this mode."""
        _all = {
            'MAS': ['config', 'hdf4', 'netcdf'],
            'ENLIL': ['config', 'netcdf'],
        }
        return _all.get(self._mode, ['config', 'netcdf'])

    def __enter__(self) -> 'EPREMBuilder':
        """Enter the building context."""
        self._logger.start
        return self

    def __exit__(self, _type, _value, _traceback) -> bool:
        """Exit the building context."""
        self._logger.stop
        return False

    def __iter__(self) -> Iterator:
        """Iterate over steps in the build process."""
        self._iteration = 0
        return self

    def __next__(self) -> BuildStep:
        """Produce the next step in the process."""
        try:
            step = self.steps[self._iteration]
            self._iteration += 1
            return BuildStep(step, **self._recipes[step])
        except IndexError:
            raise StopIteration

    def _run(self, *args, **kwargs) -> subprocess.CompletedProcess:
        """Run a process and capture output."""
        kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return subprocess.run(*args, **kwargs)

    def _prepare(self) -> bool:
        """Prepare the dependencies."""
        result = self._run('./prepare', cwd=self._build)
        self._logger.log(result)
        return result.returncode == 0

    def _configure(self) -> bool:
        """Configure the dependencies."""
        args = ['./configure', 'CC=mpicc', 'CXX=mpic++']
        if 'config' in self._deps:
            args.append(f"--with-config={self._ext_deps / 'libconfig'}")
        if 'hdf4' in self._deps:
            args.append(f"--with-hdf4={self._ext_deps / 'hdf4'}")
        if 'netcdf' in self._deps:
            args.append(f"--with-netcdf={self._ext_deps / 'netcdf'}")
        result = self._run(args, cwd=self._build)
        self._logger.log(result)
        return result.returncode == 0

    def _make(self) -> bool:
        """Create an executable from the source code."""
        reffile = Path().cwd() / '.buildref'
        reffile.touch()
        src = self._build / 'src'
        result = self._run(["make", "clean"], cwd=src)
        self._logger.log(result)
        result = self._run(["make", "-j", "4"], cwd=src)
        self._logger.log(result)
        exe = src / 'eprem'
        success = exe.stat().st_ctime > reffile.stat().st_ctime
        reffile.unlink()
        return success


class BuildRunner:
    """"""
    def __init__(self, logger: Logger) -> None:
        self._status_logger = logger

    def build(self, builder: EPREMBuilder, verbose: bool=False):
        """"""
        if verbose:
            self._status_logger.start
        with builder as b:
            for step in b:
                self._status_logger.log(
                    step.description, end=' ... ', flush=True
                )
                if step.succeeded:
                    self._status_logger.log("Succeeded")
                else:
                    self._status_logger.log("Failed")
                    raise BuildError(step.name) from None
        if self._status_logger:
            self._status_logger.stop


def main(
    target: str=None,
    mode: str=None,
    ext_deps: str=None,
    logfile: str=None,
    verbose: bool=False,
) -> None:
    """Build an EPREM distribution."""
    logger = Logger(logfile)
    eprem = EPREMBuilder(target=target, ext_deps=ext_deps, mode=mode)
    runner = BuildRunner(logger)
    runner.build(eprem, verbose=verbose)


def make_path(arg: str, default: str=None) -> Optional[Path]:
    """Convert a string into a full path, or return a default path.

    This function will attempt to create a fully-qualified path from `arg`. If
    `arg` is `None` but `default` is not `None`, this function will use
    `default`. If both are `None`, this function will return `None`.
    """
    if default:
        return Path(arg or default).expanduser().resolve()


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description=main.__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument(
        '-t',
        '--target',
        help="the EPREM distribution to build (default: current directory)"
    )
    p.add_argument(
        '-m',
        '--mode',
        help="the mode of the target EPREM distribution",
        choices=('MAS', 'ENLIL', 'default'),
        default='default',
    )
    p.add_argument(
        '-e',
        '--ext-deps',
        dest='ext_deps',
        help="path to installed external dependencies (default: user home)",
    )
    p.add_argument(
        '-l',
        '--logfile',
        help="the name of log file (default: print messages to stdout)"
    )
    p.add_argument(
        '-v',
        '--verbose',
        help="print runtime messages (default: off)",
        action='store_true',
    )
    args = p.parse_args()
    main(**vars(args))

