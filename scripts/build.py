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


class EPREMBuilder:
    """A class to manage EPREM building commands."""
    def __init__(
        self,
        target: str=None,
        ext_deps: str=None,
        logfile: str=None,
    ) -> None:
        self._target = make_path(target, '.')
        self._build = self._target / 'source' / 'eprem'
        self._ext_deps = make_path(ext_deps, '~')
        self._logger = Logger(logfile) if logfile else None

    def _run(self, *args, **kwargs) -> str:
        """Run a process and capture output."""
        kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return subprocess.run(*args, **kwargs)

    def _log(self, message: str, **kwargs) -> None:
        """Log the given message, if using a logger."""
        if self._logger:
            self._logger.log(message, **kwargs)

    def prepare(self) -> bool:
        """Prepare the dependencies."""
        try:
            result = self._run('./prepare', cwd=self._build)
            self._log(result)
            return True
        except subprocess.CalledProcessError:
            return False

    def configure(self, deps: Iterable[str]=None) -> bool:
        """Configure the dependencies."""
        args = ['./configure', 'CC=mpicc', 'CXX=mpic++']
        if 'config' in deps:
            args.append(f"--with-config={self._ext_deps / 'libconfig'}")
        if 'hdf4' in deps:
            args.append(f"--with-hdf4={self._ext_deps / 'hdf4'}")
        if 'netcdf' in deps:
            args.append(f"--with-netcdf={self._ext_deps / 'netcdf'}")
        try:
            result = self._run(args, cwd=self._build)
            self._log(result)
            return True
        except subprocess.CalledProcessError:
            return False

    def build(self) -> bool:
        """Build the source code."""
        reffile = Path().cwd() / '.buildref'
        reffile.touch()
        src = self._build / 'src'
        result = self._run(["make", "clean"], cwd=src)
        self._log(result)
        result = self._run(["make", "-j", "4"], cwd=src)
        self._log(result)
        exe = src / 'eprem'
        success = exe.stat().st_ctime > reffile.stat().st_ctime
        reffile.unlink()
        return success


class BuildError(Exception):
    """There was an error while building EPREM."""


def execute(
    step: Callable[..., bool],
    logger: Logger,
    preamble: str,
    *args,
    **kwargs,
) -> None:
    """Execute the given step in the build process."""
    logger.log(preamble, end=' ...')
    if step(*args, **kwargs):
        logger.log("Success!")
    else:
        logger.log("FAILED")
        logger.stop
        raise BuildError


def main(
    target: str=None,
    ext_deps: str=None,
    logfile: str=None,
    verbose: bool=False,
) -> None:
    """Build an EPREM distribution."""
    logger = Logger()
    eprem = EPREMBuilder(target=target, ext_deps=ext_deps, logfile=logfile)
    if verbose:
        logger.start
    execute(eprem.prepare, logger, "Preparing dependencies")
    execute(
        eprem.configure,
        logger,
        "Configuring EPREM",
        ['config', 'hdf4', 'netcdf'],
    )
    execute(eprem.build, logger, "Building EPREM")
    if logger:
        logger.stop



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

