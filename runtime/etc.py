"""Support for the EPREM runtime interface.
"""

import argparse
import functools
import os
import pathlib
import shutil
import types
import typing


PathLike = typing.Union[str, os.PathLike]


def fullpath(p: PathLike):
    """Expand and resolve the given path."""
    return pathlib.Path(p).expanduser().resolve()


def locate(
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


