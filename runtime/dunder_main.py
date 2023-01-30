"""
Create and modify runs in this EPREM project.
"""

import argparse
import typing

import project


HOME = project.fullpath(__file__).parent


cli = project.CLI(
    formatter_class=argparse.RawTextHelpFormatter,
    description=project.doc2help(__doc__),
)


@cli.include
def run(
    config: str,
    target: str=None,
    branches: typing.Union[str, typing.Iterable[str]]=None,
    nproc: int=None,
    force: bool=False,
    verbose: bool=False,
    **environment
) -> None:
    """Set up and execute a new EPREM run in all applicable branches."""
    prj = project.Project(HOME)
    prj.run(
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
cli.subcommands['run'].add_argument(
    '-v',
    '--verbose',
    help="print runtime messages",
    action='store_true',
)


@cli.include
def mv(
    source: str,
    target: str,
    branches: typing.Union[str, typing.Iterable[str]]=None,
    force: bool=False,
    verbose: bool=False,
) -> None:
    """Rename an existing EPREM run in all applicable branches."""
    prj = project.Project(HOME)
    prj.mv(
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
cli.subcommands['mv'].add_argument(
    '-v',
    '--verbose',
    help="print runtime messages",
    action='store_true',
)


@cli.include
def rm(
    target: str,
    branches: typing.Union[str, typing.Iterable[str]]=None,
    force: bool=False,
    verbose: bool=False,
) -> None:
    """Remove an existing EPREM run in all applicable branches."""
    prj = project.Project(HOME)
    prj.rm(
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
cli.subcommands['rm'].add_argument(
    '-v',
    '--verbose',
    help="print runtime messages",
    action='store_true',
)


if __name__ == "__main__":
    cli.run()
