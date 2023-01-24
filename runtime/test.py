import argparse
import datetime
import pathlib
import traceback
import typing
import sys

_RUNTIME = '~/emmrem/open/source/eprem/runtime'
sys.path.append(str(pathlib.Path(_RUNTIME).expanduser()))

import eprem


class Context: # Should this inherit from `eprem.Project`?
    """Context manager for EPREM runtime API tests."""

    def __init__(
        self,
        root: str,
        config: str,
        branches: typing.Union[str, typing.Iterable[str]]=None,
        interactive: bool=False,
        keep: bool=False,
        verbosity: int=0,
    ) -> None:
        self.root = eprem.fullpath(root)
        self.branches = branches
        self.config = eprem.fullpath(config)
        self.interactive = interactive
        self.keep = keep
        self.verbosity = verbosity
        self.project = None
        self.step = None
        self.print_status("Initialized EPREM test context")

    @property
    def silent(self):
        """Suppress runtime messages."""
        return not self.verbose

    @property
    def verbose(self):
        """Print runtime messages."""
        return self.verbosity > 0

    def __enter__(self):
        """Set up the context."""
        self.print_status("Starting EPREM runtime API tests")
        return self

    def __exit__(
        self,
        exc_type: typing.Type[Exception],
        exc_value: Exception,
        exc_tb,
    ) -> typing.Optional[typing.Literal[True]]:
        """Tear down the context."""
        self.keep or self.remove()
        return self._handle(exc_type, exc_value, exc_tb)

    def _handle(
        self,
        exc_type: typing.Type[Exception],
        exc_value: Exception,
        exc_tb,
    ) -> bool:
        """Catch and handle an exception, if necessary."""
        if exc_type:
            self.print_error(exc_type, exc_value, exc_tb)
            return True
        self.print_status("Tests finished normally")
        return False

    def print_error(
        self,
        exc_type: typing.Type[Exception],
        exc_value: Exception,
        exc_tb,
    ) -> None:
        """Print an error message, if necessary."""
        if self.verbose:
            print(self._format_exception(exc_type, exc_value, exc_tb))
            self.print_status(f"Tests ended early due to error")

    def _format_exception(
        self,
        exc_type: typing.Type[Exception],
        exc_value: Exception,
        exc_tb,
    ) -> str:
        """Prepare an exception for printing."""
        base = f"\nCaught {exc_type.__qualname__}"
        if self.verbosity <= 1 or not str(exc_value):
            # NOTE: We'll probably never get here when self.verbosity < 1 (i.e.,
            # self.verbosity == 0), but this condition provides completeness.
            return f"{base}\n"
        if self.verbosity == 2:
            return f"{base}:\n{exc_value}\n"
        # NOTE: If we're here, self.verbosity > 2.
        exc = traceback.format_exception(exc_type, exc_value, exc_tb)
        joined = ''.join(exc)
        return f"\n{joined}"

    def create(self, name: str):
        self.print_stage("create the project")
        self.project = eprem.Project(self.root / name, branches=self.branches)
        if self.verbose:
            print(f"Created project {name!r}\nin {self.root}")
        return self

    def show(self, *runs: str):
        """Show a summary of the project or the named runs."""
        target = 'run summaries' if runs else 'project summary'
        stage = f"display {target}"
        self.print_stage(stage)
        message = f"to {stage}"
        if self.prompted(message):
            self.project.show(*runs)

    def run(self, target):
        """Create a test run."""
        stage = f"create {target!r}"
        self.print_stage(stage)
        here = self.branch_string()
        base = f"to {stage}"
        prompt = f"{base} in {here}" if here else base
        if self.prompted(prompt):
            self.project.run(self.config, target, nproc=2, silent=self.silent)

    def mv(
        self,
        old: str,
        new: str,
        branches: typing.Union[str, typing.Iterable[str]]=None,
    ) -> None:
        """Rename test runs."""
        stage = f"rename {old!r} to {new!r}"
        self.print_stage(stage)
        here = self.branch_string(branches)
        base = f"to {stage}"
        prompt = f"{base} in {here}" if here else base
        if self.prompted(prompt):
            self.project.mv(old, new, branches=branches, silent=self.silent)

    def rm(
        self,
        target: str,
        branches: typing.Union[str, typing.Iterable[str]]=None,
    ) -> None:
        """Remove test runs."""
        stage = f"remove {target!r}"
        self.print_stage(stage)
        here = self.branch_string(branches)
        base = f"to {stage}"
        prompt = f"{base} from {here}" if here else base
        if self.prompted(prompt):
            self.project.rm(target, branches=branches, silent=self.silent)

    def reset(self):
        """Reset the test project."""
        stage = "reset the project"
        self.print_stage(stage)
        if self.prompted(f"to {stage}"):
            self.project.reset(silent=self.silent)

    def rename(self, target: str):
        """Rename the test project."""
        stage = "rename the project"
        self.print_stage(stage)
        if self.prompted(f"to {stage} to {target!r}"):
            self.project.rename(target, silent=self.silent)

    def remove(self):
        """Remove the test project."""
        stage = "remove the project"
        self.print_stage(stage)
        if self.prompted(f"to {stage}"):
            self.project.remove(silent=self.silent)

    def prompted(self, this: str):
        """Interact with the user, if necessary."""
        if not self.interactive:
            return True
        try:
            reply = input(f"Do you want {this}? [Y/n]: ")
        except KeyboardInterrupt as exc:
            # Add a newline to avoid cluttering output.
            end = '\n'
            raise exc
        else:
            end = ''
            return reply.lower() not in {'n', 'no'}
        finally:
            print(end)

    def print_status(self, message: str):
        """Print a status message, if necessary."""
        if self.verbose:
            print(f"\n*** {message} ***")

    def print_stage(self, message: str):
        """Print the current test stage, if necessary."""
        if self.step == None:
            self.step = 0
        if self.verbose:
            full = f"  Step {self.step}: {message}  "
            length = len(full)
            line = '=' * length
            print(f"\n{line}")
            print(f"{full}")
            print(f"{line}")
            self.step += 1

    def branch_string(self, branches=None):
        """Build a string representing the relevant branches."""
        if not self.branches:
            return
        if not branches:
            return "all branches"
        if isinstance(branches, str):
            return f"branch {branches!r}"
        return f"branches {', '.join(f'{branch!r}' for branch in branches)}"


def main(
    config: str,
    name: str=None,
    path: str=None,
    **kwargs
) -> None:
    "Test the EPREM runtime interface."
    time = datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
    context = Context(eprem.fullpath(path or '.'), config, **kwargs)
    with context.create(name or f'project_{time}') as tests:
        execute(tests)


def execute(context: Context):
    """Execute the series of tests within this context."""
    for target in ('run00', 'run01'):
        context.run(target)
    context.show()
    renamed = {
        ('run00', 'run-0'): None,
        ('run01', 'runA1'): 'A',
        ('run01', 'runB1'): 'B',
    }
    for (old, new), branches in renamed.items():
        context.mv(old, new, branches)
    context.show('*')
    removed = [
        ('runA1', None),
        ('runB1', None),
        ('run-0', 'A'),
    ]
    for (target, branches) in removed:
        context.rm(target, branches)
    context.rename(context.project.name.replace('project', 'renamed'))
    context.reset()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=eprem.doc2help(main),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        'config',
        help="path to the config file to use; may be relative",
    )
    parser.add_argument(
        '-p',
        '--path',
        help=(
            "path at which to create the test project"
            "\n(default: current directory)"
        ),
    )
    parser.add_argument(
        '-n',
        '--name',
        help=(
            "name of the test project to create"
            "\n(default: project_<date-and-time>)"
        ),
    )
    parser.add_argument(
        '-i',
        '--interactive',
        help="wait for user confirmation before each stage",
        action='store_true',
    )
    parser.add_argument(
        '-k',
        '--keep',
        help=(
            "do not remove the project"
            ", or prompt for removal in interactive mode"
        ),
        action='store_true',
    )
    parser.add_argument(
        '-v',
        '--verbose',
        dest='verbosity',
        help="print runtime messages; repeat to increase verbosity",
        action='count',
        default=0,
    )
    cliargs = vars(parser.parse_args())
    main(**cliargs)
    main(**cliargs, branches=['A', 'B'])
