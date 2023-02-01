import argparse
import datetime
import traceback
import typing

import _runtime
import etc


FILEPATH = etc.fullpath(__file__)
DIRECTORY = FILEPATH.parent


class Context:
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
        self.root = etc.fullpath(root)
        self.branches = branches
        self.config = etc.fullpath(config)
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
        self.project = _runtime.Interface(
            self.root / name,
            branches=self.branches,
        )
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

    def print_status(
        self,
        message: str,
        pad_char: str='*',
        pad_count: int=3,
        left_pad: str=None,
        right_pad: str=None,
        no_pad: bool=False,
    ) -> None:
        """Print a status message, if necessary.

        The optional parameters control the padding on either side of `message`.
        If this instance's `verbosity` is 0, this method will immediately
        return.

        Parameters
        ----------
        message : string
            The message to print.

        pad_char : string, default='*'
            The padding character to use with `pad_count`; ignored if either
            `left_pad` or `right_pad` is not ``None``.

        pad_count : int, default=3
            The number of times to repeat `pad_char` in each padding string;
            ignored if either `left_pad` or `right_pad` is not ``None``.

        left_pad : string, optional
            The string with which to pad the message on the left.

        right_pad : string, optional
            The string with which to pad the message on the right.

        no_pad : boolean, default=False
            If true, do not pad `message` before printing.

        Notes
        -----
        * When using `pad_char` and `pad_count`, this method will insert a
          single space `' '` between the left and right sets of `pad_char` and
          `message`.
        * When using `left_pad` and `right_pad`, this method will use each as
          given (i.e., without inserting `' '`).
        * When `no_pad` is true, the effect is equivalent to
          ``print_status(<message>, left_pad='', right_pad='')`` but *not*
          equivalent to ``print_status(<message>, pad_char='')`` or
          ``print_status(<message>, pad_count=0)``, which will both still
          include a single space on either side of the message.
        """
        if not self.verbose:
            return
        if no_pad:
            print(f"\n{message}\n")
            return
        if left_pad is None and right_pad is None:
            pad = pad_char * pad_count
            print(f"\n{pad} {message} {pad}")
            return
        if all(s == 'mirror' for s in (left_pad, right_pad)):
            raise ValueError(
                "Only one of left_pad or right_pad may be 'mirror'"
            ) from None
        lpad = left_pad or ''
        rpad = right_pad or ''
        if left_pad == 'mirror':
            lpad = right_pad[::-1]
        if right_pad == 'mirror':
            rpad = left_pad[::-1]
        print(f"\n{lpad}{message}{rpad}\n")

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


_TESTPRJ = 'testprj'


def main(
    directory: str=None,
    **kwargs
) -> None:
    "Test the EPREM runtime interface."
    time = datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
    context = Context(directory or '.', DIRECTORY / 'test.cfg', **kwargs)
    with context.create(f'{_TESTPRJ}_{time}') as tests:
        execute(tests)


def execute(context: Context):
    """Execute the series of tests within this context."""
    for target in ('run00', 'run01'):
        context.run(target)
    context.show()
    renamed = {
        ('run00', 'run0'): None,
        ('run01', 'run1A'): 'A',
        ('run01', 'run1B'): 'B',
    } if context.branches else {
        ('run00', 'run0'): None,
        ('run01', 'run1'): None,
    }
    for (old, new), branches in renamed.items():
        context.mv(old, new, branches)
    context.show('*')
    removed = [
        ('run1A', None),
        ('run1B', None),
        ('run0', 'A'),
    ] if context.branches else [
        ('run1', None),
    ]
    for (target, branches) in removed:
        context.rm(target, branches)
    context.rename(context.project.name.replace(_TESTPRJ, 'renamed'))
    context.reset()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=etc.doc2help(main),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        '-d',
        '--directory',
        help=(
            "directory in which to create the test project"
            "\n(default: current directory)"
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
