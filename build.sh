#!/bin/bash

# Adapted from https://stackoverflow.com/a/3352015. This will trim leading and
# trailing whitespace from a given string. For example: `trim "  a b  "` yields
# `"a b"`.
trim() {
    local var="$*"
    # leading whitespace
    var="${var#"${var%%[![:space:]]*}"}"
    # trailing whitespace
    var="${var%"${var##*[![:space:]]}"}"
    printf '%s' "$var"
}

# Set up the function that prints a section header.
DRY_RUN=
print_banner() {
    if [ $verbose == 1 ]; then
        printf "
=======================================================================
    ${DRY_RUN}$@
=======================================================================
"
    fi
}

# See https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
# - e => Exit immediately if a pipeline returns non-zero status.
# - u => Treat unset variables and parameters (with certain exceptions) as
#   errors when performing parameter expansion.
set -eu

# Define text formatting commands.
# - textbf: Use bold-face.
# - textnm: Reset to normal.
# - startul: Start underlined text.
# - endul: End underlined text.
textbf=$(tput bold)
textnm=$(tput sgr0)
startul=$(tput smul)
endul=$(tput rmul)

# Store the initial working directory.
top_dir="$(dirname "$(readlink -f "${0}")")"
logfile=${top_dir}/build.log

# Set option defaults.
verbose=0
alias=
from_clean=0

# This is the CLI's main help text.
show_help()
{
    cli_name=${0##*/}
    echo "
${textbf}NAME${textnm}
        $cli_name - Run the EPREM build process

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli_name${textnm} [${startul}OPTION${endul}]

${textbf}DESCRIPTION${textnm}
        This script will compile and link a given EPREM build. It assumes that
        you have run install.sh or taken equivalent configuration steps.

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print runtime messages.
        ${textbf}--alias=ALIAS${textnm}
                The alias of the configuration to build.
        ${textbf}--from-clean${textnm}
                Run make clean in the target src directory before building the executable.
"
}

# This will run if the CLI gets an unrecognized option.
report_bad_arg()
{
    printf "\nUnrecognized command: ${1}\n\n"
}

# Parse command-line options. See
# `/usr/share/doc/util-linux/examples/getopt-example.bash` and the `getopt`
# manual page for more information.
TEMP=$(getopt \
    -n 'setup.sh' \
    -o 'hv' \
    -l 'help,verbose,' \
    -l 'alias:' \
    -l 'from-clean' \
    -- "$@")

if [ $? -ne 0 ]; then
    echo "Failed to parse command-line arguments. Terminating."
    exit 1
fi

eval set -- "$TEMP"
unset TEMP

while [ $# -gt 0 ]; do
    case "$1" in
        '-h'|'--help')
            show_help
            exit
        ;;
        '-v'|'--verbose')
            verbose=1
            shift
            continue
        ;;
        '--alias')
            alias="${2}"
            shift 2
            continue
        ;;
        '--from-clean')
            from_clean=1
            shift
            continue
        ;;
        '--')
            shift
            break
        ;;
    esac
done

# Collect extra arguments. These may indicate that the user incorrectly entered
# an argument.
extra_args="$@"
if [ -n "${extra_args}" ]; then
    echo
    echo "Found the following unrecognized argument(s):"
    echo
    for arg in ${extra_args[@]}; do
        echo ">   $arg"
    done
    echo
    echo "Did you misspell something?"
    echo
    echo "Note that this program does not support passing arbitrary arguments"
    echo "to make. If you want (and know how to) build the code in a way that"
    echo " this program does not provide, you may directly call"
    echo "make [clean,install,...] with appropriate arguments."
    echo
    exit 1
fi

# Create a fresh build log. We do this after reading and checking command-like
# arguments because in order to avoid unnecessarily creating a new log.
> $logfile

# Declare a status flag.
status=

# Define special flag to indicate success.
success="@!SUCCESS!@"

# Define a clean-up function.
cleanup() {
    if [ "$status" != "$success" ]; then
        echo
        echo "Build failed. See $logfile for details."
        exit 1
    else
        print_banner "Done"
    fi
}

trap cleanup EXIT

# Move to the target directory.
pushd $top_dir/$alias 1> /dev/null 2>> $logfile

# Run the clean stage, if requested.
if [ ${from_clean} == 1 ]; then
    make clean &>> $logfile
fi

# Run the make stage.
make &>> $logfile

# Run the install stage.
make install &>> $logfile

# Exit the target directory.
popd 1> /dev/null 2>> $logfile

# Set the status flag to indicate success.
status=$success

