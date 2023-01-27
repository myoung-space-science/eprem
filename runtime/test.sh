#!/bin/bash

# Exit immediately if a pipeline returns non-zero status. See
# https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
set -e

# Define text formatting commands.
# - textbf: Use bold-face.
# - textnm: Reset to normal.
# - startul: Start underlined text.
# - endul: End underlined text.
textbf=$(tput bold)
textnm=$(tput sgr0)
startul=$(tput smul)
endul=$(tput rmul)

# This is the CLI's main help text.
show_help()
{
    cli_name=${0##*/}
    echo "
${textbf}NAME${textnm}
        $cli_name - Test the EPREM command-line runtime interface.

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli_name${textnm} [${startul}OPTION${endul}] ${startul}CONFIG${endul}

${textbf}DESCRIPTION${textnm}
        This script will test all subcommands of the EPREM runtime/project.py 
        command-line utility.

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-d=DIR${textnm}, ${textbf}--directory=DIR${textnm}
                The directory in which to create the test project
                (default: current directory).
        ${textbf}-n=NAME${textnm}, ${textbf}--name=NAME${textnm}
                The name of the test project to create
                (default: project_<date-and-time>).
        ${textbf}-k${textnm}, ${textbf}--keep${textnm}
                Do not automatically remove the project.
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print runtime messages; repeat to increase verbosity.
"
}

report_bad_args() {
    show_help
    for arg in $@; do
        echo "Unrecognized argument: ${arg}"
    done
    exit 1
}

# Parse command-line options. See
# `/usr/share/doc/util-linux/examples/getopt-example.bash` and the `getopt`
# manual page for more information.
TEMP=$(getopt \
    -n 'test.sh' \
    -o 'hd:n:kv' \
    -l 'help,directory:,name:,keep,verbose' \
    -- "$@")

if [ $? -ne 0 ]; then
    echo "Failed to parse command-line arguments. Terminating."
    exit 1
fi

eval set -- "$TEMP"
unset TEMP

# Set option defaults.
directory=$(pwd)
name=
keep=0
verbose=0
verbosity=0

# Parse optional command-line arguments.
while [ $# -gt 0 ]; do
    case "$1" in
        '-h'|'--help')
            show_help
            exit
        ;;
        '-d'|'--directory')
            directory="${2}"
            shift 2
            continue
        ;;
        '-n'|'--name')
            name="${2}"
            shift 2
            continue
        ;;
        '-i'|'--interactive')
            show_help
            echo "This script does not support the -i/--interactive option."
            echo "Please use ${textbf}test.py${textnm} for interactive testing."
            exit 1
        ;;
        '-k'|'--keep')
            keep=1
            shift
            continue
        ;;
        '-v'|'--verbosity')
            verbosity=$((verbosity+1))
            shift
            continue
        ;;
        '--')
            shift
            break
        ;;
    esac
done

# Set `verbose` boolean from `verbosity` integer.
if [ ${verbosity} -gt 0 ]; then
    verbose=1
fi

# Parse required command-line arguments.
config="${1}"
if [ -z "${config}" ]; then
    echo "Missing required EPREM config file"
    exit 1
fi

# Check for invalid arguments.
shift
if [ $# -gt 0 ]; then
    report_bad_args $@
fi

prog=~/emmrem/open/source/eprem/runtime/project.py
now="$(date +%Y-%m-%d)T$(date +%H:%M:%S)"
proj="project_${now}"

# TODO: Port test tasks from `test.py`.
