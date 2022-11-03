#!/bin/bash

# --- CLI template added by /home/matthew/bin/add_cli ---

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
        $cli_name - Run the full build process

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli_name${textnm} [${startul}OPTION${endul}] ${startul}ARGS${endul}

${textbf}DESCRIPTION${textnm}
        This script is designed to provide a single command for quickly (albeit prescriptively) setting up EPREM from source code. Assuming that you are running in the top-level build directory, it will execute the following steps to configure, build, and install EPREM:

        $ autoreconf --install --symlink [--verbose]
        $ ./configure ${startul}ARGS${endul} [--silent]
        $ make
        $ make install

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print runtime messages. Use of this option will enable the 
                --verbose argument to autoreconf and will disable the --silent 
                option to ./configure. By default, --verbose is disabled for 
                autoreconf and --silent is enabled for ./configure.
"
}

# This will run if the CLI gets an unrecognized option.
report_bad_arg()
{
    printf "\nUnrecognized command: ${1}\n\n"
}

# Set option defaults.
verbose=0
cf_flags=

while [ "${1}" != "" ]; do
    case "${1}" in
        -v | --verbose )
            verbose=1
            ;;
        -h | --help )
            show_help
            exit
            ;;
        * )
            cf_flags="${cf_flags} ${1}"
            ;;
    esac
    shift
done

# --- End CLI template ---


ar_flags="--install --symlink"
if [ ${verbose} == 1 ]; then
    ar_flags="${ar_flags} --verbose"
else
    cf_flags="${cf_flags} --silent"
fi

if [ ${verbose} == 1 ]; then
    echo "=================================================================="
    echo "  Setup: Running autoreconf ${ar_flags}"
    echo "=================================================================="
fi
autoreconf ${ar_flags}
if [ ${verbose} == 1 ]; then
    echo "=================================================================="
    echo "  Setup: Done"
    echo "=================================================================="
fi

if [ ${verbose} == 1 ]; then
    echo "=================================================================="
    echo "  Configure: Running ./configure ${cf_flags}"
    echo "=================================================================="
fi
./configure ${cf_flags}
if [ ${verbose} == 1 ]; then
    echo "=================================================================="
    echo "  Configure: Done"
    echo "=================================================================="
fi

if [ ${verbose} == 1 ]; then
    echo "=================================================================="
    echo "  Build: Running make"
    echo "=================================================================="
fi
make
if [ ${verbose} == 1 ]; then
    echo "=================================================================="
    echo "  Build: Done"
    echo "=================================================================="
fi

