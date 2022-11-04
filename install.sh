#!/bin/bash

print_banner() {
    if [ $verbose == 1 ]; then
        printf "
=======================================================================
    ${DRY_RUN}$@
=======================================================================
"
    fi
}

run_stage() {
    print_banner "$@"
    if [ $dry_run != 1 ]; then
        eval "$@"
        return $?
    fi
    print_banner "Done"
}

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
        ${textbf}$cli_name${textnm} [${startul}OPTION${endul}] -- ${startul}CFARGS${endul}

${textbf}DESCRIPTION${textnm}
        This script is designed to provide a single command for quickly (albeit 
        prescriptively) setting up EPREM from source code. Assuming that you 
        are running in the top-level build directory, it will execute the 
        following steps to configure, build, and (optionally) install EPREM:

        $ autoreconf --install --symlink [--verbose]
        $ ./configure ${startul}CFARGS${endul} [--silent]
        $ make
        $ [make install]

        Note that any arguments following '--' will pass to ./configure.

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print runtime messages. Use of this option will enable the 
                --verbose argument to autoreconf and will disable the --silent 
                option to ./configure. By default, --verbose is disabled for 
                autoreconf and --silent is enabled for ./configure.
        ${textbf}--install[=DIR]${textnm}
                Install the executable. If DIR is included, this will directly 
                install the executable in DIR; if not, it will install it to 
                the default location for the host system. The default action is 
                to not install the executable.
        ${textbf}--dry-run${textnm}
                Display the sequence of commands but don't run anything.
"
}

# This will run if the CLI gets an unrecognized option.
report_bad_arg()
{
    printf "\nUnrecognized command: ${1}\n\n"
}

# Set option defaults.
verbose=0
install_opt=0
dry_run=0

# Parse command-line options. See
# `/usr/share/doc/util-linux/examples/getopt-example.bash` and the `getopt`
# manual page for more information.
TEMP=$(getopt \
    -n 'install.sh' \
    -o 'hv' \
    -l 'help,verbose,install::,dry-run' \
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
        '--dry-run')
            dry_run=1
            shift
            continue
        ;;
        '--install')
            install_opt=1
            case "$2" in
                '')
                    install_dir=
                ;;
                *)
                    install_dir="${2}"
            esac
            shift 2
            continue
        ;;
        '--')
            shift
            break
        ;;
    esac
done

cf_flags="$@"

if [ ${dry_run} == 1 ]; then
    DRY_RUN="[DRY RUN] "
fi

ar_flags="--install --symlink"
if [ ${verbose} == 1 ]; then
    ar_flags="${ar_flags} --verbose"
else
    cf_flags="${cf_flags} --silent"
fi

run_stage "autoreconf ${ar_flags}"

if [ $install_opt == 1 ] && [ "x${install_dir}" != x ]; then
    cf_flags="${cf_flags} --bindir=${install_dir}"
fi

run_stage "./configure ${cf_flags}"

run_stage "make"

if [ $install_opt == 1 ]; then
    run_stage "make install"
fi
