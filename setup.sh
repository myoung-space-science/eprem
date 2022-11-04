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
        print_banner "Done"
        return $?
    fi
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
        ${textbf}$cli_name${textnm} [${startul}OPTION${endul}] -- ${startul}ARGS${endul}

${textbf}DESCRIPTION${textnm}
        This script is designed to provide a single command for quickly (albeit 
        prescriptively) setting up EPREM from source code. Assuming that you 
        are running in the top-level build directory, it will execute the 
        following steps to configure, build, and (optionally) install EPREM:

        $ autoreconf --install --symlink [--verbose]
        $ ./configure ${startul}ARGS${endul} [--silent]
        $ make
        $ [make install]

        Note that any arguments following '--' will pass to ./configure and 
        will override internally set arguments (e.g., compiler or pre-processor 
        flags set as a result of passing ${textbf}--debug${textnm} or ${textbf}--optimize${textnm}).

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
        ${textbf}--debug${textnm}
                Build EPREM for debugging. Specifically, this will pass the 
                '-g' compiler flag (by modifying the environment variable 
                ${startul}CFLAGS${endul}) and the '-DDEBUG' pre-processor 
                directive (by modifying the environment variable 
                ${startul}CPPFLAGS${endul}) to ./configure. The default 
                behavior is to not modify the environment variables.
        ${textbf}--optimize${textnm}
                Build EPREM for debugging. Specifically, this will pass the 
                '-O3' compiler flag (by modifying the environment variable 
                ${startul}CFLAGS${endul}) and the '-DNDEBUG' pre-processor 
                directive (by modifying the environment variable 
                ${startul}CPPFLAGS${endul}) to ./configure. The default 
                behavior is to not modify the environment variables.
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
debug=0
optimize=0

# Parse command-line options. See
# `/usr/share/doc/util-linux/examples/getopt-example.bash` and the `getopt`
# manual page for more information.
TEMP=$(getopt \
    -n 'install.sh' \
    -o '+hv' \
    -l 'help,verbose,install::,debug,optimize,dry-run' \
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
        '--debug')
            debug=1
            shift
            continue
        ;;
        '--optimize')
            optimize=1
            shift
            continue
        ;;
        '--')
            shift
            break
        ;;
    esac
done

# Store all un-parsed arguments to pass to ./configure.
CF_ARGS="$@"

# Update flags based on --debug or --optimize. Note that the default behavior is
# to not include environment variables in `CF_FLAGS`, rather than to include
# environment variables with existing values (e.g., CFLAGS=${CFLAGS}). This is
# intended to leave as much control as possible with the user.
if [ $debug == 1 ]; then
    SH_CFLAGS="-g ${CFLAGS}"
    SH_CPPFLAGS="-DDEBUG ${CPPFLAGS}"
    CF_FLAGS="CFLAGS=${SH_CFLAGS} CPPFLAGS=${SH_CPPFLAGS} ${CF_ARGS}"
else
    if [ $optimize == 1 ]; then
        SH_CFLAGS="-O3 ${CFLAGS}"
        SH_CPPFLAGS="-DNDEBUG ${CPPFLAGS}"
        CF_FLAGS="CFLAGS=${SH_CFLAGS} CPPFLAGS=${SH_CPPFLAGS} ${CF_ARGS}"
    else
        CF_FLAGS="${CF_ARGS}"
    fi
fi

if [ ${dry_run} == 1 ]; then
    DRY_RUN="[DRY RUN] "
fi

AR_FLAGS="--install --symlink"
if [ ${verbose} == 1 ]; then
    AR_FLAGS="${AR_FLAGS} --verbose"
else
    CF_FLAGS="${CF_FLAGS} --silent"
fi

run_stage "autoreconf ${AR_FLAGS}"

if [ $install_opt == 1 ] && [ "x${install_dir}" != x ]; then
    CF_FLAGS="${CF_FLAGS} --bindir=${install_dir}"
fi

run_stage "./configure ${CF_FLAGS}"

run_stage "make"

if [ $install_opt == 1 ]; then
    run_stage "make install"
fi
