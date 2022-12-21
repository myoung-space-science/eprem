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
        ${textbf}--with-mpi-dir=DIR${textnm}
                Look for MPI header files in DIR/include and look for MPI 
                libraries in DIR/lib.
        ${textbf}--with-ext-deps=DIR${textnm}
                Look for external dependencies in DIR/include and DIR/lib.
        ${textbf}--with-libconfig-dir=DIR${textnm}
                Look for libconfig header files in DIR/include and look for 
                libraries in DIR/lib. Overrides --with-ext-deps option.
        ${textbf}--with-netcdf-dir=DIR${textnm}
                Look for netcdf header files in DIR/include and look for 
                libraries in DIR/lib. Overrides --with-ext-deps option.
        ${textbf}--install[=DIR]${textnm}
                Install the executable. If DIR is included, this will directly 
                install the executable in DIR; if not, it will install it to 
                the default location for the host system. The default action is 
                to not install the executable. Including DIR is the equivalent 
                of passing --bindir=DIR to ./configure.
        ${textbf}--debug${textnm}
                Build EPREM for debugging. Specifically, this will pass the 
                '-g' compiler flag (by modifying the environment variable 
                ${startul}CFLAGS${endul}) and the '-DDEBUG' pre-processor 
                directive (by modifying the environment variable 
                ${startul}CPPFLAGS${endul}) to ./configure. The default 
                behavior is to not modify the environment variables.
        ${textbf}--optimize${textnm}
                Build EPREM for production. Specifically, this will pass the 
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
mpi_dir=
libconfig_dir=
netcdf_dir=
deps_dir=

# Parse command-line options. See
# `/usr/share/doc/util-linux/examples/getopt-example.bash` and the `getopt`
# manual page for more information.
TEMP=$(getopt \
    -n 'setup.sh' \
    -o '+hv' \
    -l 'help,verbose,install::,debug,optimize,dry-run' \
    -l 'with-mpi-dir:,with-libconfig-dir:,with-netcdf-dir:' \
    -l 'with-ext-deps:' \
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
        '--with-mpi-dir')
            mpi_dir="${2}"
            shift 2
            continue
        ;;
        '--with-libconfig-dir')
            libconfig_dir="${2}"
            shift 2
            continue
        ;;
        '--with-netcdf-dir')
            netcdf_dir="${2}"
            shift 2
            continue
        ;;
        '--with-ext-deps')
            deps_dir="${2}"
            shift 2
            continue
        ;;
        '--')
            shift
            break
        ;;
    esac
done

# Store all un-parsed arguments to pass to ./configure. Note that the default
# behavior when checking command-line options will be to not include environment
# variables in `CF_FLAGS`, rather than to include environment variables with
# existing values (e.g., CFLAGS=${CFLAGS}). This is intended to leave as much
# control as possible with the user.
CF_ARGS="$@"

# Update flags based on --with-ext-deps option.
if [ -n "$deps_dir" ]; then
    if [ -z "$libconfig_dir" ]; then
        libconfig_dir=${deps_dir}/libconfig
    fi
    if [ -z "$netcdf_dir" ]; then
        netcdf_dir=${deps_dir}/netcdf
    fi
fi

# Update flags based on --with-<package>-dir options.
if [ -n "$mpi_dir" ]; then
    SH_CFLAGS="  -I${mpi_dir}/include $SH_CFLAGS"
    SH_CXXFLAGS="-I${mpi_dir}/include $SH_CXXFLAGS"
    SH_CPPFLAGS="-I${mpi_dir}/include $SH_CPPFLAGS"
    SH_LDFLAGS=" -L${mpi_dir}/lib -Wl,-rpath,${mpi_dir}/lib $SH_LDFLAGS"
fi
if [ -n "$libconfig_dir" ]; then
    SH_CFLAGS="  -I${libconfig_dir}/include $SH_CFLAGS"
    SH_CXXFLAGS="-I${libconfig_dir}/include $SH_CXXFLAGS"
    SH_CPPFLAGS="-I${libconfig_dir}/include $SH_CPPFLAGS"
    SH_LDFLAGS=" -L${libconfig_dir}/lib -Wl,-rpath,${libconfig_dir}/lib $SH_LDFLAGS"
fi
if [ -n "${netcdf_dir}" ]; then
    SH_CFLAGS="  -I${netcdf_dir}/include $SH_CFLAGS"
    SH_CXXFLAGS="-I${netcdf_dir}/include $SH_CXXFLAGS"
    SH_CPPFLAGS="-I${netcdf_dir}/include $SH_CPPFLAGS"
    SH_LDFLAGS=" -L${netcdf_dir}/lib -Wl,-rpath,${netcdf_dir}/lib $SH_LDFLAGS"
fi

# Update flags based on --debug or --optimize options.
if [ $debug == 1 ]; then
    SH_CFLAGS="-g ${SH_CFLAGS}"
    SH_CXXFLAGS="-g ${SH_CXXFLAGS}"
    SH_CPPFLAGS="-DDEBUG ${SH_CPPFLAGS}"
else
    if [ $optimize == 1 ]; then
        SH_CFLAGS="-O3 ${SH_CFLAGS}"
        SH_CXXFLAGS="-O3 ${SH_CXXFLAGS}"
        SH_CPPFLAGS="-DNDEBUG ${SH_CPPFLAGS}"
    fi
fi

# Collect all non-empty flags.
SH_CFLAGS=$(trim "${SH_CFLAGS} ${CFLAGS}")
if [ -n "${SH_CFLAGS}" ]; then
    CF_FLAGS+="CFLAGS=\"${SH_CFLAGS}\" "
fi
SH_CXXFLAGS=$(trim "${SH_CXXFLAGS} ${CXXFLAGS}")
if [ -n "${SH_CXXFLAGS}" ]; then
    CF_FLAGS+="CXXFLAGS=\"${SH_CXXFLAGS}\" "
fi
SH_CPPFLAGS=$(trim "${SH_CPPFLAGS} ${CPPFLAGS}")
if [ -n "${SH_CPPFLAGS}" ]; then
    CF_FLAGS+="CPPFLAGS=\"${SH_CPPFLAGS}\" "
fi
SH_LDFLAGS=$(trim "${SH_LDFLAGS} ${LDFLAGS}")
if [ -n "${SH_LDFLAGS}" ]; then
    CF_FLAGS+="LDFLAGS=\"${SH_LDFLAGS}\" "
fi

# Check for --dry-run option.
if [ ${dry_run} == 1 ]; then
    DRY_RUN="[DRY RUN] "
fi

# Check for --verbose option.
AR_FLAGS="--install --symlink"
if [ ${verbose} == 1 ]; then
    AR_FLAGS="${AR_FLAGS} --verbose"
else
    CF_FLAGS="${CF_FLAGS} --silent"
fi

# Run Autotools setup stage, if necessary.
if [ -d "./build-aux" ]; then
    print_banner "Autotools Setup"
    echo "Found ./build-aux directory. Not running autoreconf."
else
    run_stage "autoreconf ${AR_FLAGS}"
fi

# Update configure flags based on --install option.
if [ $install_opt == 1 ] && [ "x${install_dir}" != x ]; then
    CF_FLAGS="${CF_FLAGS} --bindir=${install_dir}"
fi

# Augment local configure arguments with user-given arguments at the last
# moment, to ensure that the user's values take precedence.
CF_FLAGS=$(trim "${CF_FLAGS} ${CF_ARGS}")

# Run configure stage.
run_stage "./configure ${CF_FLAGS}"

# Run make stage.
run_stage "make"

# Run make install stage, if necessary.
if [ $install_opt == 1 ]; then
    run_stage "make install"
fi
