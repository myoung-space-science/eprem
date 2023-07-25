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

run_stage() {
    print_banner "$@"
    if [ $dry_run != 1 ]; then
        eval "$@"
        print_banner "Done"
        return $?
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

# This is the CLI's main help text.
show_help()
{
    cli_name=${0##*/}
    echo "
${textbf}NAME${textnm}
        $cli_name - Run the full build process

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli_name${textnm} [${startul}OPTION${endul}]

${textbf}DESCRIPTION${textnm}
        This script is designed to provide a single command for quickly (albeit 
        prescriptively) setting up EPREM from source code. Assuming that you 
        are running in the top-level build directory, it will execute the 
        following steps to configure, build, and (optionally) install EPREM:

        $ autoreconf --install --symlink [--verbose]
        $ ./configure [...] [--silent]
        $ make
        $ [make install]

        where the arguments to ./configure depend on the selected options.

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
        ${textbf}--download-ext-deps${textnm}
                Download external dependencies. By default, this will download
                and build packages inside ./external. However, setting 
                --with-ext-deps=DIR will cause it to download and build packages 
                inside DIR. This script will create the necessary path in either
                case.
        ${textbf}--with-libconfig-dir=DIR${textnm}
                Look for libconfig header files in DIR/include and look for 
                libraries in DIR/lib. Supersedes the --with-ext-deps option.
        ${textbf}--with-netcdf-dir=DIR${textnm}
                Look for netcdf header files in DIR/include and look for 
                libraries in DIR/lib. Supersedes the --with-ext-deps option.
        ${textbf}--install[=DIR]${textnm}
                Install the executable. If DIR is included, this will directly 
                install the executable in DIR; if not, it will install it to 
                the default location for the host system or PREFIX, if set (see 
                --prefix). The default action is to not install the executable. 
                Including DIR is the equivalent of passing --bindir=DIR to 
                ./configure. Using this in combination with --prefix=PREFIX 
                will install the executable in PREFIX/DIR.
        ${textbf}--prefix=PREFIX${textnm}
                Set the top-level installation directory (default: /usr/local). 
                This is equivalent to passing --prefix=PREFIX to ./configure. 
                Using this in combination with --install=DIR will install the 
                executable in PREFIX/DIR; using this in combination with 
                --install will install the executable in PREFIX/bin. Note that 
                using this option alone will not trigger installation.
        ${textbf}--debug${textnm}
                Build a debugging version of EPREM. Specifically, this will 
                pass the '-g' compiler flag and the '-DDEBUG' pre-processor 
                directive to ./configure. You should consider using the 
                --prefix option with this option in order to keep track of 
                different builds.
        ${textbf}--optimize${textnm}
                Build an optimized version of EPREM. Specifically, this will 
                pass the '-O3' compiler flag and the '-DNDEBUG' pre-processor 
                directive to ./configure. You should consider using the 
                --prefix option with this option in order to keep track of 
                different builds.
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
prefix=
debug=0
optimize=0
dry_run=0
mpi_dir=
libconfig_dir=
netcdf_dir=
deps_dir=
download_deps=0

# Parse command-line options. See
# `/usr/share/doc/util-linux/examples/getopt-example.bash` and the `getopt`
# manual page for more information.
TEMP=$(getopt \
    -n 'setup.sh' \
    -o 'hv' \
    -l 'help,verbose,' \
    -l 'install::,prefix:' \
    -l 'debug,optimize,dry-run' \
    -l 'with-mpi-dir:,with-libconfig-dir:,with-netcdf-dir:' \
    -l 'with-ext-deps:' \
    -l 'download-ext-deps' \
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
        '--prefix')
            prefix="${2}"
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
        '--dry-run')
            dry_run=1
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
        '--download-ext-deps')
            download_deps=1
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
    echo "to configure or make. If you want (and know how to) build the code"
    echo "in a way that this program does not provide, you may directly call"
    echo "./configure [...] && make && make install"
    echo
    exit 1
fi

# The function that will download and build external dependencies.
install_ext_deps() {
    if [ -z "${1}" ]; then
        echo "ERROR: Missing path to external-dependencies directory."
        exit 1
    fi
    local deps_dir="${1}"
    local tmp_name=tmp

    # Create and enter the top-level directory for external dependencies.
    mkdir -p ${deps_dir}
    pushd ${deps_dir} &> /dev/null

    # Create and enter the local subdirectory where this script will download
    # and build each package. Isolating the source code allows us to remove it
    # while leaving the build dependencies.
    mkdir ${tmp_name}
    pushd ${tmp_name} &> /dev/null

    # TODO: Refactor into a loop.

    # --> libconfig
    wget https://src.fedoraproject.org/repo/pkgs/libconfig/libconfig-1.5.tar.gz/a939c4990d74e6fc1ee62be05716f633/libconfig-1.5.tar.gz
    tar -xvzf libconfig-1.5.tar.gz
    pushd libconfig-1.5 &> /dev/null
    ./configure --prefix="${deps_dir}" && make && make install && make check
    popd &> /dev/null

    # --> zlib
    wget https://zlib.net/fossils/zlib-1.2.11.tar.gz
    tar -xvzf zlib-1.2.11.tar.gz
    pushd zlib-1.2.11 &> /dev/null
    ./configure --prefix="${deps_dir}" && make && make install && make check
    popd &> /dev/null

    # --> HDF5
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.17/src/hdf5-1.8.17.tar.gz
    tar -xvzf hdf5-1.8.17.tar.gz
    pushd hdf5-1.8.17 &> /dev/null
    ./configure --prefix="${deps_dir}" && make && make install && make check
    popd &> /dev/null

    # --> NetCDF4
    wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.4.1.1.tar.gz
    tar -xvzf v4.4.1.1.tar.gz
    pushd netcdf-c-4.4.1.1 &> /dev/null
    ./configure --prefix="${deps_dir}" --disable-dap-remote-tests && make && make install && make check
    popd &> /dev/null

    # Exit and remove the temporary source-code directory.
    popd &> /dev/null
    /bin/rm -rf ${tmp_name}

    # Exit the top-level external-dependencies directory.
    popd &> /dev/null
}

# Process the --download-ext-deps option.
if [ "$download_deps" == 1 ]; then
    if [ -z "$deps_dir" ]; then
        deps_dir="$(pwd)/ext-deps"
    fi
    install_ext_deps "${deps_dir}"
fi

# Update flags based on --with-ext-deps option.
if [ -n "$deps_dir" ]; then
    if [ -z "$libconfig_dir" ]; then
        libconfig_dir=${deps_dir}/libconfig
    fi
    if [ -z "$netcdf_dir" ]; then
        netcdf_dir=${deps_dir}/netcdf
    fi
fi

# Initialize temporary variables for --with-<package>-dir options.
SH_CFLAGS=
SH_CXXFLAGS=
SH_CPPFLAGS=
SH_LDFLAGS=

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

# Set any unset environment variables to null.
${CFLAGS:=}
${CXXFLAGS:=}
${CPPFLAGS:=}
${LDFLAGS:=}

# Initialize a temporary variable for ./configure options.
CF_FLAGS=

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

# Update configure flags based on installation options.
if [ $install_opt == 1 ] && [ "x${install_dir}" != x ]; then
    CF_FLAGS="${CF_FLAGS} --bindir=${install_dir}"
fi
if [ -n "$prefix" ]; then
    CF_FLAGS="${CF_FLAGS} --prefix=${prefix}"
fi

# Run configure stage.
run_stage "./configure ${CF_FLAGS}"

# Run make stage.
run_stage "make"

# Run make install stage, if necessary.
if [ $install_opt == 1 ]; then
    run_stage "make install"
fi
