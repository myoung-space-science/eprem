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

# Store the initial working directory.
script_dir="$(dirname "$(readlink -f "${0}")")"
setuplog=${script_dir}/setup.log

# Set option defaults.
verbose=0
dry_run=0
default_deps_dir=${script_dir}/deps
deps_dir=
download_deps=0

# This is the CLI's main help text.
show_help()
{
    cli_name=${0##*/}
    echo "
${textbf}NAME${textnm}
        $cli_name - Set up the EPREM build process

${textbf}SYNOPSIS${textnm}
        ${textbf}$cli_name${textnm} [${startul}OPTION${endul}]

${textbf}DESCRIPTION${textnm}
        This script will run set-up tasks necessary for building EPREM.

        ${textbf}-h${textnm}, ${textbf}--help${textnm}
                Display help and exit.
        ${textbf}-v${textnm}, ${textbf}--verbose${textnm}
                Print runtime messages.
        ${textbf}--with-deps-dir=DIR${textnm}
                Install dependencies in DIR when --download-deps is present.
        ${textbf}--download-deps${textnm}
                Download external dependencies. By default, this will download
                and build packages inside ${default_deps_dir}. However, setting 
                --with-deps-dir=DIR will cause it to download and build packages 
                inside DIR. This script will create the necessary path in either
                case.
        ${textbf}--dry-run${textnm}
                Display the sequence of commands but don't run anything.
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
    -l 'dry-run' \
    -l 'with-deps-dir:' \
    -l 'download-deps' \
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
        '--with-deps-dir')
            deps_dir="$(realpath -m "${2}")"
            shift 2
            continue
        ;;
        '--download-deps')
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
    exit 1
fi

# Create a fresh setup log. We do this after reading and checking command-like
# arguments because in order to avoid unnecessarily creating a new log.
> $setuplog

# Declare a status flag.
status=

# Define special flag to indicate success.
success="@!SUCCESS!@"

# Define a clean-up function.
cleanup() {
    if [ "$status" != "$success" ]; then
        echo
        echo "Set up failed. See $setuplog for details."
        exit 1
    else
        print_banner "Done"
    fi
}

trap cleanup EXIT

# Check for --dry-run option.
if [ ${dry_run} == 1 ]; then
    DRY_RUN="[DRY RUN] "
fi

install_libconfig() {
    local pkg_alias=libconfig
    local pkg_name=$pkg_alias-1.5
    local pkg_tar=$pkg_name.tar.gz
    local pkg_url=https://src.fedoraproject.org/repo/pkgs/libconfig/$pkg_tar/a939c4990d74e6fc1ee62be05716f633/$pkg_tar
    local pkg_args=""

    print_banner "Downloading $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        wget $pkg_url
    fi
    print_banner "Installing $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        tar -xvzf $pkg_tar &>> $setuplog
        ln -s $pkg_name $pkg_alias
        pushd $pkg_alias &> /dev/null
        ./configure $pkg_args &>> $setuplog && \
        make &>> $setuplog && \
        make check &>> $setuplog
        popd &> /dev/null
    fi
}

install_zlib() {
    local pkg_alias=zlib
    local pkg_name=$pkg_alias-1.2.11
    local pkg_tar=$pkg_name.tar.gz
    local pkg_url=https://zlib.net/fossils/$pkg_tar
    local pkg_args=""

    print_banner "Downloading $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        wget $pkg_url
    fi
    print_banner "Installing $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        tar -xvzf $pkg_tar &>> $setuplog
        ln -s $pkg_name $pkg_alias
        pushd $pkg_alias &> /dev/null
        ./configure $pkg_args &>> $setuplog && \
        make &>> $setuplog && \
        make check &>> $setuplog
        popd &> /dev/null
    fi
}

install_hdf4() {
    local pkg_alias=hdf4
    local pkg_name=hdf-4.2.16
    local pkg_tar=$pkg_name.tar.gz
    local pkg_url=https://support.hdfgroup.org/ftp/HDF/releases/HDF4.2.16/src/$pkg_tar
    local pkg_args="--disable-netcdf"

    print_banner "Downloading $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        wget $pkg_url
    fi
    print_banner "Installing $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        tar -xvzf $pkg_tar &>> $setuplog
        ln -s $pkg_name $pkg_alias
        pushd $pkg_alias &> /dev/null
        ./configure $pkg_args &>> $setuplog && \
        make &>> $setuplog && \
        make check &>> $setuplog
        popd &> /dev/null
    fi
}

install_hdf5() {
    local pkg_alias=hdf5
    local pkg_name=$pkg_alias-1.8.17
    local pkg_tar=$pkg_name.tar.gz
    local pkg_url=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.17/src/$pkg_tar
    local pkg_args=""

    print_banner "Downloading $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        wget $pkg_url
    fi
    print_banner "Installing $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        tar -xvzf $pkg_tar &>> $setuplog
        ln -s $pkg_name $pkg_alias
        pushd $pkg_alias &> /dev/null
        ./configure $pkg_args &>> $setuplog && \
        make &>> $setuplog && \
        make check &>> $setuplog
        popd &> /dev/null
    fi
}

install_netcdf() {
    local pkg_alias=netcdf
    local pkg_name=netcdf-c-4.4.1.1
    local pkg_tar=v4.4.1.1.tar.gz
    local pkg_url=https://github.com/Unidata/netcdf-c/archive/refs/tags/$pkg_tar
    local pkg_args="--disable-dap-remote-tests"

    print_banner "Downloading $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        wget $pkg_url
    fi
    print_banner "Installing $pkg_alias"
    if [ "$dry_run" == 0 ]; then
        tar -xvzf $pkg_tar &>> $setuplog
        ln -s $pkg_name $pkg_alias
        pushd $pkg_alias &> /dev/null
        ./configure $pkg_args &>> $setuplog && \
        make &>> $setuplog && \
        make check &>> $setuplog
        popd &> /dev/null
    fi
}

# The function that will download and build external dependencies.
install_ext_deps() {
    if [ -z "${1}" ]; then
        echo "ERROR: Missing path to external-dependencies directory."
        exit 1
    fi
    local top_dir="${1}"
    local tmp_dir=tmp

    # Create and enter the top-level directory for external dependencies.
    mkdir -p ${top_dir}
    pushd ${top_dir} &> /dev/null

    # Create and enter the local subdirectory where this script will download
    # and build each package. Isolating the source code allows us to remove it
    # while leaving the build dependencies.
    mkdir -p ${tmp_dir}

    # --> libconfig
    install_libconfig

    # --> zlib
    install_zlib

    # --> HDF4
    install_hdf4

    # --> HDF5
    install_hdf5

    # --> NetCDF4
    install_netcdf

    # Exit and remove the temporary source-code directory.
    /bin/rm -rf ${tmp_dir}

    # Exit the top-level external-dependencies directory.
    popd &> /dev/null
}

# Set the target directory for external dependencies.
if [ -z "$deps_dir" ]; then
    deps_dir=$default_deps_dir
fi

# Process the --download-deps option.
if [ "$download_deps" == 1 ]; then
    install_ext_deps "${deps_dir}"
fi

# Set the status flag to indicate success.
status=$success

