#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#-----------------------------------------------------------------------------#
# Initialize.
#-----------------------------------------------------------------------------#
AC_PREREQ([2.71])
AC_INIT([EPREM],[0.1.0],[Matthew.Young@unh.edu],[eprem])

#-----------------------------------------------------------------------------#
# Perform global set-up tasks.
  AS_BOX([Environment and System Setup])
#-----------------------------------------------------------------------------#

AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([src/eprem.c])
AC_CONFIG_MACRO_DIRS([m4])

#-----------------------------------------------------------------------------#
# Use `autodefs.h` instead of the default, `config.h`, because the latter is too
# close to the EPREM header file `configuration.h`.
#-----------------------------------------------------------------------------#
AC_CONFIG_HEADERS([src/autodefs.h])

#-----------------------------------------------------------------------------#
# Identify the variables containing the MPI compilers as "precious".
#-----------------------------------------------------------------------------#
AC_ARG_VAR(MPICC, [MPI C compiler command])
AC_ARG_VAR(MPICXX, [MPI C++ compiler command])

#-----------------------------------------------------------------------------#
# Determinte the canonical host-system (sets $host).
#-----------------------------------------------------------------------------#
AC_CANONICAL_HOST

#-----------------------------------------------------------------------------#
# Check for particular standard programs.
#-----------------------------------------------------------------------------#
AC_PROG_INSTALL
AC_PROG_AWK
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P
AC_PROG_RANLIB

#-----------------------------------------------------------------------------#
# Read command-line options.
#-----------------------------------------------------------------------------#
mpi_dir=
AC_ARG_WITH(
    [mpi-dir],
    AS_HELP_STRING([--with-mpi-dir@<:@=DIR@:>@], [Look for MPI header files in DIR/include and libraries in DIR/lib.]),
    [mpi_dir=$withval]
)
AS_IF(
    [test "x$mpi_dir" != x],
    [
        AM_CFLAGS="  -I$mpi_dir/include $AM_CFLAGS"
        AM_CXXFLAGS="-I$mpi_dir/include $AM_CXXFLAGS"
        AM_LDFLAGS=" -L$mpi_dir/lib -Wl,-rpath,$mpi_dir/lib $AM_LDFLAGS"
    ]
)
libconfig_dir=
AC_ARG_WITH(
    [libconfig-dir],
    AS_HELP_STRING([--with-libconfig-dir@<:@=DIR@:>@], [Look for libconfig header files in DIR/include and libraries in DIR/lib.]),
    [libconfig_dir=$withval]
)
AS_IF(
    [test "x$libconfig_dir" != x],
    [
        AM_CFLAGS="  -I$libconfig_dir/include $AM_CFLAGS"
        AM_CXXFLAGS="-I$libconfig_dir/include $AM_CXXFLAGS"
        AM_LDFLAGS=" -L$libconfig_dir/lib -Wl,-rpath,$libconfig_dir/lib $AM_LDFLAGS"
    ]
)
netcdf_dir=
AC_ARG_WITH(
    [netcdf-dir],
    AS_HELP_STRING([--with-netcdf-dir@<:@=DIR@:>@], [Look for netcdf header files in DIR/include and libraries in DIR/lib.]),
    [netcdf_dir=$withval]
)
AS_IF(
    [test "x$netcdf_dir" != x],
    [
        AM_CFLAGS="  -I$netcdf_dir/include $AM_CFLAGS"
        AM_CXXFLAGS="-I$netcdf_dir/include $AM_CXXFLAGS"
        AM_LDFLAGS=" -L$netcdf_dir/lib -Wl,-rpath,$netcdf_dir/lib $AM_LDFLAGS"
    ]
)

#-----------------------------------------------------------------------------#
# Augment environment variable values with Automake values, to ensure proper
# configuration during subsequent tests. See
# https://www.gnu.org/software/automake/manual/html_node/Program-Variables.html
# We will reset these variables to their original values at the end.
#
# Note: This is not entirely consistent with the Autotools philosophy, which
# prefers to respect the user's environment as much as possible. A more
# philosophically consistent approach would leverage `config.site` files or the
# `CONFIG_SITE` variable. However, this script assumes that the typical user
# knows very little about the build process, environment variables, etc., and
# wants to install the program as quickly as possible.
#-----------------------------------------------------------------------------#
: ${CFLAGS="${AM_CFLAGS} ${CFLAGS}"}
: ${CXXFLAGS="${AM_CXXFLAGS} ${CXXFLAGS}"}
: ${CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS}"}
: ${LDFLAGS="${AM_LDFLAGS} ${LDFLAGS}"}

#-----------------------------------------------------------------------------#
# Set up compilers.
  AS_BOX([Compiler Setup])
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Define the MPI command for MPI library tests.
#-----------------------------------------------------------------------------#
AC_MSG_CHECKING([setting MPIRUN])
AS_CASE(
    [$host],
    [*-aix*], [MPIRUN="aprun -n 2"],
    [MPIRUN="mpirun -n 2"]
)
AC_MSG_RESULT([${MPIRUN}])
AC_SUBST([MPIRUN])

#-----------------------------------------------------------------------------#
# Warn the user if CC is set to a non-MPI compiler.
#-----------------------------------------------------------------------------#
m4_define([mpi_cc_list], [mpicc mpiicc mpixlc mpipgcc])
AS_IF(
    [test x"$CC" != x],
    [
        found_mpi_cc=
        m4_define(
            [check_cc],
            [AS_IF([test "x$CC" = "x$1"], [found_mpi_cc=$CC; break;])]
        )
        m4_map([check_cc], [mpi_cc_list])
        AS_IF(
            [test "x$found_mpi_cc" = x],
            [AC_MSG_ERROR([EPREM requires an MPI C compiler (got CC=$CC)])]
        )
    ]
)

#-----------------------------------------------------------------------------#
# Set the C compiler to a suitable MPI C compiler, if necessary.
#-----------------------------------------------------------------------------#
AC_PROG_CC([mpicc mpiicc mpixlc mpipgcc])

#-----------------------------------------------------------------------------#
# Warn the user if CXX is set to a non-MPI compiler.
#-----------------------------------------------------------------------------#
m4_define([mpi_cxx_list], [mpicxx mpiCC mpic++ mpig++ mpiicpc mpipgCC mpixlC])
AS_IF(
    [test x"$CXX" != x],
    [
        found_mpi_cxx=
        m4_define(
            [check_cxx],
            [AS_IF([test "x$CXX" = "x$1"], [found_mpi_cxx=$CXX; break;])]
        )
        m4_map([check_cxx], [mpi_cxx_list])
        AS_IF(
            [test "x$found_mpi_cxx" = x],
            [AC_MSG_ERROR([EPREM requires an MPI C++ compiler (got CXX=$CXX)])]
        )
    ]
)

#-----------------------------------------------------------------------------#
# Set the C++ compiler to a suitable MPI C++ compiler, if necessary.
#-----------------------------------------------------------------------------#
AC_PROG_CXX([mpicxx mpiCC mpic++ mpig++ mpiicpc mpipgCC mpixlC])

#-----------------------------------------------------------------------------#
# Check for typedefs, structures, and compiler characteristics.
  AS_BOX([Typedefs, Structures, and Compiler Characteristics])
#
# NOTE: At this point in the script, `AC_TYPE_SIZE_T` will trigger checks for
# numerous other standard header files (e.g., stdio.h); this is preferable for
# the sake of organizing runtime messages.
#-----------------------------------------------------------------------------#
AC_TYPE_SIZE_T
AC_C_RESTRICT
AC_CHECK_TYPES([ptrdiff_t])

#-----------------------------------------------------------------------------#
# Check for required system libraries and header files.
  AS_BOX([Required System Libraries and Header Files])
#
# math.h (standard math functions): If the header file exists on the system,
# this script will additionally check that the math library defines the exp()
# function, which EPREM requires.
#-----------------------------------------------------------------------------#
AC_CHECK_HEADERS(
    math.h,
    [AC_CHECK_LIB(m, exp, [], AC_MSG_ERROR([EPREM requires the exp() function]))],
    AC_MSG_ERROR([Can't find math.h.])
)

#-----------------------------------------------------------------------------#
# Make sure we have the z compression library; it is a dependency of NetCDF.
#-----------------------------------------------------------------------------#
AC_SEARCH_LIBS(
    [deflate],
    [z],
    [],
    AC_MSG_ERROR([Can't find the z library (required for NetCDF).])
)
AC_CHECK_HEADERS(
    zlib.h,
    [],
    AC_MSG_ERROR([Can't find zlib.h.])
)

#-----------------------------------------------------------------------------#
# Search for a suitable MPI library. We currently require MPI for parallel runs
# but it should become optional if we implement a serial version.
#-----------------------------------------------------------------------------#
AC_SEARCH_LIBS(
    [MPI_Init],
    [mpi mpich],
    [],
    AC_MSG_ERROR([Can't find a suitable MPI library.])
)
AC_CHECK_HEADERS(
    [mpi.h],
    [],
    AC_MSG_ERROR([Can't find MPI headers. Try using --with-mpi-dir.])
)

#-----------------------------------------------------------------------------#
# Search for libconfig, which we require for processing input.
#-----------------------------------------------------------------------------#
AC_SEARCH_LIBS(
    [config_init],
    [config],
    [],
    AC_MSG_ERROR([Can't find a suitable version of the libconfig library.])
)
AC_CHECK_HEADERS(
    [libconfig.h],
    [],
    AC_MSG_ERROR([Can't find libconfig headers. Try using --with-libconfig-dir.])
)

#-----------------------------------------------------------------------------#
# Search for the NetCDF library, which we require for creating output.
#-----------------------------------------------------------------------------#
AC_SEARCH_LIBS(
    [nc_create],
    [netcdf],
    [],
    AC_MSG_ERROR([Can't find a suitable version of the NetCDF library.])
)
AC_CHECK_HEADERS(
    [netcdf.h],
    [],
    AC_MSG_ERROR([Can't find NetCDF4 headers. Try using --with-netcdf-dir.])
)

#-----------------------------------------------------------------------------#
# Check for library functions.
  AS_BOX([Required Library Functions])
#-----------------------------------------------------------------------------#
AC_FUNC_MALLOC
AC_CHECK_FUNCS([exp floor memset pow sqrt strstr])

#-----------------------------------------------------------------------------#
# Special tests, -lstdc++ needed, or -DMPICH_IGNORE_CXX_SEEK (Darwin issues?).
  AS_BOX([Special Tests])
#
# NOTE: I adapted these from an earlier version of this script. I don't know if
# they're strictly necessary for EPREM.
#-----------------------------------------------------------------------------#
AC_LANG_PUSH(C++)
AC_MSG_CHECKING([if empty compile works])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include <iostream>]])],
    AC_MSG_RESULT([yes (PASS)]),
    [
        AC_MSG_RESULT([no (FAIL)])
        AC_MSG_ERROR([try rerunning with LIBS=-lstdc++])
    ]
)
AC_MSG_CHECKING([if -DMPICH_IGNORE_CXX_SEEK is needed])
AC_LINK_IFELSE(
    [
        AC_LANG_PROGRAM(
            [[
            #include <iostream>
            #include <mpi.h>
            ]]
        )
    ],
    AC_MSG_RESULT([no (PASS)]),
    [
        AC_MSG_RESULT([yes (FAIL)])
        AC_MSG_ERROR([try rerunning with CPPFLAGS=-DMPICH_IGNORE_CXX_SEEK])
    ]
)
AC_LANG_POP(C++)

#-----------------------------------------------------------------------------#
# Finalize and create output.
  AS_BOX([Output])
#-----------------------------------------------------------------------------#
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
