#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#-----------------------------------------------------------------------------#
# Initialize and perform setup tasks.
#-----------------------------------------------------------------------------#
AC_PREREQ([2.71])
AC_INIT([EPREM],[0.1.0],[Matthew.Young@unh.edu],[eprem])
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
# Determinte the canonical host-system (sets $host).
#-----------------------------------------------------------------------------#
AC_CANONICAL_HOST

#-----------------------------------------------------------------------------#
# Define the MPI command for MPI library tests.
#-----------------------------------------------------------------------------#
AC_MSG_CHECKING([setting MPIRUN])
case ${host} in
    *-aix*)  MPIRUN="aprun -n 2";;
    *)       MPIRUN="mpirun -n 2";;
esac
AC_MSG_RESULT([${MPIRUN}])
AC_SUBST([MPIRUN])

#-----------------------------------------------------------------------------#
# Preserve the user's environment variables.
#-----------------------------------------------------------------------------#
user_CC="$CC"
user_CXX="$CXX"
user_CFLAGS="$CFLAGS"
user_CXXFLAGS="$CXXFLAGS"
user_CPPFLAGS="$CPPFLAGS"
user_LDFLAGS="$LDFLAGS"

#-----------------------------------------------------------------------------#
# Suppress the default compiler flags (namely, -g -O2). See documentation for
# `AC_PROG_CC` and `AC_PROG_CXX` for more information.
#-----------------------------------------------------------------------------#
: ${CFLAGS=""}

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
        AM_CPPFLAGS="-I$mpi_dir/include $AM_CPPFLAGS"
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
        AM_CPPFLAGS="-I$libconfig_dir/include $AM_CPPFLAGS"
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
        AM_CPPFLAGS="-I$netcdf_dir/include $AM_CPPFLAGS"
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
CFLAGS="${AM_CFLAGS} ${CFLAGS}"
CXXFLAGS="${AM_CXXFLAGS} ${CXXFLAGS}"
CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS}"
LDFLAGS="${AM_LDFLAGS} ${LDFLAGS}"

#-----------------------------------------------------------------------------#
# Set the C compiler to a suitable MPI C compiler, if necessary.
#-----------------------------------------------------------------------------#
MPICC=$user_CC
AS_IF(
    [test "x$MPICC" = x],
    [
        AC_CHECK_PROGS(MPICC, [mpicc mpiicc mpixlc mpipgcc], [not-found])
        AC_MSG_WARN([Setting CC=$MPICC for parallel execution.])
    ]
)
AC_ARG_VAR(MPICC, [MPI C compiler command])
AC_PROG_CC([$MPICC])

#-----------------------------------------------------------------------------#
# Set the C++ compiler to a suitable MPI C++ compiler, if necessary.
#-----------------------------------------------------------------------------#
MPICXX=$user_CXX
AS_IF(
    [test "x$MPICXX" = x],
    [
        AC_CHECK_PROGS(MPICXX, [mpicxx mpiCC mpic++ mpig++ mpiicpc mpipgCC mpixlC])
        AC_MSG_WARN([Setting CXX=$MPICXX for parallel execution.])
    ]
)
AC_ARG_VAR(MPICXX, [MPI C++ compiler command])
AC_PROG_CXX([$MPICXX])

#-----------------------------------------------------------------------------#
# Make sure we have the C math library. This explicitly checks for the existence
# of `exp`, which EPREM requires, but it also ensures that we have access to the
# full library.
#-----------------------------------------------------------------------------#
AC_SEARCH_LIBS(
    [exp],
    [m],
    [],
    AC_MSG_ERROR([Can't find the math library.])
)
AC_CHECK_HEADERS(
    math.h,
    [],
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
# Checks for typedefs, structures, and compiler characteristics.
#-----------------------------------------------------------------------------#
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

#-----------------------------------------------------------------------------#
# Checks for library functions.
#-----------------------------------------------------------------------------#
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset pow sqrt strstr])


#-----------------------------------------------------------------------------#
# Special tests, -lstdc++ needed, or -DMPICH_IGNORE_CXX_SEEK (darwin issues?).
#
# DEV: I adapted this (almost verbatim) from the `prediccs` version. There may
# be more canonical ways to perform the same tests, if they're even necessary.
#-----------------------------------------------------------------------------#
AC_LANG(C++)

AC_MSG_CHECKING([if empty compile works])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include <iostream>]])],
    [last_test=0; AC_MSG_RESULT([yes])],
    [last_test=1; AC_MSG_RESULT([no])]
)

AS_IF(
    [test "$last_test" -ne 0],
    [
        AC_SUBST(orig_LDFLAGS,[$LDFLAGS])
        AC_SUBST(LDFLAGS,["$LDFLAGS -lstdc++"])
        AC_MSG_CHECKING([again with -lstdc++])
        AC_LINK_IFELSE(
            [AC_LANG_PROGRAM([[#include <iostream>]])],
            [last_test=0; AC_MSG_RESULT([yes])],
            [last_test=1; AC_MSG_RESULT([no]); AC_SUBST(LDFLAGS,[$orig_LDFLAGS])]
        )
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
    [last_test=0],
    [last_test=1]
)

AS_IF(
    [test "$last_test" -ne 0],
    [
        AC_SUBST(orig_CPPFLAGS,[$CPPFLAGS])
        AC_SUBST(CPPFLAGS,["$CPPFLAGS -DMPICH_IGNORE_CXX_SEEK"])
        AC_LINK_IFELSE(
            [
                AC_LANG_PROGRAM(
                    [[
                    #include <iostream>
                    #include <mpi.h>
                    ]]
                )
            ],
            [AC_MSG_RESULT([yes])],
            [AC_MSG_RESULT([error]); AC_SUBST(CPPFLAGS,[$orig_CPPFLAGS])]
        )
    ],
    [AC_MSG_RESULT([no])]
)

#-----------------------------------------------------------------------------#
# Restore the user values of environment variables.
#-----------------------------------------------------------------------------#
CFLAGS="$user_CFLAGS"
CXXFLAGS="$user_CXXFLAGS"
CPPFLAGS="$user_CPPFLAGS"
LDFLAGS="$user_LDFLAGS"

#-----------------------------------------------------------------------------#
# Finalize and create output.
#-----------------------------------------------------------------------------#
AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_LDFLAGS)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
