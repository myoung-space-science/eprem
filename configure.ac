#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([EPREM], [0.1.0], [Matthew.Young@unh.edu], [eprem])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([src/eprem.c])

# Use `autodefs.h` instead of the default, `config.h`, because the latter is too
# close to the EPREM header file `configuration.h`.
AC_CONFIG_HEADERS([src/autodefs.h])

# Determinte the canonical host-system (sets $host).
AC_CANONICAL_HOST

# Define the MPI command.
#
# DEV: I adapted this from the `prediccs` version without fully understanding
# its relevance or usefulness.
AC_MSG_CHECKING([setting MPIRUN])
case ${host} in
    *-aix*)  MPIRUN="aprun -n 2";;
    *)       MPIRUN="mpirun -n 2";;
esac
AC_MSG_RESULT([${MPIRUN}])
AC_SUBST([MPIRUN])

# Save the user's environment variables.
#
# DEV: Adapted from HDF5.
user_CFLAGS="$CFLAGS"
user_CXXFLAGS="$CXXFLAGS"
user_FCFLAGS="$FCFLAGS"
user_LDFLAGS="$LDFLAGS"
user_CPPFLAGS="$CPPFLAGS"

# Suppress the default compiler flags (namely, -g -O2).
#
# DEV: I copied this from the `prediccs` version. It is also common in other
# packages (e.g., NetCDF) but I don't have any further justification for it.
: ${CFLAGS=""}

# Augment environment variable values with Automake values, to ensure proper
# configuration. See
# https://www.gnu.org/software/automake/manual/html_node/Program-Variables.html
#
# DEV: Adapted from HDF5
CFLAGS="${AM_CFLAGS} ${CFLAGS}"
FCFLAGS="${AM_FCFLAGS} ${FCFLAGS}"
CXXFLAGS="${AM_CXXFLAGS} ${CXXFLAGS}"
CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS}"
LDFLAGS="${AM_LDFLAGS} ${LDFLAGS}"

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_AWK
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P
AC_PROG_RANLIB

# Command-line option for libconfig. NOTE: These should use `@<:@=DIR@:>@`
# instead of `[=DIR]`, because `[...]` is part of the M4 quoting scheme.

# LIBS="-lconfig $LIBS"
# AC_ARG_WITH(
#     [config],
#     [
#         AS_HELP_STRING(
#             [--with-config[=DIR]],
#             [seach for libconfig include/ and lib/ in DIR, if given, or in includedir and libdir]
#         ),
#     ],
#     [CONFIG_SET=0],
#     [CONFIG_SET=1]
# )
# if test ${CONFIG_SET} = 1 ; then
#     AM_CXXFLAGS="${AM_CXXFLAGS} -I${with_config}/include"
#     AM_LDFLAGS="${AM_LDFLAGS} -L${with_config}/lib"
# fi

# Checks for libraries.

# DEV: Consider setting `CC=mpicc CXX=mpic++` by default, and implementing a
# `--disable-mpi` option for serial applications.

# Make sure we have required math functions.
#
# DEV: It may be better to combine this with `AC_CHECK_FUNCS([floor memset pow
# sqrt strstr])`, but I'm not sure about the correct order with respect to the
# following macros (e.g., `AC_CHECK_HEADER`).
AC_SEARCH_LIBS([exp], [m])

# Make sure we have zlib.
AC_SEARCH_LIBS(
    [deflate],
    [z],
    [],
    AC_MSG_ERROR([can't find the z library (required for NetCDF).])
)

# Search for a suitable MPI library. We currently require MPI for parallel runs
# but it should become optional if we implement a serial version.
AC_SEARCH_LIBS(
    [MPI_Init],
    [mpi mpich],
    [],
    AC_MSG_ERROR([can't find a suitable MPI library.])
)

# Search for libconfig, which we require for processing input.
AC_SEARCH_LIBS(
    [config_init],
    [config],
    [],
    AC_MSG_ERROR([can't find a suitable version of the libconfig library.])
)

# Search for the NetCDF library, which we require for creating output.
AC_SEARCH_LIBS(
    [nc_create],
    [netcdf],
    [],
    AC_MSG_ERROR([can't find a suitable version of the NetCDF library.])
)

# Checks for header files.

# DEV: Copied from `prediccs` version. Necessary?
# AC_CHECK_HEADERS([execinfo.h])

# --> ZLIB (dependency of NetCDF)
AC_CHECK_HEADERS(
    zlib.h,
    [],
    AC_MSG_ERROR([failed: zlib headers are required])
)

# --> MPI (required for parallel runs)
AC_CHECK_HEADERS(
    [mpi.h],
    [],
    AC_MSG_ERROR([failed: MPI headers are required])
)

# --> LIBCONFIG (unconditionally required for input)
AC_CHECK_HEADERS(
    [libconfig.h],
    [],
    AC_MSG_ERROR([failed: libconfig headers are required])
)

# --> NETCDF4 (unconditionally required for output)
AC_CHECK_HEADERS(
    [netcdf.h],
    [],
    AC_MSG_ERROR([failed: NetCDF4 headers are required])
)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
# DEV: Put 'exp' here?
AC_CHECK_FUNCS([floor memset pow sqrt strstr])


# Special tests, -lstdc++ needed, or -DMPICH_IGNORE_CXX_SEEK (darwin issues?).
#
# DEV: I adapted this (almost verbatim) from the `prediccs` version. There may
# be more canonical ways to perform the same tests, if they're even necessary.
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

# Restore the user values of environment variables.
CFLAGS="$user_CFLAGS"
CXXFLAGS="$user_CXXFLAGS"
CPPFLAGS="$user_CPPFLAGS"
LDFLAGS="$user_LDFLAGS"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
