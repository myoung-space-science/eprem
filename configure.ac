#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#-----------------------------------------------------------------------------#
# Initialization and setup.
#-----------------------------------------------------------------------------#

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

# Define the MPI command. This is relevant to MPI library tests.
AC_MSG_CHECKING([setting MPIRUN])
case ${host} in
    *-aix*)  MPIRUN="aprun -n 2";;
    *)       MPIRUN="mpirun -n 2";;
esac
AC_MSG_RESULT([${MPIRUN}])
AC_SUBST([MPIRUN])

# Save the user's environment variables. We want to preserve these in case
# subsequent tests of libraries and header files modify them.
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
# configuration during subsequent tests. See
# https://www.gnu.org/software/automake/manual/html_node/Program-Variables.html
CFLAGS="${AM_CFLAGS} ${CFLAGS}"
FCFLAGS="${AM_FCFLAGS} ${FCFLAGS}"
CXXFLAGS="${AM_CXXFLAGS} ${CXXFLAGS}"
CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS}"
LDFLAGS="${AM_LDFLAGS} ${LDFLAGS}"

#-----------------------------------------------------------------------------#
# Checks for programs.
#-----------------------------------------------------------------------------#
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_AWK
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P
AC_PROG_RANLIB

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
    [found_mpi_lib=yes; break;],
    [found_mpi_lib=no]
)
AC_ARG_WITH(
    [mpi],
    AS_HELP_STRING(
        [--with-mpi@<:@=yes,no,DIR@:>@],
        [Control use of MPI. A value of yes will cause EPREM to use the system MPI installation; a value of no will cause EPREM to build without MPI support. A value of DIR indicates that EPREM should look for MPI header files in DIR/includes and for MPI libraries in DIR/lib. Passing --with-mpi is an alias for --with-mpi=yes; passing --without-mpi is an alias for --with-mpi=no. Finally, if this option is not present, EPREM will assume a value of yes or no according to the availability of MPI.]
    ),
    [with_mpi=$withval],
    [with_mpi=$found_mpi_lib]
)

MPI_INCFLAGS=
MPI_LIBFLAGS=
AS_CASE(
    [x$with_mpi],
    dnl Enable MPI because...
    dnl The user passed --with-mpi
    dnl The user passed --with-mpi=yes
    dnl The user passed nothing and MPI is available
    [x], [use_mpi=1],
    [xyes], [use_mpi=1],
    dnl Disable MPI because...
    dnl The user passed --with-mpi=no
    dnl The user passed --without-mpi
    dnl The user passed nothing and MPI is not available
    [xno], [
        AC_MSG_ERROR([Disabling MPI (--with-mpi=no/--without-mpi) isn't supported yet.])
    ],
    dnl Default: Assume the user passed --with-mpi=DIR
    [
        use_mpi=1
        MPI_INCFLAGS="-I$with_mpi/includes"
        MPI_LIBFLAGS="L$with_mpi/lib -Wl,-rpath,$with_mpi/lib"
    ]
)
AM_CONDITIONAL([USE_MPI], [test $use_mpi = 1])
AC_MSG_CHECKING([whether to use MPI])
AM_COND_IF(
    [USE_MPI],
    dnl True
    [
        AC_MSG_RESULT([yes])
        CC=mpicc
        CXX=mpic++
        AM_CFLAGS="$MPI_INCFLAGS $AM_CFLAGS"
        AM_CPPFLAGS="$MPI_INCFLAGS $AM_CPPFLAGS"
        AM_CXXFLAGS="$MPI_INCFLAGS $AM_CXXFLAGS"
        AM_LDFLAGS="$MPI_LIBFLAGS $AM_LDFLAGS"
    ],
    dnl False
    [AC_MSG_RESULT([no])]
)
AC_CHECK_HEADERS(
    [mpi.h],
    dnl Found
    [],
    dnl Not found
    [
        AM_COND_IF(
            [USE_MPI],
            AC_MSG_ERROR([Can't find MPI headers. Use --with-mpi=no or --without-mpi to build without MPI.]),
        )
    ]
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

# libconfig: unconditionally required for input
AC_CHECK_HEADERS(
    [libconfig.h],
    [],
    AC_MSG_ERROR([failed: libconfig headers are required])
)

# libnetcdf: unconditionally required for output
AC_CHECK_HEADERS(
    [netcdf.h],
    [],
    AC_MSG_ERROR([failed: NetCDF4 headers are required])
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
# Finalization and output.
#-----------------------------------------------------------------------------#
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
