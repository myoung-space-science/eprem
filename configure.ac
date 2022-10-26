#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([eprem], [0.1.0], [Matthew.Young@unh.edu])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([src/eprem.c])

# Use `autodefs.h` instead of the default, `config.h`, because the latter is too
# close to the EPREM header file `configuration.h`.
AC_CONFIG_HEADERS([src/autodefs.h])

# Set the pre-processor, compiler, and linker flags to be empty by default.
: ${CFLAGS=""}
: ${CXXFLAGS=""}
: ${LDFLAGS=""}
: ${CPPFLAGS=""}

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

# Update compilation flags after reading command-line arguments.
CXXFLAGS=" ${CXXFLAGS} -O3 -DNDEBUG -Wall -Wextra -mtune=native"
CPPFLAGS=" ${CPPFLAGS} -O3 -DNDEBUG -Wall -Wextra -mtune=native"
CFLAGS="   ${CFLAGS}   -O3 -std=c99 -DNDEBUG -Wall -Wextra -mtune=native"

# Checks for libraries.

AC_MSG_CHECKING([setting MPIRUN])
case ${host} in
    *-aix*)  MPIRUN="aprun -n 2"
             MAKE="make";;
    *)       MPIRUN="mpirun -n 2";;
esac
AC_MSG_RESULT([${MPIRUN}])
AC_SUBST([MPIRUN])

# DEV: Copied from `prediccs` version. Necessary?
# AC_CHECK_HEADERS([execinfo.h])

# DEV: Consider setting `CC=mpicc CXX=mpic++` by default, and implementing a
# `--disable-mpi` option for serial applications.

# --> M (required math functions)
#
# DEV: It may be better to combine this with `AC_CHECK_FUNCS([floor memset pow
# sqrt strstr])`, but I'm not sure about the correct order with respect to the
# following macros (e.g., `AC_CHECK_HEADER`).
AC_SEARCH_LIBS([exp], [m])

# --> ZLIB (dependency of NetCDF)
AC_CHECK_HEADER(
    zlib.h,
    [],
    AC_MSG_ERROR([failed: zlib headers are required]),
)

# --> MPI (required for parallel runs)
AC_SEARCH_LIBS([MPI_Init], [mpi mpich])
AC_CHECK_HEADER(
    [mpi.h],
    [],
    AC_MSG_ERROR([failed: MPI headers are required]),
)

# --> LIBCONFIG (unconditionally required for input)
AC_SEARCH_LIBS([config_init], [config])
AC_CHECK_HEADER(
    [libconfig.h],
    [],
    AC_MSG_ERROR([failed: libconfig headers are required]),
)
# --> NETCDF4 (unconditionally required for output)
AC_SEARCH_LIBS([nc_create], [netcdf])
AC_CHECK_HEADER(
    [netcdf.h],
    [],
    AC_MSG_ERROR([failed: NetCDF4 headers are required])
)


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset pow sqrt strstr])


# Special tests, -lstdc++ needed, or -DMPICH_IGNORE_CXX_SEEK (darwin issues?)
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

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
