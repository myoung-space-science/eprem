#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([eprem],[0.1])
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADERS([src/epremDefs.h])
AC_CANONICAL_HOST()

AC_MSG_CHECKING([setting MPIRUN])
case $host in
    *-aix*)  MPIRUN="aprun -n 2"
             MAKE="make";;
    *)       MPIRUN="aprun -np 2";;
esac
AC_MSG_RESULT([$MPIRUN])
AC_SUBST([MPIRUN])

AS_IF([test -z "$MAKE"],[MAKE="make"])
AC_SUBST([MAKE])

#Set default flags to nothing (overwritten below).

: ${CFLAGS=" "}
: ${CXXFLAGS=" "}
: ${LDFLAGS=" "}
: ${CPPFLAGS=" "}

# Checks for programs.
AC_PROG_CXX([cc])
AC_PROG_CC([cc])
AC_PROG_F77([ftn])
AC_PROG_RANLIB

AC_ARG_WITH([netcdf], [AS_HELP_STRING([--with-netcdf],
            [Location of the netcdf package (contains /include, /lib, etc.) @<:@default=$netcdf_home@:>@])],
            [],
            [])
# default directory to look for NetCDF if the --with-netcdf-home option was used without an argument
AS_IF([test "$with_netcdf" == "yes"],[with_netcdf="/usr"])
AS_IF([test "$with_netcdf" != ""],
[
    NCINC=" -I$with_netcdf/include "
    NCLIB=" -L$with_netcdf/lib -Wl,-rpath,$with_netcdf/lib"
    CXXFLAGS=" $CXXFLAGS $NCINC "
    CPPFLAGS=" $CPPFLAGS $NCINC "
    CFLAGS="   $CFLAGS   $NCINC "
    LDFLAGS="  $LDFLAGS  $NCLIB "
])

AC_ARG_WITH([config], [AS_HELP_STRING([--with-config],
            [Location of the libconfig package (contains /include, /lib, etc.) @<:@default=$libconfig_home@:>@])],
            [],
            [])
AS_IF([test "$with_config" == "yes"],[with_config="/usr"])
AS_IF([test "$with_config" != ""],
[
   CFINC=" -I$with_config/include "
   CFLIB=" -L$with_config/lib -Wl,-rpath,$with_config/lib"
   CXXFLAGS=" $CXXFLAGS $CFINC "
   CPPFLAGS=" $CPPFLAGS $CFINC "
   CFLAGS="   $CFLAGS   $CFINC "
   LDFLAGS="  $LDFLAGS  $CFLIB "
])

# Set compiler flags.
CXXFLAGS=" $CXXFLAGS -O3 -DNDEBUG -Wall -Wextra -mtune=native"
CPPFLAGS=" $CPPFLAGS -O3 -DNDEBUG -Wall -Wextra -mtune=native"
CFLAGS="   $CFLAGS   -O3 -std=c99 -DNDEBUG -Wall -Wextra -mtune=native"

# Check for header files.
AC_CHECK_HEADERS([mpi.h],,AC_MSG_ERROR([MPI and headers required]))
AC_CHECK_HEADERS([netcdf.h],,AC_MSG_ERROR([NETCDF4 and headers required]))

AC_CHECK_HEADERS([execinfo.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_LIB([z],[zlibVersion],,AC_MSG_ERROR([libz installation required]))
AC_CHECK_LIB([jpeg],[jpeg_destroy],,AC_MSG_ERROR([libjpeg installation required]))
AC_CHECK_LIB([netcdf],[nc_inq_libvers],,AC_MSG_ERROR([NetCDF installation required]))
AC_CHECK_LIB([config],[config_init],,AC_MSG_ERROR([libconfig installation required]))

# Special tests, -lstdc++ needed, or -DMPICH_IGNORE_CXX_SEEK (darwin issues?)

AC_LANG(C++)

AC_MSG_CHECKING([if empty compile works])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
 [[
 #include <iostream>
 ]])],
 [
  last_test=0
  AC_MSG_RESULT([yes])
 ],
 [
  last_test=1
  AC_MSG_RESULT([no])
 ])

AS_IF([test "$last_test" -ne 0],
[
 AC_SUBST(orig_LDFLAGS,[$LDFLAGS])
 AC_SUBST(LDFLAGS,["$LDFLAGS -lstdc++"])
 AC_MSG_CHECKING([again with -lstdc++])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
 [[
 #include <iostream>
 ]])],
 [
  last_test=0
  AC_MSG_RESULT([yes])
 ],
 [
  last_test=1
  AC_MSG_RESULT([no])
  AC_SUBST(LDFLAGS,[$orig_LDFLAGS])
 ])
]
)

AC_MSG_CHECKING([if -DMPICH_IGNORE_CXX_SEEK is needed])
AC_LINK_IFELSE([
AC_LANG_PROGRAM(
 [[
 #include <iostream>
 #include <mpi.h>
 ]])],
 [last_test=0],
 [last_test=1])

AS_IF([test "$last_test" -ne 0],
[
AC_SUBST(orig_CPPFLAGS,[$CPPFLAGS])
AC_SUBST(CPPFLAGS,["$CPPFLAGS -DMPICH_IGNORE_CXX_SEEK"])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
 [[
 #include <iostream>
 #include <mpi.h>
 ]])],
        [AC_MSG_RESULT([yes])],
        [
   AC_MSG_RESULT([error])
   AC_SUBST(CPPFLAGS,[$orig_CPPFLAGS])
  ])
],
[
  AC_MSG_RESULT([no])
]
)

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
