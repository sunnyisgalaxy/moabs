#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([moabs], [1.3.8.1], [dsun@tamu.edu])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/bsmap/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Check home directories for external libraries
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost], [home directory for boost]), BOOST_HOME=$withval, BOOST_HOME=)
AC_ARG_WITH(zlib, AS_HELP_STRING([--with-zlib], [home directory for zlib]), ZLIB_HOME=$withval, ZLIB_HOME=)

# Assign CXXFLAGS
AS_IF([test "x$BOOST_HOME" != "x"], [AC_SUBST([CXXFLAGS], ["-I$BOOST_HOME $CXXFLAGS"]), AC_SUBST([LDFLAGS], ["-L$BOOST_HOME/lib $LDFLAGS"])], [])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h])

AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT32_T

_os=`uname`

case "$_os" in
  "Linux")
	  os_linux=yes
    ;;
  "Darwin")
	  os_darwin=yes
    ;;
  *)
    AC_MSG_ERROR([$_os operating system is not suitable to build dmake!])
    ;;
esac

AM_CONDITIONAL([OS_LINUX], [test x"$os_linux" = x"yes"])
AM_CONDITIONAL([OS_DARWIN], [test x"$os_darwin" = x"yes"])

AC_CONFIG_FILES([
	Makefile
	lib/Makefile
	src/Makefile
	src/bsmap/Makefile
	src/ExactNumCI/Makefile
	src/mcall/Makefile
	src/mcomp/Makefile
	src/bseqc2/Makefile
	src/pefilter/Makefile
	])
AC_OUTPUT
