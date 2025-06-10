AC_INIT([hdme], [0.3.0], [https://github.com/j-kieffer/hdme/issues], [hdme])
AC_PREREQ([2.62])
AC_CONFIG_FILES([Makefile hdme.pc])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([hdme_data.h])

# Required programs
AC_PROG_CC
AC_PROG_MAKE_SET
AC_CHECK_TOOL([AR], [ar])
AC_PATH_PROG([PKG_CONFIG], [pkg-config])

# Version numbers
HDME_MAJOR=0
HDME_MINOR=3
HDME_PATCH=0
HDME_VERSION="$HDME_MAJOR.$HDME_MINOR.$HDME_PATCH"
AC_SUBST([HDME_MAJOR])
AC_SUBST([HDME_MINOR]) 
AC_SUBST([HDME_PATCH])
AC_SUBST([HDME_VERSION])

# Set default AR flags
ARFLAGS="${ARFLAGS:-rcs}"
AC_SUBST([ARFLAGS])

# System detection
host_os=`uname -s 2>/dev/null | tr 'A-Z' 'a-z'` || host_os="unknown"

# Initialize variables
LIBS_EXTRA=""
INCS_EXTRA=""

# Check for required libraries
AC_CHECK_LIB([m], [sin], [], [AC_MSG_ERROR([Math library not found])])

# FLINT library (required)
AC_ARG_WITH([flint],
    [AS_HELP_STRING([--with-flint=PATH], [specify location of FLINT (required)])],
    [FLINT_DIR=$withval], [FLINT_DIR=""])

if test "x$FLINT_DIR" != "x"; then
    LIBS_EXTRA="$LIBS_EXTRA -L\"$FLINT_DIR/lib\" -lflint"
    INCS_EXTRA="$INCS_EXTRA -I\"$FLINT_DIR/include\""
    # Set flags for library test
    saved_LDFLAGS="$LDFLAGS"
    saved_CPPFLAGS="$CPPFLAGS"
    LDFLAGS="$LDFLAGS -L$FLINT_DIR/lib"
    CPPFLAGS="$CPPFLAGS -I$FLINT_DIR/include"
    AC_CHECK_LIB([flint], [flint_version], [], [
        AC_MSG_ERROR([FLINT library found in $FLINT_DIR but linking test failed])
    ])
    LDFLAGS="$saved_LDFLAGS"
    CPPFLAGS="$saved_CPPFLAGS"
elif test "x$PKG_CONFIG" != "x"; then
    PKG_CHECK_MODULES([FLINT], [flint], [
        LIBS_EXTRA="$LIBS_EXTRA $FLINT_LIBS"
        INCS_EXTRA="$INCS_EXTRA $FLINT_CFLAGS"
    ], [
        AC_CHECK_LIB([flint], [flint_version], [LIBS_EXTRA="$LIBS_EXTRA -lflint"], [
            AC_MSG_ERROR([FLINT is required but not found. Please install FLINT or specify --with-flint=PATH])
        ])
    ])
else
    AC_CHECK_LIB([flint], [flint_version], [LIBS_EXTRA="$LIBS_EXTRA -lflint"], [
        AC_MSG_ERROR([FLINT is required but not found. Please install FLINT or specify --with-flint=PATH])
    ])
fi

# MPFR library (required)
AC_ARG_WITH([mpfr],
    [AS_HELP_STRING([--with-mpfr=PATH], [specify location of MPFR (required)])],
    [MPFR_DIR=$withval], [MPFR_DIR=""])

if test "x$MPFR_DIR" != "x"; then
    LIBS_EXTRA="$LIBS_EXTRA -L\"$MPFR_DIR/lib\" -lmpfr"
    INCS_EXTRA="$INCS_EXTRA -I\"$MPFR_DIR/include\""
    # Set flags for library test
    saved_LDFLAGS="$LDFLAGS"
    saved_CPPFLAGS="$CPPFLAGS"
    LDFLAGS="$LDFLAGS -L$MPFR_DIR/lib"
    CPPFLAGS="$CPPFLAGS -I$MPFR_DIR/include"
    AC_CHECK_LIB([mpfr], [mpfr_init], [], [
        AC_MSG_ERROR([MPFR library found in $MPFR_DIR but linking test failed])
    ])
    LDFLAGS="$saved_LDFLAGS"
    CPPFLAGS="$saved_CPPFLAGS"
elif test "x$PKG_CONFIG" != "x"; then
    PKG_CHECK_MODULES([MPFR], [mpfr], [
        LIBS_EXTRA="$LIBS_EXTRA $MPFR_LIBS"
        INCS_EXTRA="$INCS_EXTRA $MPFR_CFLAGS"
    ], [
        AC_CHECK_LIB([mpfr], [mpfr_init], [LIBS_EXTRA="$LIBS_EXTRA -lmpfr"], [
            AC_MSG_ERROR([MPFR is required but not found. Please install MPFR or specify --with-mpfr=PATH])
        ])
    ])
else
    AC_CHECK_LIB([mpfr], [mpfr_init], [LIBS_EXTRA="$LIBS_EXTRA -lmpfr"], [
        AC_MSG_ERROR([MPFR is required but not found. Please install MPFR or specify --with-mpfr=PATH])
    ])
fi

# GMP library (required)
AC_ARG_WITH([gmp],
    [AS_HELP_STRING([--with-gmp=PATH], [specify location of GMP (required)])],
    [GMP_DIR=$withval], [GMP_DIR=""])

if test "x$GMP_DIR" != "x"; then
    LIBS_EXTRA="$LIBS_EXTRA -L\"$GMP_DIR/lib\" -lgmp"
    INCS_EXTRA="$INCS_EXTRA -I\"$GMP_DIR/include\""
    # Set flags for library test
    saved_LDFLAGS="$LDFLAGS"
    saved_CPPFLAGS="$CPPFLAGS"
    LDFLAGS="$LDFLAGS -L$GMP_DIR/lib"
    CPPFLAGS="$CPPFLAGS -I$GMP_DIR/include"
    AC_CHECK_LIB([gmp], [__gmpz_init], [], [
        AC_MSG_ERROR([GMP library found in $GMP_DIR but linking test failed])
    ])
    LDFLAGS="$saved_LDFLAGS"
    CPPFLAGS="$saved_CPPFLAGS"
else
    AC_CHECK_LIB([gmp], [__gmpz_init], [LIBS_EXTRA="$LIBS_EXTRA -lgmp"], [
        AC_MSG_ERROR([GMP is required but not found. Please install GMP or specify --with-gmp=PATH])
    ])
fi

# Enable/disable options
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable debug symbols])],
    [DEBUG=$enableval], [DEBUG=no])

AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--enable-openmp], [enable OpenMP support (default)])
AS_HELP_STRING([--disable-openmp], [disable OpenMP support])],
    [OPENMP=$enableval], [OPENMP=yes])

# Set up CFLAGS
BASE_CFLAGS="-Wall -Wextra -O2 -std=c99"

if test "x$DEBUG" = "xyes"; then
    BASE_CFLAGS="$BASE_CFLAGS -g"
fi

# Handle OpenMP
if test "x$OPENMP" = "xyes"; then
    case $host_os in
        darwin*)
            BASE_CFLAGS="$BASE_CFLAGS -Xpreprocessor -fopenmp"
            # Try common Homebrew locations
            if test -d "/opt/homebrew/opt/libomp/lib"; then
                LIBS_EXTRA="$LIBS_EXTRA -L\"/opt/homebrew/opt/libomp/lib\" -lomp"
                INCS_EXTRA="$INCS_EXTRA -I\"/opt/homebrew/opt/libomp/include\""
            elif test -d "/usr/local/opt/libomp/lib"; then
                LIBS_EXTRA="$LIBS_EXTRA -L\"/usr/local/opt/libomp/lib\" -lomp"
                INCS_EXTRA="$INCS_EXTRA -I\"/usr/local/opt/libomp/include\""
            else
                AC_CHECK_LIB([omp], [omp_get_num_threads], [LIBS_EXTRA="$LIBS_EXTRA -lomp"])
            fi
            ;;
        *)
            BASE_CFLAGS="$BASE_CFLAGS -fopenmp"
            ;;
    esac
fi

# Set up library names based on OS
case $host_os in
    darwin*)
        HDME_LIBNAME="libhdme.dylib"
        HDME_LIB="libhdme-$HDME_VERSION.dylib"
        EXEEXT=""
        ;;
    cygwin*|mingw*)
        HDME_LIBNAME="libhdme.dll"
        HDME_LIB="libhdme-$HDME_MAJOR.dll"
        EXEEXT=".exe"
        ;;
    *)
        HDME_LIBNAME="libhdme.so"
        HDME_LIB="libhdme.so.$HDME_VERSION"
        EXEEXT=""
        ;;
esac

# Build final flags
INCS="-I\"\$(CURDIR)\" $INCS_EXTRA"
LIBS="-L\"\$(CURDIR)\" -lm $LIBS_EXTRA"

# Set final CFLAGS if not provided by user
if test "x$CFLAGS" = "x"; then
    CFLAGS="$BASE_CFLAGS"
fi

# Substitute variables for Makefile
AC_SUBST([HDME_LIB])
AC_SUBST([HDME_LIBNAME])
AC_SUBST([INCS])
AC_SUBST([LIBS])
AC_SUBST([EXEEXT])
AC_SUBST([host_os])

AC_OUTPUT


