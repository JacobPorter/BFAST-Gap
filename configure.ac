#  See LICENSE for details
#
#  The version number in AC_INIT is for the bfast-gap package.
#  Developers - the version number should *NEVER* be changed.

AC_PREREQ(2.59);
AC_INIT([bfast-gap],[0.7.0b],[bfast-help@lists.sourceforge.net],[bfast-gap])
AC_COPYRIGHT([See LICENSE for copyright information.])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([dist-bzip2 subdir-objects])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Sets build_cpu, build_vendor and build_os variables
AC_CANONICAL_BUILD

# Check for C and C++ compilers
AC_PROG_CC
AM_PROG_CC_C_O

AC_PROG_INSTALL
AC_GNU_SOURCE

# set CFLAGS and CXXFLAGS
default_CFLAGS="-Wall -g -O2 -pthread";
extended_CFLAGS="";# "-m64 -D_FILE_OFFSET_BITS=64";

# Define some variables
GITREV="Revision: undefined$";
gitrev="`git log | grep -m 1 commit | awk '{print $2}'`";
if test -n "${gitrev}"; then
	GITREV="${gitrev}";
fi
AC_REVISION(["${GITREV}"])
AC_DEFINE_UNQUOTED(GIT_REV, "${GITREV}", [Define to be the git revision number])

AC_ARG_ENABLE(intel64, [  --enable-intel64        optimize for Intel64 CPU such as Xeon and Core2], [extended_CFLAGS="${extended_CFLAGS} -mtune=nocona"], [])
AC_ARG_ENABLE(unoptimized-sw, [  --enable-unoptimized-sw         use the unoptimized smith waterman algorithm (use for benchmarking only)], [AC_DEFINE(UNOPTIMIZED_SMITH_WATERMAN,1,[Define 1 if we want to use the unoptimized smith waterman (use for benchmarking only)])], [])
AC_ARG_ENABLE(bzlib, [  --disable-bzlib         use this option to disable bzlib support], [AC_DEFINE(DISABLE_BZLIB,1,[Define 1 if we want to disable bzlib support.])], [AC_CHECK_HEADER([bzlib.h], [ 
					LIBS="${LIBS} -lbz2";
					AC_DEFINE(HAVE_LIBBZ2, 1, [Define to 1 if you have the <bzlib.h> header file.])], 
					AC_MSG_ERROR("could not find the bzlib library.  Please use --disable-bzlib if you wish to disable bzlib support."))])

CFLAGS="${default_CFLAGS} ${extended_CFLAGS}";

# Enable large file support; disable with --disable-largefile
AC_SYS_LARGEFILE

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
# Checks for library functions.
AC_FUNC_MALLOC 
AC_FUNC_REALLOC
AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([z], [gzread], 
			 LIBS="${LIBS} -lz";
			 AC_DEFINE(HAVE_LIBZ, 1, [Define to 1 if you have the `z' library (-lz).]),
			 AC_MSG_ERROR("could not find the zlib libaray."))
# AC_CHECK_FUNCS_ONCE requires 2.63, which may not be available, so use AC_CHECH_FUNC on each 
#AC_CHECK_FUNCS_ONCE([floor pow strstr strpbrk strdup strchr strtok_r])
AC_CHECK_FUNC(floor, [AC_DEFINE(HAVE_FLOOR, 1, [Define 1 if you have the function floor.])], [])
AC_CHECK_FUNC(pow, [AC_DEFINE(HAVE_POW, 1, [Define 1 if you have the function pow.])],[])
AC_CHECK_FUNC(sqrt, [AC_DEFINE(HAVE_SQRT, 1, [Define 1 if you have the function sqrt.])],[])
AC_CHECK_FUNC(strchr, [AC_DEFINE(HAVE_STRCHR, 1, [Define 1 if you have the function strpbr.])],[])
AC_CHECK_FUNC(strdup, [AC_DEFINE(HAVE_STRDUP, 1, [Define 1 if you have the function strdup.])],[])
AC_CHECK_FUNC(strpbrk, [AC_DEFINE(HAVE_STRPBRK, 1, [Define 1 if you have the function strpbrk.])],[])
AC_CHECK_FUNC(strstr, [AC_DEFINE(HAVE_STRSTR, 1, [Define 1 if you have the function strstr.])],[])
AC_CHECK_FUNC(strtok_r, [AC_DEFINE(HAVE_STRTOK_R, 1, [Define 1 if you have the function strpbrk.])],[])

AC_FUNC_FSEEKO

# Check types
AC_CHECK_TYPE(int8_t)
AC_CHECK_TYPE(int32_t)
AC_CHECK_TYPE(int64_t)
AC_CHECK_TYPE(uint8_t)
AC_CHECK_TYPE(uint32_t)
AC_CHECK_TYPE(uint64_t)

# Set proprocessor SIZEOF_short_int, SIZEOF_int, SIZEOF_long_int
AC_CHECK_SIZEOF([short int])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stdint.h stdlib.h string.h sys/time.h unistd.h float.h zlib.h fcntl.h])

AC_C_INLINE

AC_CONFIG_FILES([Makefile
				 bfast-gap/Makefile
				 butil/Makefile
				 scripts/Makefile
				 tests/Makefile
				 ])

AC_OUTPUT
