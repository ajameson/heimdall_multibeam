dnl @synopsis SWIN_LIB_THRUST
dnl 
AC_DEFUN([SWIN_LIB_THRUST],
[
  AC_PROVIDE([SWIN_LIB_THRUST])

  AC_REQUIRE([SWIN_PACKAGE_OPTIONS])
	AC_REQUIRE([SWIN_LIB_CUDA])

  SWIN_PACKAGE_OPTIONS([thrust])

  AC_MSG_CHECKING([for THRUST Library installation])

  if test x"$THRUST" == x; then
    THRUST=thrust
  fi

  if test "$have_thrust" != "user disabled"; then

		# require version >= 1.6
		AC_LANG_PUSH(C++)
		CXXFLAGS="$CXXFLAGS $CUDA_CFLAGS"

	 	SWIN_PACKAGE_FIND([thrust],[version.h])

		SWIN_PACKAGE_TRY_COMPILE([thrust], [#include <thrust/version.h>],
                             [#if (THRUST_VERSION < 100600)
this_function_should_fail();
#endif])

		AC_LANG_POP

  fi
  
  AC_MSG_RESULT([$have_thrust])

  if test x"$have_thrust" = xyes; then

    AC_DEFINE([HAVE_THRUST],[1],
              [Define if the THRUST Library is present])
    [$1]

  else
    :
    [$2]
  fi

  THRUST_LIBS="$thrust_LIBS"
  THRUST_CFLAGS="$thrust_CFLAGS"

  AC_SUBST(THRUST_LIBS)
  AC_SUBST(THRUST_CFLAGS)
  AM_CONDITIONAL(HAVE_THRUST,[test "$have_thrust" = yes])

])

