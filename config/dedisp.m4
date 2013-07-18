dnl @synopsis SWIN_LIB_DEDISP
dnl 
AC_DEFUN([SWIN_LIB_DEDISP],
[
  AC_PROVIDE([SWIN_LIB_DEDISP])

  AC_REQUIRE([SWIN_PACKAGE_OPTIONS])
  SWIN_PACKAGE_OPTIONS([dedisp])

  AC_MSG_CHECKING([for DEDISP Library installation])

  if test "$have_dedisp" != "user disabled"; then

    SWIN_PACKAGE_FIND([dedisp],[dedisp.h])
    SWIN_PACKAGE_TRY_COMPILE([dedisp],[#include <dedisp.h>],
													   [dedisp_get_error_string (DEDISP_NO_ERROR);])

    SWIN_PACKAGE_FIND([dedisp],[libdedisp.*])
    SWIN_PACKAGE_TRY_LINK([dedisp],[#include <dedisp.h>],
                          [dedisp_get_error_string (DEDISP_NO_ERROR);],
                          [-ldedisp])
  fi

  AC_MSG_RESULT([$have_dedisp])

  if test x"$have_dedisp" = xyes; then

    AC_DEFINE([HAVE_DEDISP],[1],
              [Define if the DEDISP Library is present])
    [$1]

  else
    :
    [$2]
  fi

  DEDISP_LIBS="$dedisp_LIBS"
  DEDISP_CFLAGS="$dedisp_CFLAGS"

  AC_SUBST(DEDISP_LIBS)
  AC_SUBST(DEDISP_CFLAGS)
  AM_CONDITIONAL(HAVE_DEDISP,[test "$have_dedisp" = yes])

])

