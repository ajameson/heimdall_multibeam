dnl @synopsis SWIN_LIB_PSRDADA
dnl 
AC_DEFUN([SWIN_LIB_PSRDADA],
[
  AC_PROVIDE([SWIN_LIB_PSRDADA])

  AC_REQUIRE([SWIN_PACKAGE_OPTIONS])
  SWIN_PACKAGE_OPTIONS([psrdada])

  AC_MSG_CHECKING([for PSRDADA Library installation])

  if test x"$PSRDADA" == x; then
    PSRDADA=psrdada
  fi

  if test "$have_psrdada" != "user disabled"; then

    SWIN_PACKAGE_FIND([psrdada],[dada_def.h])
    SWIN_PACKAGE_TRY_COMPILE([psrdada],[#include <dada_def.h>])

    SWIN_PACKAGE_FIND([psrdada],[lib$PSRDADA.*])
    SWIN_PACKAGE_TRY_LINK([psrdada],[#include <dada_hdu.h>],
                          [ dada_hdu_create (0);],
                          [-l$PSRDADA])

  fi

  AC_MSG_RESULT([$have_psrdada])

  if test x"$have_psrdada" = xyes; then

    AC_DEFINE([HAVE_PSRDADA],[1],
              [Define if the PSRDADA Library is present])
    [$1]

  else
    :
    [$2]
  fi

  PSRDADA_LIBS="$psrdada_LIBS"
  PSRDADA_CFLAGS="$psrdada_CFLAGS"

  AC_SUBST(PSRDADA_LIBS)
  AC_SUBST(PSRDADA_CFLAGS)
  AM_CONDITIONAL(HAVE_PSRDADA,[test "$have_psrdada" = yes])

])

