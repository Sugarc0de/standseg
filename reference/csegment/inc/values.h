/*
 * values.h -- minimal replacement for the glibc <values.h> that segment.h
 * expects.  macOS does not ship one.  Values match glibc's definitions
 * exactly (glibc simply forwards to <limits.h>/<float.h>).
 */
#ifndef VALUES_H_SHIM
#define VALUES_H_SHIM

#include <limits.h>
#include <float.h>

#ifndef MAXSHORT
#define MAXSHORT  SHRT_MAX
#endif
#ifndef MAXINT
#define MAXINT    INT_MAX
#endif
#ifndef MAXLONG
#define MAXLONG   LONG_MAX
#endif
#ifndef MAXFLOAT
#define MAXFLOAT  FLT_MAX
#endif
#ifndef MAXDOUBLE
#define MAXDOUBLE DBL_MAX
#endif
#ifndef MINFLOAT
#define MINFLOAT  FLT_MIN
#endif
#ifndef MINDOUBLE
#define MINDOUBLE DBL_MIN
#endif

#endif
