#ifndef __GETOPT_H_
#define __GETOPT_H_

#ifdef _GETOPT_API
#undef _GETOPT_API
#endif

//#include <tchar.h>
#define TCHAR char
#define _T(x) (x)

// Standard GNU options
#define	null_argument		0 /*Argument Null*/
#define	no_argument		0 /*Argument Switch Only*/
#define required_argument	1 /*Argument Required*/
#define optional_argument	2 /*Argument Optional*/

// Shorter Versions of options
#define ARG_NULL 0 /*Argument Null*/
#define ARG_NONE 0 /*Argument Switch Only*/
#define ARG_REQ 1  /*Argument Required*/
#define ARG_OPT 2  /*Argument Optional*/

// Change behavior for C\C++
#ifdef __cplusplus
#define _BEGIN_EXTERN_C extern "C" {
#define _END_EXTERN_C }
#define _GETOPT_THROW throw()
#else
#define _BEGIN_EXTERN_C
#define _END_EXTERN_C
#define _GETOPT_THROW
#endif

_BEGIN_EXTERN_C

extern TCHAR *optarg;
extern int optind;
extern int opterr;
extern int optopt;

struct option {
	const TCHAR* name;
	int has_arg;
	int *flag;
	TCHAR val;
};

extern int getopt(int argc, TCHAR *const *argv, const TCHAR *optstring) _GETOPT_THROW;
extern int getopt_long(int ___argc, TCHAR *const *___argv, const TCHAR *__shortopts, const struct option *__longopts, int *__longind) _GETOPT_THROW;
extern int getopt_long_only(int ___argc, TCHAR *const *___argv, const TCHAR *__shortopts, const struct option *__longopts, int *__longind) _GETOPT_THROW;
_END_EXTERN_C

// Undefine so the macros are not included
#undef _BEGIN_EXTERN_C
#undef _END_EXTERN_C
#undef _GETOPT_THROW

#endif  // __GETOPT_H_
