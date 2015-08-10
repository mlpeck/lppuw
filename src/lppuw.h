#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "lp_lib.h"

/* some possible failure modes for LP setup */

#define MAKELP_FAIL -1
#define SETOBJ_FAIL -2
#define RESIZELP_FAIL -3
#define ADDCONSTR_FAIL -4

