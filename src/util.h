#ifndef ZJF_UTIL_H
#define ZJF_UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

void suicidef(const char* format, ...);
FILE* fopen_or_die(const char* name,  const char*  mode);

#endif
