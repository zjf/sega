#include "util.h"

void suicidef(const char* format, ...)
{
  va_list args;
  va_start (args, format);
  fflush  (stdout);
  fprintf (stderr, "FAILURE: ");
  
  if (format != NULL)
  {
    vfprintf (stderr, format, args); 
    fprintf  (stderr, "\n");
  }
  va_end (args);

  exit (0);  
}


FILE* fopen_or_die(const char* name,  const char*  mode)
{
  FILE* f;
  f = fopen (name, mode);
  if (f == NULL)
    suicidef ("fopen_or_die failed to open \"%s\" for \"%s\"", name, mode);

  return f;
}

