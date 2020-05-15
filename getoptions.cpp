/**
 * Copyright (c) 1996 Massachusetts Institute of Technology 
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 **/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "getoptions.h"

/* Used by example programs to evaluate command line options */

void get_options(int argc, char *argv[], const char *specs[], int *types,...)
{
     va_list ap;
     int type, i;
     int *intval;
     double *doubleval;
     long *longval;
     char * stringval; 

     va_start(ap, types);

     while (((type = *types++) != 0) && (specs != 0)) {
	  switch (type) {
	      case INTARG:
		   intval = (int *) va_arg(ap, int *);
		   for (i = 1; i < (argc - 1); i++)
			if (!(strcmp(argv[i], specs[0]))) {
			     *intval = atoi(argv[i + 1]);
			     argv[i][0] = 0;
			     argv[i + 1][0] = 0;
			}
		   break;
	      case DOUBLEARG:
		   doubleval = (double *) va_arg(ap, double *);
		   for (i = 1; i < (argc - 1); i++)
			if (!(strcmp(argv[i], specs[0]))) {
			     *doubleval = atof(argv[i + 1]);
			     argv[i][0] = 0;
			     argv[i + 1][0] = 0;
			}
		   break;
	      case LONGARG:
		   longval = (long *) va_arg(ap, long *);
		   for (i = 1; i < (argc - 1); i++)
			if (!(strcmp(argv[i], specs[0]))) {
			     *longval = atol(argv[i + 1]);
			     argv[i][0] = 0;
			     argv[i + 1][0] = 0;
			}
		   break;
	      case BOOLARG:
		   intval = (int *) va_arg(ap, int *);
		   *intval = 0;
		   for (i = 1; i < argc; i++)
			if (!(strcmp(argv[i], specs[0]))) {
			     *intval = 1;
			     argv[i][0] = 0;
			}
		   break;
	      case STRINGARG:
		  stringval = (char *) va_arg(ap, char *);
		   for (i = 1; i < (argc - 1); i++)
			if (!(strcmp(argv[i], specs[0]))) {
			     strcpy(stringval, (char *)argv[i + 1]);
			     argv[i][0] = 0;
			     argv[i + 1][0] = 0;
			}
		   break;
	      case BENCHMARK:
		   intval = (int *) va_arg(ap, int *);
		   *intval = 0;
		   for (i = 1; i < argc; i++) {
			if (!(strcmp(argv[i], specs[0]))) {
			     *intval = 2;
			     if ((i + 1) < argc) {
				  if (!(strcmp(argv[i + 1], "short")))
				       *intval = 1;
				  if (!(strcmp(argv[i + 1], "medium")))
				       *intval = 2;
				  if (!(strcmp(argv[i + 1], "long")))
				       *intval = 3;
				  argv[i + 1][0] = 0;
			     }
			     argv[i][0] = 0;
			}
		   }
		   break;
	  }
	  specs++;
     }
     va_end(ap);

     for (i = 1; i < argc; i++)
	  if (argv[i][0] != 0)
	       printf("\nInvalid option: %s\n", argv[i]);

}
