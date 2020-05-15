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

/*
 * Function to evaluate argv[]. specs is a 0 terminated array of command 
 * line options and types an array that stores the type as one of 
 */
#define INTARG 1
#define DOUBLEARG 2
#define LONGARG 3
#define BOOLARG 4
#define STRINGARG 5
#define BENCHMARK 6
/*
 * for each specifier. Benchmark is specific for cilk samples. 
 * -benchmark or -benchmark medium sets integer to 2
 * -benchmark short returns 1
 * -benchmark long returns 2
 * a boolarg is set to 1 if the specifier appears in the option list ow. 0
 * The variables must be given in the same order as specified in specs. 
 */

void get_options(int argc, char *argv[], const char *specs[], int *types,...);
