/**
 * Copyright (c) 2018 I-Ting Angelina Lee  
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
 **/

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#ifdef INSTRUMENT_RTS
#include <cilk/common.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "getoptions.h"
#include "ktiming.h"
#include "papi.h"

#ifndef RAND_MAX
#define RAND_MAX 32767
#endif

#define REAL double 

#define EPSILON (1.0E-6)

#define THRESH 16


unsigned long rand_nxt = 0;

int cilk_rand(void) {
    int result;
    rand_nxt = rand_nxt * 1103515245 + 12345;
    result = (rand_nxt >> 16) % ((unsigned int) RAND_MAX + 1);
    return result;
}

/*
 * Naive sequential algorithm, for comparison purposes
 */
void matrixmul(REAL *C, REAL *A, REAL *B, int n) {

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      REAL s = (REAL)0;
      for(int k = 0; k < n; ++k) {
        s += A[i*n+k] * B[k*n+j];
      }
      C[i*n+j] = s;
    }
  }
}

/*
 * Compare two matrices.  Print an error message if they differ by
 * more than EPSILON.
 * return 0 if no error, and return non-zero if error.
 */
static int compare_matrix(REAL *A, REAL *B, int n) {

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      // compute the relative error c 
      REAL c = A[n*i + j] - B[n*i + j];
      if (c < 0.0) c = -c;

      c = c / A[n*i + j];
      if (c > EPSILON) { 
        return -1; 
      }
    }
  }
  return 0;
}

//init the matrix to all zeros 
void zero(REAL *M, int n){
    int i;
    for(i = 0; i < n * n; i++) {
        M[i] = 0.0;
    }
}

//init the matrix to random numbers
void init(REAL *M, int n){
    int i,j;
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++){
            M[i*n+j] = (REAL) cilk_rand();
        }
    }
}

void mm_base(REAL *A, REAL *B, REAL *C, int size){
    for( int i = 0; i < size; ++i ){
        for( int j = 0; j < size; ++j ){
            REAL s = (REAL)0;
            for( int k = 0; k < size; ++k ){
                s += A[ i * size + k ] * B[ j * size + k ];
            }
            C[i * size + j] += s;
        }
    }
}

//recursive parallel solution to matrix multiplication - row major order
void mat_mul_par(REAL *A, REAL *B, REAL *C, int n, int orig_n) {

    if(n <= THRESH) {
        mm_base(A, B, C, n);
        return;
    }

    //partition each matrix into 4 sub matrices
    //each sub-matrix points to the start of the z pattern
    REAL *A1 = &A[0];
    REAL *A2 = &A[(n*n) >> 2]; //bit shift to divide by 2
    REAL *A3 = &A[(n*n) >> 1];
    REAL *A4 = &A[((n*n) >> 2) + ((n*n) >>1)];

    REAL *B1 = &B[0];
    REAL *B2 = &B[(n*n) >> 1];
    REAL *B3 = &B[(n*n) >> 2];
    REAL *B4 = &B[((n*n) >> 2) + ((n*n) >>1)];
    
    REAL *C1 = &C[0];
    REAL *C2 = &C[(n*n) >> 2];
    REAL *C3 = &C[(n*n) >> 1];
    REAL *C4 = &C[((n*n) >> 2) + ((n*n) >>1)];

    //recrusively call the sub-matrices for evaluation in parallel
    cilk_spawn mat_mul_par(A1, B1, C1, n >> 1, orig_n);
    cilk_spawn mat_mul_par(A1, B2, C2, n >> 1, orig_n);
    cilk_spawn mat_mul_par(A3, B1, C3, n >> 1, orig_n);
    mat_mul_par(A3, B2, C4, n >> 1, orig_n);
    cilk_sync; //wait here for first round to finish

    cilk_spawn mat_mul_par(A2, B3, C1, n >> 1, orig_n);
    cilk_spawn mat_mul_par(A2, B4, C2, n >> 1, orig_n);
    cilk_spawn mat_mul_par(A4, B3, C3, n >> 1, orig_n);
    mat_mul_par(A4, B4, C4, n >> 1, orig_n);
    cilk_sync; //wait here for all second round to finish

}

void transformA(REAL *src, REAL *des, int ori_row, int ori_col, int cur_n, int ori_n, int cur_size){
    int des_index = 0;
    if(cur_n == THRESH){
        //transform matrix to z layout
        for(int row = 0; row < THRESH; row++){
            int row_index = ori_row + row;
            for(int col = 0; col < THRESH; col++){
                int col_index = ori_col + col;
                des[des_index] = src[row_index*ori_n + col_index];
                des_index++;
            }
        }
        return;
    }
    cilk_spawn transformA(src, des, ori_row, ori_col, cur_n/2, ori_n, cur_size/4);
    cilk_spawn transformA(src, des + cur_size/4, ori_row, ori_col+cur_n/2, cur_n/2, ori_n, cur_size/4);
    cilk_spawn transformA(src, des + cur_size/2, ori_row+cur_n/2, ori_col, cur_n/2, ori_n, cur_size/4);
    transformA(src, des + cur_size/4 + cur_size/2, ori_row + cur_n/2, ori_col + cur_n/2, cur_n/2, ori_n, cur_size/4);
	cilk_sync;
}

void transformB(REAL *src, REAL *des, int ori_row, int ori_col, int cur_n, int ori_n, int cur_size){
    int des_index = 0;
    if(cur_n == THRESH){
        //transform matrix to z layout
        for(int col = 0; col < THRESH; col++){
            int col_index = ori_col + col;
            for(int row = 0; row < THRESH; row++){
                int row_index = ori_row + row;
                des[des_index] = src[row_index*ori_n + col_index];
                des_index++;
            }
        }
        return;
    }
    cilk_spawn transformB(src, des, ori_row, ori_col, cur_n/2, ori_n, cur_size/4);
    cilk_spawn transformB(src, des + cur_size/4, ori_row+cur_n/2, ori_col, cur_n/2, ori_n, cur_size/4);
    cilk_spawn transformB(src, des + cur_size/2, ori_row, ori_col+cur_n/2, cur_n/2, ori_n, cur_size/4);
    transformB(src, des + cur_size/4 + cur_size/2, ori_row + cur_n/2, ori_col + cur_n/2, cur_n/2, ori_n, cur_size/4);
	cilk_sync;
}

void transformResult(REAL *src, REAL *des, int ori_row, int ori_col, int cur_n, int ori_n, int cur_size){
    int des_index = 0;
    if(cur_n == THRESH){
        //transform matrix to z layout
        for(int row = 0; row < THRESH; row++){
            int row_index = ori_row + row;
            for(int col = 0; col < THRESH; col++){
                int col_index = ori_col + col;
                des[row_index*ori_n + col_index] = src[des_index];
                des_index++;
            }
        }
        return;
    }
    cilk_spawn transformResult(src, des, ori_row, ori_col, cur_n/2, ori_n, cur_size/4);
    cilk_spawn transformResult(src + cur_size/4, des, ori_row, ori_col+cur_n/2, cur_n/2, ori_n, cur_size/4);
    cilk_spawn transformResult(src + cur_size/2, des, ori_row+cur_n/2, ori_col, cur_n/2, ori_n, cur_size/4);
    transformResult(src + cur_size/4 + cur_size/2, des, ori_row + cur_n/2, ori_col + cur_n/2, cur_n/2, ori_n, cur_size/4);
	cilk_sync;
}

const char *specifiers[] = {"-n", "-c", "-h", 0};
int opt_types[] = {INTARG, BOOLARG, BOOLARG, 0};

int usage(void) {
  fprintf(stderr, 
      "\nUsage: mm_dac [-n #] [-c]\n\n"
      "Multiplies two randomly generated n x n matrices. To check for\n"
      "correctness use -c\n");
  return 1;
}

int main(int argc, char *argv[]) {

    int n = 2048;  
    int verify = 0;  
    int help = 0;

    get_options(argc, argv, specifiers, opt_types, &n, &verify, &help);
    if (help || argc == 1) return usage();

    REAL *A, *B, *C, *A_T, *B_T, *C_T;

    A = (REAL *) malloc(n * n * sizeof(REAL)); //source matrix 
    B = (REAL *) malloc(n * n * sizeof(REAL)); //source matrix
    C = (REAL *) malloc(n * n * sizeof(REAL)); //result matrix
	A_T = (REAL *) malloc(n * n * sizeof(REAL)); //source matrix 
    B_T = (REAL *) malloc(n * n * sizeof(REAL)); //source matrix
    C_T = (REAL *) malloc(n * n * sizeof(REAL)); //result matrix
    
    init(A, n);
    init(B, n);
    zero(C, n);
	zero(C_T, n);
	
	transformA( A, A_T, 0, 0, n, n,n*n );
	transformB( B, B_T, 0, 0, n, n,n*n );

#ifdef INSTRUMENT_RTS
    __cilkrts_init();
    __cilkrts_reset_timing();
#endif
    clockmark_t begin_rm = ktiming_getmark(); 
    mat_mul_par(A_T, B_T, C_T, n, n);
    clockmark_t end_rm = ktiming_getmark();
#ifdef INSTRUMENT_RTS
    __cilkrts_accum_timing();
#endif

    printf("Elapsed time in seconds: %f\n", ktiming_diff_sec(&begin_rm, &end_rm));
	transformResult(C_T, C, 0, 0, n, n,n*n );

    if(verify) {
        printf("Checking results ... \n");
        REAL *C2 = (REAL *) malloc(n * n * sizeof(REAL));
        matrixmul(C2, A, B, n);
        verify = compare_matrix(C, C2, n);
        free(C2);
    }

    if(verify) {
        printf("WRONG RESULT!\n");
    } else {
        printf("\nCilk Example: matrix multiplication\n");
        printf("Options: n = %d\n\n", n);
    }

    //clean up memory
    delete [] A;
    delete [] B;
    delete [] C;
	
    return 0;
}
