#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "assert.h"
#include "pngio.h"

#define PNG_DEBUG 3

#include "mat.h"

using namespace VEC;

PngIO::PngIO(bool bw): _bw(bw) {
}

void PngIO::write(char *file, MatI &mat) {
    // This could be made faster (faster printing of matrix [access pointer
    // directly] and input to matrix2png stdin
    char *_tmp = (char *)"tmp.tmp.tmp";
    FILE *pOUT = fopen(_tmp, "w");
    assert(pOUT);
    
    // Print data to file:
    fputs("CORNERLABEL", pOUT);
    for (int i = 0; i < mat.cols(); ++i) {
        fprintf(pOUT, "\t%d", i);
    }
    fputs("\n", pOUT);

    for (int m = 0; m < mat.rows(); ++m) {
        fprintf(pOUT, "%d", m);
        for (int n = 0; n < mat.cols(); ++n) {
            fprintf(pOUT, "\t%d", mat(m,n));
            //printf("%d ", mat(m,n));
        }
        //printf("\n");
        fputs("\n", pOUT);
    }
    fclose(pOUT);
    
    // CREATE system call to matrix2png:
    char str1[1000]; 
    strcpy (str1, "matrix2png -data ");
    strcat(str1, _tmp);
    strcat(str1, " ");
    if (_bw) {
        strcat(str1, "-mincolor white -maxcolor black");
    }
    else {
        strcat(str1, "-mincolor green -maxcolor red");
    }
    strcat(str1, " >");
    strcat(str1, file);
    printf("*****************************************************\n");
    printf("Calling: %s\n", str1);
    int ret = system(str1);
    printf("SYSTEM RETURNED %d\n", ret);
    printf("*****************************************************\n");

    // CLEANUP:
    remove(_tmp);
}

void PngIO::write(char *file, MatF &mat) {
    // This could be made faster (faster printing of matrix [access pointer
    // directly] and input to matrix2png stdin
    char *_tmp = (char *)"tmp.tmp.tmp";
    FILE *pOUT = fopen(_tmp, "w");
    assert(pOUT);
    
    // Print data to file:
    fputs("CORNERLABEL", pOUT);
    for (int i = 0; i < mat.cols(); ++i) {
        fprintf(pOUT, "\t%d", i);
    }
    fputs("\n", pOUT);

    for (int m = 0; m < mat.rows(); ++m) {
        fprintf(pOUT, "%d", m);
        for (int n = 0; n < mat.cols(); ++n) {
            fprintf(pOUT, "\t%f", mat(m,n));
            //printf("%d ", mat(m,n));
        }
        //printf("\n");
        fputs("\n", pOUT);
    }
    fclose(pOUT);
    
    // CREATE system call to matrix2png:
    char str1[1000]; 
    strcpy (str1, (char *)"matrix2png -data ");
    strcat(str1, _tmp);
    strcat(str1, " ");
    if (_bw) {
        strcat(str1, "-mincolor white -maxcolor black");
    }
    else {
        strcat(str1, "-mincolor green -maxcolor red");
    }
    strcat(str1, " >");
    strcat(str1, file);
    printf("*****************************************************\n");
    printf("Calling: %s\n", str1);
    int ret = system(str1);
    printf("SYSTEM RETURNED %d\n", ret);
    printf("*****************************************************\n");

    // CLEANUP:
    remove(_tmp);
}

