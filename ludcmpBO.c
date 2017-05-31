#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <assert.h>
#include <cmath>

#define TINY 1.0e-10;
#define FREE_ARG char*


void ludcmpBO( double** AA, int nIn, int* indx, double* dd)
{
    double dum, sum, temp;
    int imax = 0;

    double* vv = (double *) malloc( (size_t) (nIn * sizeof(double)) );
    if (!vv) fprintf(stdout,"Allocation failure in ludcmpBO\n");

    *dd = 1.0;
    for( int i = 0; i < nIn; i++ ) {
        double big = 0.0;
        for( int j = 0; j < nIn; j++ ) {
            if ((temp = fabs(AA[i][j])) > big) big = temp;
        }
        if (big == 0.0) fprintf(stdout,"Singular matrix in routine ludcmp!\n");
        vv[i] = 1.0 / big;
    }

    for( int j = 0; j < nIn; j++ ) {
        for( int i = 0; i < j; i++ ) {
            sum = AA[i][j];
            for( int k = 0; k < i; k++) sum -= AA[i][k]*AA[k][j];
            AA[i][j] = sum;
        }
        double big = 0.0;
        for( int i = j; i < nIn; i++ ) {
            sum = AA[i][j];
            for( int k = 0; k < j; k++ ) {
                sum -= AA[i][k]*AA[k][j];
            }
            AA[i][j] = sum;
            if( (dum = vv[i]*fabs(sum)) >= big ) {
                big = dum;
                imax = i;
            }
        }
        if( j != imax ) {
            for( int k = 0; k < nIn; k++ ) {
                dum = AA[imax][k];
                AA[imax][k] = AA[j][k];
                AA[j][k] = dum;
            }
            *dd = -(*dd);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if( AA[j][j] == 0.0 ) AA[j][j] = TINY;
        if( j != nIn ) {
            dum = 1.0 / AA[j][j];
            for( int i = j+1; i < nIn; i++ ) AA[i][j] *= dum;
        }
    }

    free( (FREE_ARG) (vv) );

}
