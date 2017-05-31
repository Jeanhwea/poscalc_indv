#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <assert.h>
#include <cmath>

void lubksbBO( double ** AA, int nIn, int * indx, double bb[])
{
    int ii = -1;

    for( int i = 0; i < nIn; i++ ) {
        int ip = indx[i];
        double sum = bb[ip];
        bb[ip] = bb[i];
        if( ii != -1 ) {
            for( int j = ii; j <= i-1; j++ ) sum -= AA[i][j]*bb[j];
        }
        else if( sum ) {
            ii = i;
        }
        bb[i] = sum;
    }
    for( int i = nIn-1; i >= 0; i-- ) {
        double sum = bb[i];
        for( int j = i+1; j < nIn; j++ ) sum -= AA[i][j]*bb[j];
        bb[i] = sum/AA[i][i];
    }
}
