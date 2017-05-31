/*
 * Program poscalc_indv.cpp
 *
 * Author: Bjorn Olofsson
 * Compile: mload lang/6.1
 *          CC poscalc_indv.cc lubksbBO.c ludcmpBO.c -o poscalc_indv
 * Date of last update: 27.09.2001
 *
 *
 * Calculates new positions for either receivers or shots using
 * first break picks
 *
 * The position of each individual receiver is calculated independently.
 *
 *
 * Methode: "Ausgleichung nach vermittelnden Beobachtungen"
 *
 ****  L_i + v_i = f_i(xHat_j)
 ****            = f_i(xStar_j + delta_xHat_j)
 ****
 **** Taylor expansion, linearising:
 ****  v_i = (df_i / dxHat_j)_0 * delta_xHat_j - l_i  /  l_i = L_i - f_i(x_Hat_j)
 *
 * L_i       : Observations (supposed to be correct in a least-square's sense)
 * v_i       : Corrections (will be minimised)
 *  (f_i       : functional relationship)
 * xStar_j   : 'Estimated' parameters (original values for the parameters)
 * xHat_j    : 'Compensated' parameters (newly calculated values)
 * delta_xHat_j : 'Reduced' vectors (corrections to the original parameter values)
 *
 *
 *----------------------------------------------------------------------*
 * Implementation of the above problem
 *
 * L_i     : First break pick time  --> fbpTime[i]
 * xStar_j : Original receiver (x,y,z) position and time delay --> xStar[j], (j=1..4)
 * xHat_j  : Newly calculated receiver position and time delay --> xHat[j],  (j=1..4)
 * l_i --> ll[i]
 * (df_i / dxHat_j)_0 --> AA[i][j]
 *
 *
 *----------------------------------------------------------------------*
 *
 * Input file format:

#  1117                               (Unique receiver station number)
360.0 1485.0 8.0  546980.3 6887101.0  (rcv depth, water velocity, src depth, rcvx, rcvy)
546703.6    6887542.0    459.500      (srcx, srcy, first break time) (1)
546703.7    6887266.5    353.500      (srcx, srcy, first break time) (2)
546703.7    6887291.5    360.500      (srcx, srcy, first break time) (3)
546704.1    6887492.0    437.000      (srcx, srcy, first break time) (4)
...
#  1118                               (Unique receiver station number)
359.0 1485.0 8.0  547005.3 6887101.0  (rcv depth, water velocity, src depth, rcvx, rcvy)
...

 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <assert.h>
#include <cmath>

using namespace std;

void ludcmpBO( double ** AA, int nIn, int * indx, double * dd ); // From numerical recipes in C
void lubksbBO( double ** AA, int nIn, int * indx, double bb[] ); // From numerical recipes in C


void sub_poscalc_indv(
    double timeDelay, // Time delay to apply to FB picks in advance, default=0.0
    bool verbose,
    int nObs, // Number of observations
    double velocity,
    double rcvx,
    double rcvy,
    double rcvdep,
    double** xSrc,
    double* fbpTime,
    int nUnknowns,
    double* xDelta,
    double* xStddev
)
{

    //  const unsigned short int nUnknowns = 5; // Unknowns : rcvx, rcvy, rcvdep, timeDelay/velocity
    const unsigned short int nStations = 1; // Stations : 1 receiver only

    double ** AA; // Matrix A
    double * ll; // Vector l_i
    // double * ssStar; // Vectors S*i, approximation / distances between shots and receiver, by positions
    // double ssx, ssy, ssz; // Help variables...
    // double x, y, z, t; // Help variables
    // double distance; // Maximum distance to process, distance calc. from approximate positions
    double ** NN; // Coefficient matrix N
    double ** NNinv; // Inverse of NN
    double * nn; // normal vector...
    double * xStar; // Approximate unknowns / original receiver position and time delay
    double * xHat; // Compensated position / new receiver position
    double * delta_xHat; // 'Reduced' positions
    double * sigma_xx; // Covariance matrix --> only diagonal elements!
    double sigma2; // Variance
    double * vv; // Correction vector

    double * srcRcvDist; // Distance between rcv & src (offset and depth accounted for)
    double left_side, right_side, ss;

    int * indxTmp; // Only needed for LU programs
    double dTmp, * colTmp; // Only needed for LU programs
    double ** LUmatrix;

    AA = NULL;
    ll = NULL;
    NN = NNinv = NULL;

    //--------------------------------------------------------------------------------
    // Command line entries
    xStar = new double [nUnknowns];
    assert(xStar != NULL);
    srcRcvDist = new double[nObs];
    assert( srcRcvDist != NULL );

    //  velocity / 1000.0 in order to convert from [m/s] into [m/ms]
    velocity /= 1000.0;

    xStar[0] = rcvx; // X
    xStar[1] = rcvy; // Y
    if( nUnknowns > 2 ) {
        xStar[2] = rcvdep; // Z, dept
        if( nUnknowns > 3 ) {
            xStar[3] = timeDelay;  // use variable 'timeDelay' instead of xStar[3] where possible
            if( nUnknowns > 4 ) {
                xStar[4] = velocity; // [m/ms]]
            }
        }
    }

    double sum;
    for( int iobs = 0; iobs < nObs; iobs++ ) {
        sum = 0.0;
        for( int j = 0; j < nUnknowns; j++ ) {
            ss = xStar[j] - xSrc[iobs][j];
            sum += ss * ss;
        }
        srcRcvDist[iobs] = sqrt(sum); // Approximate distance
    }

    //--------------------------------------------------------------------------------
    // Initialize some fields

    nn = new double [nUnknowns];
    assert(nn != NULL);
    delta_xHat = new double [nUnknowns];
    assert(delta_xHat != NULL);
    xHat = new double [nUnknowns];
    assert(xHat != NULL);
    sigma_xx = new double [nUnknowns];
    assert(sigma_xx != NULL);

    colTmp = new double [nUnknowns];
    assert(colTmp != NULL);
    indxTmp = new int [nUnknowns];
    assert(indxTmp != NULL);

    //--------------------------------------------------
    // Initialize vectors and matrizes:

    ll = new double [nObs];
    assert(ll != NULL);
    vv = new double [nObs];
    assert(vv != NULL);

    AA = new double * [nObs];
    assert(AA != NULL);
    for( int i = 0; i < nObs; i++ ) {
        AA[i] = new double [nUnknowns]; // 4 unknowns: rcvx, rcvy, rcvdep, timeDelay
        assert(AA[i] != NULL);
    }

    NN = new double * [nUnknowns];
    assert(NN != NULL);
    NNinv = new double * [nUnknowns];
    assert(NNinv != NULL);
    LUmatrix = new double * [nUnknowns];
    assert(LUmatrix != NULL);
    for( int j = 0; j < nUnknowns; j++ ) {
        NN[j] = new double [nUnknowns]; // 3 unknowns: rcvx, rcvy, rcvdep
        assert(NN[j] != NULL);
        NNinv[j] = new double [nUnknowns];
        assert(NNinv[j] != NULL);
        LUmatrix[j] = new double [nUnknowns];
        assert(LUmatrix[j] != NULL);
    }

    //--------------------------------------------------
    // Step 1: Setting initial matrizes and vectors

    for( int i = 0; i < nObs; i++ ) {
        ll[i] = fbpTime[i] - (srcRcvDist[i]/velocity + timeDelay);
        AA[i][0] = (xStar[0] - xSrc[i][0]) / (srcRcvDist[i] * velocity);
        AA[i][1] = (xStar[1] - xSrc[i][1]) / (srcRcvDist[i] * velocity);
        if( nUnknowns > 2 )
            AA[i][2] = (xStar[2] - xSrc[i][2]) / (srcRcvDist[i] * velocity);
        if( nUnknowns > 3 )
            AA[i][3] = 1.0;
        if( nUnknowns > 4 )
            AA[i][4] = srcRcvDist[i] / (velocity*velocity);

        if (verbose) {
            fprintf(stdout,"ll[%d] = %f srcrcvdist= %f\n",i, ll[i], srcRcvDist[i]);
            fprintf(stdout,"  [%d] = xStar=%f yStar=%f zStar=%f\n",i, xStar[0], xStar[1], xStar[2]);
            fprintf(stdout,"  [%d] = xSrc[0]=%f xSrc[1]=%f zSrc=%f\n",i, xSrc[i][0], xSrc[i][1], xSrc[i][2]);
        }
    }

    //  if (1) { exit(-1); }

    //--------------------------------------------------
    // Step 2: Calculate coefficient matrix NN and normal vector nn
    // NN = AA(T) * AA    / AA(T): Transposed AA
    // nn = AA(T) * ll

    for( int j = 0; j < nUnknowns; j++) {
        for( int k = 0; k < nUnknowns; k++) {
            sum = 0;
            for( int i = 0; i < nObs; i++) {
                sum += AA[i][j] * AA[i][k];
            }
            NN[j][k] = sum;
        }
    }

    for( int j = 0; j < nUnknowns; j++) {
        sum = 0;
        for( int i = 0; i < nObs; i++) {
            sum += AA[i][j] * ll[i];
        }
        nn[j] = sum;
    }

    //--------------------------------------------------
    // Step 3: Invert matrix NN --> NNinv
    //

    // Copy matrix NN into temporary matrix LUmatrix
    for( int j = 0; j < nUnknowns; j++) {
        for( int k = 0; k < nUnknowns; k++) {
            LUmatrix[j][k] = NN[j][k];
        }
    }

    ludcmpBO(LUmatrix,nUnknowns,indxTmp,&dTmp);

    for( int j = 0; j < nUnknowns; j++) {
        for( int k = 0; k < nUnknowns; k++) colTmp[k] = 0.0;
        colTmp[j] = 1.0;
        lubksbBO(LUmatrix,nUnknowns,indxTmp,colTmp);
        for( int k = 0; k < nUnknowns; k++) NNinv[k][j] = colTmp[k];
    }

    //--------------------------------------------------
    // Test paragraph only:
    // Print NN and NNinv matrix

    if( verbose ) {
        for( int i = 0; i < nUnknowns; i++) {
            for( int j = 0; j < nUnknowns; j++) {
                fprintf(stdout,"%-10.4f ",NN[i][j]);
            }
            cout << endl;
        }
        cout << endl;
        for( int i = 0; i < nUnknowns; i++) {
            for( int j = 0; j < nUnknowns; j++) {
                fprintf(stdout,"%-10.4f ",NNinv[i][j]);
            }
            cout << endl;
        }
        cout << endl;

        fprintf(stdout,"Check matrix inversion: NN * NNinv =\n");
        // Check matrix inversion
        for( int i = 0; i < nUnknowns; i++) {
            for( int j = 0; j < nUnknowns; j++) {
                sum = 0.0;
                for( int k = 0; k < nUnknowns; k++) {
                    sum += NN[i][k]*NNinv[k][j];
                }
                fprintf(stdout,"%-10.4f ",sum);
                //          LUmatrix[i][j] = sum;
            }
            cout << endl;
        }
    }

    //--------------------------------------------------
    // Step 4: Calculate reduced vector delta_xHat and
    // compensated parameters (vector xHat)
    //

    for( int j = 0; j < nUnknowns; j++) {
        sum = 0.0;
        for( int k = 0; k < nUnknowns; k++) {
            sum += NNinv[j][k]*nn[k];
        }
        delta_xHat[j] = sum;
        xHat[j] = xStar[j] + delta_xHat[j];
    }

    //--------------------------------------------------
    // Step 5: Calculate correction vector vv
    //

    for( int i = 0; i < nObs; i++) {
        sum = 0.0;
        for( int k = 0; k < nUnknowns; k++) {
            sum += AA[i][k]*delta_xHat[k];
        }
        vv[i] = sum - ll[i];
    }

    //--------------------------------------------------
    // Step 6: Calculate variance matrix sigma2 and
    // covariance matrix sigma_xx
    // Don't need full matrix sigma_xx --> only diagonal elements!

    sum = 0.0;
    for( int i = 0; i < nObs; i++) {
        sum += vv[i]*vv[i];
    }
    sigma2 = sum / (nUnknowns-nStations);
    for( int k = 0; k < nUnknowns; k++) {
        sigma_xx[k] = sigma2 * NNinv[k][k];
    }

    //--------------------------------------------------
    // Step 7: Cross check results
    //  Method 1: v(T)*v = -v(T)*ll
    //  Method 2: A(T)*v = 0

    left_side = sigma2 * (nUnknowns-nStations);
    sum = 0.0;
    for( int i = 0; i < nObs; i++) {
        sum += vv[i]*ll[i];
    }
    right_side = -sum;

    if( verbose ) {
        fprintf(stdout,"Error check 1:  %-8.4f-%-8.4f = %-14.10f\n",left_side,right_side,left_side-right_side);
        for( int i = 0; i < nObs; i++) {
            sum = 0.0;
            for( int k = 0; k < nUnknowns; k++) {
                sum += AA[i][k]*vv[k];
            }
            fprintf(stdout,"Error check 2(%-2d):  0 = %-14.10f\n",i,sum);
        }
    }

    //--------------------------------------------------
    // Output solutions and standard deviation for each Unknown

    for( int k = 0; k < nUnknowns; k++) {
        xDelta[k]  = delta_xHat[k];
        xStddev[k] = sqrt(fabs(sigma_xx[k]));
    }
    // Convert velocity back to [m/s], and reverse polarity (why? only when reversed the sign is OK)
    if( nUnknowns > 4 ) {
        xDelta[4]  *= -1000;
        xStddev[4] *= 1000;
    }
}

// void sub_poscalc_indv(
//     double timeDelay, // Time delay to apply to FB picks in advance, default=0.0
//     bool verbose,
//     int nObs, // Number of observations
//     double velocity,
//     double rcvx,
//     double rcvy,
//     double rcvdep,
//     double** xSrc,
//     double* fbpTime,
//     int nUnknowns,
//     double* xDelta,
//     double* xStddev
// )
int main(int argc, char *argv[])
{

    double timeDelay = 0.0;
    bool verbose = 1;
    int nObs = 1;
    // sub_poscalc_indv(timeDelay, verbose, nObs, velocity, rcvx, rcvy, rcvdep, xSrc, fbpTime, nUnknowns, xDelta, xStddev);
    return 0;
}
