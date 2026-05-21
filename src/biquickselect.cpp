#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <functional>

#include <execution>

#include "mex.h"


void biquickselect(double *array, int n, int k1, int k2) {
    std::nth_element(std::execution::par, array, array + k2 - 1, array + n, std::greater<double>());
    std::nth_element(std::execution::par, array, array + k1 - 1, array + k2 - 1, std::greater<double>());
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // inputs -- 3: array, interval [k1, k2]
    // outputs -- 1: partially sorted array

    // input check
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("nrhs: ", " Input Variables Needed");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("nlhs: ", " Output Variables Needed");
    }

    if (!mxIsDouble(prhs[0])) {
        mexErrMsgIdAndTxt("VarType:", "Input array must be double in type");
    }
    if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2])) {
        mexErrMsgIdAndTxt("VarType:", "Input interval nodes must be double in type");
    }
    if (mxGetNumberOfElements(prhs[1]) != 1 || mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("VarType:", "Input interval nodes must be scalars");
    }

    // dimension check
    size_t rows = mxGetM(prhs[0]);
    size_t cols = mxGetN(prhs[0]);
    if (rows != 1 && cols != 1) {
        mexErrMsgIdAndTxt("ArrCheck:", "Input must be an array");
    }
    // get length and pointer
    mwSize n = (mwSize)mxGetNumberOfElements(prhs[0]);
    double *arr = mxGetPr(prhs[0]);

    // obtain scalars
    double k = mxGetScalar(prhs[1]);
    if (std::fabs(k - std::round(k)) > 1e-8) {
        mexErrMsgIdAndTxt("IntCheck:", "Input interval nodes must be integers");
    }
    int k1 = static_cast<int>(k);

    k = mxGetScalar(prhs[2]);
    if (std::fabs(k - std::round(k)) > 1e-8) {
        mexErrMsgIdAndTxt("IntCheck:", "Input interval nodes must be integers");
    }
    int k2 = static_cast<int>(k);

    // outputs
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double *output_arr = mxGetPr(plhs[0]);

    // copy array
    memcpy(output_arr, arr, n * sizeof(double));

    // partial sort
    biquickselect(output_arr, n, k1, k2);

    return;
}

// void biquickselect(double *array, int n, int k) {
//     if (k > n / 2) {
//         std::nth_element(std::execution::par, array, array + k - 2, array + n, std::greater<double>());
        
//         auto maxIt = std::max_element(std::execution::par, array + k - 1, array + n - 1);
//         if (maxIt != (array + k - 1)) {
//             std::swap(*maxIt, array[k - 1]);
//         }
//     } else {
//         std::nth_element(std::execution::par, array, array + k - 1, array + n, std::greater<double>());

//         auto minIt = std::min_element(std::execution::par, array, array + k - 2);
//         if (minIt != (array + k - 2)) {
//             std::swap(*minIt, array[k - 2]);
//         }
//     }
// }

// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//     // inputs -- 2: array, index k
//     // outputs -- 1: sorted array

//     // input check
//     if (nrhs != 2) {
//         mexErrMsgIdAndTxt("nrhs: ", " Input Variables Needed");
//     }
//     if (nlhs != 1) {
//         mexErrMsgIdAndTxt("nlhs: ", " Output Variables Needed");
//     }

//     if (!mxIsDouble(prhs[0])) {
//         mexErrMsgIdAndTxt("VarType:", "Input array must be double in type");
//     }
//     if (!mxIsDouble(prhs[1])) {
//         mexErrMsgIdAndTxt("VarType:", "Input index must be double in type");
//     }
//     if (mxGetNumberOfElements(prhs[1]) != 1) {
//         mexErrMsgIdAndTxt("VarType:", "Input index must be a scalar");
//     }

//     // dimension check
//     size_t rows = mxGetM(prhs[0]);
//     size_t cols = mxGetN(prhs[0]);
//     if (rows != 1 && cols != 1) {
//         mexErrMsgIdAndTxt("ArrCheck:", "Input must be an array");
//     }
//     // get length and pointer
//     mwSize n = (mwSize)mxGetNumberOfElements(prhs[0]);
//     double *arr = mxGetPr(prhs[0]);

//     // obtain scalars
//     double k0 = mxGetScalar(prhs[1]);
//     if (std::fabs(k0 - std::round(k0)) > 1e-8) {
//         mexErrMsgIdAndTxt("IntCheck:", "Input index must be an integer");
//     }
//     int k = static_cast<int>(k0);

//     // outputs
//     plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
//     double *output_arr = mxGetPr(plhs[0]);

//     // copy array
//     memcpy(output_arr, arr, n * sizeof(double));

//     // partial sort
//     biquickselect(output_arr, n, k);

//     return;
// }