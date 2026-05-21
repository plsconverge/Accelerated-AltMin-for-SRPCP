#include <cmath>
#include <cstring>
#include <algorithm>
#include <functional>
// #include <parallel/algorithm>

#include <execution>

#include "mex.h"

void quickselect(double *array, int n, int k) {
    std::nth_element(std::execution::par, array, array + k - 1, array + n, std::greater<double>());
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // inputs -- 2: array, index k
    // outputs -- 1: sorted array

    // input check
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("nrhs: ", " Input Variables Needed");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("nlhs: ", " Output Variables Needed");
    }

    if (!mxIsDouble(prhs[0])) {
        mexErrMsgIdAndTxt("VarType:", "Input array must be double in type");
    }
    if (!mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt("VarType:", "Input index must be double in type");
    }
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("VarType:", "Input index must be a scalar");
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
    double k0 = mxGetScalar(prhs[1]);
    if (std::fabs(k0 - std::round(k0)) > 1e-8) {
        mexErrMsgIdAndTxt("IntCheck:", "Input index must be an integer");
    }
    int k = static_cast<int>(k0);

    // outputs
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double *output_arr = mxGetPr(plhs[0]);

    // copy array
    memcpy(output_arr, arr, n * sizeof(double));

    // partial sort
    quickselect(output_arr, n, k);

    return;
}
