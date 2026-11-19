#include "mex.h"
#include "coneik_core.hpp"
#include <cstddef>
#include <exception>
#include <limits>

// [Q,codes]=coneik_inverse_mex(gamma_Nx2,alpha_1x2)
// Q is N-by-4 [q11 q21 q12 q22], codes is N-by-1 real double:
// 0 fast solution, -1 safely unreachable, -2 precision/identity fallback,
// -3 invalid row. All nonsuccess rows of Q are NaN. No implicit clamping.

namespace {
bool full_real_double_matrix(const mxArray* value) noexcept {
    return mxIsDouble(value) && !mxIsComplex(value) && !mxIsSparse(value) &&
           mxGetNumberOfDimensions(value)==2;
}
}

#if defined(_WIN32)
#define CONEIK_MEX_EXPORT __declspec(dllexport)
#else
#define CONEIK_MEX_EXPORT
#endif

extern "C" CONEIK_MEX_EXPORT void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    if (nrhs!=2 || nlhs>2) {
        mexErrMsgIdAndTxt("coneik_inverse_mex:ArgumentCount",
            "Use [Q,codes]=coneik_inverse_mex(gamma_Nx2,alpha_1x2), with at most two outputs.");
    }
    if (!full_real_double_matrix(prhs[0]) || mxGetN(prhs[0])!=2) {
        mexErrMsgIdAndTxt("coneik_inverse_mex:InvalidInput",
            "gamma must be a full real double N-by-2 matrix; nonfinite rows return code -3.");
    }
    if (!full_real_double_matrix(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=2) {
        mexErrMsgIdAndTxt("coneik_inverse_mex:InvalidModel",
            "alpha must be a full real double 1-by-2 row of cone angles in degrees.");
    }
    const double* alpha=mxGetPr(prhs[1]);
    if (!std::isfinite(alpha[0]) || !std::isfinite(alpha[1]) ||
        !(alpha[0]>0.0 && alpha[0]<90.0 && alpha[1]>0.0 && alpha[1]<90.0)) {
        mexErrMsgIdAndTxt("coneik_inverse_mex:InvalidModel",
            "Both cone angles must be finite and strictly between 0 and 90 degrees.");
    }
    const mwSize n=mxGetM(prhs[0]);
    if (n>std::numeric_limits<std::size_t>::max()/(4*sizeof(double))) {
        mexErrMsgIdAndTxt("coneik_inverse_mex:InputTooLarge", "N-by-4 output size exceeds addressable memory.");
    }
    // A call used only for input validation must not dereference plhs[0].
    if (nlhs==0) return;
    try {
        const coneik_native::Model model(alpha[0],alpha[1]);
        plhs[0]=mxCreateDoubleMatrix(n,4,mxREAL);
        double* q=mxGetPr(plhs[0]);
        double* codes=nullptr;
        if (nlhs>1) {
            plhs[1]=mxCreateDoubleMatrix(n,1,mxREAL);
            codes=mxGetPr(plhs[1]);
        }
        const double* input=mxGetPr(prhs[0]);
        for (mwSize i=0;i<n;++i) {
            const auto result=coneik_native::solve_position(model,input[i],input[i+n]);
            for (mwSize j=0;j<4;++j) q[i+j*n]=result.q[j];
            if (codes) codes[i]=static_cast<int>(result.status);
        }
    } catch (const std::exception& error) {
        mexErrMsgIdAndTxt("coneik_inverse_mex:NativeError", "%s",error.what());
    } catch (...) {
        mexErrMsgIdAndTxt("coneik_inverse_mex:NativeError", "Unexpected native exception.");
    }
}
