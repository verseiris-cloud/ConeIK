#include "mex.h"
#include "coneik_precision.hpp"
#include <cstddef>
#include <exception>
#include <limits>

// [Q,codes,dps]=coneik_precision_mex(original_gamma_Nx2,original_alpha_1x2)
// 0: agreed 50/80- or 100/200-digit regular solution, -1: agreed unreachable,
// -2: undecided or exact identity (caller must use the full MATLAB path),
// -3: invalid input row. Nonsuccess Q rows are NaN. dps is 80 or 200 only when
// classification is agreed; it is NaN for undecided/invalid rows.
namespace {
bool real_matrix(const mxArray* value) noexcept {
    return mxIsDouble(value)&&!mxIsComplex(value)&&!mxIsSparse(value)&&
        mxGetNumberOfDimensions(value)==2;
}
// One exact-geometry entry: retain expensive multiprecision constants across
// single-point calls, without retaining targets, answers, or an unbounded map.
// A geometry change reconstructs the model. Unloading this MEX clears it.
coneik_native::PrecisionModel& cached_model(double a,double b) {
    static std::optional<coneik_native::PrecisionModel> model;
    static double previous_a=0,previous_b=0;
    if (!model || a!=previous_a || b!=previous_b) {
        model.emplace(a,b);previous_a=a;previous_b=b;
    }
    return *model;
}
}
#if defined(_WIN32)
#define CONEIK_EXPORT __declspec(dllexport)
#else
#define CONEIK_EXPORT
#endif
extern "C" CONEIK_EXPORT void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[]) {
    if (nrhs!=2||nlhs>3) mexErrMsgIdAndTxt("coneik_precision_mex:ArgumentCount",
        "Use [Q,codes,dps]=coneik_precision_mex(gamma_Nx2,alpha_1x2), with at most three outputs.");
    if (!real_matrix(prhs[0])||mxGetN(prhs[0])!=2) mexErrMsgIdAndTxt("coneik_precision_mex:InvalidInput",
        "gamma must be a full real double N-by-2 matrix.");
    if (!real_matrix(prhs[1])||mxGetM(prhs[1])!=1||mxGetN(prhs[1])!=2) mexErrMsgIdAndTxt("coneik_precision_mex:InvalidModel",
        "alpha must be a full real double 1-by-2 row.");
    const double* alpha=mxGetPr(prhs[1]);
    if (!std::isfinite(alpha[0])||!std::isfinite(alpha[1])||
        !(alpha[0]>0&&alpha[0]<90&&alpha[1]>0&&alpha[1]<90)) mexErrMsgIdAndTxt("coneik_precision_mex:InvalidModel",
        "Both cone angles must be finite and strictly between 0 and 90 degrees.");
    const mwSize n=mxGetM(prhs[0]);
    if (n>std::numeric_limits<std::size_t>::max()/(4*sizeof(double))) mexErrMsgIdAndTxt("coneik_precision_mex:InputTooLarge",
        "N-by-4 output size exceeds addressable memory.");
    if (nlhs==0) return;
    try {
        auto& model=cached_model(alpha[0],alpha[1]);
        plhs[0]=mxCreateDoubleMatrix(n,4,mxREAL);double* q=mxGetPr(plhs[0]);
        double* codes=nullptr;double* dps=nullptr;
        if (nlhs>1) {plhs[1]=mxCreateDoubleMatrix(n,1,mxREAL);codes=mxGetPr(plhs[1]);}
        if (nlhs>2) {plhs[2]=mxCreateDoubleMatrix(n,1,mxREAL);dps=mxGetPr(plhs[2]);}
        const double* input=mxGetPr(prhs[0]);
        for (mwSize i=0;i<n;++i) {
            const auto out=model.solve(input[i],input[i+n]);
            for (mwSize j=0;j<4;++j) q[i+j*n]=out.q[j];
            if (codes) codes[i]=static_cast<int>(out.status);
            if (dps) dps[i]=model.precision_dps();
        }
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("coneik_precision_mex:NativeError","%s",e.what());
    } catch (...) {
        mexErrMsgIdAndTxt("coneik_precision_mex:NativeError","Unexpected native exception.");
    }
}
