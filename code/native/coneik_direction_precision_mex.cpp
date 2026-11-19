#include "mex.h"
#include "coneik_direction_precision.hpp"
#include <cstddef>
#include <exception>
#include <limits>

// [Q,codes,dps]=coneik_direction_precision_mex(original_XYZ_Nx3,alpha_1x2)
// 0: agreed regular solution, -1: agreed unreachable, -2: delegate,
// -3: invalid row; 1: exact outer-45 boundary; 2: exact inner-45 boundary.
// Negative-code Q is NaN. Confirmed dps is 80, 200 or 400; identity dps=80
// describes finite-precision azimuth recovery, not the exact dyadic predicate.
// Delegated/invalid dps is NaN. Old meanings and three-output ABI are intact.
namespace {
bool real_matrix(const mxArray* value) noexcept {
    return mxIsDouble(value)&&!mxIsComplex(value)&&!mxIsSparse(value)&&
        mxGetNumberOfDimensions(value)==2;
}
// One exact-geometry entry: retain expensive multiprecision constants across
// single-point calls, without retaining targets, answers, or an unbounded map.
// A geometry change reconstructs the model. Unloading this MEX clears it.
coneik_native::DirectionPrecisionModel& cached_model(double a,double b) {
    static std::optional<coneik_native::DirectionPrecisionModel> model;
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
    if (nrhs!=2||nlhs>3) mexErrMsgIdAndTxt("coneik_direction_precision_mex:ArgumentCount",
        "Use [Q,codes,dps]=coneik_direction_precision_mex(R_Nx3,alpha_1x2), with at most three outputs.");
    if (!real_matrix(prhs[0])||mxGetN(prhs[0])!=3) mexErrMsgIdAndTxt("coneik_direction_precision_mex:InvalidInput",
        "R must be a full real double N-by-3 matrix.");
    if (!real_matrix(prhs[1])||mxGetM(prhs[1])!=1||mxGetN(prhs[1])!=2) mexErrMsgIdAndTxt("coneik_direction_precision_mex:InvalidModel",
        "alpha must be a full real double 1-by-2 row.");
    const double* alpha=mxGetPr(prhs[1]);
    if (!std::isfinite(alpha[0])||!std::isfinite(alpha[1])||
        !(alpha[0]>0&&alpha[0]<90&&alpha[1]>0&&alpha[1]<90)) mexErrMsgIdAndTxt("coneik_direction_precision_mex:InvalidModel",
        "Both cone angles must be finite and strictly between 0 and 90 degrees.");
    const mwSize n=mxGetM(prhs[0]);
    if (n>std::numeric_limits<std::size_t>::max()/(4*sizeof(double))) mexErrMsgIdAndTxt("coneik_direction_precision_mex:InputTooLarge",
        "N-by-4 output size exceeds addressable memory.");
    if (nlhs==0) return;
    try {
        const coneik_native::direction_precision_detail::Exact45Geometry exact(alpha[0],alpha[1]);
        plhs[0]=mxCreateDoubleMatrix(n,4,mxREAL);double* q=mxGetPr(plhs[0]);
        double* codes=nullptr;double* dps=nullptr;
        if (nlhs>1) {plhs[1]=mxCreateDoubleMatrix(n,1,mxREAL);codes=mxGetPr(plhs[1]);}
        if (nlhs>2) {plhs[2]=mxCreateDoubleMatrix(n,1,mxREAL);dps=mxGetPr(plhs[2]);}
        const double* input=mxGetPr(prhs[0]);
        for (mwSize i=0;i<n;++i) {
            const auto boundary=exact.solve(input[i],input[i+n],input[i+2*n]);
            if (boundary.code>0) {
                for (mwSize j=0;j<4;++j) q[i+j*n]=boundary.q[j];
                if (codes) codes[i]=boundary.code;
                if (dps) dps[i]=boundary.dps;
                continue;
            }
            // Exact boundary rows need no adaptive geometry construction.
            auto& model=cached_model(alpha[0],alpha[1]);
            const auto out=model.solve(input[i],input[i+n],input[i+2*n]);
            for (mwSize j=0;j<4;++j) q[i+j*n]=out.q[j];
            if (codes) codes[i]=static_cast<int>(out.status);
            if (dps) dps[i]=model.precision_dps();
        }
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("coneik_direction_precision_mex:NativeError","%s",e.what());
    } catch (...) {
        mexErrMsgIdAndTxt("coneik_direction_precision_mex:NativeError","Unexpected native exception.");
    }
}
