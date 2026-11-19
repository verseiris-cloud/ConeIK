#include "mex.h"
#include "coneik_core.hpp"
#include <algorithm>
#include <optional>

// Optional success-only scalar raw-XYZ kernel. No reachability/boundary
// decision is exported. Decline leaves Q/residual NaN and accepted=false.
// [Q,residual,accepted]=coneik_direction_double_mex(R_1x3,alpha_1x2,tolerance)
// Compile C++17 without fast-math or floating-point contraction.
namespace {
bool matrix(const mxArray* x,mwSize rows,mwSize cols) noexcept {
    return mxIsDouble(x)&&!mxIsComplex(x)&&!mxIsSparse(x)&&
        mxGetNumberOfDimensions(x)==2&&mxGetM(x)==rows&&mxGetN(x)==cols;
}
struct Geometry {
    double a,b,outer,inner,sa,sb,ca,cb,product,difference;
    Geometry(double aa,double bb) noexcept:a(aa),b(bb) {
        using namespace coneik_native;
        const auto s=detail::two_sum(a,b),d=detail::two_sum(a,-b);
        outer=s.hi*pi/180;inner=std::abs(d.hi)*pi/180;
        const auto p=detail::sincos_degrees(a),q=detail::sincos_degrees(b);
        sa=p.s;ca=p.c;sb=q.s;cb=q.c;product=sa*sb;
        const auto hi=detail::sincos_degrees(d.hi),lo=detail::sincos_degrees(d.lo);
        difference=hi.s*lo.c+hi.c*lo.s;
    }
};
const Geometry& geometry(double a,double b) {
    static std::optional<Geometry> cached;
    if (!cached||cached->a!=a||cached->b!=b) cached.emplace(a,b);
    return *cached;
}
double degrees_atan2(double y,double x) noexcept {
    return std::atan2(y,x)/coneik_native::rad;
}
bool solve(const double* original,const Geometry& m,double tol,
           std::array<double,4>& q,std::array<double,2>& residual) noexcept {
    using namespace coneik_native;
    if (tol<1e-10) return false;
    double scale=0;
    for (int j=0;j<3;++j) {
        if (!std::isfinite(original[j])) return false;
        scale=std::max(scale,std::abs(original[j]));
    }
    if (scale==0) return false;
    double r[3];for (int j=0;j<3;++j) r[j]=original[j]/scale;
    const double norm=std::sqrt((r[0]*r[0]+r[1]*r[1])+r[2]*r[2]);
    for (int j=0;j<3;++j) {
        r[j]/=norm;
        if (original[j]!=0&&(r[j]==0||detail::subnormal_nonzero(r[j]))) return false;
    }
    const double rho=std::hypot(r[1],r[2]),beta=std::atan2(rho,r[0]);
    const double go=m.outer-beta,gi=beta-m.inner;
    if (rho==0||m.product==0||go<=4*guard||gi<=4*guard) return false;
    const double U2=std::sin((m.outer+beta)/2)*std::sin(go/2);
    const double V2=std::sin((beta+m.inner)/2)*std::sin(gi/2);
    if (!(U2>0&&V2>0&&std::isfinite(U2)&&std::isfinite(V2))) return false;
    const double U=std::sqrt(U2),V=std::sqrt(V2);
    if ((guard*180/pi)*(1/U+1/V+1/rho)>2.5e-9) return false;
    const double eta=2*degrees_atan2(U,V);
    const double delta=degrees_atan2(2*U*V,m.sa*m.difference+2*m.ca*V2);
    const auto result=detail::finish(degrees_atan2(r[2],r[1]),delta,eta);
    if (result.status!=Status::success) return false;
    q=result.q;
    for (int k=0;k<2;++k) {
        const auto u=detail::sincos_degrees(q[2*k]),v=detail::sincos_degrees(q[2*k+1]);
        const double y=m.sa*m.cb+m.ca*m.sb*v.c,z=m.sb*v.s;
        const double fx=m.ca*m.cb-m.sa*m.sb*v.c,fy=u.c*y-u.s*z,fz=u.s*y+u.c*z;
        const double cx=fy*r[2]-fz*r[1],cy=fz*r[0]-fx*r[2],cz=fx*r[1]-fy*r[0];
        residual[k]=degrees_atan2(std::sqrt((cx*cx+cy*cy)+cz*cz),(fx*r[0]+fy*r[1])+fz*r[2]);
        if (!std::isfinite(residual[k])||residual[k]>std::min(tol/4,2.5e-11)) return false;
    }
    return true;
}
}
#if defined(_WIN32)
#define CONEIK_EXPORT __declspec(dllexport)
#else
#define CONEIK_EXPORT
#endif
// Batch gateway. Same scalar acceptance and two-branch residual checks.
// [Q,accepted,residual]=coneik_direction_double_batch_mex(R_Nx3,alpha_1x2,tolerance)
extern "C" CONEIK_EXPORT void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    if (nrhs!=3||nlhs>3)
        mexErrMsgIdAndTxt("coneik_direction_double_batch_mex:ArgumentCount","Three inputs and at most three outputs are required.");
    const mxArray* input=prhs[0];
    if (!mxIsDouble(input)||mxIsComplex(input)||mxIsSparse(input)||mxGetNumberOfDimensions(input)!=2||mxGetN(input)!=3)
        mexErrMsgIdAndTxt("coneik_direction_double_batch_mex:InvalidInput","R must be a full real double N-by-3 matrix.");
    if (!matrix(prhs[1],1,2))
        mexErrMsgIdAndTxt("coneik_direction_double_batch_mex:InvalidModel","alpha must be a full real double 1-by-2 row.");
    const double* alpha=mxGetPr(prhs[1]);
    for (int j=0;j<2;++j) if (!(std::isfinite(alpha[j])&&alpha[j]>0&&alpha[j]<90))
        mexErrMsgIdAndTxt("coneik_direction_double_batch_mex:InvalidModel","Cone angles must be finite and strictly between 0 and 90 degrees.");
    if (!matrix(prhs[2],1,1)||!std::isfinite(mxGetScalar(prhs[2]))||mxGetScalar(prhs[2])<=0||mxGetScalar(prhs[2])>=180)
        mexErrMsgIdAndTxt("coneik_direction_double_batch_mex:InvalidTolerance","Tolerance must be a finite real double scalar in (0,180).");
    if(nlhs==0)return;
    const mwSize n=mxGetM(input);const double* in=mxGetPr(input);const double nan=std::numeric_limits<double>::quiet_NaN();
    plhs[0]=mxCreateDoubleMatrix(n,4,mxREAL);double* Q=mxGetPr(plhs[0]);
    mxLogical* accepted=nullptr;double* errors=nullptr;
    if(nlhs>1){plhs[1]=mxCreateLogicalMatrix(n,1);accepted=mxGetLogicals(plhs[1]);}
    if(nlhs>2){plhs[2]=mxCreateDoubleMatrix(n,2,mxREAL);errors=mxGetPr(plhs[2]);}
    if(n==0)return; // Avoid null-pointer arithmetic on empty MATLAB storage.
    std::fill(Q,Q+4*n,nan);if(errors)std::fill(errors,errors+2*n,nan);
    const auto& m=geometry(alpha[0],alpha[1]);const double tol=mxGetScalar(prhs[2]);
    for(mwSize k=0;k<n;++k){
        const double r[3]={in[k],in[k+n],in[k+2*n]};std::array<double,4> q;std::array<double,2> residual;
        if(solve(r,m,tol,q,residual)){
            for(int j=0;j<4;++j)Q[k+j*n]=q[j];
            if(accepted)accepted[k]=true;
            if(errors)for(int j=0;j<2;++j)errors[k+j*n]=residual[j];
        }
    }
}
