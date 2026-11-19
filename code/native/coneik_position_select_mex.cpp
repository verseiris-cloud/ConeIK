#include "mex.h"
#include "coneik_core.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>

// Optional post-selector ONLY. The caller supplies two verified ordinary
// analytic branches. NaN means decline, never unreachable or no solution.
// Strict C++17 FP; no contraction/fast-math. No target or result cache.
namespace {
constexpr double eps=std::numeric_limits<double>::epsilon();
constexpr double reliable_turn=1048576; // conservative subset of MATLAB's 2^40
bool matrix(const mxArray* x,mwSize r,mwSize c) {
    return mxIsDouble(x)&&!mxIsComplex(x)&&!mxIsSparse(x)&&
        mxGetNumberOfDimensions(x)==2&&mxGetM(x)==r&&mxGetN(x)==c;
}
double wrap(double x) {
    x=std::fmod(x,360.0);
    if (x>=180) x-=360;
    if (x< -180) x+=360;
    return x;
}
double guard(double a,double b) {
    return 128*eps*std::max(1.0,std::max(std::abs(a),std::abs(b)));
}
// 1 selected, 0 no representative, -1 numerical ambiguity -> whole decline.
int lift(double base,double lo,double hi,double reference,double& value) {
    double low=-std::numeric_limits<double>::infinity(),high=-low;
    if (std::isfinite(lo)) low=std::ceil((lo-base)/360);
    if (std::isfinite(hi)) high=std::floor((hi-base)/360);
    double near=std::round((reference-base)/360);
    near=std::min(std::max(near,low),high);
    if (!std::isfinite(near)||std::abs(near)>reliable_turn) return -1;
    // Preserve the middle turn itself (including -0), rather than near+0.
    const double turns[3]={near-1,near,near+1};
    double distance=std::numeric_limits<double>::infinity();int found=0;
    for (double turn:turns) {
        if (std::abs(turn)>reliable_turn) continue;
        const double candidate=base+360*turn;
        if (!std::isfinite(candidate)) return -1;
        if ((std::isfinite(lo)&&std::abs(candidate-lo)<=guard(candidate,lo))||
            (std::isfinite(hi)&&std::abs(candidate-hi)<=guard(candidate,hi))) return -1;
        if (candidate<lo||candidate>hi) continue;
        const double d=std::abs(candidate-reference);
        if (!std::isfinite(d)) return -1;
        if (found&&std::abs(d-distance)<=guard(d,distance)) return -1;
        if (!found||d<distance) {value=candidate;distance=d;found=1;}
    }
    return found;
}
struct Geometry {
    double a,b,sa,sb,ca,cb;
    Geometry(double aa,double bb):a(aa),b(bb) {
        const auto x=coneik_native::detail::sincos_degrees(a);
        const auto y=coneik_native::detail::sincos_degrees(b);
        sa=x.s;ca=x.c;sb=y.s;cb=y.c;
    }
};
const Geometry& geometry(double a,double b) {
    static std::optional<Geometry> g;
    if (!g||g->a!=a||g->b!=b) g.emplace(a,b);
    return *g;
}
bool select(const double* base,const double* r,const double* ref,const double* bounds,
            const double* weights,const Geometry& m,double tolerance,double* selected,int& branch) {
    if (tolerance<1e-10) return false;
    for (int j=0;j<4;++j) if (!std::isfinite(base[j])) return false;
    for (int j=0;j<3;++j) if (!std::isfinite(r[j])) return false;
    for (int j=0;j<2;++j) {
        if (!std::isfinite(ref[j])||!std::isfinite(weights[j])||weights[j]<0||weights[j]>1) return false;
        if (std::isnan(bounds[j])||std::isnan(bounds[j+2])||bounds[j]>bounds[j+2]||
            bounds[j]==INFINITY||bounds[j+2]==-INFINITY) return false;
    }
    if (weights[0]==0&&weights[1]==0) return false;
    const double norm=std::hypot(std::hypot(r[0],r[1]),r[2]);
    if (!std::isfinite(norm)||std::abs(norm-1)>64*eps) return false;
    double candidates[2][2],cost[2];bool feasible[2]={false,false};
    for (int k=0;k<2;++k) {
        bool okay=true;
        for (int j=0;j<2;++j) {
            const double b=wrap(base[2*k+j]);
            const int state=lift(b,bounds[j],bounds[j+2],ref[j],candidates[k][j]);
            if (state<0) return false;
            if (state==0) {okay=false;continue;}
            const double phase_loss=std::abs(wrap(wrap(candidates[k][j])-b));
            if (phase_loss>std::min(1e-10,tolerance/32)) return false;
        }
        if (!okay) continue;
        const double d0=candidates[k][0]-ref[0],d1=candidates[k][1]-ref[1];
        cost[k]=std::hypot(std::sqrt(weights[0])*d0,std::sqrt(weights[1])*d1);
        if (!std::isfinite(cost[k])||!std::isfinite(std::hypot(d0,d1))) return false;
        feasible[k]=true;
    }
    if (!feasible[0]&&!feasible[1]) return false;
    int best=feasible[0]?0:1;
    if (feasible[0]&&feasible[1]) {
        if (std::abs(cost[0]-cost[1])<=guard(cost[0],cost[1])) return false;
        if (cost[1]<cost[0]) best=1;
    }
    const auto u=coneik_native::detail::sincos_degrees(candidates[best][0]);
    const auto v=coneik_native::detail::sincos_degrees(candidates[best][1]);
    const double y=m.sa*m.cb+m.ca*m.sb*v.c,z=m.sb*v.s;
    const double x=m.ca*m.cb-m.sa*m.sb*v.c,yy=u.c*y-u.s*z,zz=u.s*y+u.c*z;
    const double cx=yy*r[2]-zz*r[1],cy=zz*r[0]-x*r[2],cz=x*r[1]-yy*r[0];
    const double residual=std::atan2(std::hypot(std::hypot(cx,cy),cz),(x*r[0]+yy*r[1])+zz*r[2])/coneik_native::rad;
    if (!std::isfinite(residual)||residual>std::min(tolerance/4,2.5e-11)) return false;
    selected[0]=candidates[best][0];selected[1]=candidates[best][1];branch=best+1;
    return true;
}
}
#if defined(_WIN32)
#define CONEIK_EXPORT __declspec(dllexport)
#else
#define CONEIK_EXPORT
#endif
extern "C" CONEIK_EXPORT void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    if (nrhs!=7||nlhs<1||nlhs>3||!matrix(prhs[0],1,4)||!matrix(prhs[1],1,3)||
        !matrix(prhs[2],2,1)||!matrix(prhs[3],2,2)||!matrix(prhs[4],2,1)||
        !matrix(prhs[5],1,2)||!matrix(prhs[6],1,1))
        mexErrMsgIdAndTxt("coneik_position_select_mex:Arguments","Internal post-selector expects seven validated double arrays and one to three outputs.");
    const double* alpha=mxGetPr(prhs[5]);const double tolerance=mxGetScalar(prhs[6]);
    if (!(std::isfinite(alpha[0])&&std::isfinite(alpha[1])&&alpha[0]>0&&alpha[0]<90&&alpha[1]>0&&alpha[1]<90&&
          std::isfinite(tolerance)&&tolerance>0&&tolerance<180))
        mexErrMsgIdAndTxt("coneik_position_select_mex:Model","Invalid geometry or tolerance.");
    plhs[0]=mxCreateDoubleMatrix(2,1,mxREAL);double* q=mxGetPr(plhs[0]);
    q[0]=q[1]=std::numeric_limits<double>::quiet_NaN();int branch=0;
    const bool accepted=select(mxGetPr(prhs[0]),mxGetPr(prhs[1]),mxGetPr(prhs[2]),mxGetPr(prhs[3]),
        mxGetPr(prhs[4]),geometry(alpha[0],alpha[1]),tolerance,q,branch);
    if (nlhs>1) plhs[1]=mxCreateLogicalScalar(accepted);
    if (nlhs>2) plhs[2]=mxCreateDoubleScalar(accepted?branch:mxGetNaN());
}
