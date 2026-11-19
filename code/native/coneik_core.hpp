#ifndef CONEIK_CORE_HPP
#define CONEIK_CORE_HPP

// Portable binary64 ConeIK position kernel, with factored axial SP4/SP1 (C++17).
// Angles are degrees. Model geometry is r=Rx(q1)Rz(a1)Rx(q2)Rz(a2)ex.
// Status::needs_precision is NOT a solution: the caller must pass the
// original target and geometry to a precision/identity-capable solver.
// Compile without fast-math or floating-point contraction. Dispatch guards
// are conservative engineering rules, not certified floating-point bounds.

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace coneik_native {

inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double rad = pi / 180.0;
inline constexpr double guard = 32.0 * std::numeric_limits<double>::epsilon();
inline constexpr double joint_budget_deg = 1e-8;

enum class Status : int {
    success = 0,
    unreachable = -1,
    needs_precision = -2,
    invalid_input = -3
};

struct Solution {
    // q11,q21,q12,q22: negative q2 branch first before periodic adjustment.
    std::array<double,4> q;
    Status status;
    explicit Solution(Status code = Status::needs_precision) noexcept
        : q{{std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN()}}, status(code) {}
};

namespace detail {
inline double wrap_degrees(double value) noexcept {
    if (value >= -180.0 && value < 180.0) return value == 0.0 ? 0.0 : value;
    value = std::fmod(value,360.0);
    if (value >= 180.0) value -= 360.0;
    else if (value < -180.0) value += 360.0;
    return value == 0.0 ? 0.0 : value;
}

struct SinCos { double s, c; };
inline SinCos sincos_degrees(double value) noexcept {
    value = wrap_degrees(value);
    const int quadrant = static_cast<int>(std::round(value/90.0));
    const double reduced = (value - 90.0*quadrant)*rad;
    const double s = std::sin(reduced), c = std::cos(reduced);
    switch (quadrant) {
        case 1: return {c,-s};
        case -1: return {-c,s};
        case 2: case -2: return {-s,-c};
        default: return {s,c};
    }
}

struct Expansion { double hi,lo; };
inline Expansion two_sum(double a,double b) noexcept {
    const double hi=a+b, z=hi-a;
    return {hi,(a-(hi-z))+(b-z)};
}

inline bool subnormal_nonzero(double value) noexcept {
    return value != 0.0 && std::abs(value) < std::numeric_limits<double>::min();
}

inline bool degree_underflow(double value) noexcept {
    return value != 0.0 && std::abs(value) < std::numeric_limits<double>::min()/rad;
}

inline Solution finish(double psi,double delta,double eta) noexcept {
    Solution result(Status::success);
    result.q={{wrap_degrees(psi+delta),-eta,wrap_degrees(psi-delta),eta}};
    for (int branch=0;branch<2;++branch) {
        const int j=2*branch;
        const double actual=result.q[j]+result.q[j+1];
        if (actual>180.0) result.q[j+1]-=360.0;
        else if (actual< -180.0) result.q[j+1]+=360.0;
    }
    for (double value:result.q) {
        if (!std::isfinite(value)) return Solution(Status::needs_precision);
    }
    return result;
}

// Unequal-cone axial positive-gap path. Restrict to exactly representable
// degree sums/differences and well-scaled geometry. The near-boundary gap
// is an exact Sterbenz subtraction; do not form it from atan2(beta).
// The acceptance budget below is an engineering guard, not an interval proof.
inline Solution solve_unequal_axis(double g1,double g2,double a,double b) noexcept {
    if (a==b || a<1 || a>80 || b<1 || b>80 || std::abs(a-b)<0.25) return Solution{};
    const auto sum=two_sum(a,b),diff=two_sum(std::max(a,b),-std::min(a,b));
    if (sum.lo!=0 || diff.lo!=0) return Solution{};
    double beta,psi;
    if (g1==0) {beta=std::abs(g2);psi=g2<0 ? -180.0 : 0.0;}
    else if (g2==0) {beta=std::abs(g1);psi=g1>0 ? -90.0 : 90.0;}
    else if (g1==-180 && std::abs(g2)>=90) {beta=180-std::abs(g2);psi=g2>0 ? -180.0 : 0.0;}
    else if (g2==-180 && std::abs(g1)>=90) {beta=180-std::abs(g1);psi=g1>0 ? -90.0 : 90.0;}
    else return Solution{};
    const double outer=sum.hi,inner=diff.hi;
    if (!((beta>=0.9*outer && beta<=outer)||(beta>=inner && beta<=1.1*inner))) return Solution{};
    const double go=outer-beta,gi=beta-inner;
    if (!(go>0 && gi>0)) return Solution{}; // Keep identities/outside semantics.
    const double U2=sincos_degrees((outer+beta)/2).s*sincos_degrees(go/2).s;
    const double V2=sincos_degrees((beta+inner)/2).s*sincos_degrees(gi/2).s;
    if (!(U2>0 && V2>0)) return Solution{};
    const double U=std::sqrt(U2),V=std::sqrt(V2);
    const auto ac=sincos_degrees(a);
    const double t1=ac.s*sincos_degrees(a-b).s,t2=2*ac.c*V2;
    const double A=t1+t2,B=2*U*V,radius=std::hypot(A,B);
    const double budget=128*std::numeric_limits<double>::epsilon()/rad*
        (1+(std::abs(t1)+std::abs(t2)+B)/radius);
    if (!(radius>0 && budget<1e-10)) return Solution{};
    const double eta=2*std::atan2(U,V)/rad,delta=std::atan2(B,A)/rad;
    constexpr double cut=16384*std::numeric_limits<double>::epsilon();
    if (std::abs(std::abs(psi+delta-eta)-180)<=cut ||
        std::abs(std::abs(psi-delta+eta)-180)<=cut) return Solution{};
    return finish(psi,delta,eta);
}

// Extended equal-cone positive-gap coverage. Exact binary64 degree identities
// retain the boundary gap. Supplementary targets require Sterbenz subtraction
// in 180-|g|; |g|<90, non-axes and unequal cones are deliberately not covered.
inline Solution solve_equal_axis(double g1,double g2,double a) noexcept {
    if (a<1.0 || a>80.0) return Solution{};
    double beta,psi;
    if (g1==0.0) {beta=std::abs(g2);psi=g2<0.0 ? -180.0 : 0.0;}
    else if (g2==0.0) {beta=std::abs(g1);psi=g1>0.0 ? -90.0 : 90.0;}
    else if (g1==-180.0 && std::abs(g2)>=90.0) {
        beta=180.0-std::abs(g2);psi=g2>0.0 ? -180.0 : 0.0;
    } else if (g2==-180.0 && std::abs(g1)>=90.0) {
        beta=180.0-std::abs(g1);psi=g1>0.0 ? -90.0 : 90.0;
    } else return Solution{};
    const double outer=2.0*a;
    if (beta<0.9*outer || beta>1.1*outer) return Solution{};
    const double gap=outer-beta; // Exact in the local window (Sterbenz).
    if (!(gap>0.0)) return Solution{}; // Identity and negative-gap routes stay unchanged.
    const double U=std::sqrt(sincos_degrees((outer+beta)/2.0).s*sincos_degrees(gap/2.0).s);
    const double V=sincos_degrees(beta/2.0).s;
    if (!(U>0.0 && V>0.0)) return Solution{};
    const double eta=2.0*std::atan2(U,V)/rad;
    const double delta=std::atan2(U,sincos_degrees(a).c*V)/rad;
    // Avoid changing a periodic representative through a few rounding ulps.
    constexpr double cut_guard=16384.0*std::numeric_limits<double>::epsilon();
    if (std::abs(std::abs((psi+delta)-eta)-180.0)<=cut_guard ||
        std::abs(std::abs((psi-delta)+eta)-180.0)<=cut_guard) return Solution{};
    return finish(psi,delta,eta);
}

// IK-Geo SP4 -> SP1 specialized to the default model's exact axial target.
// Factoring the SP4 discriminant in original degrees retains the outer gap.
// This is not a generic SP4 least-squares solution or a boundary tolerance.
inline Solution solve_default_axis(double g1,double g2) noexcept {
    double beta,psi;
    if (g1==0.0) {beta=std::abs(g2);psi=g2<0.0 ? -180.0 : 0.0;}
    else if (g2==0.0) {beta=std::abs(g1);psi=g1>0.0 ? -90.0 : 90.0;}
    else return solve_equal_axis(g1,g2,20.0);
    if (beta<36.0 || beta>44.0) return Solution(Status::needs_precision);
    const double gap=40.0-beta; // Exact by Sterbenz in this restricted range.
    if (gap<0.0) return Solution(Status::unreachable);
    if (gap==0.0) return Solution(Status::needs_precision); // Existing identity handler.
    const double U=std::sqrt(sincos_degrees((40.0+beta)/2.0).s*sincos_degrees(gap/2.0).s);
    const double V=sincos_degrees(beta/2.0).s;
    if (!(U>0.0 && V>0.0)) return Solution(Status::needs_precision);
    return finish(psi,std::atan2(U,sincos_degrees(20.0).c*V)/rad,2.0*std::atan2(U,V)/rad);
}
} // namespace detail

class Model {
public:
    // Throws std::invalid_argument for mathematical domain violations.
    // Positive angles that underflow are retained and request precision.
    explicit Model(double alpha1_deg=20.0,double alpha2_deg=20.0)
        : a_(alpha1_deg),b_(alpha2_deg) {
        if (!std::isfinite(a_) || !std::isfinite(b_) ||
            !(a_>0.0 && a_<90.0 && b_>0.0 && b_<90.0)) {
            throw std::invalid_argument("Cone angles must be finite and strictly between 0 and 90 degrees.");
        }
        default_=(a_==20.0 && b_==20.0);
        const auto ac=detail::sincos_degrees(a_),bc=detail::sincos_degrees(b_);
        sa_=ac.s; sb_=bc.s; ca_=ac.c;
        product_=sa_*sb_;
        outer_=(a_+b_)*rad; inner_=std::abs(a_-b_)*rad;
        const auto difference=detail::two_sum(a_,-b_);
        const auto dh=detail::sincos_degrees(difference.hi);
        const auto dl=detail::sincos_degrees(difference.lo);
        sin_difference_=dh.s*dl.c+dh.c*dl.s;
        if (default_) { sa_=std::sin(pi/9.0); ca_=std::cos(pi/9.0); }
    }
    double alpha1_deg() const noexcept {return a_;}
    double alpha2_deg() const noexcept {return b_;}

private:
    double a_,b_,sa_,sb_,ca_,product_,outer_,inner_,sin_difference_;
    bool default_;
    friend Solution solve_position(const Model&,double,double) noexcept;
};

inline Solution solve_position(const Model& model,double pitch_deg,double yaw_deg) noexcept {
    if (!std::isfinite(pitch_deg) || !std::isfinite(yaw_deg)) {
        return Solution(Status::invalid_input);
    }
    const double g1=detail::wrap_degrees(pitch_deg),g2=detail::wrap_degrees(yaw_deg);
    if (detail::degree_underflow(g1) || detail::degree_underflow(g2) ||
        detail::degree_underflow(model.a_) || detail::degree_underflow(model.b_)) {
        return Solution(Status::needs_precision);
    }
    double x,ry,rz;
    if (std::abs(g1)<=90.0 && std::abs(g2)<=90.0) {
        const double cp=std::cos(g1*rad);
        x=cp*std::cos(g2*rad); ry=cp*std::sin(g2*rad); rz=-std::sin(g1*rad);
    } else {
        const auto p=detail::sincos_degrees(g1),y=detail::sincos_degrees(g2);
        x=p.c*y.c; ry=p.c*y.s; rz=-p.s;
    }
    const double rho=std::hypot(ry,rz),beta=std::atan2(rho,x);
    if (rho==0.0 || detail::subnormal_nonzero(ry) || detail::subnormal_nonzero(rz) ||
        detail::subnormal_nonzero(rho) || detail::subnormal_nonzero(beta)) {
        return Solution(Status::needs_precision);
    }
    const double psi=std::atan2(rz,ry)/rad;
    if (model.default_) {
        const double a=pi/9.0, gap=2.0*a-beta;
        if (std::abs(gap)<=guard) return detail::solve_default_axis(g1,g2);
        if (gap< -guard) return Solution(Status::unreachable);
        const double product=std::sin(a+beta/2.0)*std::sin(a-beta/2.0);
        if (!(product>0.0) || !std::isfinite(product) ||
            guard*std::cos(beta/2.0)/std::sqrt(product)/rad>joint_budget_deg) {
            return detail::solve_default_axis(g1,g2);
        }
        const double u=std::sin(beta/2.0)/model.sa_;
        const double v=std::sqrt(product)/model.sa_;
        const double eta=2.0*std::atan2(v,u)/rad;
        const double delta=std::atan2(v,model.ca_*u)/rad;
        return detail::finish(psi,delta,eta);
    }
    const double go=model.outer_-beta,gi=beta-model.inner_;
    if (model.product_==0.0 || std::abs(go)<=guard || std::abs(gi)<=guard) {
        if (model.a_==model.b_) return detail::solve_equal_axis(g1,g2,model.a_);
        return detail::solve_unequal_axis(g1,g2,model.a_,model.b_);
    }
    if (go< -guard || gi< -guard) return Solution(Status::unreachable);
    const double U2=std::sin((model.outer_+beta)/2.0)*std::sin(go/2.0);
    const double V2=std::sin((beta+model.inner_)/2.0)*std::sin(gi/2.0);
    if (!(U2>0.0 && V2>0.0) || !std::isfinite(U2) || !std::isfinite(V2)) {
        if (model.a_==model.b_) return detail::solve_equal_axis(g1,g2,model.a_);
        return detail::solve_unequal_axis(g1,g2,model.a_,model.b_);
    }
    const double U=std::sqrt(U2),V=std::sqrt(V2);
    if ((guard/rad)*(1.0/U+1.0/V+1.0/rho)>joint_budget_deg) {
        if (model.a_==model.b_) return detail::solve_equal_axis(g1,g2,model.a_);
        return detail::solve_unequal_axis(g1,g2,model.a_,model.b_);
    }
    const double A=model.sa_*model.sin_difference_+2.0*model.ca_*V2;
    const double B=2.0*U*V;
    const double eta=2.0*std::atan2(U,V)/rad;
    const double delta=std::atan2(B,A)/rad;
    return detail::finish(psi,delta,eta);
}

} // namespace coneik_native
#endif // CONEIK_CORE_HPP
