#ifndef CONEIK_PRECISION_HPP
#define CONEIK_PRECISION_HPP
// Optional 50/80 then 100/200 decimal-digit binary multiprecision layer. This header
// depends on Boost, while the ordinary coneik_core.hpp remains standalone.
#include "coneik_core.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <optional>

namespace coneik_native {
namespace precision_detail {
template<unsigned D> using Float = boost::multiprecision::number<
    boost::multiprecision::cpp_bin_float<D,boost::multiprecision::digit_base_10>>;

struct Trial {
    Status status=Status::needs_precision;
    std::array<double,4> q{}; // Raw representatives, before binary64 wrapping.
};

// Exact identities continue to use the established MATLAB handler. Test
// exact two-term degree sums, never a rounded a+b or a boundary tolerance.
inline bool identity(double p,double y,double a,double b) noexcept {
    if ((p==0 && y==0)||(p==-180 && y==-180)) return true;
    detail::Expansion beta;
    if (std::abs(p)==90 || std::abs(y)==90) beta={90,0};
    else if (p==0) beta={std::abs(y),0};
    else if (y==0) beta={std::abs(p),0};
    else if (p==-180) beta=detail::two_sum(180,-std::abs(y));
    else if (y==-180) beta=detail::two_sum(180,-std::abs(p));
    else return false;
    const auto outer=detail::two_sum(a,b);
    const auto inner=detail::two_sum(std::max(a,b),-std::min(a,b));
    return (beta.hi==outer.hi && beta.lo==outer.lo)||
           (beta.hi==inner.hi && beta.lo==inner.lo);
}

template<unsigned D> struct Geometry {
    using T=Float<D>;
    T rad,a,b,outer,inner,sa,ca,sd,margin;
    bool is_default;
    Geometry(double aa,double bb):
        rad(boost::math::constants::pi<T>()/180),
        a(T(aa)*rad),b(T(bb)*rad),outer(a+b),inner(abs(a-b)),
        sa(sin(a)),ca(cos(a)),sd(sin(a-b)),
        margin(pow(T(10),-static_cast<int>(D)+20)),
        is_default(aa==20 && bb==20) {}

    Trial solve(double p,double y) const {
        // Direct binary64 -> higher precision binary construction is exact.
        // In particular do not initialize through a rounded direction or text.
        const T pp=T(p)*rad,yy=T(y)*rad;
        const T cp=cos(pp),x=cp*cos(yy),ry=cp*sin(yy),rz=-sin(pp);
        const T beta=atan2(sqrt(ry*ry+rz*rz),x);
        const T go=outer-beta,gi=beta-inner;
        Trial out;
        if (abs(go)<=margin || (is_default ? beta==0 : abs(gi)<=margin)) return out;
        if (go<0 || (!is_default && gi<0)) {
            out.status=Status::unreachable;return out;
        }
        T eta,delta;
        if (is_default) {
            // Retain the default symbolic solver's expression sequence.
            const T h=beta/2,product=sin(a+h)*sin(a-h);
            if (!(product>0)) return out;
            const T u=sin(h)/sa,v=sqrt(product)/sa;
            eta=2*atan2(v,u)/rad;delta=atan2(v,ca*u)/rad;
        } else {
            const T U2=sin((outer+beta)/2)*sin(go/2);
            const T V2=sin((beta+inner)/2)*sin(gi/2);
            if (!(U2>0 && V2>0)) return out;
            const T U=sqrt(U2),V=sqrt(V2);
            eta=2*atan2(U,V)/rad;
            delta=atan2(2*U*V,sa*sd+2*ca*V2)/rad;
        }
        const T psi=atan2(rz,ry)/rad;
        out.q={{static_cast<double>(psi+delta),static_cast<double>(-eta),
                static_cast<double>(psi-delta),static_cast<double>(eta)}};
        for (double q:out.q) if (!std::isfinite(q)) return Trial{};
        out.status=Status::success;return out;
    }
};
} // namespace precision_detail

class PrecisionModel {
    double a_,b_;
    double last_dps_=std::numeric_limits<double>::quiet_NaN();
    std::optional<precision_detail::Geometry<50>> medium_low_;
    std::optional<precision_detail::Geometry<80>> medium_high_;
    std::optional<precision_detail::Geometry<100>> low_;
    std::optional<precision_detail::Geometry<200>> high_;
public:
    PrecisionModel(double a,double b):a_(a),b_(b) {
        // Reuse the public mathematical domain validation.
        const Model checked(a,b);(void)checked;
    }
    Solution solve(double pitch,double yaw) {
        last_dps_=std::numeric_limits<double>::quiet_NaN();
        if (!std::isfinite(pitch)||!std::isfinite(yaw)) return Solution(Status::invalid_input);
        const double p=detail::wrap_degrees(pitch),y=detail::wrap_degrees(yaw);
        if (precision_detail::identity(p,y,a_,b_)) return Solution{};
        // Bounded geometry and ordinary-magnitude input only. Each trial
        // retains original doubles, rejects uncertain boundary signs, and
        // both levels must agree on classification AND all four angles.
        // This is validated adaptive precision, not a formal interval bound.
        const double scale=std::hypot(p,y);
        if (a_>=1 && a_<=80 && b_>=1 && b_<=80 && scale>=1e-18) {
            if (!medium_low_) medium_low_.emplace(a_,b_);
            const auto lo=medium_low_->solve(p,y);
            if (lo.status!=Status::needs_precision) {
                if (!medium_high_) medium_high_.emplace(a_,b_);
                const auto hi=medium_high_->solve(p,y);
                auto result=agree(lo,hi);
                if (result.status==Status::success || result.status==Status::unreachable) {
                    last_dps_=80;return result;
                }
            }
        }
        if (!low_) low_.emplace(a_,b_);
        const auto lo=low_->solve(p,y);
        // Without a trustworthy 100-digit trial, this two-level stage cannot
        // establish agreement. Leave the full symbolic ladder to its caller.
        if (lo.status==Status::needs_precision) return Solution{};
        if (!high_) high_.emplace(a_,b_);
        const auto hi=high_->solve(p,y);
        auto result=agree(lo,hi);
        if (result.status==Status::success || result.status==Status::unreachable) last_dps_=200;
        return result;
    }
    double precision_dps() const noexcept {return last_dps_;}
private:
    static Solution agree(const precision_detail::Trial& lo,const precision_detail::Trial& hi) {
        if (lo.status!=hi.status) return Solution{};
        if (hi.status==Status::unreachable) return Solution(Status::unreachable);
        if (hi.status!=Status::success) return Solution{};
        for (std::size_t j=0;j<4;++j) {
            if (std::abs(detail::wrap_degrees(hi.q[j]-lo.q[j]))>1e-12) return Solution{};
        }
        Solution out(Status::success);out.q=hi.q;
        for (int j=0;j<4;j+=2) {
            out.q[j]=detail::wrap_degrees(out.q[j]);
            const double actual=out.q[j]+out.q[j+1];
            if (actual>180) out.q[j+1]-=360;
            else if (actual< -180) out.q[j+1]+=360;
        }
        return out;
    }
};
} // namespace coneik_native
#endif
