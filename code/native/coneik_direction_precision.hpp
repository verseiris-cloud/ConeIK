#ifndef CONEIK_DIRECTION_PRECISION_HPP
#define CONEIK_DIRECTION_PRECISION_HPP
// Optional raw-XYZ 50/80 then 100/200/400 decimal-digit acceleration. Identities
// and unconfirmed results belong to the caller's full symbolic ladder.
#include "coneik_core.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/constants/constants.hpp>
#include <cstdint>
#include <cstring>
#include <optional>

namespace coneik_native {
namespace direction_precision_detail {
template<unsigned D> using Float=boost::multiprecision::number<
    boost::multiprecision::cpp_bin_float<D,boost::multiprecision::digit_base_10>>;

// A binary64 value is exactly significand*2^exponent. No floating square,
// normalization, or tolerance participates in the equality predicate.
struct Dyadic { std::uint64_t significand; int exponent; };
inline Dyadic unsigned_dyadic(double value) noexcept {
    static_assert(sizeof(double)==sizeof(std::uint64_t) &&
        std::numeric_limits<double>::is_iec559,"IEEE binary64 required");
    std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));
    const unsigned exponent=static_cast<unsigned>((bits>>52)&0x7ffu);
    const std::uint64_t fraction=bits&((std::uint64_t(1)<<52)-1);
    return {exponent ? fraction|(std::uint64_t(1)<<52) : fraction,
            exponent ? static_cast<int>(exponent)-1075 : -1074};
}
inline bool exact_pythagorean(double x,double y,double z) {
    using boost::multiprecision::cpp_int;
    const Dyadic xx=unsigned_dyadic(x),yy=unsigned_dyadic(y),zz=unsigned_dyadic(z);
    int exponent=xx.exponent;
    if (yy.significand) exponent=std::min(exponent,yy.exponent);
    if (zz.significand) exponent=std::min(exponent,zz.exponent);
    const auto square=[exponent](Dyadic v) -> cpp_int {
        if (!v.significand) return cpp_int(0);
        cpp_int result=cpp_int(v.significand)*v.significand;
        result<<=static_cast<unsigned>(2*(v.exponent-exponent));
        return result;
    };
    return square(xx)==square(yy)+square(zz);
}
template<unsigned D> inline double azimuth(double y,double z) {
    using T=Float<D>;
    static const T degrees=T(180)/boost::math::constants::pi<T>();
    return static_cast<double>(atan2(T(z),T(y))*degrees);
}
struct ExactBoundary {
    int code=-2; // 1 outer, 2 inner; -2 is not an accepted exact identity.
    std::array<double,4> q{};
    double dps=std::numeric_limits<double>::quiet_NaN();
};
class Exact45Geometry {
    int code_=0;bool reversed_=false;
public:
    Exact45Geometry(double a,double b) noexcept:reversed_(a<b) {
        const auto outer=detail::two_sum(a,b),inner=detail::two_sum(a,-b);
        if (outer.hi==45.0 && outer.lo==0.0) code_=1;
        else if (std::abs(inner.hi)==45.0 && inner.lo==0.0) code_=2;
    }
    ExactBoundary solve(double x,double y,double z) const {
        ExactBoundary out;
        if (!code_ || !(x>0) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return out;
        if (!exact_pythagorean(x,y,z)) return out;
        // Only classification is exact. Azimuth is finite precision, checked
        // at two levels and by the caller's unchanged independent FK residual.
        const double lo=azimuth<50>(y,z),hi=azimuth<80>(y,z);
        if (!std::isfinite(lo) || !std::isfinite(hi) ||
            std::abs(detail::wrap_degrees(hi-lo))>1e-12) return out;
        double psi=hi;
        if (code_==2 && reversed_) psi-=180;
        out.q={{psi,code_==1?0.0:-180.0,psi,code_==1?0.0:180.0}};
        for (int j=0;j<4;j+=2) {
            out.q[j]=detail::wrap_degrees(out.q[j]);
            const double actual=out.q[j]+out.q[j+1];
            if (actual>180) out.q[j+1]-=360;
            else if (actual< -180) out.q[j+1]+=360;
        }
        out.code=code_;out.dps=80;return out;
    }
};
struct Trial {
    Status status=Status::needs_precision;
    std::array<double,4> q{}; // Raw degrees, before representative adjustment.
};
template<unsigned D> struct Geometry {
    using T=Float<D>;
    T rad,a,b,outer,inner,sa,ca,sd,margin,outer_s2,inner_s2;
    Geometry(double aa,double bb):rad(boost::math::constants::pi<T>()/180),
        a(T(aa)*rad),b(T(bb)*rad),outer(a+b),inner(abs(a-b)),
        sa(sin(a)),ca(cos(a)),sd(sin(a-b)),
        margin(pow(T(10),-static_cast<int>(D)+20)),
        outer_s2(pow(sin(outer/2),2)),inner_s2(pow(sin(inner/2),2)) {}
    Trial solve(double xx,double yy,double zz) const {
        // Promote original binary64 values exactly BEFORE all arithmetic.
        // cpp_bin_float's exponent range accommodates squared double extrema.
        const T x(xx),y(yy),z(zz);
        const T beta=atan2(sqrt(y*y+z*z),x),go=outer-beta,gi=beta-inner;
        Trial out;
        if (abs(go)<=margin || abs(gi)<=margin) return out;
        if (go<0 || gi<0) {out.status=Status::unreachable;return out;}
        // sin((A+B)/2)*sin((A-B)/2) = sin^2(A/2)-sin^2(B/2).
        // Geometry-only sine values are cached. Rationalize the positive-x
        // half-angle expression to retain arbitrarily small transverse values;
        // for negative x, n-x has no subtractive cancellation.
        // Angular gap guards above still own boundary uncertainty. Nonpositive
        // factors below and disagreement between precision levels fail closed.
        const T rho2=y*y+z*z,n=sqrt(x*x+rho2);
        const T half_s2=x>=0 ? rho2/(2*n*(n+x)) : (n-x)/(2*n);
        const T U2=outer_s2-half_s2;
        const T V2=half_s2-inner_s2;
        if (!(U2>0 && V2>0)) return out;
        const T U=sqrt(U2),V=sqrt(V2),eta=2*atan2(U,V)/rad;
        const T delta=atan2(2*U*V,sa*sd+2*ca*V2)/rad,psi=atan2(z,y)/rad;
        out.q={{static_cast<double>(psi+delta),static_cast<double>(-eta),
                static_cast<double>(psi-delta),static_cast<double>(eta)}};
        for (double q:out.q) if (!std::isfinite(q)) return Trial{};
        out.status=Status::success;return out;
    }
};
} // namespace direction_precision_detail

class DirectionPrecisionModel {
    double a_,b_,last_dps_=std::numeric_limits<double>::quiet_NaN();
    std::optional<direction_precision_detail::Geometry<50>> fast_low_;
    std::optional<direction_precision_detail::Geometry<80>> fast_high_;
    std::optional<direction_precision_detail::Geometry<100>> low_;
    std::optional<direction_precision_detail::Geometry<200>> middle_;
    std::optional<direction_precision_detail::Geometry<400>> high_;
public:
    DirectionPrecisionModel(double a,double b):a_(a),b_(b) {const Model checked(a,b);(void)checked;}
    Solution solve(double x,double y,double z) {
        last_dps_=std::numeric_limits<double>::quiet_NaN();
        if (!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z)||(x==0&&y==0&&z==0))
            return Solution(Status::invalid_input);
        // Structural axis identities retain the caller's established handler.
        if (y==0 && z==0) return Solution{};
        // This performance gate uses only finite, nonzero component exponents;
        // it is not an accuracy certificate or a geometry restriction. Avoid
        // low-tier work for extreme transverse ratios, without dividing doubles
        // or changing the original values passed to multiprecision arithmetic.
        const double scale=std::max(std::abs(x),std::max(std::abs(y),std::abs(z)));
        const double transverse=std::max(std::abs(y),std::abs(z));
        if (std::ilogb(scale)-std::ilogb(transverse)<=80) {
            if (!fast_low_) fast_low_.emplace(a_,b_);
            const auto fast_lo=fast_low_->solve(x,y,z);
            if (fast_lo.status!=Status::needs_precision) {
                if (!fast_high_) fast_high_.emplace(a_,b_);
                const auto fast_hi=fast_high_->solve(x,y,z);
                auto fast_result=agree(fast_lo,fast_hi);
                if (confirmed(fast_result)) {last_dps_=80;return fast_result;}
            }
        }
        // A skipped or unconfirmed low pair continues the existing ladder.
        if (!low_) low_.emplace(a_,b_);
        const auto lo=low_->solve(x,y,z);
        if (!middle_) middle_.emplace(a_,b_);
        const auto mid=middle_->solve(x,y,z);
        auto result=agree(lo,mid);
        if (confirmed(result)) {last_dps_=200;return result;}
        // No trustworthy 200-digit result => a 400-only result cannot agree.
        if (mid.status==Status::needs_precision) return Solution{};
        if (!high_) high_.emplace(a_,b_);
        const auto hi=high_->solve(x,y,z);
        result=agree(mid,hi);
        if (confirmed(result)) last_dps_=400;
        return result;
    }
    double precision_dps() const noexcept {return last_dps_;}
private:
    static bool confirmed(const Solution& s) noexcept {return s.status==Status::success||s.status==Status::unreachable;}
    static Solution agree(const direction_precision_detail::Trial& lo,const direction_precision_detail::Trial& hi) {
        if (lo.status!=hi.status) return Solution{};
        if (hi.status==Status::unreachable) return Solution(Status::unreachable);
        if (hi.status!=Status::success) return Solution{};
        // Adaptive agreement is empirically validated, not an interval proof.
        for (std::size_t j=0;j<4;++j)
            if (std::abs(detail::wrap_degrees(hi.q[j]-lo.q[j]))>1e-12) return Solution{};
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
