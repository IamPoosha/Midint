#pragma once

#include <x86intrin.h>

#include <array>
#include <iostream>

// TODO: seperate into several files

#if __has_include("boost/multiprecision/cpp_int.hpp")
#include <boost/multiprecision/cpp_int.hpp>
#define YAY_WE_HAVE_BOOST
static_assert(sizeof(boost::multiprecision::limb_type) * __CHAR_BIT__ == 64);
#endif

namespace {
template <size_t sz>
struct midint;

template <size_t num_bits, class = std::enable_if_t<num_bits % 64 == 0>>
using Midint = midint<num_bits / 64>;

#ifdef YAY_WE_HAVE_BOOST
template <size_t num_bits>
using boost_uint = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<num_bits, num_bits, boost::multiprecision::unsigned_magnitude>>;

template <size_t num_bits>
using boost_int = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<num_bits, num_bits, boost::multiprecision::signed_magnitude>>;
#endif

template <typename boost_int_type>
void trim_boost(boost_int_type &var) {
    if (var.backend().size() > 1 and var.backend().limbs()[var.backend().size() - 1] == 0) {
        size_t new_size = var.backend().size() - 1;
        while (new_size > 1 and var.backend().limbs()[new_size - 1] == 0) {
            --new_size;
        }
        var.backend().resize(new_size, 0);
    }
}

// Number is stored as Big Endian.
// sz = number of limbs. Everything is modulo 2^(64*sz), i.e. 2's complement.
// The sign is implicit. By default, functionality is unisgned.
template <size_t sz>
struct midint {
    static_assert(0 < sz and sz < 128, "Notice midint expects to receive number of 64-limbs, not bits");
    static constexpr size_t num_bits = 64 * sz;
    static constexpr size_t num_limbs = sz;
    using uchar_t = unsigned char;
    using i64 = long long unsigned int;

    // uninitialized
    midint() {}
    midint(std::array<i64, sz> lo_to_hi) : L(std::move(lo_to_hi)) {}
    explicit midint(int64_t num) {
        L[0] = num;
        std::fill_n(L.data() + 1, sz - 1, num >> 63);
    }
    template <size_t sz2>
    explicit midint(midint<sz2> other) {
        static_assert(sz <= sz2);
        *this = other.template low_part<sz>();
    }
#ifdef YAY_WE_HAVE_BOOST
   private:
    explicit midint(i64 *boost_limbs, size_t size) {
        size = std::min(size, sz);
        std::copy_n(boost_limbs, size, L.begin());
        std::fill_n(L.begin() + size, sz - size, 0);
    }

   public:
    explicit midint(const std::string &str) : midint(boost_int<64 * sz>(str)) {}

    explicit midint(const boost_uint<64 * sz> &num) : midint((i64 *)num.backend().limbs(), (size_t)num.backend().size()) {}
    explicit midint(const boost_int<64 * sz> &num) : midint((i64 *)num.backend().limbs(), (size_t)num.backend().size()) {
        if (num.sign() == -1) {
            negate_inplace();
        }
    }
    explicit operator boost_uint<64 * sz>() const {
        boost_uint<64 * sz> res;
        // can we do without mallocs?
        res.backend().resize(sz, 0);
        std::copy_n(L.begin(), sz, (i64 *)res.backend().limbs());
        trim_boost(res);
        return res;
    }
    explicit operator boost_int<64 * sz>() const {
        boost_int<64 * sz> res;
        bool minus_sign = (int64_t(L[sz - 1]) >> 63);
        res.backend().resize(sz, 0);
        if (minus_sign) {
            midint copy = *this;
            copy.negate_inplace();
            std::copy_n(copy.L.begin(), sz, (i64 *)res.backend().limbs());
            res.backend().sign(-1);
        } else {
            std::copy_n(L.begin(), sz, (i64 *)res.backend().limbs());
        }
        trim_boost(res);
        return res;
    }
    // Notice this gives unsigned string.
    explicit operator std::string() const { return to_unsigned_string(); }

    std::string to_unsigned_string() const {
        std::stringstream ss;
        ss << boost_uint<64 * sz>(*this);
        return ss.str();
    }

    std::string to_signed_string() const {
        std::stringstream ss;
        ss << boost_int<64 * sz>(*this);
        return ss.str();
    }

    std::string to_limbs_string() const {
        std::stringstream ss;
        ss << "midint<" << sz << ">{";
        ss << L[0];
        for (size_t i = 1; i < sz; ++i) {
            ss << ", " << L[i];
        }
        ss << "}";
        return ss.str();
    }
#endif
    explicit operator int64_t() const { return L[0]; }

    template <size_t sz2>
    midint<sz2> high_part() const {
        static_assert(sz2 <= sz);
        midint<sz2> R;
        for (size_t i = 0; i < sz2; ++i) {
            R.L[i] = L[i + sz - sz2];
        }
        return R;
    }

    template <size_t sz2>
    midint<sz2> low_part() const {
        static_assert(sz2 <= sz);
        midint<sz2> R;
        for (size_t i = 0; i < sz2; ++i) {
            R.L[i] = L[i];
        }
        return R;
    }

    template <size_t extend_to>
    midint<extend_to> extend_unsigned() const {
        static_assert(extend_to >= sz);
        midint<extend_to> R;
        for (size_t i = 0; i < sz; ++i) {
            R.L[i] = L[i];
        }
        for (size_t i = sz; i < extend_to; ++i) {
            R.L[i] = 0;
        }
        return R;
    }
    template <size_t extend_to>
    midint<extend_to> extend_signed() const {
        static_assert(extend_to >= sz);
        midint<extend_to> R;
        for (size_t i = 0; i < sz; ++i) {
            R.L[i] = L[i];
        }
        if (int64_t(L[sz - 1]) >= 0) {
            for (size_t i = sz; i < extend_to; ++i) {
                R.L[i] = 0;
            }
        } else {
            for (size_t i = sz; i < extend_to; ++i) {
                R.L[i] = -1;
            }
        }
        return R;
    }

    template <size_t num_words64_to_shift_by>
    midint<sz + num_words64_to_shift_by> shift_left() const {
        midint<sz + num_words64_to_shift_by> R;
        for (size_t i = 0; i < sz; ++i) {
            R.L[i + num_words64_to_shift_by] = L[i];
        }
        for (size_t i = 0; i < num_words64_to_shift_by; ++i) {
            R.L[i] = 0;
        }
        return R;
    }

    uchar_t add_inplace_return_carry(const midint o) {
        uchar_t carry = 0;
        for (size_t i = 0; i < sz; ++i) {
            carry = _addcarry_u64(carry, L[i], o.L[i], &(L[i]));
        }
        return carry;
    }

    template <uchar_t initial_carry = 0>
    uchar_t sub_inplace_return_carry(const midint o) {
        uchar_t carry = initial_carry;
        for (size_t i = 0; i < sz; ++i) {
            carry = _subborrow_u64(carry, L[i], o.L[i], &(L[i]));
        }
        return carry;
    }

    midint &operator+=(const midint o) {
        if constexpr (sz == 1) {
            L[0] += o.L[0];
        } else if constexpr(sz == 2) {
            L[0] += o.L[0];
            L[1] += o.L[1] + (L[0] < o.L[0]);
        } else {
            this->add_inplace_return_carry(o);
        }
        return *this;
    }

    midint &operator-=(const midint o) {
        if constexpr (sz == 1) {
            L[0] -= o.L[0];
        } else if constexpr (sz == 2) {
            L[1] -= o.L[1] + (L[0] < o.L[0]);
            L[0] -= o.L[0];
        } else {
            this->sub_inplace_return_carry(o);
        }
        return *this;
    }

    midint &negate_inplace() {
        if constexpr (sz == 1) {
            L[0] = -L[0];
        } else if constexpr (sz == 2) {
            L[1] = (~L[1]) + (L[0] == 0);
            L[0] = -L[0];
        } else {
            uchar_t carry = 0;
            for (size_t i = 0; i < sz; ++i) {
                carry = _subborrow_u64(carry, 0, L[i], &(L[i]));
            }
        }
        return *this;
    }
#define BITWISE_OPERATOR(op)                      \
    midint &operator op(const midint o) {         \
        if constexpr (sz <= 8 and sz % 2 == 0) {  \
            auto LL = (__int128_t *)L.data();     \
            auto oLL = (__int128_t *)o.L.data();  \
            for (size_t i = 0; i < sz / 2; ++i) { \
                LL[i] op oLL[i];                  \
            }                                     \
        } else {                                  \
            for (size_t i = 0; i < sz; ++i) {     \
                L[i] op o.L[i];                   \
            }                                     \
        }                                         \
        return *this;                             \
    }
    BITWISE_OPERATOR(^=)
    BITWISE_OPERATOR(|=)
    BITWISE_OPERATOR(&=)
#undef BITWISE_OPERATOR

// disables compiler from loop unrolling.
// Might be beneficial in case where loop breaks early.
// From: Google Benchmark Framework
#define DoNotOptimize(var) asm volatile(""            \
                                        : "+r,m"(var) \
                                        :             \
                                        : "memory");

    bool operator<(const midint o) const {
        if constexpr (sz == 1) {
            return L[0] < o.L[0];
        } else if constexpr (sz == 2) {
            return *(__int128_t *)L.data() < *(__int128_t *)o.L.data();
        }
        for (size_t i = sz; i--;) {
            // DoNotOprimize(i)
            if (L[i] != o.L[i]) {
                return (L[i] < o.L[i]);
            }
        }
        return false;
    }
    bool operator>(const midint o) const { return o < *this; }
    bool operator<=(const midint o) const {
        if constexpr (sz == 1) {
            return L[0] <= o.L[0];
        } else if constexpr (sz == 2) {
            return *(__int128_t *)L.data() <= *(__int128_t *)o.L.data();
        } else if constexpr (sz <= 8) {
            return midint(*this).sub_inplace_return_carry<1>(o);
        }
        for (size_t i = sz; i--;) {
            if (L[i] != o.L[i]) {
                return (L[i] < o.L[i]);
            }
        }
        return true;
    }
    bool operator>=(const midint o) const { return o <= *this; }
    bool operator==(const midint o) const {
        if constexpr (sz == 1) {
            return L[0] == o.L[0];
        } else if constexpr (sz == 2) {
            return *(__int128_t *)L.data() == *(__int128_t *)o.L.data();
        }
        for (size_t i = 0; i < sz; ++i) {
            if (L[i] != o.L[i]) {
                return false;
            }
        }
        return true;
    }
    bool operator!=(const midint o) const {
        return not(*this == o);
    }

    template <size_t sz2>
    midint operator*=(const midint<sz2> B) {
        static_assert(sz2 <= sz);
        return (*this = this->multiply<sz, sz2>(B));
    }

    // performs a full multiplication.
    template <size_t sz2>
    midint<sz + sz2> full_multiply(const midint<sz2> B) const {
        return this->multiply<sz + sz2, sz2>(B);
    }

    template <size_t out_size, size_t sz2>
    midint<out_size> multiply(const midint<sz2> B) const {
        // TODO: experiment with karatsuba
        if constexpr (sz2 < sz) {
            return B.template multiply<out_size>(*this);
        } else {
            static_assert(out_size <= sz + sz2 and out_size >= sz2);
            if constexpr (sz > 8) {
                constexpr size_t _8 = 8;
                auto R1 = this->low_part<_8>().template multiply<std::min(out_size, _8 + sz2)>(B);
                auto R2 = this->high_part<sz - _8>().template multiply<out_size - _8>(
                    B.template low_part<std::min(out_size - _8, sz2)>());
                return R1.template extend_unsigned<out_size>() += R2.template shift_left<_8>();
            } else if constexpr (sz2 > 12) {
                constexpr size_t _8 = 8;
                auto R1 = B.template low_part<_8>().template multiply<std::min(out_size, _8 + sz)>(*this);
                auto R2 = B.template high_part<sz2 - _8>().template multiply<out_size - _8>(
                    this->template low_part<std::min(out_size - _8, sz)>());
                return R1.template extend_unsigned<out_size>() += R2.template shift_left<_8>();
            } else {
                return this->school_multiply<out_size, sz2>(B);
            }
        }
    }

    template <size_t out_size, size_t sz2>
    midint<out_size> polynomial_multiply(const midint<sz2> B) const {
        // TODO: write with intrinsics, compare to school.
        // From boost's multiplication.
        // Comba Multiplier - based on Paul Comba's Exponentiation cryptosystems on the IBM PC, 1990
        midint<out_size> R;

        static_assert(out_size <= sz + sz2 and out_size >= std::max(sz, sz2));
        constexpr size_t as = sz, bs = sz2, rs = out_size;

        __uint128_t carry = 0;
        __uint128_t temp = 0;
        constexpr size_t limb_bits = sizeof(i64) * 8;
        constexpr size_t lim = std::min<size_t>(rs, as + bs - 1);
        for (int r = 0; r < lim; ++r) {
            i64 overflow = 0;
            size_t i = std::min<size_t>(r, as - 1);
            size_t j = r - i;
            size_t k = std::min<size_t>(i + 1, bs - j);

            temp = carry;
            carry += static_cast<__uint128_t>(L[i]) * (B.L[j]);
            overflow += (carry < temp);
            while (--k) {
                temp = carry;
                carry += static_cast<__uint128_t>(L[--i]) * (B.L[++j]);
                overflow += (carry < temp);
            }
            R.L[r] = static_cast<i64>(carry);
            carry = (static_cast<__uint128_t>(overflow) << limb_bits) | (carry >> limb_bits);
        }
        if constexpr (rs == as + bs) {
            R.L[out_size - 1] = carry;
        }
        return R;
    }

    template <size_t out_size, size_t sz2>
    midint<out_size> school_multiply(const midint<sz2> B) const {
        static_assert(out_size <= sz + sz2 and out_size >= sz2 and sz2 >= sz);
        midint<out_size> R;
        if constexpr (sz == 1) {
            // simpler logic in this case
            if constexpr (out_size == 1) {
                R.L[0] = L[0] * B.L[0];
                return R;
            }
            R.L[0] = _mulx_u64(L[0], B.L[0], &R.L[1]);
            uchar_t c = 0;
            for (size_t i = 1; i < sz2; ++i) {
                i64 temp;
                c = _addcarry_u64(c, R.L[i], _mulx_u64(L[0], B.L[i], &temp), &R.L[i]);
                if(i+1 < out_size) {
                    R.L[i+1] = temp;
                }
            }
            if constexpr (sz2 < out_size) {
                R.L[sz2] += c;
            }
        } else {
            R.L[0] = _mulx_u64(L[0], B.L[0], &R.L[1]);
            uchar_t c = 0;
            for (size_t i = 1; i < sz2; ++i) {
                i64 temp;
                c = _addcarry_u64(c, R.L[i], _mulx_u64(L[0], B.L[i], &temp), &R.L[i]);
                if (i + 1 < out_size) {
                    R.L[i + 1] = temp;
                }
            }
            if (sz2 < out_size) {
                R.L[sz2] += c;
            }
            for (size_t j = 1; j < std::min(sz, out_size - 1); ++j) {
                size_t lim = std::min(j + sz2, out_size);
                i64 d;
                uchar_t c1 = 0;
                c = _addcarry_u64(0, R.L[j], _mulx_u64(L[j], B.L[0], &d), &R.L[j]);
                for (size_t i = j + 1; i + 1 < lim; ++i) {
                    c1 = _addcarry_u64(c1, R.L[i], d, &R.L[i]);
                    c = _addcarry_u64(c, R.L[i], _mulx_u64(L[j], B.L[i - j], &d), &R.L[i]);
                }
                c1 = _addcarry_u64(c1, R.L[lim - 1], d, &R.L[lim - 1]);
                if (lim == out_size) {
                    _addcarry_u64(c, R.L[lim - 1], L[j] * B.L[lim - j - 1], &R.L[lim - 1]);
                } else {
                    c = _addcarry_u64(c, R.L[lim - 1], _mulx_u64(L[j], B.L[lim - j - 1], &d), &R.L[lim - 1]);
                    if (lim == j + sz2) {
                        _addcarry_u64(c, c1, d, &R.L[lim]);
                    } else {
                        _addcarry_u64(c, R.L[lim], d, &R.L[lim]);
                        R.L[lim] += c1;
                    }
                }
            }
            // the loop above is not designed for the case j+1 == out_size, hence this case
            if constexpr (out_size == sz) {
                R.L[sz - 1] += L[sz - 1] * B.L[0];
            }
        }
        return R;
    }

    // non-full squaring
    // inspired by paper "Speeding Up Big-Numbers Squaring" due to Gueron, Krasnov.
    // Equivalent to: return this->half_multiply(*this);
    midint<sz> half_square() const {
        return this->template square<sz>();
    }

    template <size_t out_size = 2 * sz>
    midint<out_size> square() const {
        if constexpr(out_size < sz) {
            return this->low_part<out_size>().template square<out_size>();
        } else if constexpr (out_size > 8) {
            constexpr size_t h = sz / 5 * 2;
            auto lo = this->low_part<h>();
            auto hi = this->high_part<sz - h>();
            auto R1 = lo.template square<std::min(out_size, 2*h)>();
            auto res = R1.template extend_unsigned<out_size>();
            auto R2 = hi.template multiply<std::min(out_size - h, sz)>(lo);
            constexpr size_t R2_limbs = R2.num_limbs;
            uchar_t c = R2.add_inplace_return_carry(R2);
            auto _2R2 = R2.template extend_unsigned<out_size - h>();
            if constexpr(out_size - h > R2_limbs) {
                _2R2.L[R2_limbs] = c;
            }
            res += _2R2.template shift_left<h>();
            if constexpr(out_size > 2*h) {
                auto R3 = hi.template square<out_size - 2*h>();
                res += R3.template shift_left<2*h>();
            }
            return res;
        } else {
            return this->template two_step_square<out_size>();
        }
    }
    
    std::array<i64, sz> L;

   private:
template <size_t out_size = 2 * sz>
    midint<out_size> two_step_square() const {
        static_assert(out_size <= 2 * sz and out_size >= sz);

        if constexpr (sz == 1) {
            midint<out_size> R;
            if constexpr (out_size == 2) {
                R.L[0] = _mulx_u64(L[0], L[0], &R.L[1]);
            } else {
                R.L[0] = L[0] * L[0];
            }
            return R;
        } else {
            // L[i] x L[j], j > i
            midint<out_size> R = this->square_aux<out_size>();

            // multiplying by 2 and adding L[i] x L[i]
            i64 d;
            uchar_t c = 0;
            uchar_t c1 = 0;
            R.L[0] = _mulx_u64(L[0], L[0], &d);
            for (size_t i = 1; i < out_size; ++i) {
                c1 = _addcarry_u64(c1, R.L[i], R.L[i], &R.L[i]);
                if (i % 2 == 1) {
                    c = _addcarry_u64(c, R.L[i], d, &R.L[i]);
                } else {
                    c = _addcarry_u64(c, R.L[i], _mulx_u64(L[i / 2], L[i / 2], &d), &R.L[i]);
                }
            }
            return R;
        }
    }


    template <size_t out_size>
    midint<out_size> square_aux() const {
        // TODO: write naively with int128_t and see if improves...
        // L[i] x L[j], j > i
        static_assert(out_size <= 2 * sz and out_size >= sz);
        midint<out_size> R;
        R.L[0] = R.L[out_size - 1] = 0;
        if constexpr (sz > 1) {
            {
                // L[0] x L[i], i > 0
                R.L[1] = _mulx_u64(L[0], L[1], &R.L[2]);
                uchar_t c = 0;
                for (size_t i = 2; i < sz; ++i) {
                    i64 temp;
                    c = _addcarry_u64(c, R.L[i], _mulx_u64(L[0], L[i], &temp), &R.L[i]);
                    if (i + 1 < out_size) {
                        R.L[i + 1] = temp;
                    }
                }
                if constexpr (sz < out_size) {
                    R.L[sz] += c;
                }
            }
            for (size_t i = 1; i + 1 < sz; ++i) {
                // L[i] x L[j], j > i
                i64 d = 0;
                uchar_t c = 0;
                uchar_t c1 = 0;
                size_t lim = std::min(sz, out_size - i);
                for (size_t j = i + 1; j < lim; ++j) {
                    c1 = _addcarry_u64(c1, R.L[i + j], d, &R.L[i + j]);
                    c = _addcarry_u64(c, R.L[i + j], _mulx_u64(L[i], L[j], &d), &R.L[i + j]);
                }
                if (i + lim < out_size) {
                    _addcarry_u64(c, c1, d, &R.L[i + lim]);
                }
            }
        }
        return R;
    }
};

#ifdef YAY_WE_HAVE_BOOST
template <size_t sz>
std::ostream &operator<<(std::ostream &os, const midint<sz> &num) {
    return os << num.to_unsigned_string();
}
#endif

}  // namespace

#undef YAY_WE_HAVE_BOOST