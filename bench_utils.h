#pragma once

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <type_traits>

// in seconds since epoch
double get_time() {
    return double(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count()) /
           1E9;
}

template <typename... Args>
void println(Args&&... args) {
    int A[] = {0, ((std::cout << std::forward<Args>(args) << " "), 0)...};
    std::cout << std::endl;
}

template <typename Container>
auto sum(const Container& container) {
    typename std::remove_const<typename Container::value_type>::type sm(0);
    for (const auto& val : container) {
        sm += val;
    }
    return sm;
}

std::mt19937_64& get_rng(size_t seed = 0) {
    static std::mt19937_64 rng(0xDEFAUL);
    if (seed) {
        rng.seed(seed);
    }
    return rng;
}

#define BENCH_BASE(func, num_iters, name)                                                                                    \
    ({                                                                                                                       \
        double start_time = get_time();                                                                                      \
        size_t res = func(num_iters);                                                                                        \
        double end_time = get_time();                                                                                        \
        double total_time = end_time - start_time;                                                                           \
        double op_time = (total_time + 1E-50 * res) / num_iters;                                                             \
        if (std::string(name).size() > 0) {                                                                                  \
            println("ops per second: ", name, ": 2^", -std::log2(op_time), "   ||    time:", total_time, "  ||   iters: 2^", \
                    std::log2(num_iters));                                                                                   \
        }                                                                                                                    \
        op_time;                                                                                                             \
    })
#define BENCH(func, num_iters, name) BENCH_BASE(func, (num_iters), name)
#define BENCH_QUIET(func, num_iters) BENCH_BASE(func, num_iters, "")

#define COMPARE(func1, func2, inp)                                                                           \
    ({                                                                                                       \
        size_t input = (inp);                                                                                \
        get_rng(0xBAD6EED);                                                                                  \
        auto res1 = func1(input);                                                                            \
        get_rng(0xBAD6EED);                                                                                  \
        auto res2 = func2(input);                                                                            \
        if (res1 == res2) {                                                                                  \
            println("CHECK SUCCESS: ", #func1, "==", #func2, "  ||   input=", input, " output=", res1);      \
        } else {                                                                                             \
            println("CHECK FAIL: ", #func1, "!=", #func2, "on input=", input, "   ||  values=", res1, res2); \
        }                                                                                                    \
    })

#define COMPARE_EQ(val1, val2)                                                                         \
    ({                                                                                                 \
        if (val1 == val2) {                                                                            \
            println("CHECK SUCCESS: ", #val1, "==", #val2, "  ||   value=", val1);                     \
        } else {                                                                                       \
            println("CHECK FAIL: ", #val1, "!=", #val2, "where first=", val1, "   ||  second=", val2); \
        }                                                                                              \
    })

uint64_t perform_int64_muls(size_t num_iters) {
    size_t x = 1;
    for (size_t i = 1; i <= num_iters; ++i) {
        x *= i;
    }
    return x;
}

double get_num_cycles_per_second(size_t num_iters) {
    double op_time = BENCH_QUIET(perform_int64_muls, num_iters);
    return 1. / op_time;
}

//
//
//
//
// specific things for our benchmarks:
//
//
//
//
#include "external/ttmath/ttmath.h"
#include "midint.h"

template <size_t num_bits, class = std::enable_if_t<num_bits % 64 == 0>>
using ttmath_uint = ttmath::UInt<num_bits / 64>;

template <typename T>
void random_init(T& var_to_init) = delete;

template <size_t sz>
void random_init(midint<sz>& var_to_init) {
    for (size_t i = sz; i != 0; --i) {
        var_to_init.L[i - 1] = get_rng()();
    }
    // So that things will not zero out deliberately.
    // Boost exploits things like that.
    var_to_init.L[0] |= 1;
}

template <size_t sz>
void random_init(ttmath::UInt<sz>& var_to_init) {
    midint<sz> temp;
    random_init(temp);
    for (size_t i = 0; i < sz; ++i) {
        var_to_init.table[i] = temp.L[i];
    }
}

template <unsigned int num_bits>
void random_init(boost_uint<num_bits>& var_to_init) {
    // synchronized with midint's random
    constexpr size_t num_limbs = num_bits / 64;
    midint<num_limbs> temp;
    random_init(temp);
    var_to_init = boost_uint<num_bits>(temp);
}

template <typename T>
size_t hashit(const T& x) = delete;

template <size_t size, size_t num_limbs>
size_t hashit(const std::array<midint<num_limbs>, size>& x) {
    return size_t(boost_uint<midint<num_limbs>::num_bits>(sum(x)) % ~size_t(0));
}

template <size_t size, unsigned int num_bits>
size_t hashit(const std::array<boost_uint<num_bits>, size>& x) {
    return size_t(sum(x) % ~size_t(0));
}

template <size_t size, size_t num_limbs>
size_t hashit(const std::array<ttmath::UInt<num_limbs>, size>& x) {
    std::array<midint<num_limbs>, size> temp;
    for (size_t j = 0; j < size; ++j) {
        for (size_t i = 0; i < num_limbs; ++i) {
            temp[j].L[i] = x[j].table[i];
        }
    }
    return hashit(temp);
}

template <typename T>
struct helper {};

template <size_t sz>
struct helper<midint<sz>> {
    static constexpr size_t num_bits = midint<sz>::num_bits;
    static constexpr size_t num_limbs = num_bits / 64;
    static constexpr bool is_boost_int = false;
    static constexpr bool is_midint = true;
};

template <>
struct helper<boost_uint<64>> {
    // we use this specialization due to a bug of g++
    static constexpr size_t num_bits = 64;
    static constexpr size_t num_limbs = 64 / 64;
    static constexpr bool is_boost_int = true;
    static constexpr bool is_midint = false;
};

template <>
struct helper<boost_uint<128>> {
    // again, g++ bug!
    static constexpr size_t num_bits = 128;
    static constexpr size_t num_limbs = 128 / 64;
    static constexpr bool is_boost_int = true;
    static constexpr bool is_midint = false;
};

template <>
struct helper<boost_uint<192>> {
    // again, g++ bug!
    static constexpr size_t num_bits = 192;
    static constexpr size_t num_limbs = 192 / 64;
    static constexpr bool is_boost_int = true;
    static constexpr bool is_midint = false;
};

template <>
struct helper<boost_uint<256>> {
    // again, g++ bug!
    static constexpr size_t num_bits = 256;
    static constexpr size_t num_limbs = 256 / 64;
    static constexpr bool is_boost_int = true;
    static constexpr bool is_midint = false;
};

template <>
struct helper<boost_uint<384>> {
    // again, g++ bug!
    static constexpr size_t num_bits = 384;
    static constexpr size_t num_limbs = 384 / 64;
    static constexpr bool is_boost_int = true;
    static constexpr bool is_midint = false;
};

template <>
struct helper<boost_uint<512>> {
    // again, g++ bug!
    static constexpr size_t num_bits = 512;
    static constexpr size_t num_limbs = 512 / 64;
    static constexpr bool is_boost_int = true;
    static constexpr bool is_midint = false;
};

template <>
struct helper<boost_uint<768>> {
    // again, g++ bug!
    static constexpr size_t num_bits = 768;
    static constexpr size_t num_limbs = 768 / 64;
    static constexpr bool is_boost_int = true;
    static constexpr bool is_midint = false;
};

template <>
struct helper<boost_uint<1024>> {
    // again, g++ bug!
    static constexpr size_t num_bits = 1024;
    static constexpr size_t num_limbs = 1024 / 64;
    static constexpr bool is_boost_int = true;
    static constexpr bool is_midint = false;
};

// template <size_t bit_count>
// struct helper<boost_uint<bit_count>> {
//     static constexpr size_t num_bits = bit_count;
//     static constexpr size_t num_limbs = bit_count / 64;
//     static constexpr bool is_boost_int = true;
//     static constexpr bool is_midint = false;
// };

template <size_t sz>
struct helper<ttmath::UInt<sz>> {
    static constexpr size_t num_bits = sz * 64;
    static constexpr size_t num_limbs = num_bits / 64;
    static constexpr bool is_boost_int = false;
    static constexpr bool is_midint = false;
};

#define preamble()                                         \
    size_t operator()(size_t iters) const {                \
        constexpr size_t reps = 4;                         \
        constexpr size_t num_bits = helper<T>::num_bits;   \
        constexpr size_t num_limbs = helper<T>::num_limbs; \
        std::array<T, 2 * reps> regs;                      \
        size_t num_iters = iters / reps;                   \
        for (size_t j = 0; j < regs.size(); ++j) {         \
            random_init(regs[j]);                          \
        }  //

#define epilogue()       \
    return hashit(regs); \
    }

template <size_t num_bits>
using Midint_ = Midint<num_bits>;

template <size_t num_bits>
using ttmath_uint_ = ttmath_uint<num_bits>;