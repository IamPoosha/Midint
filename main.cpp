// runme:
//    clang++-14 main.cpp -O3 -march=native -o wide -std=c++17 && ./wide 1
// assembly view:
//    clang++-14 main.cpp -O3 -march=native -S -o wide.s -std=c++17

#include "bench_utils.h"

template <typename T>
struct add {
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            regs[j] += regs[j + reps];
        }
    }
    epilogue();
};

template <typename T>
struct half_mult {
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            regs[j] *= regs[j + reps];
        }
    }
    epilogue();
};

template <typename T, typename = void>
struct full_mult {};

template <typename T>
using full_mult_ = full_mult<T>;

template <typename T>
struct full_mult<T, std::enable_if_t<helper<T>::is_midint>> {
    // midint version
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            Midint<2 * num_bits> res = regs[j].template full_multiply(regs[j + reps]);
            regs[j] = res.template high_part<num_bits>();
            regs[j] ^= res.template low_part<num_bits>();
        }
    }
    epilogue();
};

template <typename T>
struct full_mult<T, std::enable_if_t<helper<T>::is_boost_int>> {
    // boost_uint version
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            boost_uint<2 * num_bits> res = boost_uint<2 * num_bits>(regs[j]) * regs[j + reps];
            res ^= res >> num_bits;
            regs[j] = boost_uint<num_bits>(res);
        }
    }
    epilogue();
};

template <typename T>
struct xorit {
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            regs[j] ^= regs[j + reps];
        }
    }
    epilogue();
};

template <typename T>
struct less_than {
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            // usually compares only top limb
            if (regs[j] < regs[j + reps]) {
                regs[j] ^= regs[j + reps];
            }
        }
    }
    epilogue();
};

template <typename T>
struct equals_to {
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            // true only O(1) times. Usually compares all limbs
            if (regs[j] != regs[j + reps]) {
                regs[j] |= regs[j + reps];
                regs[j + reps] += (T(regs[j]) *= regs[j + reps]);
                regs[j + reps] += regs[j];
            }
        }
    }
    epilogue();
};

template <typename T, typename = void>
struct half_square {
    // midint version
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            // apparently for 64 x 64 we lose the ability to use ymm/zmm
            // thus giving some slowdown.
            regs[j] *= regs[j];
            // fast and prevents regs[j] from being 1.
            regs[j] ^= regs[j + reps];
        }
    }
    epilogue();
};

template <typename T>
using half_square_ = half_square<T>;

template <typename T>
struct half_square<T, std::enable_if_t<helper<T>::is_midint>> {
    // midint version
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            // apparently for 64 x 64 we lose the ability to use ymm/zmm
            // thus giving some slowdown.
            regs[j] = regs[j].half_square();
            // fast and prevents regs[j] from being 1.
            regs[j] ^= regs[j + reps];
        }
    }
    epilogue();
};

template <typename T, typename = void>
struct full_square {};

template <typename T>
using full_square_ = full_square<T>;

template <typename T>
struct full_square<T, std::enable_if_t<helper<T>::is_midint>> {
    // midint version
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            Midint<2 * num_bits> res = regs[j].template square();
            regs[j] = res.template high_part<num_bits>();
            regs[j] ^= res.template low_part<num_bits>();
            regs[j] ^= regs[j + reps];  // avoid zeroing out.
        }
    }
    epilogue();
};

template <typename T>
struct full_square<T, std::enable_if_t<helper<T>::is_boost_int>> {
    // boost_uint version
    preamble();
    for (size_t i = 0; i < num_iters; ++i) {
        for (size_t j = 0; j < reps; ++j) {
            boost_uint<2 * num_bits> res = boost_uint<2 * num_bits>(regs[j]) * regs[j];
            res ^= res >> num_bits;
            regs[j] = boost_uint<num_bits>(res);
            regs[j] ^= regs[j + reps];  // avoid zeroing out.
        }
    }
    epilogue();
};

template <template <class> typename tester, template <size_t> typename type>
void run_benchmark(size_t cycles_per_test, std::array<double, 8> rates, std::string operation_name, std::string type_name) {
    println();
    BENCH(tester<type<64>>(), cycles_per_test * std::exp2(rates[0]), operation_name + "~" + type_name + "~64");
    BENCH(tester<type<128>>(), cycles_per_test * std::exp2(rates[1]), operation_name + "~" + type_name + "~128");
    BENCH(tester<type<192>>(), cycles_per_test * std::exp2(rates[2]), operation_name + "~" + type_name + "~192");
    BENCH(tester<type<256>>(), cycles_per_test * std::exp2(rates[3]), operation_name + "~" + type_name + "~256");
    BENCH(tester<type<384>>(), cycles_per_test * std::exp2(rates[4]), operation_name + "~" + type_name + "~384");
    BENCH(tester<type<512>>(), cycles_per_test * std::exp2(rates[5]), operation_name + "~" + type_name + "~512");
    BENCH(tester<type<768>>(), cycles_per_test * std::exp2(rates[6]), operation_name + "~" + type_name + "~768");
    BENCH(tester<type<1024>>(), cycles_per_test * std::exp2(rates[7]), operation_name + "~" + type_name + "~1024");
}

int main(int argc, char** argv) {
    if (argc != 2) {
        println("Error! enter approximate time in seconds to run each benchmark! Number of arguments should be 1, but is", argc-1);
        return 1;
    }
    double time_each_test = atof(argv[1]);
    if (time_each_test < 0.01 or time_each_test > 100) {
        println("Error! strange arguments", argv[1]);
        return 2;
    }
    double cycles_per_sec = get_num_cycles_per_second(std::exp2(30));
    size_t num_cycles_per_test = time_each_test * cycles_per_sec;

    std::cout.precision(3);
    size_t num_tests = 180;
    println("Expected to run for", time_each_test * num_tests, "seconds    ||    cycles_per_sec=2^", std::log2(cycles_per_sec));
    println();

    // Casting comparisons
    std::string unsigned_minus_2 = "6277101735386680763835789423207666416102355444464034512894";
    std::string two_power_192 = "6277101735386680763835789423207666416102355444464034512896";
    std::string unsigned_minus_2_128 = "340282366920938463463374607431768211454";
    COMPARE_EQ(std::string(midint<3>("2")), std::string("2"));
    COMPARE_EQ(std::string(midint<3>("-2")), std::string(unsigned_minus_2));
    COMPARE_EQ(std::string(midint<3>(unsigned_minus_2)), std::string(unsigned_minus_2));
    COMPARE_EQ(std::string(midint<3>(boost_int<192>("-2"))), std::string(unsigned_minus_2));
    COMPARE_EQ(boost_uint<192>(midint<3>("-2")), boost_uint<192>("-2"));
    COMPARE_EQ(boost_uint<192>(midint<3>("2")), boost_uint<192>("2"));
    COMPARE_EQ(boost_int<192>(midint<3>("2")), boost_int<192>("2"));
    COMPARE_EQ(boost_int<192>(midint<3>("-2")), boost_int<192>("-2"));
    COMPARE_EQ(midint<3>("2").to_signed_string(), std::string("2"));
    COMPARE_EQ(midint<3>("-2").to_signed_string(), std::string("-2"));
    COMPARE_EQ(midint<3>(unsigned_minus_2).to_signed_string(), std::string("-2"));
    COMPARE_EQ(midint<3>("-2").to_unsigned_string(), std::string(unsigned_minus_2));
    COMPARE_EQ(midint<3>(two_power_192), midint<3>(0));

    COMPARE_EQ(std::string(midint<2>(boost_int<128>("-2"))), std::string(unsigned_minus_2_128));
    COMPARE_EQ(boost_uint<128>(midint<2>("-2")), boost_uint<128>("-2"));
    COMPARE_EQ(boost_uint<128>(midint<2>("2")), boost_uint<128>("2"));
    COMPARE_EQ(boost_int<128>(midint<2>("2")), boost_int<128>("2"));
    COMPARE_EQ(boost_int<128>(midint<2>("-2")), boost_int<128>("-2"));

    // A random selection of tests. More tests for mults due to complexity.
    COMPARE(add<Midint<64>>(), add<boost_uint<64>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(add<Midint<512>>(), add<boost_uint<512>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_square<Midint<128>>(), half_square<boost_uint<128>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_square<Midint<192>>(), half_square<boost_uint<192>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_square<Midint<256>>(), half_square<boost_uint<256>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_square<Midint<512>>(), half_square<boost_uint<512>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_square<Midint<1024>>(), half_square<boost_uint<1024>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_mult<Midint<64>>(), half_mult<boost_uint<64>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_mult<Midint<128>>(), half_mult<boost_uint<128>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_mult<Midint<192>>(), half_mult<boost_uint<192>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(half_mult<Midint<512>>(), half_mult<boost_uint<512>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_mult<Midint<64>>(), full_mult<boost_uint<64>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_mult<Midint<128>>(), full_mult<boost_uint<128>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_mult<Midint<192>>(), full_mult<boost_uint<192>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_mult<Midint<256>>(), full_mult<boost_uint<256>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_mult<Midint<512>>(), full_mult<boost_uint<512>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_square<Midint<64>>(), full_square<boost_uint<64>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_square<Midint<128>>(), full_square<boost_uint<128>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_square<Midint<192>>(), full_square<boost_uint<192>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_square<Midint<256>>(), full_square<boost_uint<256>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_square<Midint<512>>(), full_square<boost_uint<512>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(full_square<Midint<1024>>(), full_square<boost_uint<1024>>(), num_cycles_per_test * std::exp2(-10));
    COMPARE(xorit<Midint<512>>(), xorit<boost_uint<512>>(), num_cycles_per_test * std::exp2(-10));

    run_benchmark<add, Midint_>(num_cycles_per_test, {4.2, -0.2, -0.8, -1.8, -1.9, -2.6, -4.3, -4.7}, "add", "Midint");
    run_benchmark<add, boost_uint>(num_cycles_per_test, {4.2, 0.4, -4.5, -4.9, -5.3, -5.3, -5.6, -6.0}, "add", "boost_uint");
    run_benchmark<add, ttmath_uint_>(num_cycles_per_test, {-1.8, -2.4, -2.8, -3.4, -3.8, -4.4, -4.8, -5.2}, "add", "ttmath_uint");

    // TODO: partial half_mult

    run_benchmark<half_mult, Midint_>(num_cycles_per_test, {0.3, -2.0, -3.0, -3.7, -4.8, -5.6, -7.1, -7.8}, "half_mult", "Midint");
    run_benchmark<half_mult, boost_uint>(num_cycles_per_test, {0.3, -2.0, -6.0, -6.3, -7.0, -7.6, -8.6, -9.2}, "half_mult", "boost_uint");
    run_benchmark<half_mult, ttmath_uint_>(num_cycles_per_test, {-2.6, -4.8, -7.2, -7.7, -9.0, -9.4, -10.7, -11.1}, "half_mult", "ttmath_uint");

    // TODO: partial full_mult

    run_benchmark<full_mult_, Midint_>(num_cycles_per_test, {-0.1, -2.4, -3.5, -4.4, -5.5, -6.9, -8.0, -9.0}, "full_mult", "Midint");
    run_benchmark<full_mult_, boost_uint>(num_cycles_per_test, {-0.7, -6.3, -7.3, -7.3, -8.0, -8.5, -9.5, -10.2}, "full_mult", "boost_uint");

    //

    run_benchmark<xorit, Midint_>(num_cycles_per_test, {6.6, 4.6, -0.2, 3.6, 2.6, 2.6, -1.3, -2.3}, "xor", "Midint");
    run_benchmark<xorit, boost_uint>(num_cycles_per_test, {6.6, 4.6, -3.4, -4.0, -3.9, -4.1, -4.6, -3.6}, "xor", "boost_uint");
    run_benchmark<xorit, ttmath_uint_>(num_cycles_per_test, {6.6, 0.6, 0.1, -1.4, -2.3, -2.8, -1.5, -2.4}, "xor", "ttmath_uint");

    //

    run_benchmark<less_than, Midint_>(num_cycles_per_test, {0.3, -1.3, -2.3, -1.2, -1.3, -1.5, -2.5, -2.1}, "less_than", "Midint");
    run_benchmark<less_than, boost_uint>(num_cycles_per_test, {0.0, -1.6, -1.5, -2.9, -3.1, -3.3, -1.5, -1.5}, "less_than", "boost_uint");
    run_benchmark<less_than, ttmath_uint_>(num_cycles_per_test, {0.2, -1.2, -1.7, -2.0, -1.8, -0.6, -2.5, -3.0}, "less_than", "ttmath_uint");

    //

    run_benchmark<equals_to, Midint_>(num_cycles_per_test, {-1.5, -1.4, -2.0, -1.9, -2.5, -2.8, -3.4, -4.0}, "equals_to", "Midint");
    run_benchmark<equals_to, boost_uint>(num_cycles_per_test, {-1.5, -1.5, -3.2, -3.8, -3.9, -4.8, -5.3, -5.7}, "equals_to", "boost_uint");
    run_benchmark<equals_to, ttmath_uint_>(num_cycles_per_test, {0.2, -1.8, -2.5, -3.1, -3.3, -3.8, -4.7, -4.8}, "equals_to", "ttmath_uint");

    //

    run_benchmark<half_square_, Midint_>(num_cycles_per_test, {-0.4, -1.4, -2.4, -3.3, -4.3, -5.1, -6.2, -7.2}, "half_square", "Midint");
    run_benchmark<half_square_, boost_uint>(num_cycles_per_test, {-0.4, -1.4, -6.5, -6.7, -7.2, -7.7, -8.6, -9.3}, "half_square", "boost_uint");
    run_benchmark<half_square_, ttmath_uint_>(num_cycles_per_test, {-2.6, -4.5, -7.2, -7.8, -9.0, -9.4, -10.6, -11.0}, "half_square", "ttmath_uint");

    //

    run_benchmark<full_square_, Midint_>(num_cycles_per_test, {-1.0, -2.9, -3.6, -4.5, -5.6, -6.7, -7.6, -8.4}, "full_square", "Midint");
    run_benchmark<full_square_, boost_uint>(num_cycles_per_test, {-1.0, -6.3, -7.5, -7.5, -8.1, -8.7, -9.5, -10.2}, "full_square", "boost_uint");

    return 0;
}
