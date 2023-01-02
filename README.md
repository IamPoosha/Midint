# Midint
A medium-sized integer arithmetic library. Lightweight design for efficiency and flexibility.

## Purpose
This tiny library defines a `Midint<n>` structure which represents an (unsigned) integer with `n` bits (supports only `n % 64 == 0`). It allows simple arithmetic and logic of these numbers, with the exception of divison which is currently not implemented[^1].

The library is designed for simplicity -- you may seamlessly convert between `Midint` and boost's fixed size `cpp_int`. Also, the following plots demonstrate the performance gain relative to `cpp_int`.

[^1]: Since `Midint` is easily convertible to `cpp_int`, you may still divide and gain performance in case your application is not divisons-intensive.

## Code samples
Python's `0xFEDCBA9876543210 ** 8 // 2**(512-64)` may be written as
```
auto a = Midint<64>(0xFEDCBA9876543210);
auto b = a.square();
Midint<256> c = b.square();
Midint<512> d = c.square();
std::cout << d.high_part<64>();
```

## How to compile
Just `#include "midint.h"` wherever you use the library.

Feel free to leave comments about this project :)

### Alternatives
GMP also has a set of fixed-size arithmetic functions.

For documentation, see: https://gmplib.org/manual/Low_002dlevel-Functions

For implementation, see: https://github.com/alisw/GMP/tree/master/mpn/generic

Due to the complexity of the integration with that code, it was not included in the benchmark. 

## Performance comparison

<picture>
<img src="https://github.com/ohadkel/Midint/blob/main/figures_clang/add.png" alt="Addition efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/ohadkel/Midint/blob/main/figures_clang/half_mult.png" alt="Multiplication efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/ohadkel/Midint/blob/main/figures_clang/full_mult.png" alt="Full multiplication efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/ohadkel/Midint/blob/main/figures_clang/half_square.png" alt="Squaring efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/ohadkel/Midint/blob/main/figures_clang/less_than.png" alt="Comparison (usually determined by top limb) efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/ohadkel/Midint/blob/main/figures_clang/xor.png" alt="Xor efficiency comparison"/>
</picture>
