# Midint
A medium-sized integer arithmetic library. Lightweight design for efficiency and flexibility.

## Purpose
This tiny library defines a `Midint<n>` structure which represents an (unsigned) integer with `n` bits (supports only `n % 64 == 0`). It allows simple arithmetic and logic of these numbers, with the exception of divison which is currently not implemented.

The library is designed for simplicy -- you may seamlessly convert between boost's fixed size `cpp_int`. Also, the following plots demonstrate the performance gain relative to `cpp_int`.

## Code samples
Python's `0xFEDCBA9876543210 ** 8 // 2**(512-64)` may be written as
```
auto a = Midint<64>(0xFEDCBA9876543210);
Midint<128> b = a.square();
Midint<256> c = b.square();
Midint<512> d = c.square();
std::cout << d.high_part<64>();
```

## How to compile
Just `#include "miding.h"` wherever you use the library.

Feel free to leave comments about this project :)

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
