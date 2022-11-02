# Midint
A medium-sized integer arithmetic library. Lightweight Design for efficiency and flexibility.

## Purpose
This tiny library defines a `Midint<n>` structure which represents an (unsigned) integer with `n` bits (supports only `n % 64 == 0`). It allows simple arithmetic and logic of these numbers, with the exception of divison which is currently not implemented.

The library is designed for simplicy -- you may seamlessly convert between boost's fixed size `cpp_int`. Also, the following plots demonstrate the performance gain relative to `cpp_int`.

## Code samples
TODO

## How to compile
TODO

Feel free to leave comments about this project :)

## Performance comparison

<picture>
<img src="https://github.com/IamPoosha/Midint/blob/main/figures_clang/add.png" alt="Addition efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/IamPoosha/Midint/blob/main/figures_clang/half_mult.png" alt="Multiplication efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/IamPoosha/Midint/blob/main/figures_clang/full_mult.png" alt="Full multiplication efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/IamPoosha/Midint/blob/main/figures_clang/half_square.png" alt="Squaring efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/IamPoosha/Midint/blob/main/figures_clang/less_than.png" alt="Comparison (usually determined by top limb) efficiency comparison"/>
</picture>

<picture>
<img src="https://github.com/IamPoosha/Midint/blob/main/figures_clang/xor.png" alt="Xor efficiency comparison"/>
</picture>
