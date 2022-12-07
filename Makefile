default: main

main: main.cpp bench_utils.h midint.h
	clang++-14 main.cpp -O3 -march=native -o main -std=c++17

clean:
	rm main
