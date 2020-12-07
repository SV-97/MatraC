# MatraC
Small C++ Matrix library

I created this because I needed some matrices in C++ and wanted to try a few things with templates / find out what's possible. One should be able to clean the code up quite a bit once C++20 is widely available.

* `util.h` contains a few useful utilities like a `PANIC` macro, a string formatting function and a bit of cursed template hackery
* `DynMax.h` contains dynamically sized, dense and sparse matrices (easily extendable to other data representations)
* `Matrix.h` contains statically sized, fully stack allocatable matrices
* `Range.h` contains what the name says. Ranges

One could probably deduplicate a bit of code between dynamic and static matrices and the template stuff definitely isn't nice to read as it is, but it's quite nice to work with.
