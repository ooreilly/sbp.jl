# sbp.jl

[![Build Status](https://travis-ci.org/ooreilly/sbp.jl.svg?branch=master)](https://travis-ci.org/ooreilly/sbp.jl)
[![Codecov](https://codecov.io/gh/ooreilly/sbp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ooreilly/sbp.jl)

This Julia package can be used for discretizing partial differential
equations using the summation-by-parts approach. All schemes are discretized by
sparse matrices and dense vectors. This approach is computationally inefficient
and should therefore only be used for educational purposes. 
In addition, this package is necessary for reproducing some of the published
numerical examples presented in the repository
https://github.com/ooreilly/sbp. 

## Installation
Open Julia, and press `]` to go to package manner. Type
```
pkg> add https://github.com/ooreilly/sbp.jl/.

```
Type `ctrl` + `d` to exit.

## Usage
Take a look at this [example](examples/wave_equation_1d.jl) to get an idea of
how to use this package.


## Tests
Tests for this package can be run by using Julia's package manner.
```
pkg> test sbp
```
The tests themselves can be found in the directory [test](test/) and are worth
looking at for discovering many of the features that this package has to offer.
