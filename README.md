# DAQP.jl

[![](https://github.com/darnstrom/DAQP.jl/workflows/CI/badge.svg)](https://github.com/darnstrom/DAQP.jl/actions)
[![](http://codecov.io/gh/darnstrom/DAQP.jl/graphs/badge.svg)](http://codecov.io/github/darnstrom/DAQP.jl)

[DAQP.jl](https://github.com/darnstrom/DAQP.jl) is a Julia wrapper for the
Quadratic Programming solver [DAQP](https://github.com/darnstrom/daqp).

## License

DAQP.jl is licensed under the [MIT license](https://github.com/darnstrom/DAQP.jl/blob/main/LICENSE).

The underlying solver, [darnstrom/daqp](https://github.com/darnstrom/daqp) is
licensed under the [MIT license](https://github.com/darnstrom/daqp/blob/master/LICENSE).

## Installation

Install DAQP.jl usinng the Julia package manager:

```julia
import Pkg
Pkg.add("DAQP")
```

## Use with JuMP

To use DAQP with JuMP, do:
```
using JuMP, DAQP
model = Model(DAQP.Optimizer)
```

## Documentation

General information about the solver is available at [https://darnstrom.github.io/daqp/](https://darnstrom.github.io/daqp/),
and specifics for the Julia interface are available at
[https://darnstrom.github.io/daqp/start/julia](https://darnstrom.github.io/daqp/start/julia). 
