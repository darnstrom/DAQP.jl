module TestMOIDAQP

import DAQP 
using MathOptInterface
using Test

const MOI = MathOptInterface

const OPTIMIZER = MOI.instantiate(
    MOI.OptimizerWithAttributes(DAQP.Optimizer),
)

const BRIDGED = MOI.instantiate(
    MOI.OptimizerWithAttributes(DAQP.Optimizer),
    with_bridge_type = Float64,
)

# See the docstring of MOI.Test.Config for other arguments.
const CONFIG = MOI.Test.Config(
    # Modify tolerances as necessary.
    atol = 1e-6,
    rtol = 1e-6,
    # Use MOI.LOCALLY_SOLVED for local solvers.
    optimal_status = MOI.OPTIMAL,
    # Pass attributes or MOI functions to `exclude` to skip tests that
    # rely on this functionality.
    exclude = Any[MOI.VariableName, MOI.delete],
)

"""
    runtests()

This function runs all functions in the this Module starting with `test_`.
"""
function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

"""
    test_runtests()

This function runs all the tests in MathOptInterface.Test.

Pass arguments to `exclude` to skip tests for functionality that is not
implemented or that your solver doesn't support.
"""
function test_runtests()
    MOI.Test.runtests(
        BRIDGED,
        CONFIG,
        exclude = [ # Bin only supported for strictly convex objective in DAQP
            "test_constraint_ZeroOne_",
            "test_variable_solve_ZeroOne_",
            "test_cpsat_",
            "test_linear_Indicator_ON_ONE",
            "test_solve_ObjectiveBound_MAX_SENSE_IP",
            "test_solve_ObjectiveBound_MIN_SENSE_IP",
            "test_solve_SOS2_",
            "test_variable_solve_Integer_with_",
        ],
        exclude_tests_after = v"1.22.0",
    )
    return
end

"""
    test_SolverName()

You can also write new tests for solver-specific functionality. Write each new
test as a function with a name beginning with `test_`.
"""
function test_SolverName()
    @test MOI.get(DAQP.Optimizer(), MOI.SolverName()) == "DAQP"
    return
end

"""
    test_row_index_integer_width()

Row indices must not be dispatched on `Int`. On 32-bit platforms `Int === Int32`,
while MOI index values and attributes such as `MOI.NumberOfVariables` are `Int64`
on every platform, so mixing the two (e.g. `row - n`) promotes the row index to
`Int64` and no `::Int` method matches. See DAQP.jl issue #25.
"""
function test_row_index_integer_width()
    @test fieldtype(DAQP.Optimizer, :rows) == Dict{Int64,Int64}
    idxmap = typeof(MOI.Utilities.IndexMap())
    affine = MOI.ScalarAffineFunction{Cdouble}
    interval = MOI.Interval{Cdouble}
    for T in (Int32, Int64)
        @test hasmethod(DAQP.extract_offset, Tuple{Vector{Cdouble},T,affine})
        @test hasmethod(DAQP.extract_offset, Tuple{Vector{Cdouble},T,MOI.VariableIndex})
        @test hasmethod(DAQP.extract_A, Tuple{Matrix{Cdouble},affine,T,idxmap})
        @test hasmethod(DAQP.extract_A, Tuple{Matrix{Cdouble},MOI.VariableIndex,T,idxmap})
        @test hasmethod(DAQP.extract_b,
                        Tuple{Vector{Cdouble},Vector{Cdouble},Vector{Cint},T,affine,interval})
        @test hasmethod(DAQP.extract_b,
                        Tuple{Vector{Cdouble},Vector{Cdouble},Vector{Cint},T,MOI.VariableIndex,interval})
        @test hasmethod(DAQP.extract_b,
                        Tuple{Vector{Cdouble},Vector{Cdouble},Vector{Cint},T,interval})
    end
    return
end

end # module TestMOIDAQP 

# This line at tne end of the file runs all the tests!
TestMOIDAQP.runtests()
