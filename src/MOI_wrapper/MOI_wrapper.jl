using MathOptInterface
export Optimizer

## Const definitions

const MOI = MathOptInterface
const MOIU = MOI.Utilities

const Interval = MOI.Interval{Float64}

const Affine = MOI.ScalarAffineFunction{Cdouble}
const Quadratic = MOI.ScalarQuadraticFunction{Cdouble}
const VectorAffine = MOI.VectorAffineFunction{Cdouble}

const SupportedVectorSets = Union{MOI.Zeros,MOI.Nonnegatives,MOI.Nonpositives}

const Interval = MOI.Interval{Cdouble}
const LessThan = MOI.LessThan{Cdouble}
const GreaterThan = MOI.GreaterThan{Cdouble}
const EqualTo = MOI.EqualTo{Cdouble}
const SupportedSets= Union{GreaterThan, LessThan, EqualTo, Interval, MOI.ZeroOne}
# mappings between MOI and internal definitions

# TODO
const DAQPtoMOITerminationStatus = Dict([
    DAQP.SOFT_OPTIMAL    =>  MOI.OPTIMAL,
    DAQP.OPTIMAL         =>  MOI.OPTIMAL,
    DAQP.INFEASIBLE      =>  MOI.INFEASIBLE,
    DAQP.CYCLE           =>  MOI.SLOW_PROGRESS,
    DAQP.UNBOUNDED       =>  MOI.DUAL_INFEASIBLE,
    DAQP.ITERLIMIT       =>  MOI.ITERATION_LIMIT,
    DAQP.NONCONVEX       =>  MOI.INVALID_MODEL,
    DAQP.OVERDETERMINED  =>  MOI.INVALID_OPTION,
])

const DAQPtoMOIPrimalStatus = Dict([
    DAQP.SOFT_OPTIMAL    =>  MOI.NEARLY_FEASIBLE_POINT,
    DAQP.OPTIMAL         =>  MOI.FEASIBLE_POINT,
    DAQP.INFEASIBLE      =>  MOI.INFEASIBLE_POINT,
    DAQP.CYCLE           =>  MOI.UNKNOWN_RESULT_STATUS,
    DAQP.UNBOUNDED       =>  MOI.UNKNOWN_RESULT_STATUS,
    DAQP.ITERLIMIT       =>  MOI.UNKNOWN_RESULT_STATUS,
    DAQP.NONCONVEX       =>  MOI.UNKNOWN_RESULT_STATUS,
    DAQP.OVERDETERMINED  =>  MOI.UNKNOWN_RESULT_STATUS,
])

const DAQPtoMOIDualStatus = Dict([
    DAQP.SOFT_OPTIMAL    =>  MOI.FEASIBLE_POINT,
    DAQP.OPTIMAL         =>  MOI.FEASIBLE_POINT,
    DAQP.INFEASIBLE      =>  MOI.INFEASIBLE_POINT,
    DAQP.CYCLE           =>  MOI.UNKNOWN_RESULT_STATUS,
    DAQP.UNBOUNDED       =>  MOI.UNKNOWN_RESULT_STATUS,
    DAQP.ITERLIMIT       =>  MOI.UNKNOWN_RESULT_STATUS,
    DAQP.NONCONVEX       =>  MOI.UNKNOWN_RESULT_STATUS,
    DAQP.OVERDETERMINED  =>  MOI.UNKNOWN_RESULT_STATUS,
])


## Main interface struct

mutable struct Optimizer <: MOI.AbstractOptimizer
    model::DAQP.Model
    has_results::Bool
    is_empty::Bool
    sense::MOI.OptimizationSense
    objconstant::Cdouble
    rowranges::Dict{Int, UnitRange{Int}}
    setup_time::Cdouble
    silent::Bool
    info #TODO: specify type

    function Optimizer(; user_settings...)
        model = DAQP.Model()
        has_results = false
        is_empty = true
        sense = MOI.MIN_SENSE
        objconstant = 0.0 
        rowranges = Dict{Int, UnitRange{Int}}()
        setup_time = 0.0
        silent=true
        optimizer = new(model,has_results,is_empty,sense,
                        objconstant,rowranges,setup_time,silent,nothing)
        for (key, value) in user_settings
            MOI.set(optimizer, MOI.RawOptimizerAttribute(string(key)), value)
        end
        return optimizer
    end
end

## Required basic methods

# reset the optimizer
function MOI.empty!(optimizer::Optimizer)
    #just make a new model, keeping current settings
    tmp_settings = settings(optimizer.model)
    optimizer.model = DAQP.Model()
    isnothing(tmp_settings) || settings(optimizer.model,tmp_settings)

    optimizer.has_results = false
    optimizer.is_empty = true
    optimizer.sense = MOI.MIN_SENSE # model parameter, so needs to be reset
    optimizer.objconstant = 0.0 
    optimizer.rowranges = Dict{Int, UnitRange{Int}}()
end

MOI.is_empty(optimizer::Optimizer) = optimizer.is_empty

function MOI.optimize!(optimizer::Optimizer)
    if(!optimizer.is_empty)
        ~,~,~,optimizer.info=DAQP.solve(optimizer.model)
        optimizer.has_results = true
    end
    return
end

function Base.show(io::IO, optimizer::Optimizer)

    myname = MOI.get(optimizer, MOI.SolverName())
    if optimizer.is_empty
        print(io,"Empty $(myname) - Optimizer")

    else
        println(io, "$(myname) - Optimizer")
        println(io, " : Has results: $(optimizer.has_results)")
        println(io, " : Objective constant: $(optimizer.objconstant)")
        println(io, " : Sense: $(optimizer.sense)")

        if optimizer.has_results
            println(io, " : Problem status: $(MOI.get(optimizer,MOI.RawStatusString()))")
            value = round(MOI.get(optimizer,MOI.ObjectiveValue()),digits=3)
            println(io, " : Optimal objective: $(value)")
            println(io, " : Iterations: $(MOI.get(optimizer,MOI.SimplexIterations()))")
            solvetime = round.(optimizer.model.info.solve_time*1000,digits=2)
            println(io, " : Solve time: $(solvetime)ms")
        end
    end
end


## Solver Attributes, get/set

MOI.get(opt::Optimizer, ::MOI.SolverName)        = "DAQP" 
MOI.get(opt::Optimizer, ::MOI.SolverVersion)     = "0.3.0" 
MOI.get(opt::Optimizer, ::MOI.RawSolver)         = opt.model
MOI.get(opt::Optimizer, ::MOI.ResultCount)       = opt.has_results ? 1 : 0
MOI.get(opt::Optimizer, ::MOI.NumberOfVariables) = opt.model.n
MOI.get(opt::Optimizer, ::MOI.SolveTimeSec)      = opt.info.solve_time+opt.setup_time
MOI.get(opt::Optimizer, ::MOI.RawStatusString)   = string(opt.info.status)
MOI.get(opt::Optimizer, ::MOI.SimplexIterations) = Int64(opt.info.iterations)

function MOI.get(opt::Optimizer, a::O) where {O<:Union{ MOI.DualObjectiveValue, 
                                                        MOI.ObjectiveValue}}
    MOI.check_result_index_bounds(opt, a)
    rawobj = opt.info.fval + opt.objconstant
    return opt.sense == MOI.MIN_SENSE ? rawobj : -rawobj
end


MOI.supports(::Optimizer, ::MOI.TerminationStatus) = true
function MOI.get(opt::Optimizer, ::MOI.TerminationStatus)
    opt.has_results || return MOI.OPTIMIZE_NOT_CALLED
    return DAQPtoMOITerminationStatus[opt.info.exitflag]
end

MOI.supports(::Optimizer, ::MOI.PrimalStatus) = true
function MOI.get(opt::Optimizer, attr::MOI.PrimalStatus)
    if !opt.has_results || attr.result_index != 1
        return MOI.NO_SOLUTION
    else
        return DAQPtoMOIPrimalStatus[opt.info.exitflag]
    end
end

MOI.supports(::Optimizer, a::MOI.DualStatus) = true
function MOI.get(opt::Optimizer, attr::MOI.DualStatus)
    if !opt.has_results || attr.result_index != 1
        return MOI.NO_SOLUTION
    end
    return DAQPtoMOIDualStatus[opt.info.exitflag]
end


MOI.supports(::Optimizer, ::MOI.VariablePrimal) = true
function MOI.get(opt::Optimizer, a::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(opt, a)
    return opt.info.x[vi.value]
end

MOI.supports(::Optimizer, ::MOI.ConstraintDual) = true
function MOI.get(
    opt::Optimizer,
      a::MOI.ConstraintDual,
     ci::MOI.ConstraintIndex{F, S}
) where {F, S <: MOI.AbstractSet}

    MOI.check_result_index_bounds(opt, a)
    rows = constraint_rows(opt.rowranges, ci)

    if(opt.sense != MOI.FEASIBILITY_SENSE)
        λ=abs.(opt.info.λ[rows]) # λ for lower bounds ≥ 0 in MOI
    else
        λ = (length(rows)>1) ? zeros(Cdouble,length(rows)) : 0.0
    end
    return λ
end

MOI.supports(::Optimizer, ::MOI.ConstraintPrimal) = true
function MOI.get(
    opt::Optimizer,
      a::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{F, S}
) where {F, S <: MOI.AbstractSet}

     MOI.check_result_index_bounds(opt, a)
     rows = constraint_rows(opt.rowranges, ci)
     n = opt.model.qpj.n
     Ax = (last(rows)<=n) ? opt.info.x[rows] : opt.model.qpj.A[:,rows.-n]'*opt.info.x
     return min.(opt.model.qpj.bupper[rows]-Ax,Ax-opt.model.qpj.blower[rows])
end

#Currently there is no internal printing in DAQP
MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.get(opt::Optimizer, ::MOI.Silent) = opt.silent
MOI.set(opt::Optimizer, ::MOI.Silent, v::Bool) = (opt.silent=v)

MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true
MOI.get(opt::Optimizer, param::MOI.RawOptimizerAttribute) =
    getproperty(settings(opt.model), Symbol(param.name))
MOI.set(opt::Optimizer, param::MOI.RawOptimizerAttribute, value) =
    settings(opt.model, Dict(Symbol(param.name)=>value))


# not currently supported
MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = false
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = false

## Supported constraint types

MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:MOI.VectorAffineFunction},
    ::Type{<:SupportedVectorSets}
) = true
MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:MOI.ScalarAffineFunction},
    ::Type{<:SupportedSets}
) = true

MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:MOI.VariableIndex},
    ::Type{<:SupportedSets}
) = true

## Supported objective functions

MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{<:Union{
       #MOI.ScalarAffineFunction,
       MOI.ScalarQuadraticFunction,
    }}
) = true


## copy_to interface

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    idxmap = MOIU.IndexMap(dest, src)

    copy_to_check_attributes(dest,src)

    #assemble the constraints data
    assign_constraint_row_ranges!(dest.rowranges, idxmap, src)
    A, bupper, blower, sense = process_constraints(dest, src, idxmap)

    #assemble the objective data
    dest.sense = MOI.get(src, MOI.ObjectiveSense())
    H, f, dest.objconstant = process_objective(dest, src, idxmap)

    println(H)
    println(f)
    display(A)
    println(bupper)
    println(blower)
    println(sense)
    # setup solver
    exitflag, dest.setup_time = DAQP.setup(dest.model,H,f,A,bupper, blower, sense;A_rowmaj=true)
    if(exitflag >= 0)
        dest.is_empty = false
    else
        error("DAQP currently only supports strictly convex objectives")
    end

    return idxmap
end

function copy_to_check_attributes(dest, src)

    #allowable model attributes
    for attr in MOI.get(src, MOI.ListOfModelAttributesSet())
        if attr == MOI.Name()           ||
           attr == MOI.ObjectiveSense() ||
           attr isa MOI.ObjectiveFunction
            continue
        end
        throw(MOI.UnsupportedAttribute(attr))
    end

    #allowable variable attributes
    for attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        if attr == MOI.VariableName()
            continue
        end
        throw(MOI.UnsupportedAttribute(attr))
    end

    #allowable constraint types and attributes
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        if !MOI.supports_constraint(dest, F, S)
            throw(MOI.UnsupportedConstraint{F, S}())
        end
        for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F, S}())
            if attr == MOI.ConstraintName()
                continue
            end
            throw(MOI.UnsupportedAttribute(attr))
        end
    end

    return nothing
end

#Set up index map from `src` variables and constraints to `dest` variables and constraints.

function MOIU.IndexMap(dest::Optimizer, src::MOI.ModelLike)

    idxmap = MOIU.IndexMap()

    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = MOI.VariableIndex(i)
    end
    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        MOI.supports_constraint(dest, F, S) || throw(MOI.UnsupportedConstraint{F, S}())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci in cis_src
            i += 1
            idxmap[ci] = MOI.ConstraintIndex{F, S}(i)
        end
    end

    return idxmap
end


function assign_constraint_row_ranges!(
    rowranges::Dict{Int, UnitRange{Int}},
    idxmap::MOIU.IndexMap,
    src::MOI.ModelLike
)

    startrow = 1+MOI.get(src, MOI.NumberOfVariables())
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
        for ci_src in cis_src
            set = MOI.get(src, MOI.ConstraintSet(), ci_src)
            ci_dest = idxmap[ci_src]
            if(F!=MOI.VariableIndex)
                endrow = startrow + MOI.dimension(set) - 1
                rowranges[ci_dest.value] = startrow : endrow
                startrow = endrow + 1
            else
                var_id = MOI.get(src, MOI.ConstraintFunction(), ci_src).value
                rowranges[ci_dest.value] = var_id:var_id
            end
        end
    end

    return nothing
end

function constraint_rows(
    rowranges::Dict{Int, UnitRange{Int}},
    ci::MOI.ConstraintIndex{<:Any, <:MOI.AbstractScalarSet}
)
    rowrange = rowranges[ci.value]
    length(rowrange) == 1 || error()
    return first(rowrange)
end

constraint_rows(
    rowranges::Dict{Int, UnitRange{Int}},
    ci::MOI.ConstraintIndex{<:Any, <:MOI.AbstractVectorSet}
) = rowranges[ci.value]

constraint_rows(
    optimizer::Optimizer,
    ci::MOI.ConstraintIndex
) = constraint_rows(optimizer.rowranges, ci)


## Extract Constraints
#constraints: 
# blower[1:n]     ≤  x ≤  bupper[1:n] 
# blower[n+1:end] ≤ Ax ≤  bupper[n+1:end]
function process_constraints(
    dest::Optimizer,
    src::MOI.ModelLike,
    idxmap
)

    rowranges = dest.rowranges
    n = MOI.get(src, MOI.NumberOfVariables())
    m = max(n,mapreduce(last, max, values(rowranges), init=0))

    A = zeros(Cdouble,n,m-n)
    bupper = Vector{Cdouble}(undef, m)
    blower = Vector{Cdouble}(undef, m)
    offset = Vector{Cdouble}(undef, m)
    sense  = Vector{Cint}(undef,m)

    bupper[1:n].=1e30;
    blower[1:n].=-1e30;
    offset[1:n].=0;
    sense[1:n].= IMMUTABLE


    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        process_constraints!(A, bupper, blower, sense, offset,
                            src, idxmap, rowranges, 
                            F, S)
    end
    bupper .-= offset 
    blower .-= offset 

    return (A, bupper, blower, sense)

end

function process_constraints!(
    A::Matrix{Cdouble},
    bupper::Vector{Cdouble},
    blower::Vector{Cdouble},
    sense::Vector{Cint},
    offset::Vector{Cdouble},
    src::MOI.ModelLike,
    idxmap,
    rowranges::Dict{Int,UnitRange{Int}},
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    n = MOI.get(src, MOI.NumberOfVariables())
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
    for ci in cis_src
        s = MOI.get(src, MOI.ConstraintSet(), ci)
        f = MOI.get(src, MOI.ConstraintFunction(), ci)
        rows = constraint_rows(rowranges, idxmap[ci])
        extract_offset(offset, rows, f)
        extract_A(A, f, rows.-n, idxmap)
        extract_b(bupper, blower, sense, rows, f, s)
    end
    return
end

extract_offset(::Vector{Cdouble}, ::Int, ::MOI.VariableIndex) = nothing
extract_A(::Matrix{Cdouble},::MOI.VariableIndex,::Int, Any) = nothing

function extract_offset(offset::Vector{Cdouble}, row::Int, f::Affine)
    offset[row] = MOI.constant(f, Cdouble)
    return
end

function extract_offset(
    offset::Vector{Float64},
    rows::UnitRange{Int},
    f::VectorAffine,
)
    for (i, row) in enumerate(rows)
        offset[row] = f.constants[i]
    end
end


function extract_A(
    A::Matrix{Cdouble},
    f::Affine,
    row::Int,
    idxmap,
)
    for term in f.terms
        var = term.variable
        col = idxmap[var].value
        A[col,row] = term.coefficient # colmaj -> rowmaj
    end
end

function extract_A(
    A::Matrix{Cdouble},
    f::VectorAffine,
    rows::UnitRange{Int},
    idxmap,
)
    for term in f.terms
        row = rows[term.output_index]
        var = term.scalar_term.variable
        col = idxmap[var].value
        A[col,row] = term.scalar_term.coefficient # colmaj -> rowmaj
    end
end

function extract_b(
    bupper::Vector{Cdouble},
    blower::Vector{Cdouble},
    sense::Vector{Cint},
    row::Int,
    f::Affine,
    s::SupportedSets,
)
    extract_b(bupper,blower,sense, row, MOI.Interval(s))
    return
end
function extract_b(
    bu::Vector{Cdouble},
    bl::Vector{Cdouble},
    sense::Vector{Cint},
    row::Int,
    f::MOI.VariableIndex,
    s::SupportedSets,
)

    println(s)
    i = (s==MOI.ZeroOne()) ? MOI.Interval(0,1) : MOI.Interval(s)
    extract_b(bu,bl,sense, row, MOI.Interval(max(i.lower,bl[row]),min(i.upper,bu[row])))
    if(s==MOI.ZeroOne()) sense[row] = BINARY end # Mark binary constraints
end

function extract_b(
    bupper::Vector{Cdouble},
    blower::Vector{Cdouble},
    sense::Vector{Cint},
    row::Int,
    interval::Interval,
)
    bupper[row] = interval.upper #TODO: add min/max for simple...
    blower[row] = interval.lower
    sense[row] = (interval.lower == interval.upper) ?  DAQP.EQUALITY : 0
    return
end

lower(::MOI.Zeros, ::Int) = 0.0
lower(::MOI.Nonnegatives, ::Int) = 0.0
lower(::MOI.Nonpositives, ::Int) = -Inf
upper(::MOI.Zeros, ::Int) = 0.0
upper(::MOI.Nonnegatives, ::Int) = Inf
upper(::MOI.Nonpositives, ::Int) = 0.0

function extract_b(
    bu::Vector{Cdouble},
    bl::Vector{Cdouble},
    sense::Vector{Cint},
    rows::UnitRange{Int},
    f,
    s::S,
) where {S<:SupportedVectorSets}
    for (i, row) in enumerate(rows)
        bu[row] = upper(s, i)
        bl[row] = lower(s, i)
        sense[row] = (bu[row]== bl[row]) ? DAQP.ACTIVE+DAQP.IMMUTABLE : 0
    end
end


## Extract Objective
# Construct the objective minimize `0.5 x' H x + f' x + c

function process_objective(dest::Optimizer, src::MOI.ModelLike, idxmap)
    sense = dest.sense
    n = MOI.get(src, MOI.NumberOfVariables())

    if sense == MOI.FEASIBILITY_SENSE
        H = Matrix{Cdouble}(I(n))# TODO: use nothing instead 
        f = zeros(n); 
        c = 0.0
    else
        function_type = MOI.get(src, MOI.ObjectiveFunctionType())
        f = zeros(Cdouble,n)

        if function_type == Affine  #TODO: DAQP supports pure LPs, but still experimental
            faffine = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction}())
            H = nothing 
            process_objective_linearterm!(f, faffine.terms, idxmap)
            c = faffine.constant

        elseif function_type == Quadratic 
            H = zeros(Cdouble,n,n)
            fquadratic = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction}())
            for term in fquadratic.quadratic_terms
                i,j = Int(idxmap[term.variable_1].value),Int(idxmap[term.variable_2].value)
                println("element $i,$j  is $(term.coefficient)")
                H[i,j] += term.coefficient
                if(i!=j)
                    H[j,i] += term.coefficient
                end
            end
            process_objective_linearterm!(f, fquadratic.affine_terms, idxmap)
            c = fquadratic.constant

        else
            throw(MOI.UnsupportedAttribute(MOI.ObjectiveFunction{function_type}()))
        end

        if sense == MOI.MAX_SENSE
            H = -H
            f = -f
            c = -c
        end

    end
    return (H, f, c)
end


function process_objective_linearterm!(
    f::Vector{Cdouble},
    terms::Vector{<:MOI.ScalarAffineTerm},
    idxmapfun::Function = identity
) 
    f .= 0
    for term in terms
        var = term.variable
        coeff = term.coefficient
        f[idxmapfun(var).value] += coeff
    end
    return nothing
end

function process_objective_linearterm!(
    f::Vector{Cdouble},
    terms::Vector{<:MOI.ScalarAffineTerm},
    idxmap::MOIU.IndexMap
)
    process_objective_linearterm!(f, terms, var -> idxmap[var])
end
