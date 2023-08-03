"""
# Example calls
    xstar, fval, exitflag, info = DAQP.quadprog(H,f,A,bupper)
    xstar, fval, exitflag, info = DAQP.quadprog(H,f,A,bupper,blower,sense) 

finds the solution `xstar` to the quadratic program

```
 minimize	   0.5 x' H x + f' x
subject to   blower <= A x <= bupper
```
If `bupper` and `blower` have more elements than rows of `A`, the first
elements are interpreted as simple bounds. For example:

```
    A = [7.0 8.0]
    blower = [-4.0; -5.0; -6.0]
    bupper = [ 1.0;  2.0;  3.0]
```
is interpreted as

```
        -4.0 <= x₁ <= 1.0
        -5.0 <= x₂ <= 2.0
    -6.0 <= 7 x₁ + 8 x₂ <= 3.0

```

# Input 
* `H`           - cost matrix
* `f`           - cost vector
* `A`           - linear constraint matrix
* `buppe`       - upper bounds for constraints
* `blower`      - lower bounds for constraints (default: -Inf)
* `sense`       - constraint types, as a vector of Cints (default: 0). Example types:
  * `0 ` : inequality
  * `1 ` : active inequality (used as warm start)
  * `5 ` : equality
  * `8 ` : soft (allowed to be violated if necessary)
  * `16` : binary (either upper or lower bound should hold with equality)

# Output
* `xstar`       - solution
* `fval`        - objective function value for `xstar`. 
* `exitflag`    - flag from solver (>0 success, <0 failure) 
* `info`        - tuple containing profiling information from the solver. 

"""
function quadprog(H::Matrix{Float64},f::Vector{Float64}, 
        A::Matrix{Float64},bupper::Vector{Float64},blower::Vector{Float64}=Float64[],sense::Vector{Cint}=Cint[];A_rowmaj=false,settings=nothing)
    return quadprog(QPj(H,f,A,bupper,blower,sense;A_rowmaj);settings)
end
function quadprog(qpj::QPj;settings=nothing)
    # TODO: check validity of dimensions
    # Setup QP
    qp = QPc(qpj);

    # Setup output struct
    xstar = zeros(Float64,qp.n); 
    lam= zeros(Float64,qp.m); 
    result= Ref(DAQPResult(xstar,lam));
    ptr_settings = settings isa DAQPSettings ? Ref(settings) : Ptr{DAQP.DAQPSettings}(C_NULL)


    ccall((:daqp_quadprog, DAQP.libdaqp), Nothing,
          (Ref{DAQP.DAQPResult},Ref{DAQP.QPc},Ref{DAQP.DAQPSettings}), 
          result,Ref(qp),ptr_settings)

    info = (x = xstar, λ=lam, fval=result[].fval,
            exitflag=result[].exitflag,
            status = DAQP.flag2status[result[].exitflag],
            solve_time = result[].solve_time,
            setup_time = result[].setup_time,
            iterations= result[].iter, nodes = result[].nodes)
    return xstar,result[].fval,result[].exitflag,info
end

"""
    xstar, fval, exitflag, info = DAQP.linprog(f,A,bupper)
    xstar, fval, exitflag, info = DAQP.linprog(f,A,bupper,blower,sense)

finds the solution `xstar` to the linear program

```
min_x	f' x
subject to 
    blower[1:ms]	<= x[1:ms] <= bupper[1:ms]
    blower[ms+1:m]  <= A*x 	   <= bupper[ms+1:m],
```
where `m = length(bupper)` and `ms = m-size(A,2)`.

# Input 
* `f`  			- linear term in objective function, n-vector 
* `A`  			- constraint normals, (`(m-ms) x n`)-matrix
* `bupper` 		- upper bounds for constraints, `m`-vector
* `blower` 		- lower bounds for constraints, `m`-vector (default: -Inf)
* `sense` 		- constraint types,  `m`-vector of Cints (default: 0). Example types:
  * `0` : inequality
  * `5` : equality

# Output
* `xstar` 		- solution provided by solver
* `fval` 		- objective function value for `xstar`. 
* `exitflag` 	- flag from solver (>0 success, <0 failure) 
* `info` 		- tuple containing profiling information from the solver. 

"""
function linprog(f::Vector{Float64}, 
        A::Matrix{Float64},bupper::Vector{Float64},blower::Vector{Float64}=Float64[],sense::Vector{Cint}=Cint[];A_rowmaj=false)
    d = DAQP.Model() 
    DAQP.setup(d,QPj(zeros(0,0),f,A,bupper,blower,sense;A_rowmaj))
    return DAQP.solve(d);
end
"""
    d = DAQP.Model() 
creates an empty optimization model `d`. 


The following functions acts on such models: 

* `DAQP.setup(d,H,f,A,bupper,blower,sense)`: setup a QP problem (see `DAQP.quadprog` for problem details)
* `DAQP.solve(d)`: solve a populated model
* `DAQP.update(d,H,f,A,bupper,blower,sense)`: update an existing model 
* `dict = DAQP.settings(d)`: return a Dictionary with the current settings for the model `d` 
* `DAQP.settings(d,dict)`: update the settings for `d` with the Dictionary `dict` 
"""
mutable struct Model 
    work::Ptr{DAQP.Workspace}
    qpj::QPj
    qpc::QPc
    qpc_ptr::Ptr{DAQP.QPc}
    has_model::Bool
    function Model()
        # Setup initial model
        work = Libc.calloc(1,sizeof(DAQP.Workspace))
        daqp= new(Ptr{DAQP.Workspace}(work))
        daqp.qpc_ptr = Libc.calloc(1,sizeof(DAQP.QPc))
        ccall((:allocate_daqp_settings,DAQP.libdaqp),Nothing,(Ptr{DAQP.Workspace},),work)
        finalizer(DAQP.delete!, daqp)
        daqp.has_model=false
        return daqp 
    end
end

function delete!(daqp::DAQP.Model)
    if(daqp.work != C_NULL)
        ccall((:free_daqp_workspace,DAQP.libdaqp),Nothing,(Ptr{DAQP.Workspace},),daqp.work)
        Libc.free(daqp.work);
        daqp.work = C_NULL
        Libc.free(daqp.qpc_ptr);
        daqp.qpc_ptr = C_NULL
    end
end

function setup(daqp::DAQP.Model, qp::DAQP.QPj)
    daqp.qpj = qp
    daqp.qpc = DAQP.QPc(daqp.qpj)
    old_settings = settings(daqp); # in case setup fails
    unsafe_store!(daqp.qpc_ptr,daqp.qpc)
    setup_time = Cdouble(0);

    if(isempty(qp.H) && !isempty(qp.f))# LP
        # ensure their is no binary constraint
        @assert(!any((qp.sense.&BINARY).==BINARY),
                "DAQP requires the objective to be strictly convex to support binary variables")
        # ensure proximal-point iterations are used for LPs
        (old_settings.eps_prox == 0) && settings(daqp,Dict(:eps_prox=>1))
    end

    exitflag = ccall((:setup_daqp,DAQP.libdaqp),Cint,(Ptr{DAQP.QPc}, Ptr{DAQP.Workspace}, Ptr{Cdouble}), daqp.qpc_ptr, daqp.work, Ref{Cdouble}(setup_time))
    if(exitflag < 0)
        # XXX: if setup fails DAQP currently clears settings
        ccall((:allocate_daqp_settings,DAQP.libdaqp),Nothing,(Ptr{DAQP.Workspace},),daqp.work)
        settings(daqp,old_settings)
    else
        daqp.has_model = true
        # Quick fix to initialize x for proximal
        # will be fixed in DAQP_jll 0.3.2
        if((isempty(qp.H)&&!isempty(qp.f)) || old_settings.eps_prox != 0)
            workspace = unsafe_load(daqp.work);
            xinit = zeros(workspace.n)
            unsafe_copyto!(workspace.x,pointer(xinit),workspace.n);
        end
        # should be unsafe_copy_to
    end

    return exitflag, setup_time
end

function setup(daqp::DAQP.Model, H::Matrix{Cdouble},f::Vector{Cdouble},
        A::Matrix{Cdouble},bupper::Vector{Cdouble},blower::Vector{Cdouble}=Cdouble[],
        sense::Vector{Cint}=Cint[];A_rowmaj=false)
    return setup(daqp,QPj(H,f,A,bupper,blower,sense;A_rowmaj))
end

function solve(daqp::DAQP.Model)
    if(!daqp.has_model) return  zeros(0), NaN, -10, [] end
    xstar = zeros(Float64,daqp.qpc.n); 
    lam = zeros(Float64,daqp.qpc.m); 
    result= Ref(DAQPResult(xstar,lam));

    exitflag=ccall((:daqp_solve, DAQP.libdaqp), Cint,
                   (Ref{DAQP.DAQPResult},Ref{DAQP.Workspace}), 
                   result,daqp.work)

    info = (x = xstar, λ=lam, fval=result[].fval,
            exitflag=result[].exitflag,
            status = DAQP.flag2status[result[].exitflag],
            solve_time = result[].solve_time,
            setup_time = result[].setup_time,
            iterations= result[].iter, nodes = result[].nodes)
    return xstar,result[].fval,result[].exitflag,info
end

function settings(daqp::DAQP.Model)
    workspace = unsafe_load(daqp.work);
    if(workspace.settings != C_NULL)
        return unsafe_load(workspace.settings)
    end
end
function settings(daqp::DAQP.Model, new_settings::DAQP.DAQPSettings)
    workspace = unsafe_load(daqp.work);
    if(workspace.settings != C_NULL)
        unsafe_store!(workspace.settings,new_settings)
    end
    return new_settings
end
function settings(daqp::DAQP.Model,changes::Dict{Symbol,<:Any})
    workspace = unsafe_load(daqp.work);
    if(workspace.settings == C_NULL) return end
    settings = unsafe_load(workspace.settings)
    new = [haskey(changes,f) ? changes[f] : getfield(settings,f)
           for f in fieldnames(DAQP.DAQPSettings)];
    new_settings = DAQP.DAQPSettings(new...)
    unsafe_store!(workspace.settings,new_settings);
    return new_settings;
end

function update(daqp::DAQP.Model, H,f,A,bupper,blower,sense) 
    update_mask = Cint(0);
    work = unsafe_load(daqp.work);
    if(!isnothing(H) && work.n == size(H,1) && work.n == size(H,2))
        daqp.qpj.H[:].=H[:]
        update_mask +=1
    end
    if(!isnothing(A) && size(A,1)==(work.m-work.ms) && size(A,2)==work.n)
        daqp.qpj.A[:].=A'[:]
        update_mask+=2
    end

    if(!isnothing(f) && length(f)==work.n)
        daqp.qpj.f[:].=f[:]
        update_mask+=4
    end

    if(!isnothing(bupper) && !isnothing(blower) &&
       length(bupper)==work.m && length(blower)==work.m)
        daqp.qpj.bupper[:].=bupper[:]
        daqp.qpj.blower[:].=blower[:]
        update_mask+=8
    end

    if(!isnothing(sense) && length(sense)== work.m)
        daqp.qpj.sense[:] .= sense[:]
        update_mask+=16
    end
    daqp.qpc = QPc(daqp.qpj);
    unsafe_store!(work.qp,daqp.qpc);

    exitflag = ccall((:update_ldp,DAQP.libdaqp),Cint,(Cint,Ptr{DAQP.Workspace},), update_mask, daqp.work);
end

function codegen(d::DAQP.Model; fname="daqp_workspace", dir="",)
    @assert(d.has_model, "setup the model before code generation")
    exitflag = ccall((:render_daqp_workspace, libdaqp),Cvoid,
                     (Ptr{DAQP.Workspace},Cstring,Cstring,), d.work,fname,dir);
end

function setup_c_workspace(n)::Ptr{Cvoid}
    p = Libc.calloc(1,sizeof(DAQP.Workspace));
    ccall((:allocate_daqp_workspace,libdaqp), Cvoid, (Ptr{Cvoid},Cint,Cint),p, n, 0);
    ccall((:allocate_daqp_settings,libdaqp), Cvoid, (Ptr{Cvoid},),p);
    return p
end

function free_c_workspace(p::Ptr{Cvoid})
    ccall((:free_daqp_workspace,libdaqp), Cvoid, (Ptr{Cvoid},),p)
    Libc.free(p)
end

function init_c_workspace_ldp(p::Ptr{Cvoid},A::Matrix{Cdouble},bupper::Vector{Cdouble},blower::Vector{Cdouble},sense::Vector{Cint}; max_radius=nothing) 
    # Set fval_bound to maximal radius for early termination
    if(!isnothing(max_radius))
        d_work = unsafe_load(Ptr{DAQP.Workspace}(p));
        unsafe_store!(Ptr{Cdouble}(d_work.settings+fieldoffset(DAQPSettings,8)),max_radius);
    end

    # Pass pointers for LDP to DAQP
    unsafe_store!(Ptr{Ptr{Cdouble}}(p+fieldoffset(DAQP.Workspace,5)),pointer(A))
    unsafe_store!(Ptr{Ptr{Cdouble}}(p+fieldoffset(DAQP.Workspace,6)),pointer(bupper))
    unsafe_store!(Ptr{Ptr{Cdouble}}(p+fieldoffset(DAQP.Workspace,7)),pointer(blower))
    unsafe_store!(Ptr{Ptr{Cint}}(p+fieldoffset(DAQP.Workspace,10)),pointer(sense))
end

function isfeasible(p::Ptr{Cvoid}, m=nothing, ms=nothing ;validate=false)::Bool
    # Update A and bupper/blower dimensions 
    !isnothing(m)  && unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,3)),m)
    !isnothing(ms) && unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,4)),ms)

    exitflag =ccall((:daqp_ldp,DAQP.libdaqp), Int32, (Ptr{Cvoid},),p);

    if(validate && exitflag == -1)
        daqp_ws = unsafe_load(Ptr{DAQP.Workspace}(p))
        m,ms,n = daqp_ws.m, daqp_ws.ms, daqp_ws.n
        A = unsafe_wrap(Matrix{Cdouble}, daqp_ws.M, (n,m), own=false)
        b= unsafe_wrap(Vector{Cdouble}, daqp_ws.dupper, m+ms, own=false)
        AS = copy(unsafe_wrap(Vector{Cint}, daqp_ws.WS, daqp_ws.n_active, own=false))
        AS .+= 1 # Offset for C -> Julia
        lam_star = unsafe_wrap(Vector{Cdouble}, daqp_ws.lam_star, daqp_ws.n_active, own=false)
        err  = dot(b[AS],lam_star)+norm(A[:,AS]*lam_star)
        if(err>0)
            @warn "Couldn't validate infeas. with Frakas (err:$(err), fval=$(daqp_ws.fval))"
        end
    end

    # Reset the workspace 
    ccall((:deactivate_constraints,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);
    ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);
    #if(exitflag != 1 && exitflag != -1)
    #    @warn "exitflag: $exitflag"
    #end
    return exitflag == 1
end
