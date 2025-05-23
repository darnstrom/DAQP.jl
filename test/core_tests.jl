using LinearAlgebra
using Random
using Test
using DAQP
include(joinpath(dirname(@__FILE__), "utils.jl"))

# API Tests
nQPs,nLPs = 100,100;
n = 100; m = 500; ms = 50;
nAct = 80
kappa = 1e2
tol = 1e-4
@testset "Quadprog (C )" begin
  for nQP in 1:nQPs
	xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa);
	x,fval,exitflag,info = DAQP.quadprog(H,f,A,bupper,blower,sense);
	@test norm(xref-x) < tol;
  end
  # Test quadprog interface by passing settings 
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa);
  s = DAQP.settings(DAQP.Model(),Dict(:iter_limit=>1))
  x,fval,exitflag,info = DAQP.quadprog(H,f,A,bupper,blower,sense;settings=s)
  @test exitflag == -4
end

@testset "Quadprog (one-sided)" begin
  for nQP in 1:10
	_,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa);
    blower = fill(-1e30,length(bupper))
	xref,fval,exitflag,info = DAQP.quadprog(H,f,A,bupper,blower,sense);
	x,fval,exitflag,info = DAQP.quadprog(H,f,A,bupper);
	@test norm(xref-x) < tol;
  end
end

@testset "Linprog (C )" begin
    for nQP in 1:nQPs
        xref,f,A,bupper,blower,sense = generate_test_LP(n,m,ms);
        x,fval,exitflag,info = DAQP.linprog(f,A,bupper,blower,sense);
        @test abs(f'*(xref-x)) < tol;
    end
end

@testset "Linprog (one-sided)" begin
  for nQP in 1:10
	_,f,A,bupper,blower,sense = generate_test_LP(n,m,ms);
    blower = fill(-1e30,length(bupper))
	xref,fval,exitflag,info = DAQP.linprog(f,A,bupper,blower,sense);
	x,fval,exitflag,info = DAQP.linprog(f,A,bupper);
    @test abs(f'*(xref-x)) < tol;
  end
end

@testset "BnB" begin
    nb = 10
    ϵb= 1e-5
    nMIQPs = 5
    for nMIQP = 1:nMIQPs
        M = randn(n,n)
        H = M'*M
        A = randn(m-ms,n);
        bupper = 20*rand(m); blower = -20*rand(m) # ensure that origin is feasible
        f = 100*randn(n); f[1:nb].=-abs.(f[1:nb]) # make it lucrative to avoid origin
        # Make first nb variables binary
        bupper[1:nb].=1
        blower[1:nb].=0
        sense = zeros(Cint,m)
        sense[1:nb].=DAQP.BINARY
        x,fval,exitflag,info = DAQP.quadprog(H,f,A,bupper,blower,sense);
        @test exitflag == 1 # was able to solve problem
        @test all((abs.(x[1:nb].-1.0).<ϵb) .| (abs.(x[1:nb]).<ϵb)) # is binary feasible
    end
end

@testset "Model interface" begin
  # Setup model and solve problem
  d = DAQP.Model()
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
  DAQP.setup(d,H,f,A,bupper,blower,sense)
  x,fval,exitflag,info = DAQP.solve(d)
  @test norm(xref-x) < tol
  @test norm(H*x+[I(n)[1:ms,:];A]'*info.λ+f) < tol

  # Test access settings
  s = DAQP.settings(d)
  @test s.primal_tol==1e-6
  DAQP.settings(d,Dict(:primal_tol=>1e-5))
  s = DAQP.settings(d)
  @test s.primal_tol==1e-5
  
  # Update existing model with new problem
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
  DAQP.update(d,H,f,A,bupper,blower,sense)
  x,fval,exitflag,info = DAQP.solve(d)
  @test norm(xref-x) < tol
  DAQP.delete!(d); # Test manually running destructor
end

@testset "C LDP interface" begin
  # Setup model and solve problem
  n = 10; m = 50; ms = 5; nAct =0;
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
  p=DAQP.setup_c_workspace(n)
  A = A'[:,:] # since row major...
  DAQP.init_c_workspace_ldp(p,A,bupper,blower,sense;max_radius=1e30) 
  @test isfeasible(p,m,ms)
  bupper[1] = -1e30 #Make trivially infeasible
  @test !isfeasible(p,m,ms;validate=true)
  work = unsafe_load(Ptr{DAQP.Workspace}(p));
  @test work.n == n
  DAQP.free_c_workspace(p)
end

@testset "Code generation" begin
  n = 5; m = 5; ms = 5; nAct =2;
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
  sense[1] = DAQP.BINARY
  d = DAQP.Model()
  DAQP.setup(d,H,f,A,bupper,blower,sense)
  DAQP.codegen(d,dir="codegen",src=true)
  rm("codegen",recursive=true)
end
