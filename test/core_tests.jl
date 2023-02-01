using LinearAlgebra
using Random
using Test
using DAQP
include(joinpath(dirname(@__FILE__), "utils.jl"))

# API Tests
nQPs = 100;
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
end

@testset "Quadprog (JL)" begin
  qpj = DAQP.QPj()
  for nQP in 1:22
	xref,H,f,A,bupper,blower,sense = generate_test_QP(20,100,0,16,1e2);
	x,lam,AS,J,iter= DAQP.daqp_jl(H,f,[A;-A],[bupper;-blower],[sense;sense],Int64[]);
	@test norm(xref-x) < tol;
  end
  # Test warm start
  xref,H,f,A,bupper,blower,sense = generate_test_QP(20,100,0,16,1e2);
  x,lam,AS,J,iter= DAQP.daqp_jl(H,f,[A;-A],[bupper;-blower],[sense;sense],Int64[1]);
  @test norm(xref-x) < tol;
  # Test infeasible problem
  x,lam,AS,J,iter= DAQP.daqp_jl([1.0 0; 0 1],zeros(2),[1.0 0;-1 0],[1;-2],zeros(Cint,2),Int64[]);
  @test isinf(J)
  # Test unconstrained problem
  x,lam,AS,J,iter= DAQP.daqp_jl((H=[1.0 0; 0 1], f=zeros(2),
								 A=[1.0 0;0 1], b=[1;1],
								 senses=zeros(Cint,2)),Int64[]);
  @test isempty(AS)
  # Test with pivoting  TODO: add informative example
  xref,H,f,A,bupper,blower,sense = generate_test_QP(20,100,0,16,1e2);
  x,lam,AS,J,iter= DAQP.daqp_jl(H,f,[A;-A],[bupper;-blower],[sense;sense],Int64[]; pivoting=true);
  @test norm(xref-x) < tol;

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
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
  p=DAQP.setup_c_workspace(n)
  DAQP.init_c_workspace_ldp(p,A,bupper,blower,sense;max_radius=1.0) 
  work = unsafe_load(Ptr{DAQP.Workspace}(p));
  @test work.n == n
  DAQP.free_c_workspace(p)
end
