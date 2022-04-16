using LinearAlgebra
using Random
using Test
using DAQP
include(joinpath(dirname(@__FILE__), "utils.jl"))

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
  for nQP in 1:25
	xref,H,f,A,bupper,blower,sense = generate_test_QP(20,100,0,16,1e2);
	x,fval,exitflag,info = DAQP.daqp_jl(H,f,[A;-A],[bupper;-blower],[sense;sense],Int64[]);
	@test norm(xref-x) < tol;
  end
end

@testset "Model interface" begin
  # Setup model and solve problem
  d = DAQP.Model()
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
  DAQP.setup(d,H,f,A,bupper,blower,sense)
  x,fval,exitflag,info = DAQP.solve(d)
  @test norm(xref-x) < tol

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
end

@testset "Model interface" begin
  # Setup model and solve problem
  d = DAQP.Model()
  xref,H,f,A,bupper,blower,sense = generate_test_QP(n,m,ms,nAct,kappa)
  DAQP.setup(d,H,f,A,bupper,blower,sense)
  x,fval,exitflag,info = DAQP.solve(d)
  @test norm(xref-x) < tol

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
