using LinearAlgebra 

function daqp_ldp_jl(M,d,AS0,senses;settings=DAQPSettings())
  # Initial AS is union AS0 and equality constraints  
  AS = findall((senses.&IMMUTABLE).!=0); # Maybe this could be done better (i.e., input the indices..)
  AS = AS ∪ AS0;
  IS = setdiff(1:length(d),AS);
  lambda=ones(length(AS));
  L = zeros(0,0);
  D = zeros(0);
  for (k,ind) in enumerate(AS) # Setup LDL' for AS
	m = M[ind,:];
	L,D=updateLDLadd(L,D,M[AS[1:k-1],:]*m,m'*m);
  end

  iter = 1;
  fin = false;
  ASs = falses(length(d),0);
  inequality_mask = ((senses.&IMMUTABLE).==0)
  while(~fin)
	ASs = [ASs falses(length(d),1)];
	ASs[AS,end].=true;
	singular_ind= findfirst(D.<=settings.zero_tol); # Zero element in D signifies singularity
	if(isnothing(singular_ind)) # "Normal iteration 
	  if(!isempty(AS))
		y= -forward_L(L,d[AS]);
		z= y./D;
		lambda_star= backward_L(L,z);
	  else
		lambda_star=Float64[]
	  end
	  if(isempty(AS)||all(lambda_star[inequality_mask[AS]].>= 0)) # Local optimum
		if(isempty(AS))
		  mu = d[IS];
		else
		  mu= M[IS,:]*(M[AS,:]'*lambda_star)+d[IS];
		end
		lambda = lambda_star;
		if(all(mu.>= -settings.primal_tol)) #Global optimum
		  fin = true;
		else # Negative mu
		  _,free_ind = findmin(mu[:,1]);# Decide which component to free
		  m = M[IS[free_ind],:];
		  L,D = updateLDLadd(L,D,M[AS,:]*m,m'*m);
		  #Update working set
		  push!(AS,IS[free_ind]);
		  deleteat!(IS,free_ind);
		  push!(lambda,0);
		end
	  else #Constrained stationary point not primal feasible 
		p = lambda_star-lambda;
		block_inds = inequality_mask[AS] .& (lambda_star.<0)
		lambda,AS,IS,L,D = remove_constraint(lambda,AS,IS,L,D,p,block_inds);
	  end
	else #Singular 
	  p = compute_singulardirection(L,D,singular_ind);
	  block_inds = inequality_mask[AS] .& (p.<-settings.zero_tol)
	  if(sum(block_inds)==0) # Infeasible  
		return zeros(0),zeros(0),AS,Inf,iter,ASs;
	  end
	  lambda,AS,IS,L,D = remove_constraint(lambda,AS,IS,L,D,p,block_inds);
	end
	if(!fin)
	  iter = iter + 1;
	end
  end
  # Optimum found 
  lam_opt = zeros(length(d),1);
  if(isempty(AS))
	u=zeros(size(M,2)); # Unconstrained solution
  else
	lam_opt[AS] = lambda;
	u=(M[AS,:]'*lambda);
  end
  
  J = dot(u,u);
  return u,lam_opt,AS,J,iter,ASs;
end

function remove_constraint(lambda,AS,IS,L,D,p,block_inds) 
  alpha_cands = -lambda[block_inds]./p[block_inds];
  alpha,alpha_ind = findmin(alpha_cands);#Find first blocking constraint
  AS_tmp = AS[block_inds];
  fix_ind = findfirst(x->x==AS_tmp[alpha_ind],AS);
  lambda = lambda + alpha*p; #Update iterates and working set
  push!(IS,AS[fix_ind]);
  deleteat!(AS,fix_ind);
  deleteat!(lambda,fix_ind);
  L,D = updateLDLremove(L,D,fix_ind); 
  return lambda,AS,IS,L,D
end


function daqp_jl(H,f,A,b,sense,AS;settings=DAQPSettings())
  R = cholesky((H+H')/2);
  M = A/R.U;
  v = (R.L)\f;
  d = b+M*v;;
  # normalize
  norm_factor = 0
  for i in 1:size(M,1)
	norm_factor = norm(M[i,:],2); 
	M[i,:]./=norm_factor;
	d[i,:]./=norm_factor;
  end

  u,lam_opt,AS,Ju,iter,ASs = daqp_ldp_jl(M,d,AS,sense;settings)
  # Form solution to nomial problem
  if(!isempty(u))
	x_opt=R.U\(-u-v);
	J = 0.5*(Ju - dot(v,v));
  else # Infeasible
	x_opt = zeros(0);
	J = Inf;
  end
  
  return x_opt,lam_opt,AS,J,iter,ASs;
end

function daqp_jl(QP,AS;settings=DAQPSettings())
  return daqp_jl(QP.H,QP.f,QP.A,QP.b,QP.senses,AS;settings);
end


# ========= Factorization code ============
function forward_L(L,b)
  # Solve L x = b
  n = size(b,1);
  x = deepcopy(b);
  for i in 1:n
	for j in 1:(i-1)
	  x[i,:] -= L[i,j]*x[j,:];
	end
  end
  return x
end

function backward_L(L,b)
  # Solve L'x = b
  n = size(b,1);
  x = deepcopy(b); 
  for i = n:-1:1
	for j = i+1:n
	  x[i,:] -= L[j,i]*x[j,:];
	end
  end
  return x
end

function updateLDLadd(L,D,b,bet;precision=Float64)
  # Return Lup, Dup
  # compute Lup Dup s.t. Lup Dup Lup' = [A;a'] [A;a']'
  # when A A' = LDL'
  #b = A*a, bet = a'*a
  if(isempty(D))
	return ones(1,1),[D;bet] 
  end
  k = size(L,1);
  #c = L\b; # Implement own backwards solve?
  c = forward_L(L,b); # Implement own backwards solve?
  l = zeros(precision,k);
  for i = 1:k
	if(D[i]>1e-9)
	  l[i] = c[i]/D[i];
	end
  end
  d = bet-(l.*D)'*l;
  if(d<1e-9)
	d = 0;
  end
  Lup = [L zeros(precision,k,1);l' 1];
  return Lup,[D;d]
end

function updateLDLremove(L,D,i)
  # Output Lup, Dup
  # ith row removed
  # TODO no need to extract L1 and L3
  k = size(L,1);
  L1 = L[1:i-1,1:i-1];
  D1 = D[1:i-1];
  L3 = L[i+1:k,1:i-1];
 
  d = D[i];
  l =L[i+1:k,i];
  L2 =L[i+1:k,i+1:k];
  D2 = D[i+1:k];
  
  L2_tilde,D2_tilde =rankone_add_LDL(L2,D2,l,d);
  Lup = [L1 zeros(i-1,k-i);L3 L2_tilde];
  Dup = [D1;D2_tilde];
  return Lup, Dup
end

function rankone_add_LDL(L,D,l,delta;precision=Float64)
  # Out put Lup, Dup
  # Algorithm C1 Methods for modidfying matrix factorizations - Gill et al 1974
  # Compute Lup Dup such that Lup Dup Lup' = LDL'+delta*l l'
  n = length(l);
  Dup = zeros(precision,n) ;
  Lup = Matrix{precision}(I, n, n);

  a = delta;
  w = l ;
  for j = 1:n
	p = w[j] ;
	Dup[j] = D[j] + a*p^2;
	b = p*a/Dup[j];
	a = D[j]*a/Dup[j];
	for r = j+1:n
	  w[r] = w[r] - p*L[r,j] ;
	  Lup[r,j] = L[r,j] + b*w[r];
	end
  end
  return Lup,Dup
end

function compute_singulardirection(L,D,i;precision=Float64)
  k = length(D);p = zeros(precision,k);
  L1 = L[1:i-1,1:i-1]; l = L[i,1:i-1];
  p[1:i-1] = -backward_L(L1,l);
  p[i]= 1;
  return p
end
