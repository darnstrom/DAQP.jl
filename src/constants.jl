# Constraint types
const ACTIVE = 1
const LOWER = 2 
const IMMUTABLE= 4
const SOFT= 8 
const BINARY= 16 

# Exit Flags
const flag2status= Dict{Int,Symbol}(2 => :Soft_Optimal,
                                    1 => :Optimal,
								   -1 => :Primal_Infeasible,
								   -2 => :Cycling,
								   -3 => :Unbounded,
								   -4 => :Iteration_Limit,
                                   -5 => :Nonconvex_Problem,
                                   -6 => :Initial_Overdetermined
								   )
