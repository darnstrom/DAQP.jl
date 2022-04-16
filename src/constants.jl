# Constraint types
const ACTIVE = 1
const LOWER = 2 
const IMMUTABLE= 4
const SOFT= 8 

# Exit Flags
const flag2status= Dict{Int,Symbol}(1 => :Optimal,
								   -1 => :Primal_Infeasible,
								   -2 => :Cycling,
								   -3 => :Unbounded,
								   -4 => :Iteration_Limit
								   )
