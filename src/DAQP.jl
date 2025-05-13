module DAQP

using Reexport
@reexport using DAQPBase
@reexport using DAQPBase: DAQPSettings, DAQPResult
Model = DAQPBase.Model
QPj,QPc,Workspace = DAQPBase.QPj,DAQPBase.QPc, DAQPBase.Workspace

ACTIVE = DAQPBase.ACTIVE 
LOWER = DAQPBase.LOWER
IMMUTABLE = DAQPBase.IMMUTABLE
EQUALITY = DAQPBase.EQUALITY 
SOFT = DAQPBase.SOFT
BINARY = DAQPBase.BINARY 

SOFT_OPTIMAL = DAQPBase.SOFT_OPTIMAL
OPTIMAL = DAQPBase.OPTIMAL
INFEASIBLE = DAQPBase.INFEASIBLE
CYCLE = DAQPBase.CYCLE
UNBOUNDED = DAQPBase.UNBOUNDED
ITERLIMIT = DAQPBase.ITERLIMIT
NONCONVEX = DAQPBase.NONCONVEX
OVERDETERMINED = DAQPBase.OVERDETERMINED

flag2status = DAQPBase.flag2status

codegen = DAQPBase.codegen
setup_c_workspace = DAQPBase.setup_c_workspace
free_c_workspace = DAQPBase.free_c_workspace 
init_c_workspace_ldp = DAQPBase.init_c_workspace_ldp
delete! = DAQPBase.delete!

include("MOI_wrapper/MOI_wrapper.jl")
end
