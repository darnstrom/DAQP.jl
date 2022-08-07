using Test
@testset "Core" begin
    include("core_tests.jl")
end

@testset "MOI" begin
    include("MOI_wrapper.jl")
end


