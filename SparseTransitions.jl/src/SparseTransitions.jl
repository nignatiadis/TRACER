module SparseTransitions

using Gurobi
using JuMP
using LinearAlgebra:tr,I

include("main_func.jl")

export fast_minim

end # module
