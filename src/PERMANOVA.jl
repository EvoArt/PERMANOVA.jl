module PERMANOVA

using Random, LinearAlgebra, Distances, StatsModels,Statistics,TexTables,LoopVectorization, NamedArrays,DataFrames
include("perm2.jl")
include("output.jl")
if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end
export permanova, 
    hydra,
    permute, 
    Euclidean,
    SqEuclidean,
    PeriodicEuclidean,
    Cityblock,
    TotalVariation,
    Chebyshev,
    Minkowski,
    Jaccard,
    BrayCurtis,
    RogersTanimoto,
    @formula

precompile(A_mul_B!,(Array{Float64},Array{Float64},Array{Float64}))
precompile(permanova,(DataFrame ,Array{Float64},FormulaTerm ,Int64))
end
