module Hydra

using Random, LinearAlgebra, Distances, NamedArrays, StatsModels,Statistics,TexTables,LoopVectorization
include("hydra2.jl")
include("output.jl")

export hydra2, 
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


end
