"""
    PERMANOVA

Permutational multivariate analysis of variance on a dissimilarity matrix
(Anderson 2001; McArdle & Anderson 2001), with the same statistics as
`vegan::adonis2`.

```julia
permanova(df, D, @formula(1 ~ treatment + site))
permanova(df, Y, BrayCurtis(), @formula(1 ~ treatment))
pairwise_permanova(df, D, :treatment)
```
"""
module PERMANOVA

using LinearAlgebra, Random, Statistics
using Distances
using StatsModels
using Tables
import Printf
using Printf: @sprintf

export permanova, hydra, pairwise_permanova, PermanovaResult, PairwisePermanova, adjust_p, @formula
# metrics re-exported since 0.1
export Euclidean, SqEuclidean, PeriodicEuclidean, Cityblock, TotalVariation, Chebyshev,
       Minkowski, Jaccard, BrayCurtis, RogersTanimoto

include("design.jl")
include("fit.jl")
include("pairwise.jl")
include("show.jl")

end
