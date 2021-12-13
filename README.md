# PERMANOVA

[![Build Status](https://github.com/EvoArt/Hydra.jl/workflows/CI/badge.svg)](https://github.com/EvoArt/Hydra.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/EvoArt/Hydra.jl?svg=true)](https://ci.appveyor.com/project/EvoArt/Hydra-jl)
[![Coverage](https://codecov.io/gh/EvoArt/Hydra.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/EvoArt/Hydra.jl)

PERMANOVA implementation based on the work of [McArdle and Anderson](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/0012-9658%282001%29082%5B0290%3AFMMTCD%5D2.0.CO%3B2). This package aims to provide similar functionality to the `adonis2` in the R package [`vegan`](https://cran.r-project.org/web/packages/vegan/index.html). Thus, implementation details are more similar to `adonis2` than to the original work by McArdle and Anderson. P-values are calculated via permuting of data (not residuals) and calculating sequential sums of squares. In keeping with names from mythology, we provide an alias to the `permanova` function :`hydra` the many headed Lernaean Hydra represents the multivariate response data we aim to tackle here. Though perhaps Heracles or Iolis (the eventual slayers of the Hydra) would be more apt.

<img src="https://github.com/EvoArt/Hydra/blob/master/docs/Sargent_Hercules.jpg" alt="drawing" width="400"/>

The function `permanova`/`hydra2` expects:

*   data: a table/dataframe with a column containing the independent variables. 
*   y: the dependent variables, where each row is an observation
*   metric: distance metric to be used.
* formula: a [StatsModels.jl](https://juliastats.org/StatsModels.jl/stable/formula/) formula 

Alternatively, instead of y and metric, pass in a distance matrix D.
The function retruns a `PSummary` struct containing `table` - a [TexTables.jl](https://jacobadenbaum.github.io/TexTables.jl/stable/) formatted ANOVA table for display purposes and `results` - a NamedArray for easier access to specific results.

## Example
```julia
using PERMANOVA, DataFrames,Distances
x = rand(1:4,100)
y = rand(100,5)

preds = [[rand(["a","b","c","d"]) for i in 1:100] for j in 1:3]
df = DataFrame(preds,[:X,:Y,:Z])
permanova(df,y,BrayCurtis,@formula(1~X+Y))
```

## TODO
*   Pairwise PerMANOVA

## Experimental
*   `blocks` (strata) e.g. `permanova(df,y,BrayCurtis,@formula(1~X+Y), blocks = @formula(block ~ Z))`
    
