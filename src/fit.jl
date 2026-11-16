"""
    PermanovaResult

An ANOVA-style table from [`permanova`](@ref). Rows are the tested terms
(or `"Model"`), then `"Residual"` and `"Total"`. Columns: `terms`, `df`,
`sumsq`, `r2`, `F`, `p` (`NaN` where not defined). It is a Tables.jl table,
so `DataFrame(result)` works.
"""
struct PermanovaResult
    terms::Vector{String}
    df::Vector{Int}
    sumsq::Vector{Float64}
    r2::Vector{Float64}
    F::Vector{Float64}
    p::Vector{Float64}
    permutations::Int
    by::Symbol
    n::Int
    formula::String
end

Tables.istable(::Type{PermanovaResult}) = true
Tables.columnaccess(::Type{PermanovaResult}) = true
Tables.columns(r::PermanovaResult) =
    (term=r.terms, Df=r.df, SumOfSqs=r.sumsq, R2=r.r2, F=r.F, P=r.p)

"""
    permanova(data, D, formula; permutations=999, by=:terms, strata=nothing,
              rng=Random.default_rng(), contrasts=Dict())
    permanova(data, Y, metric, formula; kw...)

PERMANOVA of the dissimilarity matrix `D` (or of the rows of `Y` under a
Distances.jl `metric`, e.g. `BrayCurtis()`) on the variables of `data` (any
table, e.g. a `DataFrame`) given by the right-hand side of `formula`. The
left-hand side is ignored: write `@formula(1 ~ a + b)`.

The sums of squares, R², and pseudo-F statistics are those of
`vegan::adonis2`. P-values are `(1 + #{F_perm ≥ F}) / (1 + permutations)`,
from permutations of the observations (of rows and columns of `D`).

- `by=:terms`: terms added sequentially, first to last (order matters in
  unbalanced designs); `:margin`: each term given all others (terms within
  interactions are not tested); `:model`: one test of the whole model.
- `strata`: a vector of block labels; permutations stay within blocks.
- `permutations`: a number, or a matrix with one permutation of `1:n` per row.
- `contrasts`: StatsModels contrasts or hints, e.g. `Dict(:dose => DummyCoding())`
  to treat a numeric column as categorical.

Rows with a missing value in a formula variable are dropped (with a warning).
"""
function permanova(data, D::AbstractMatrix, formula::FormulaTerm; permutations=999, by::Symbol=:terms,
                   strata=nothing, rng::AbstractRNG=Random.default_rng(), contrasts=Dict{Symbol,Any}())
    by in (:terms, :margin, :model) || throw(ArgumentError("by must be :terms, :margin or :model, got :$by"))
    des = design(data, formula; contrasts=contrasts)
    N = size(D, 1)
    size(D, 2) == N || throw(ArgumentError("D must be a square dissimilarity matrix, got $(size(D))"))
    nrow = Tables.rowcount(Tables.columntable(data))
    nrow == N || throw(DimensionMismatch("data has $nrow rows but D has $N"))
    any(<(-sqrt(eps())), D) && throw(ArgumentError("dissimilarities must be non-negative"))
    keep = des.keep
    n = length(keep)
    G = gower(D[keep, keep])
    trG = tr(G)
    blocks = des.blocks
    nt = length(blocks)
    Qfull, ranks = basis(blocks, 1:nt, n)
    K = size(Qfull, 2)
    K > 0 || throw(ArgumentError("the terms have no variation (model rank 0)"))
    dfres = n - K - 1
    dfres > 0 || throw(ArgumentError("no residual degrees of freedom: $n observations, model rank $K"))

    # each test: a list of bases, and a map from their per-column fitted
    # sums of squares to the numerator sums of squares
    if by === :terms
        tested = findall(>(0), ranks)
        labels = des.labels[tested]
        dfs = ranks[tested]
        ends = cumsum(ranks)[tested]
        bases = [Qfull]
        numer = c -> (cs = cumsum(c[1]); [cs[e] - (e > d ? cs[e-d] : 0.0) for (e, d) in zip(ends, dfs)])
    elseif by === :model
        labels, dfs = ["Model"], [K]
        bases = [Qfull]
        numer = c -> [sum(c[1])]
    else
        scope = _margin_scope(des)
        reduced = [basis(blocks, setdiff(1:nt, [t]), n)[1] for t in scope]
        dfs = [K - size(Q, 2) for Q in reduced]
        ok = dfs .> 0
        scope, reduced, dfs = scope[ok], reduced[ok], dfs[ok]
        isempty(scope) && throw(ArgumentError("no term can be tested by margin"))
        labels = des.labels[scope]
        bases = vcat([Qfull], reduced)
        numer = c -> (full = sum(c[1]); [full - sum(c[j+1]) for j in eachindex(dfs)])
    end
    stat = c -> (ss = numer(c); res = trG - sum(c[1]); (ss ./ dfs) ./ (res / dfres))

    c0 = [_colss(G, Q) for Q in bases]
    ss0 = numer(c0)
    res0 = trG - sum(c0[1])
    F0 = stat(c0)

    pmat = permutations isa AbstractMatrix ? _check_perms(permutations, N, keep) : nothing
    nperm = pmat === nothing ? Int(permutations) : size(pmat, 1)
    nperm >= 0 || throw(ArgumentError("permutations must be ≥ 0"))
    blocksof = strata === nothing ? [collect(1:n)] : _strata_blocks(strata, keep)
    hits = _permute(G, bases, stat, F0, nperm, pmat, blocksof, rng)
    p = (hits .+ 1) ./ (nperm + 1)

    terms = vcat(labels, "Residual", "Total")
    df = vcat(dfs, dfres, n - 1)
    sumsq = vcat(ss0, res0, trG)
    return PermanovaResult(terms, df, sumsq, sumsq ./ trG, vcat(F0, NaN, NaN), vcat(p, NaN, NaN),
                           nperm, by, n, string(formula.rhs))
end

permanova(data, Y::AbstractMatrix, metric::PreMetric, formula::FormulaTerm; kw...) =
    permanova(data, pairwise(metric, Y; dims=1), formula; kw...)
permanova(data, Y::AbstractMatrix, metric::Type{<:PreMetric}, formula::FormulaTerm, args...; kw...) =
    permanova(data, Y, metric(), formula, args...; kw...)
# positional permutation count, as in 0.1
permanova(data, D::AbstractMatrix, formula::FormulaTerm, nperm::Integer; kw...) =
    permanova(data, D, formula; permutations=nperm, kw...)
permanova(data, Y::AbstractMatrix, metric::PreMetric, formula::FormulaTerm, nperm::Integer; kw...) =
    permanova(data, Y, metric, formula; permutations=nperm, kw...)

"`hydra` is another name for [`permanova`](@ref)."
const hydra = permanova

# Gower's centred matrix -1/2 (I - 11'/n) D² (I - 11'/n)
function gower(D)
    A = -0.5 .* Matrix{Float64}(D) .^ 2
    A .-= mean(A; dims=1)
    A .-= mean(A; dims=2)
    return Symmetric(A)
end

# per-column fitted sum of squares q_jᵀ G q_j
function _colss(G, Q)
    W = G * Q
    return [dot(view(W, :, j), view(Q, :, j)) for j in axes(Q, 2)]
end

# Count permutations with F ≥ observed (minus vegan's tolerance). Permuting
# the observations of G is the same as permuting the rows of each basis:
# tr(Qᵀ G[p,p] Q) = tr(Q̃ᵀ G Q̃) with Q̃ = Q[invperm(p), :].
function _permute(G, bases, stat, F0, nperm, pmat, blocksof, rng)
    nt = length(F0)
    hits = zeros(Int, nt)
    nperm == 0 && return hits
    n = size(G, 1)
    nchunk = min(nperm, Threads.nthreads())
    seeds = rand(rng, UInt64, nchunk)
    counts = [zeros(Int, nt) for _ in 1:nchunk]
    EPS = sqrt(eps())
    Gm = Matrix(G)
    Threads.@threads for c in 1:nchunk
        crng = Xoshiro(seeds[c])
        Qp = [similar(Q) for Q in bases]
        W = [similar(Q) for Q in bases]
        perm = collect(1:n)
        ip = similar(perm)
        cs = [zeros(size(Q, 2)) for Q in bases]
        for b in c:nchunk:nperm
            if pmat === nothing
                _shuffle_blocks!(perm, blocksof, crng)
            else
                perm .= view(pmat, b, :)
            end
            for i in 1:n
                ip[perm[i]] = i
            end
            for (k, Q) in enumerate(bases)
                @inbounds for j in axes(Q, 2), i in 1:n
                    Qp[k][i, j] = Q[ip[i], j]
                end
                mul!(W[k], Gm, Qp[k])
                for j in axes(Q, 2)
                    cs[k][j] = dot(view(W[k], :, j), view(Qp[k], :, j))
                end
            end
            Fb = stat(cs)
            for t in 1:nt
                counts[c][t] += Fb[t] >= F0[t] - EPS
            end
        end
    end
    foreach(c -> hits .+= c, counts)
    return hits
end

function _shuffle_blocks!(perm, blocksof, rng)
    for idx in blocksof
        vals = idx[randperm(rng, length(idx))]
        for (i, v) in zip(idx, vals)
            perm[i] = v
        end
    end
    return perm
end

function _strata_blocks(strata, keep)
    s = collect(strata)[keep]
    any(ismissing, s) && throw(ArgumentError("strata has missing values"))
    return [findall(==(u), s) for u in unique(s)]
end

function _check_perms(P, N, keep)
    size(P, 2) == N || throw(DimensionMismatch("permutation matrix has $(size(P, 2)) columns, need $N"))
    length(keep) == N || throw(ArgumentError("a permutation matrix cannot be used when rows are dropped for missing values"))
    for r in eachrow(P)
        isperm(r) || throw(ArgumentError("each row of the permutation matrix must be a permutation of 1:$N"))
    end
    return P
end

# terms not contained in a higher-order term (R's drop.scope)
function _margin_scope(des)
    nt = length(des.vars)
    return [t for t in 1:nt if !any(u -> u != t && length(des.vars[u]) > length(des.vars[t]) &&
                                           issubset(des.vars[t], des.vars[u]), 1:nt)]
end
