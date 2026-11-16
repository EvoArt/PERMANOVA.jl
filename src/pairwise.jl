"""
    PairwisePermanova

Result of [`pairwise_permanova`](@ref): one row per pair of groups, with
`group1`, `group2`, `n`, `F`, `r2`, `p` and `p_adjusted`. A Tables.jl table.
"""
struct PairwisePermanova
    variable::String
    group1::Vector{String}
    group2::Vector{String}
    n::Vector{Int}
    F::Vector{Float64}
    r2::Vector{Float64}
    p::Vector{Float64}
    p_adjusted::Vector{Float64}
    adjust::Symbol
    permutations::Int
end

Tables.istable(::Type{PairwisePermanova}) = true
Tables.columnaccess(::Type{PairwisePermanova}) = true
Tables.columns(r::PairwisePermanova) =
    (group1=r.group1, group2=r.group2, n=r.n, F=r.F, R2=r.r2, P=r.p, P_adjusted=r.p_adjusted)

"""
    pairwise_permanova(data, D, variable::Symbol; permutations=999, adjust=:holm,
                       strata=nothing, rng=Random.default_rng())
    pairwise_permanova(data, Y, metric, variable; kw...)

A PERMANOVA of `D ~ variable` for every pair of levels of the categorical
column `variable`, on the observations of those two levels only. The
p-values are corrected for multiple testing with `adjust` (`:holm`,
`:bh`, `:bonferroni` or `:none`).
"""
function pairwise_permanova(data, D::AbstractMatrix, variable::Symbol; permutations=999,
                            adjust::Symbol=:holm, strata=nothing, rng::AbstractRNG=Random.default_rng())
    tbl = Tables.columntable(data)
    haskey(tbl, variable) || throw(ArgumentError("`$variable` is not a column of the data"))
    x = tbl[variable]
    lv = sort(unique(skipmissing(x)))
    length(lv) >= 2 || throw(ArgumentError("`$variable` has fewer than two levels"))
    f = FormulaTerm(ConstantTerm(1), Term(variable))
    hint = Dict{Symbol,Any}(variable => StatsModels.DummyCoding())
    g1, g2, ns, Fs, r2s, ps = String[], String[], Int[], Float64[], Float64[], Float64[]
    for a in 1:length(lv)-1, b in a+1:length(lv)
        idx = findall(v -> !ismissing(v) && (v == lv[a] || v == lv[b]), x)
        sub = (; variable => [x[i] for i in idx])
        st = strata === nothing ? nothing : collect(strata)[idx]
        r = permanova(sub, D[idx, idx], f; permutations=permutations, strata=st, rng=rng, contrasts=hint)
        push!(g1, string(lv[a])); push!(g2, string(lv[b])); push!(ns, length(idx))
        push!(Fs, r.F[1]); push!(r2s, r.r2[1]); push!(ps, r.p[1])
    end
    nperm = permutations isa Integer ? Int(permutations) : 0
    return PairwisePermanova(string(variable), g1, g2, ns, Fs, r2s, ps, adjust_p(ps, adjust), adjust, nperm)
end

pairwise_permanova(data, Y::AbstractMatrix, metric::PreMetric, variable::Symbol; kw...) =
    pairwise_permanova(data, pairwise(metric, Y; dims=1), variable; kw...)
pairwise_permanova(data, Y::AbstractMatrix, metric::Type{<:PreMetric}, variable::Symbol; kw...) =
    pairwise_permanova(data, Y, metric(), variable; kw...)

"""
    adjust_p(p, method=:holm)

Adjust p-values for multiple testing: `:holm`, `:bh` (Benjamini-Hochberg),
`:bonferroni` or `:none`, as R's `p.adjust`.
"""
function adjust_p(p, method::Symbol=:holm)
    m = length(p)
    m == 0 && return Float64[]
    method === :none && return collect(float.(p))
    method === :bonferroni && return min.(1.0, m .* p)
    o = sortperm(p)
    q = similar(p, Float64)
    if method === :holm
        run = 0.0
        for (r, i) in enumerate(o)
            run = max(run, min(1.0, (m - r + 1) * p[i]))
            q[i] = run
        end
    elseif method === :bh
        run = 1.0
        for r in m:-1:1
            i = o[r]
            run = min(run, m / r * p[i])
            q[i] = min(run, 1.0)
        end
    else
        throw(ArgumentError("adjust must be :holm, :bh, :bonferroni or :none, got :$method"))
    end
    return q
end
