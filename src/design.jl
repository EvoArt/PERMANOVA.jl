# The design: one block of model-matrix columns per term, in formula order,
# via StatsModels (treatment contrasts by default, as R).

struct Design
    labels::Vector{String}             # term labels as R writes them ("a", "a:b")
    blocks::Vector{Matrix{Float64}}    # model-matrix columns of each term
    vars::Vector{Vector{Symbol}}       # variables of each term
    keep::Vector{Int}                  # rows used (complete cases)
end

_label(t) = replace(string(t), " & " => ":")

function design(data, formula::FormulaTerm; contrasts=Dict{Symbol,Any}())
    tbl = Tables.columntable(data)
    f = FormulaTerm(ConstantTerm(1), formula.rhs)       # the left side is not used
    vars = collect(StatsModels.termvars(f.rhs))
    for v in vars
        haskey(tbl, v) || throw(ArgumentError("variable `$v` in the formula is not a column of the data"))
    end
    n = Tables.rowcount(tbl)
    keep = [i for i in 1:n if all(v -> !ismissing(tbl[v][i]), vars)]
    if length(keep) < n
        @warn "$(n - length(keep)) row(s) dropped because of missing values"
    end
    sub = NamedTuple{Tuple(vars)}(Tuple(_nomissing(tbl[v][keep]) for v in vars))
    sch = StatsModels.schema(f, sub, contrasts)
    f = StatsModels.apply_schema(f, sch, StatsModels.StatisticalModel)
    rhs = f.rhs isa MatrixTerm ? f.rhs.terms : (f.rhs,)
    terms = [t for t in rhs if !(t isa InterceptTerm) && !(t isa ConstantTerm)]
    isempty(terms) && throw(ArgumentError("the formula has no terms on its right-hand side"))
    m = length(keep)
    blocks = [Matrix{Float64}(reshape(modelcols(t, sub), m, :)) for t in terms]
    return Design(_label.(terms), blocks, [collect(StatsModels.termvars(t)) for t in terms], keep)
end

_nomissing(v) = [x for x in v]       # narrows Union{Missing,T} to T

# Orthonormal basis of the blocks `which`, built term by term after the
# intercept; returns the basis and the rank added by each block. A column
# set whose residual (after earlier terms) is below 1e-7 of its norm adds no
# rank, as R's QR in lm/adonis2 does for aliased terms.
function basis(blocks, which, n)
    Q = ones(n, 1) ./ sqrt(n)
    ranks = Int[]
    for b in which
        X = blocks[b]
        R = X .- Q * (Q' * X)
        R .-= Q * (Q' * R)                       # second pass for stability
        F = svd(R)
        tol = 1e-7 * max(maximum(norm, eachcol(X); init=0.0), eps())
        r = count(>(tol), F.S)
        push!(ranks, r)
        r > 0 && (Q = hcat(Q, F.U[:, 1:r]))
    end
    return Q[:, 2:end], ranks
end
