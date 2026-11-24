using PERMANOVA, Test, Distances, DataFrames, Random, Statistics, LinearAlgebra, StableRNGs

# Anderson (2001): SS_total = Σ d²/N, SS_within = Σ_groups Σ d²/n_g
function brute_oneway(D, g)
    N = length(g)
    sst = sum(D[i, j]^2 for i in 1:N for j in i+1:N) / N
    ssw = 0.0
    for l in unique(g)
        ii = findall(==(l), g)
        ssw += sum((D[i, j]^2 for i in ii for j in ii if i < j); init=0.0) / length(ii)
    end
    a = length(unique(g))
    return sst, sst - ssw, ((sst - ssw) / (a - 1)) / (ssw / (N - a))
end

# Euclidean case: sums of squares of a multivariate linear model
function lm_ss(X, Y)
    Yc = Y .- mean(Y; dims=1)
    H = X * pinv(X)
    return sum(abs2, H * Yc)
end

@testset "PERMANOVA" begin
    rng = StableRNG(1)
    n = 60
    df = DataFrame(a=rand(rng, ["x", "y", "z"], n), b=rand(rng, ["u", "v"], n), c=randn(rng, n))
    Y = rand(rng, n, 6) .+ (df.a .== "x") .* 0.3

    @testset "one-way: brute-force formula" begin
        D = pairwise(BrayCurtis(), Y; dims=1)
        r = permanova(df, D, @formula(1 ~ a); permutations=199, rng=StableRNG(2))
        sst, ssa, F = brute_oneway(D, df.a)
        @test r.sumsq[end] ≈ sst
        @test r.sumsq[1] ≈ ssa
        @test r.F[1] ≈ F
        @test r.df == [2, n - 3, n - 1]
        @test r.r2[1] ≈ ssa / sst
        @test r.terms == ["a", "Residual", "Total"]
        @test 0 < r.p[1] <= 1
    end

    @testset "Euclidean: sequential SS of a linear model" begin
        D = pairwise(Euclidean(), Y; dims=1)
        r = permanova(df, D, @formula(1 ~ a + b + c + a & b); permutations=0)
        one = ones(n)
        A = Float64.(hcat(df.a .== "y", df.a .== "z"))
        B = Float64.(df.b .== "v")
        AB = A .* B
        X1 = hcat(one, A); X2 = hcat(X1, B); X3 = hcat(X2, df.c); X4 = hcat(X3, AB)
        ss = [lm_ss(X1, Y), lm_ss(X2, Y), lm_ss(X3, Y), lm_ss(X4, Y)]
        @test r.sumsq[1:4] ≈ [ss[1]; diff(ss)]
        @test r.terms[1:4] == ["a", "b", "c", "a:b"]
        @test r.df[1:4] == [2, 1, 1, 2]
        @test sum(r.sumsq[1:5]) ≈ r.sumsq[6]
        # whole-model and marginal tests
        m = permanova(df, D, @formula(1 ~ a + b + c); by=:model, permutations=0)
        @test m.sumsq[1] ≈ lm_ss(X3, Y)
        g = permanova(df, D, @formula(1 ~ a + b + c); by=:margin, permutations=0)
        @test g.sumsq[3] ≈ lm_ss(X3, Y) - lm_ss(X2, Y)
        @test g.sumsq[1] ≈ lm_ss(X3, Y) - lm_ss(hcat(one, B, df.c), Y)
        # a and b are inside a:b, so only a:b and c are marginal
        gi = permanova(df, D, @formula(1 ~ a * b + c); by=:margin, permutations=0)
        @test Set(gi.terms[1:end-2]) == Set(["c", "a:b"])
    end

    @testset "data + metric forms, positional count, alias" begin
        r1 = permanova(df, Y, BrayCurtis(), @formula(1 ~ a); permutations=99, rng=StableRNG(3))
        r2 = permanova(df, Y, BrayCurtis, @formula(1 ~ a), 99; rng=StableRNG(3))
        r3 = hydra(df, pairwise(BrayCurtis(), Y; dims=1), @formula(1 ~ a), 99; rng=StableRNG(3))
        @test isequal(r1.F, r2.F) && isequal(r2.F, r3.F)
        @test isequal(r1.p, r2.p) && isequal(r2.p, r3.p)
        @test DataFrame(r1).term == ["a", "Residual", "Total"]
    end

    @testset "power and calibration" begin
        D = pairwise(BrayCurtis(), Y; dims=1)
        @test permanova(df, D, @formula(1 ~ a); permutations=499, rng=StableRNG(4)).p[1] < 0.01
        # unrelated variable: p-values roughly uniform
        ps = [permanova(DataFrame(z=rand(StableRNG(10 + s), ["p", "q"], n)), D, @formula(1 ~ z);
                        permutations=99, rng=StableRNG(s)).p[1] for s in 1:40]
        @test 0.25 < mean(ps) < 0.75
    end

    @testset "permutation matrix, strata, missing values, numeric as factor" begin
        D = pairwise(BrayCurtis(), Y; dims=1)
        P = reduce(vcat, (randperm(StableRNG(i), n)' for i in 1:50))
        p1 = permanova(df, D, @formula(1 ~ a); permutations=P)
        p2 = permanova(df, D, @formula(1 ~ a); permutations=P)
        @test isequal(p1.p, p2.p) && p1.permutations == 50
        # with strata, every permutation keeps each observation in its block
        s = permanova(df, D, @formula(1 ~ a); strata=df.b, permutations=99, rng=StableRNG(5))
        @test 0 < s.p[1] <= 1
        dm = allowmissing(df)
        dm.a[[1, 5]] .= missing
        r = @test_logs (:warn, r"2 row") permanova(dm, D, @formula(1 ~ a); permutations=0)
        @test r.df[end] == n - 3
        dose = rand(StableRNG(6), 1:3, n)
        dd = DataFrame(dose=dose)
        rc = permanova(dd, D, @formula(1 ~ dose); permutations=0)
        rf = permanova(dd, D, @formula(1 ~ dose); permutations=0,
                       contrasts=Dict(:dose => PERMANOVA.StatsModels.DummyCoding()))
        @test rc.df[1] == 1 && rf.df[1] == 2
    end

    @testset "pairwise" begin
        D = pairwise(BrayCurtis(), Y; dims=1)
        pw = pairwise_permanova(df, D, :a; permutations=199, rng=StableRNG(7))
        @test pw.group1 == ["x", "x", "y"] && pw.group2 == ["y", "z", "z"]
        ii = findall(in(("y", "z")), df.a)
        one = permanova(df[ii, :], D[ii, ii], @formula(1 ~ a); permutations=0)
        @test pw.F[3] ≈ one.F[1]
        @test all(pw.p_adjusted .>= pw.p)
        @test pw.p[1] < 0.05 && pw.p[2] < 0.05
        @test nrow(DataFrame(pw)) == 3
    end

    @testset "errors and display" begin
        D = pairwise(BrayCurtis(), Y; dims=1)
        @test_throws ArgumentError permanova(df, D, @formula(1 ~ nothere))
        @test_throws DimensionMismatch permanova(df[1:10, :], D, @formula(1 ~ a))
        @test_throws ArgumentError permanova(df, D, @formula(1 ~ a); by=:foo)
        r = permanova(df, D, @formula(1 ~ a + b); permutations=99)
        s = sprint(show, MIME("text/plain"), r)
        @test occursin("Residual", s) && occursin("Pr(>F)", s)
        @test occursin("pair", sprint(show, MIME("text/plain"), pairwise_permanova(df, D, :b; permutations=9)))
        @test adjust_p([0.01, 0.04, 0.03, 0.2], :holm) ≈ [0.04, 0.09, 0.09, 0.2]
        @test adjust_p([0.01, 0.04, 0.03, 0.2], :bh) ≈ [0.04, 0.16 / 3, 0.16 / 3, 0.2]
        @test adjust_p([0.01, 0.04, 0.03, 0.2], :bonferroni) ≈ [0.04, 0.16, 0.12, 0.8]
        @test adjust_p([0.01, 0.5], :none) == [0.01, 0.5]
        @test adjust_p(Float64[]) == Float64[]
        @test_throws ArgumentError adjust_p([0.01, 0.5], :fdr)
        @test sprint(show, r) == "PermanovaResult(a, b)"
        pwb = pairwise_permanova(df, D, :b; permutations=9)
        @test sprint(show, pwb) == "PairwisePermanova(b, 1 pairs)"
        @test_throws ArgumentError pairwise_permanova(df, D, :nothere)
        @test_throws ArgumentError pairwise_permanova(DataFrame(a=fill("x", n)), D, :a)
    end

    @testset "pairwise: data + metric forms, Tables.jl" begin
        D = pairwise(BrayCurtis(), Y; dims=1)
        ref = pairwise_permanova(df, D, :a; permutations=0)
        @test pairwise_permanova(df, Y, BrayCurtis(), :a; permutations=0).F ≈ ref.F
        @test pairwise_permanova(df, Y, BrayCurtis, :a; permutations=0, adjust=:bh).F ≈ ref.F
        @test PERMANOVA.Tables.columnaccess(typeof(ref))
        r = permanova(df, D, @formula(1 ~ a); permutations=0)
        @test PERMANOVA.Tables.columnaccess(typeof(r))
        @test DataFrame(r).term == r.terms
    end
end
