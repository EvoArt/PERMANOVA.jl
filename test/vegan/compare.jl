# Compare PERMANOVA.jl with vegan::adonis2 on reference.R's output.
# Usage: julia --project=. test/vegan/compare.jl OUTDIR   (with DataFrames available)
using PERMANOVA, DelimitedFiles, Test
dir = ARGS[1]
readcsv(f) = (m = readdlm(joinpath(dir, f), ','; header=true); (m[1], vec(String.(strip.(m[2], '"')))))
dat, hdr = readcsv("data.csv")
data = (a=String.(strip.(string.(dat[:, 1]), '"')), b=String.(strip.(string.(dat[:, 2]), '"')), c=Float64.(dat[:, 3]))
Y = Float64.(readcsv("Y.csv")[1])
P = Int.(readcsv("perm.csv")[1])
function check(file, f, by)
    ref, _ = readcsv(file)
    r = permanova(data, Y, BrayCurtis(), f; permutations=P, by=by)
    k = size(ref, 1)
    @testset "$by" begin
        @test r.terms == String.(strip.(string.(ref[:, 1]), '"'))
        @test r.df == Int.(ref[:, 2])
        @test r.sumsq ≈ Float64.(ref[:, 3]) rtol = 1e-10
        @test r.r2 ≈ Float64.(ref[:, 4]) rtol = 1e-10
        nt = k - 2
        @test r.F[1:nt] ≈ Float64.(ref[1:nt, 5]) rtol = 1e-10
        @test r.p[1:nt] == Float64.(ref[1:nt, 6])          # same permutations, same p
    end
end
@testset "vegan adonis2" begin
    check("terms.csv", @formula(1 ~ a + b + c + a & b), :terms)
    check("margin.csv", @formula(1 ~ a + b + c), :margin)
    check("model.csv", @formula(1 ~ a + b + c), :model)
end
