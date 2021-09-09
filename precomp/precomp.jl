
using Pkg
Pkg.activate(".")
using DataFrames,Distances,  BenchmarkTools
x = rand(1:4,100)
y = rand(100,5) .+ (x .* (x .==3))

preds = [["a","b","c","d"][i] for i in x]
preds2 = [rand(["a","b","c","d"]) for i in x]
preds3 = [rand(["a","b","c","d"]) for i in x]
df = DataFrame([preds,preds2,preds3],[:X,:Z,:W])

using SnoopCompile,StatsModels,ProfileView
tinf1 =@snoopi_deep begin
    using PERMANOVA
end
ProfileView.view(flamegraph(tinf1))
tinf2 =@snoopi_deep begin
    Hyde = hydra(df,y,BrayCurtis,@formula(1~X+Z+W),10000)
end
ProfileView.view(flamegraph(tinf2))
ttot, pcs = SnoopCompile.parcel(tinf2)
SnoopCompile.write("./precomp/precompiles", pcs)