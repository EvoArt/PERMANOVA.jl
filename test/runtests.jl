
using PerMANOVA, DataFrames,Distances, RCall, BenchmarkTools
using Test
R"library(vegan)"

x = rand(1:4,100)
z = rand(1:4,100)
y = rand(100,5) .+ (x .* (x .==3) .*(z .>2))

preds = [["a","b","c","d"][i] for i in x]
preds2 = [["a","b","c","d"][i] for i in z]
preds3 = [rand(["a","b","c","d"]) for i in x]
df = DataFrame([preds,preds2,preds3],[:X,:Z,:W])
Hyde = hydra(df,y,BrayCurtis,@formula(1~X+Z+W),1000)
R"ad =adonis2($y ~X+Z+W,$df,1000)"
n_terms = 3

@test sum(R"ad$Df" .== Hyde.results[:,"Df"]) == n_terms +2
@test sum(R"ad$SumOfSqs".≈ Hyde.results[:,"SumOfSqs"]) == n_terms +2
@test sum(R"ad$R2".≈ Hyde.results[:,"R2"]) == n_terms +2
@test sum(R"ad$F".≈ Hyde.results[:,"F"]) == n_terms


x = rand(1:4,1000)
z = rand(1:4,1000)
y = rand(1000,5) .+ (x .* (x .==3) .*(z .>2))

preds = [["a","b","c","d"][i] for i in x]
preds2 = [["a","b","c","d"][i] for i in z]
preds3 = [rand(["a","b","c","d"]) for i in x]
df = DataFrame([preds,preds2,preds3],[:X,:Z,:W])
Hyde = hydra(df,y,BrayCurtis,@formula(1~X+Z+W),10000)
R"ad =adonis2($y ~X+Z+W,$df,10000)"
n_terms = 3

@test sum(R"ad$Df" .== Hyde.results[:,"Df"]) == n_terms +2
@test sum(R"ad$SumOfSqs".≈ Hyde.results[:,"SumOfSqs"]) == n_terms +2
@test sum(R"ad$R2".≈ Hyde.results[:,"R2"]) == n_terms +2
@test sum(R"ad$F".≈ Hyde.results[:,"F"]) == n_terms
