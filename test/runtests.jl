
using PERMANOVA, DataFrames,Distances, RCall, BenchmarkTools
using Test
R"library(vegan)"

x = rand(1:4,100)
y = rand(100,5) .+ (x .* (x .==3))

preds = [["a","b","c","d"][i] for i in x]
preds2 = [rand(["a","b","c","d"]) for i in x]
preds3 = [rand(["a","b","c","d"]) for i in x]
df = DataFrame([preds,preds2,preds3],[:X,:Z,:W])
Hyde = hydra(df,y,BrayCurtis,@formula(1~X+Z+W),100)
R"adonis2($y ~X+Z+W,$df,1000)"

J = @benchmark permanova($df,$y,BrayCurtis,@formula(1~X+Z+W),999)
R = @benchmark R"adonis2($y ~X+Z+W,$df,999)"
@test mean(J.times) < mean(R.times)

J = @benchmark permanova($df,$y,BrayCurtis,@formula(1~X+Z+W),9999)
R = @benchmark R"adonis2($y ~X+Z+W,$df,9999)"
@test mean(J.times) < mean(R.times)
