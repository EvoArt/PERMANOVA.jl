
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
for nperm in 1000:1000:20000
    J = @benchmark permanova(df,y,BrayCurtis,@formula(1~X+Z+W),nperm)
    R = @benchmark R"adonis2($y ~X+Z+W,$df,nperm)"
end


