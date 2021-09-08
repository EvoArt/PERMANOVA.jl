
function get_output(coefs,Df,sq,r2,f_terms,Gres,Tot,p,n)
n_term = length(coefs)
DF = RegCol("Df")
R2 = RegCol("R²")
sumsq = RegCol("SumOfSqs")
F = RegCol("F")
P = RegCol("P")


for i in 1:n_term
setcoef!(R2, coefs[i]=>r2[i])
setcoef!(DF, coefs[i]=>Df[i])
setcoef!(sumsq, coefs[i]=>sq[i])
setcoef!(F, coefs[i]=>f_terms[i])
setcoef!(P, coefs[i]=>p[i])
end
coefArray = hcat(Df,sq,r2,f_terms,p)
coefArray= vcat(coefArray,[n-1-sum(Df) 1-sum(r2) Gres NaN NaN ])
coefArray= vcat(coefArray,[n-1 1 Tot NaN NaN ])


setcoef!(R2, "Residual"=>1-sum(r2))
setcoef!(DF, "Residual"=>n-1-sum(Df))
setcoef!(sumsq, "Residual"=>Gres)
setcoef!(R2, "Total"=>1)
setcoef!(DF, "Total"=>n-1)
setcoef!(sumsq, "Total"=>Tot)
return NamedArray(coefArray,(vcat(coefs...,"Residual","Total"),["Df","SumOfSqs","R2","F","P"])), hcat(DF,sumsq,R2,F,P)
end

