function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(get_fitted),Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Int64})   # time: 17.953636
    Base.precompile(Tuple{typeof(permanova),DataFrame,Matrix{Float64},DataType,FormulaTerm{ConstantTerm{Int64}, Tuple{Term, Term, Term}},Int64})   # time: 5.8466363
    Base.precompile(Tuple{typeof(get_fitted_residuals),Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},Int64})   # time: 0.6654432
    Base.precompile(Tuple{typeof(permanova),DataFrame,Matrix{Float64},FormulaTerm{ConstantTerm{Int64}, Tuple{Term, Term, Term}},Int64})   # time: 0.523199
    Base.precompile(Tuple{typeof(QRfit),Matrix{Float64},Hermitian{Float64, Matrix{Float64}}})   # time: 0.0766731
    Base.precompile(Tuple{typeof(term_fit),Matrix{Float64},Hermitian{Float64, Matrix{Float64}}})   # time: 0.0535826
end
