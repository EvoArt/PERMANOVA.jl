function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(pairwise),BrayCurtis,Adjoint{Float64, Matrix{Float64}},Adjoint{Float64, Matrix{Float64}}})   # time: 0.241804
end
