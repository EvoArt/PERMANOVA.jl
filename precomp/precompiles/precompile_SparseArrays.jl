function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(hcat),Matrix{Int64},Matrix{Float64}})   # time: 0.0145908
end
