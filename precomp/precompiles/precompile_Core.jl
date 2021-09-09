function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Type{NamedTuple{(:X, :Z, :W), T} where T<:Tuple},Core.Tuple{Base.AbstractVector{T} where T, Base.AbstractVector{T} where T, Base.AbstractVector{T} where T}})   # time: 0.0056522
end
