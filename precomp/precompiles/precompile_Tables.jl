function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{Type{Schema},Vector{Symbol},Vector{DataType}})   # time: 0.0035495
    Base.precompile(Tuple{typeof(which(_columntable,(Tables.Schema{names, types},Any,)).generator.gen),Any,Any,Any,Any,Any})   # time: 0.0033783
end
