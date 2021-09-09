function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(setindex!),RegCol{1},FNum{Int64},TableIndex{3}})   # time: 0.0012955
end
