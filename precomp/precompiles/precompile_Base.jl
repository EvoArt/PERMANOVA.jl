function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(deepcopy_internal),Tuple{Symbol, Symbol, Symbol},IdDict{Any, Any}})   # time: 0.0203449
    Base.precompile(Tuple{typeof(deepcopy_internal),Tuple{Int64},IdDict{Any, Any}})   # time: 0.0182394
    Base.precompile(Tuple{typeof(deepcopy_internal),Tuple{Int64, Int64, Int64},IdDict{Any, Any}})   # time: 0.0170211
    Base.precompile(Tuple{typeof(deepcopy_internal),Tuple{Symbol},IdDict{Any, Any}})   # time: 0.0164719
    Base.precompile(Tuple{typeof(symdiff),Vector{String},Vector{String}})   # time: 0.0073157
    Base.precompile(Tuple{typeof(_shrink),Function,Vector{String},Tuple{Vector{String}}})   # time: 0.0035465
    Base.precompile(Tuple{typeof(deepcopy_internal),Vector{Int32},IdDict{Any, Any}})   # time: 0.0031766
    Base.precompile(Tuple{typeof(copyto!),Matrix{Float64},Int64,Vector{Float64},Int64,Int64})   # time: 0.0019458
    Base.precompile(Tuple{typeof(setindex_widen_up_to),Vector{DataType},Any,Int64})   # time: 0.001777
    Base.precompile(Tuple{typeof(fieldcount),Type{Tuple{String, String, String}}})   # time: 0.0013452
    Base.precompile(Tuple{typeof(vectorfilter),Function,Vector{String}})   # time: 0.0012713
end
