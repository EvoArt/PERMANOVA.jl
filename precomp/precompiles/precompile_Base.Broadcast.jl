function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(materialize),Broadcasted{DefaultArrayStyle{1}, Nothing, typeof(/), Tuple{Broadcasted{DefaultArrayStyle{1}, Nothing, typeof(/), Tuple{Vector{Float64}, Vector{Int64}}}, Float64}}})   # time: 0.0245907
    Base.precompile(Tuple{typeof(materialize),Broadcasted{DefaultArrayStyle{1}, Nothing, typeof(/), Tuple{Vector{Float64}, Float64}}})   # time: 0.0153087
    Base.precompile(Tuple{typeof(materialize!),Vector{Float64},Broadcasted{DefaultArrayStyle{1}, Nothing, typeof(/), Tuple{Vector{Float64}, Float64}}})   # time: 0.0049648
    Base.precompile(Tuple{typeof(broadcasted),Function,Vector{Float64},Float64})   # time: 0.0036101
    Base.precompile(Tuple{typeof(broadcasted),Function,Broadcasted{DefaultArrayStyle{1}, Nothing, typeof(/), Tuple{Vector{Float64}, Vector{Int64}}},Float64})   # time: 0.0034639
end
