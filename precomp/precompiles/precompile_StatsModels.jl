function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(concrete_term),Term,Vector{String},Nothing})   # time: 0.4439514
    Base.precompile(Tuple{typeof(modelcols),MatrixTerm{Tuple{CategoricalTerm{DummyCoding, String, 3}}},NamedTuple{(:X, :Z, :W), Tuple{Vector{String}, Vector{String}, Vector{String}}}})   # time: 0.0605238
    Base.precompile(Tuple{typeof(+),Tuple{Term, Term},Term})   # time: 0.0308432
    Base.precompile(Tuple{Type{ContrastsMatrix},Matrix{Float64},Vector{String},Vector{String},DummyCoding})   # time: 0.0124063
    Base.precompile(Tuple{typeof(schema),Vector{Term},NamedTuple{(:X, :Z, :W), Tuple{Vector{String}, Vector{String}, Vector{String}}},Dict{Symbol, Any}})   # time: 0.0070532
    Base.precompile(Tuple{typeof(getindex),ContrastsMatrix{DummyCoding, String, String},Vector{String},Function})   # time: 0.0040177
    Base.precompile(Tuple{typeof(contrasts_matrix),DummyCoding,Int64,Int64})   # time: 0.0032101
    Base.precompile(Tuple{Type{CategoricalTerm},Symbol,ContrastsMatrix{DummyCoding, String, String}})   # time: 0.0029878
    Base.precompile(Tuple{typeof(~),ConstantTerm{Int64},Tuple{Term, Term, Term}})   # time: 0.0017909
    Base.precompile(Tuple{typeof(termnames),DummyCoding,Vector{String},Int64})   # time: 0.0012678
    Base.precompile(Tuple{typeof(collect_matrix_terms),CategoricalTerm{DummyCoding, String, 3}})   # time: 0.0012
end
