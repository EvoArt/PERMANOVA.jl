struct PSummary
    results
    table
end

function Base.display(x::PSummary) 
     println("");println(x.table)
end

# swap columns i and j of a, in-place
function swaprows!(a::AbstractMatrix, i, j)
    i == j && return
    rows = axes(a,1)
    @boundscheck i in rows || throw(BoundsError(a, (i,:)))
    @boundscheck j in rows || throw(BoundsError(a, (j,:)))
    for k in axes(a,1)
        @inbounds a[i,k],a[j,k] = a[j,k],a[i,k]
    end
end
# like permute!! applied to each row of a, in-place in a (overwriting p).
function permuterows!!(a::AbstractMatrix, p::AbstractVector{<:Integer})
    #require_one_based_indexing(a, p)
    count = 0
    start = 0
    while count < length(p)
        ptr = start = findnext(!iszero, p, start+1)::Int
        next = p[start]
        count += 1
        while next != start
            swaprows!(a, ptr, next)
            p[ptr] = 0
            ptr = next
            next = p[next]
            count += 1
        end
        p[ptr] = 0
    end
    a
end



function A_mul_B!(C::Array{Float64}, A::Array{Float64}, B::Array{Float64})
    @turbo for n ∈ indices((C,B), 2), m ∈ indices((C,A), 1)
        Cmn = zero(eltype(C))
        for k ∈ indices((A,B), (2,1))
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] = Cmn
    end
end
function sumdiag(x::Array{Float64},n::Int64)
    s = 0
    @turbo for i in 1:n
    s += x[i,i]
    end
    return s
    end
function sumdiag(x::Array{Float64},y::Array{Float64},n::Int64)
    s = 0
    r = 0
    @turbo for i in 1:n
        s += x[i,i]
        r += y[i,i] - x[i,i]
    end
    return s,r
end

function dblcen(D)
    A = D .- mean(D,dims = 1)
    return A .- mean(A,dims = 2) 
end

function QRfit(x,y)
    Q, R = LinearAlgebra.qr(x)
    R_inv_Q_trans = R^(-1)*Q'
    β = R^(-1)*Q' *y
    fitted =   x * β
    residuals = y .-fitted
    return sum(Diagonal(fitted) ), sum(Diagonal(residuals) ),R_inv_Q_trans
end

function get_fitted_residuals(R_inv_Q_trans::Array{Float64},x::Array{Float64},y::Array{Float64},c::Array{Float64},n ::Int64)
    β = R_inv_Q_trans*y
    A_mul_B!(c, x, β)

    return sumdiag(c,y,n)
end
function get_fitted(R_inv_Q_trans::Array{Float64},x::Array{Float64},y::Array{Float64},c::Array{Float64},n ::Int64)
    β = R_inv_Q_trans *y
    A_mul_B!(c, x, β)
    return sumdiag(c,n)
end

function term_fit(x,y)
    Q, R = LinearAlgebra.qr(x)
    R_inv_Q_trans = R^(-1)*Q' 
    β = R_inv_Q_trans*y
    fitted =   x * β
    return sum(Diagonal(fitted)), R_inv_Q_trans
end

function unpack_formula(form)
    rhs = form.rhs
    lhs = form.lhs
    if typeof(rhs) == Term # if theres only a single rhs term
        return [form], [string(rhs)]
    else
        return     FormulaTerm.([lhs],rhs), string.(rhs)
    end
end

function permanova(data::DataFrame ,D ::Array{Float64}, formula::FormulaTerm = @formula(1~1), n_perm ::Int64 = 999)
    
    n = size(D)[1]
    G = Hermitian(-0.5 * dblcen(D .^2))
    fitted, residual = 0.0,0.0
    formulae,terms = unpack_formula(formula) #list of model formulae for analysis and names for coef table
    n_terms = length(terms) 
    mod_mats = Vector{Matrix}(undef,n_terms)
    R_inv_Q_trans = Vector{Matrix}(undef,n_terms)
    Df = zeros(Int64,n_terms)
    sumsq = zeros(Float64,n_terms)
    term_matrix = ones(Int64,n,1) # just the intercept (not run, but added to later matrices)
    @inbounds for i in 1:n_terms
        mat_size = size(term_matrix)[2]
        formula = formulae[i]
        term_matrix = hcat(term_matrix,StatsModels.modelmatrix(formula,data))
        mod_mats[i] =term_matrix # store model matrices for permutation
        Df[i] = size(term_matrix)[2] -mat_size
        if i == n_terms # if full model
            fitted, residual,R_inv_Q_trans[i] = QRfit(term_matrix,G)
            sumsq[i] =fitted-sum(sumsq)
        else
        tmp,R_inv_Q_trans[i] = term_fit(term_matrix,G)
        sumsq[i] =tmp-sum(sumsq)
        end
    end
    Tot = fitted + residual
    Resid_Df = n - 1- sum(Df)
    f_terms = (sumsq ./Df) ./(residual /Resid_Df)
    r2 = sumsq ./Tot
    C = Vector{Array{Float64}}(undef,n_terms)
    for i in 1:n_terms
        C[i] = Array{Float64}(undef,size(mod_mats[i])[1],size(G)[2])
    end

    p = permute(G, n, n_terms, mod_mats,R_inv_Q_trans,n_perm,C)  

    regtab = get_output(terms,Df,sumsq,r2,f_terms,residual,Tot,p,n) 
    return PSummary(regtab[1],regtab[2])
    

end

function permanova(data::DataFrame,M::Array{Float64},metric ::DataType, formula::FormulaTerm = @formula(1~1), n_perm ::Int64 = 999)
    D = pairwise(metric(),M',M')
return  permanova(data,D, formula,n_perm )
end

hydra = permanova
function permute(G ::Hermitian, n ::Int64, n_terms ::Int64, mod_mats ::Vector{Matrix},R_inv_Q_trans::Vector{Matrix},n_perm ::Int64, C::Vector)  
   
    inds = collect(1:n)
    indscopy = copy(inds)
    perms = Array{Float64}(undef,n_terms,n_perm+1)
    fit = zeros(n_terms)
    f_terms = zeros(n_terms)
    Gres = 0.0
    g = Array(G)
    @inbounds for j in 1:n_perm +1
        permutecols!!(g, inds), setup=(copyto!(indscopy, inds))
        permuterows!!(g, inds), setup=(copyto!(indscopy, inds))
        #g .=view(G,inds,inds)
      # Gres = Resfit(Qs[end],Rs[end],mod_mats[end],g)
       fit .= zeros(n_terms)
       @inbounds for i in 1:n_terms
        if i == n_terms
            prevsum = sum(fit)
        fit[i], Gres = get_fitted_residuals(R_inv_Q_trans[i],mod_mats[i],g,C[i],n)
        fit[i] -=prevsum
        else
        fit[i] = get_fitted(R_inv_Q_trans[i],mod_mats[i],g,C[i],n) - sum(fit)
        end
    end
    f_terms .= (fit) ./(Gres)
    perms[:,j] .= f_terms
    shuffle!(inds)
    
end
p = sum( perms[:,1] .<perms[:,2:end],dims = 2) ./n_perm
return p
end
