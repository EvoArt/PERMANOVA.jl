struct PSummary
    results
    table
end

function Base.display(x::PSummary) 
     println("");println(x.table)
end

function A_mul_B!(C, A, B)
    @turbo for n ∈ indices((C,B), 2), m ∈ indices((C,A), 1)
        Cmn = zero(eltype(C))
        for k ∈ indices((A,B), (2,1))
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] = Cmn
    end
end
function sumdiag(x,n)
    s = 0
    @turbo for i in 1:n
    s += x[i,i]
    end
    return s
    end
function sumdiag(x,y,n)
    s = 0
    r = 0
    @turbo for i in 1:n
        s += x[i,i]
        r += y[i,i] - x[i,i]
    end
    return s,r
end

function doube_center(D)
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

function get_fitted_residuals(R_inv_Q_trans,x,y,c,n)
    β = R_inv_Q_trans*y
    A_mul_B!(c, x, β)

    return sumdiag(c,y,n)
end
function get_fitted(R_inv_Q_trans,x,y,c,n)
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

function permanova(data,D, formula = @formula(1~1), n_perm = 999;blocks = false)

    n = size(D)[1]
    formulae,terms = unpack_formula(formula) #list of model formulae for analysis and names for coef table
    try
        inds = vec(sum(ismissing.(Array(data[!,terms])),dims = 2) .==0)
        cc = sum(inds) # no. complete cases
        if cc < n
            data = data[(1:n)[inds],:]
            D = D[inds,inds]
            @warn "$(n-cc) data row(s) dropped due to missing values."
            n = cc
        end
    catch e
        @warn "Unable to check for missing values. Please make sure there are no missing values in relevant variables."
    end

    data = float.(data)

    G = Hermitian(-0.5 * doube_center(float.(D) .^2))
    fitted, residual = 0.0,0.0
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
    if blocks == false
        p = permute(G, n, n_terms, mod_mats,R_inv_Q_trans,n_perm,C)  
    else
        block_mat = StatsModels.modelmatrix(blocks,data)
        p = permute(G, n, n_terms, mod_mats,R_inv_Q_trans,n_perm,C,block_mat) 
    end 



    regtab = get_output(terms,Df,sumsq,r2,f_terms,residual,Tot,p,n) 
    return PSummary(regtab[1],regtab[2])
    

end

function permanova(data,M,metric ::DataType, formula= @formula(1~1), n_perm = 999;blocks = false)
    D = pairwise(metric(),M',M')
return  permanova(data,D, formula,n_perm , blocks = blocks)
end

hydra = permanova
function permute(G ::Hermitian, n, n_terms , mod_mats ,R_inv_Q_trans,n_perm , C)    
    inds = collect(1:n)
    indscopy = copy(inds)
    perms = Array{Float64}(undef,n_terms,n_perm+1)
    fit = zeros(n_terms)
    f_terms = zeros(n_terms)
    res = 0.0
    g = Array(G)
    @inbounds for j in 1:n_perm +1
        g .= g[inds,inds]
       fit .= zeros(n_terms)
       @inbounds for i in 1:n_terms
            if i == n_terms
                prevsum = sum(fit)
                fit[i], res = get_fitted_residuals(R_inv_Q_trans[i],mod_mats[i],g,C[i],n)
                fit[i] -=prevsum
            else
                fit[i] = get_fitted(R_inv_Q_trans[i],mod_mats[i],g,C[i],n) - sum(fit)
            end
    end
    f_terms .= (fit) ./(res)
    perms[:,j] .= f_terms
    shuffle!(inds)
end
p = sum( perms[:,1] .<=perms[:,1:end],dims = 2) ./(n_perm+1)
return p
end

function permute(G ::Hermitian, n , n_terms , mod_mats ,R_inv_Q_trans,n_perm , C,block_mat)  
    blockviews = []
    inds = collect(1:n)
    for col in eachcol(block_mat)
        push!(blockviews,view(inds,Bool.(collect(col))))
    end
    perms = Array{Float64}(undef,n_terms,n_perm+1)
    fit = zeros(n_terms)
    f_terms = zeros(n_terms)
    res = 0.0
    g = Array(G)
    @inbounds for j in 1:n_perm +1
        g .= g[inds,inds]
       fit .= zeros(n_terms)
       @inbounds for i in 1:n_terms
        if i == n_terms
            prevsum = sum(fit)
            fit[i], res = get_fitted_residuals(R_inv_Q_trans[i],mod_mats[i],g,C[i],n)
            fit[i] -=prevsum
        else
            fit[i] = get_fitted(R_inv_Q_trans[i],mod_mats[i],g,C[i],n) - sum(fit)
        end
    end
    f_terms .= (fit) ./(res)
    perms[:,j] .= f_terms
    shuffle!.(blockviews)
    
end
p = sum( perms[:,1] .<=perms[:,1:end],dims = 2) ./(n_perm +1)
return p
end
