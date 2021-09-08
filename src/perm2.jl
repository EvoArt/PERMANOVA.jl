

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
function dblcen(D)
    A = D .- mean(D,dims = 1)
    return A .- mean(A,dims = 2) 
end

function QRfit(x,y)
    Q, R = LinearAlgebra.qr(x)
    β = (R^(-1)*Q' *y)
    fitted =   x * β#sum(x .* β, dims = 2)
    residuals = y .-fitted
    return sum(Diagonal(fitted) ), sum(Diagonal(residuals) ),Q,R
end

function QRfit(Q,R,x,y)
    β = (R^(-1)*Q' *y)
    fitted =   x * β#sum(x .* β, dims = 2)
    residuals = y .-fitted
    return sum(Diagonal(fitted) )/ sum(Diagonal(residuals))
end

function Resfit(Q,R,x,y)
    β = (R^(-1)*Q' *y)
    fitted =   x * β#sum(x .* β, dims = 2)
    residuals = y .-fitted
    return sum(Diagonal(fitted) ),sum(Diagonal(residuals))
end
function Fitfit(Q,R,x,y)
    β = (R^(-1)*Q' *y)
    fitted =   x * β#sum(x .* β, dims = 2)
    return sum(Diagonal(fitted) )
end
function Fitfit(Q,R,x,y,n ::Int64)
    β = (R^(-1)*Q' *y)
    fitted =   x * β#sum(x .* β, dims = 2)
    return sumdiag(fitted,n)
end
function Resfit(Q,R,x,y,c,n)
    β = (R^(-1)*Q' *y)
    A_mul_B!(c, x, β)#sum(x .* β, dims = 2)

    return sumdiag(c,y,n)#sum(Diagonal(c) ),sum(Diagonal(resids))
end
function Fitfit(Q,R,x,y,c,n ::Int64)
    β = (R^(-1)*Q' *y)
    A_mul_B!(c, x, β)#sum(x .* β, dims = 2)
    return sumdiag(c,n)#sum(Diagonal(c) )
end



function term_fit(x,y)
    Q, R = LinearAlgebra.qr(x)
    β = (R^(-1)*Q' *y)
    fitted =   x * β
    return sum(Diagonal(fitted)), Q, R
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

function permanova(data,D, formula = @formula(1~1), n_perm ::Int64 = 999; pairs = false)
    
    n = size(D)[1]
    G = Hermitian(-0.5 * dblcen(D .^2))
    fitted, residual = 0.0,0.0
    formulae,terms = unpack_formula(formula) #list of model formulae for analysis and names for coef table
    n_terms = length(terms) 
    mod_mats = Vector{Matrix}(undef,n_terms)
    Qs = Vector{Matrix}(undef,n_terms)
    Rs = Vector{Matrix}(undef,n_terms)
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
            fitted, residual,Qs[i], Rs[i] = QRfit(term_matrix,G)
            sumsq[i] =fitted-sum(sumsq)
        else
        tmp,Qs[i], Rs[i] = term_fit(term_matrix,G)
        sumsq[i] =tmp-sum(sumsq)
        end
    end
    Tot = fitted + residual
    Resid_Df = n - 1- sum(Df)
    f_terms = (sumsq ./Df) ./(residual /Resid_Df)
    r2 = sumsq ./Tot
    p = permute(G, n, n_terms, mod_mats,Qs,Rs,n_perm)  

    regtab = get_output(terms,Df,sumsq,r2,f_terms,residual,Tot,p,n) 
    return regtab

end


function permanova(data,M,metric ::DataType, formula = @formula(1~1), n_perm ::Int64 = 999; pairs = false)
    D = pairwise(metric(),M',M')
return  permanova(data,D, formula,n_perm ; pairs = pairs)
end


function permute(G ::Hermitian, n ::Int64, n_terms ::Int64, mod_mats ::Vector{Matrix},Qs ::Vector{Matrix},Rs ::Vector{Matrix},n_perm ::Int64)  
    C = [Array{Float64}(undef,size(mod_mats[i])[1],size(G)[2]) for i in 1:n_terms]
    
    inds = collect(1:n)
    perms = Array{Float64}(undef,n_terms,n_perm+1)
    fit = zeros(n_terms)
    f_terms = zeros(n_terms)
    Gres = 0.0
    g =Array{Float64}(undef,n,n)
    @inbounds for j in 1:n_perm +1
        g .=view(G,inds,inds)
      # Gres = Resfit(Qs[end],Rs[end],mod_mats[end],g)
       fit .= zeros(n_terms)
       @inbounds for i in 1:n_terms
        if i == n_terms
            prevsum = sum(fit)
        fit[i], Gres = Resfit(Qs[i],Rs[i],mod_mats[i],g,C[i],n)
        fit[i] -=prevsum
        else
        fit[i] = Fitfit(Qs[i],Rs[i],mod_mats[i],g,C[i],n) - sum(fit)
        end
    end
    f_terms .= (fit) ./(Gres)
    perms[:,j] .= f_terms
    shuffle!(inds)
    
end
p = sum( perms[:,1] .<perms[:,2:end],dims = 2) ./n_perm
return p
end

