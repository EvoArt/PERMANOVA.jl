#The key to the non-parametric method describedhere 
#is that the sum of squared distances between points and 
#their centroid is equal to (and can be calculated directly from) 
#the sum of squared interpoint distancesdivided by the number of points
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1442-9993.2001.01070.pp.x

SS(D,N) = sum(D .^2)/N 

function F(D,W,N,n,a)
    SST = SS(D,N)
    SSW = SS(W,n)
    SSA = SST - SSW
    return (SSA/(a-1))/(SSW/(N-a))
end

function F(SST ::Real,W,N,n,a)

    SSW = SS(W,n)
    SSA = SST - SSW
    return (SSA/(a-1))/(SSW/(N-a))
end

function getinds(group,N)
    W =Vector{CartesianIndex}(undef,0)
    notW = Vector{CartesianIndex}(undef,0)
    for i in 1:(N-1)
        for j in i+1:N 
            if i !==1j
                if group[i] == group[j]
                    push!(W,(CartesianIndex(i,j)))
                else push!(notW,(CartesianIndex(i,j)))
                end
            end
        end
    end
    return W, notW
end

function perm1(D,group, n_perm = 1000)
    N =  length(group)
    a = length(unique(group))
    n = N/a
    Winds ,notWinds =getinds(group,N)

    @inbounds W = D[Winds]
    @inbounds notW= D[notWinds]
    
    Dvec = vcat(W,notW)
    
    SST = SS(Dvec,N) 
   
    f = F(SST,W,N,n,a)
    Fs = Vector{Float64}(undef,n_perm)
    inds = 1:length(W)
    
    @inbounds for i in 1:n_perm
        shuffle!(Dvec)
        Fs[i] = F(SST,view(Dvec,inds),N,n,a)
    end
    P = sum(Fs .>= f)/n_perm

    g = unique(group)
    if a>2
        ppairs = NamedArray(zeros(a,a),( g,g ), ("group","group"))
        tpairs = NamedArray(zeros(a,a),( g,g ), ("group","group"))
        Threads.@threads for i in 1:a-1
            for j in i+1:a
                boolmask = (group .== g[i]) .| (group .== g[j])
                d = D[boolmask,boolmask]
                fstat,pstat =permutest(d,group[boolmask], n_perm)
                ppairs[j,i] = pstat
                tpairs[j,i] = sqrt(fstat)
            end
        end
    return (F = f,P = P, fpairs = tpairs,ppairs = ppairs)
    else
        return (F = f,P = P)
    end

end
