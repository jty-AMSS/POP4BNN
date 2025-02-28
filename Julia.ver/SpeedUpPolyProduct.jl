import DynamicPolynomials

function polydot(A,x)
    # speed up compute A*X
    # where x is vectors of monomals, A is float Matrix
    n=length(x)
    m=size(A,1)
    F=Vector{Polynomial{true, Float64}}(undef,m)
    for i=1:m
        F[i]=polynomial(A[i,:],x)
    end
return F
end

function AplusBdotF(a,b,f)
    
    # compute (a.+b).*f
    # where a,b,f are vectors 
    # b can be a scalar
    m=length(f)
    F=Vector{Polynomial{true, Float64}}(undef,m)
if length(b)>1
    for i=1:m

        F[i]=a[i]*f[i]+b[i]*f[i]
    end
end
if length(b)==1
    for i=1:m
        F[i]=a[i]*f[i]+b*f[i]
    end
end
return F
    
end


function ComputeYMinusKTimesAXPlusB(y,k,a,x,b)
    # compute (y+k)*(<a,x>+b)
    # length(y)==length(k)=length(b)=1
    # length(a)=length(x)=n
    a=vec(a)
    mono=[y.*x;y;x;1]
    f=[a;b;k*a;k*b]
    F=polynomial(f,mono)
    return F
end


function ComputeYMinusKTimesAXPlusB_vec(y,k,A,x,b)
    # compute (y+k)*(A*x+b)
    # length(y)==length(k)=length(b)=m
    # length(a)=length(x)=n
    m=length(y)
    F=Vector{Polynomial{true, Float64}}(undef,m)
    b=vec(b)
    if length(k)==1
    for i=1:m
        F[i]=ComputeYMinusKTimesAXPlusB(y[i],k,A[i,:],x,b[i])
    end
    return F
    end     
    if length(k)>1
        for i=1:m
            F[i]=ComputeYMinusKTimesAXPlusB(y[i],k[i],A[i,:],x,b[i])
        end
        return F
        end  
    
end

