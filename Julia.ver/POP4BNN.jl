

import LinearAlgebra
import TSSOS
import DynamicPolynomials
import DelimitedFiles
import JuMP

include("SpeedUpPolyProduct.jl")

function vecnorm(A, p, dim)
    if (p == 1) && (dim == 2)
        At = abs.(A)
        o = ones(size(At, 2), 1)
        return At * o
    end
end

function POP4BNNstd(A, b, f) 
    # This function is for L_inf verify
    # Problem statement:
    # min f(x_0,x_1,...,x_L)
    # s.t.  x_{k}=sign(A_{k-1}x_{k-1}+b_{k-1}),k=1,...,L
    # x_0 \in [-1,1]^n
    #  type of A:  Vector{Matrix{Float64}}
    #  type of b:  Vector{Matrix{Float64}}
    #  type of f:  Vector{{Float64}} 
    L = length(A) + 1;
    N = Vector{Int}(undef, L);
    N[1] = size(A[1], 2)
    for i = 1:length(A)
        N[i+1] = size(A[i], 1)
    end
    
    
    Layer_X = Vector{Vector{PolyVar{true}}}(undef, 0)
    
    n = N[1]
    # x=nothing;
    nt=sum(N);

    @polyvar y[1:nt];
    x=y[end:-1:1];
    n_start=1;
    for i = 1:L
        n_end=n_start+N[i]-1;
        tempX=x[n_start:n_end];
        push!(Layer_X, tempX);
        n_start=n_end+1;
    end
    x0=x[1:N[1]];
    # Obj_function
    G = Vector{Polynomial{true,Float64}}(undef, 0)
    
    push!(G,dot(f,Layer_X[end]));
    # L_inf
    
    x_0 = zeros(N[1], 1);
    vareps = 1.0;
    
    G = [G; (x0 - x_0 .+ vareps) .* (x_0 .+ vareps - x0)];
    
    
    for i = 1:(L-1)
        At = A[i]
        bt = b[i]
        xt = Layer_X[i]
        yt = Layer_X[i+1]
        println("Generate equs starting")
        F1=ComputeYMinusKTimesAXPlusB_vec(yt,1,At,xt,bt)
        # G = [G; F1]
        println("Generate equs stage 2")
        # G = [G; (yt .- 1) .* (At * xt + bt)]
        # F0=polydot(At,xt).+bt;F1=AplusBdotF(yt,-1,F0)
        F2=ComputeYMinusKTimesAXPlusB_vec(yt,-1,At,xt,bt)
        # G = [G; F1]
        At_norm=vecnorm(At,1,2);
        # F0=polydot(At,xt).+At_norm;F1=AplusBdotF(yt,-1,F0);F1=(-1).*F1;
         # G=[G;-1*(yt.-1).*(At_norm.+At*xt)];
         println("Generate equs stage 3")
        F3=ComputeYMinusKTimesAXPlusB_vec(yt,-1,At,xt,At_norm);
        println("Generate equs stage--3.1")
        F3=-F3;
        # G = [G; F1];
        # G=[G;(yt.+1).*(At_norm.-At*xt)]; 
        println("Generate equs stage 4")
        F4=ComputeYMinusKTimesAXPlusB_vec(yt,1,-At,xt,At_norm);
        println("Now marge the vector")
        # G = [G; F1];
        G=[G;F1;F2;F3;F4]
       
        
        println("Generate equs ended")
    end
    
    
    var = Vector{PolyVar{true}}(undef, 0);
    NumBVar=0;
    for i = 2:L
        append!(var, Layer_X[i])
    end
    NumBVar=length(var);
    append!(var, Layer_X[1]);
    d=1;
    G=vec(G);
    NN=length(var)
    @polyvar x_Vars[1:NN]

    opt,sol,data = cs_tssos_first(G,y, d, TS=false,nb=NumBVar,solver="Mosek");
    return opt,sol,data
end

