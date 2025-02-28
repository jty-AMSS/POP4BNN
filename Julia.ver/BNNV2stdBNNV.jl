# This code can transform an L_inf BNN verifying problem to a Std BNN verifying problem, i.e.
# 
# min <f,x[end]>
# s.t. L<= x[0]<=U 
#      x is computed by BNN with wrights W[i] and bias b[i]

# to a Std BNNV problem 
# min <f,x[end]>
# s.t. -1<= x[0]<=1
#      x is computed by BNN with wrights W_new[i] and bias b_new[i]
using LinearAlgebra
using TSSOS
using DynamicPolynomials
using DelimitedFiles
using JuMP
using GLPK
using LinearAlgebra

include("POP4BNN.jl")

function runBNN(A,b,t)
    y=t;
    o=t
    y0=[]
    for i=1:length(A)
        y0=A[i]*y+b[i]
        y=sign.(y0)
        o=[o;y0];
    end
    return y,y0,o
end

function BNNV2stdBNNV_inf(A,b,L,U)
    k=(U-L)./2;
    c=(L+U)./2;
    A1=deepcopy(A)
    b1=deepcopy(b)
    b1[1]=A[1]*c+b1[1];
    A1[1]=A[1]*diagm(vec(k))

    Norm_A_inf=maximum(abs.(vec(A1[1])));
    A1[1]=1/Norm_A_inf *A1[1];
    b1[1]=1/Norm_A_inf *b1[1];


    # Need: adding remove not active nodes:
    A2,b2=ReduceStdBNN(A1,b1)
        # normalize:A1[1],b1[1]

    return A2,b2,L,U
end

function nv(A)
        At = abs.(A)
        o = ones(size(At, 2), 1)
        return At * o
end

function find(v)
    index=[]
    for i=1:length(v)
        if v[i]==1
            push!(index,i)
            end
        end
        return index
end
function  ReduceStdBNN(A1,b1)
# for BNN with weight and bias A,B, and input region \in [-1,1]^n
# generate a small BNN which reomve all of the node that are always  disactive or active

A=deepcopy(A1)
b=deepcopy(b1)
Index_pos=[];
Index_neg=[];
for i=1:(length(A)-1)
z=nv(A[i])
push!(Index_pos,[])
push!(Index_neg,[])
Index_pos[i]=find(b[i].>z);
Index_neg[i]=find(b[i].<-z);
num_pos=length(Index_pos[i]);
num_neg=length(Index_neg[i]);
b[i+1]= b[i+1]+ A[i+1][:,Index_pos[i]]*ones(num_pos,1);
b[i+1]= b[i+1]- A[i+1][:,Index_neg[i]]*ones(num_neg,1);

# A{i+1}(:,[Index_neg{i};Index_pos{i}])=[];
A[i+1]=A[i+1][:,setdiff(1:end, [Index_neg[i];Index_pos[i]])];

# A{i}([Index_neg{i};Index_pos{i}],:)=[];
A[i]=A[i][setdiff(1:end, [Index_neg[i];Index_pos[i]]),:];

# b{i}([Index_pos{i};Index_neg{i}],:)=[];
b[i]=b[i][setdiff(1:end, [Index_neg[i];Index_pos[i]]),:];
end

return A,b
end


function preprocessingBNNForLastLayer(A1,b1,f1)
    # function [A,b,f,c]=BNN_verify_preprocessing_last_layer(A,b,f)

    #     % input:
    #     % A,b: std BNN verificiation prob with input -1<=x<=1
    #     % f: obj function (linear)
        
    #     % output:
    #     % A,b: std BNN verificiation prob with input -1<=x<=1
    #     % new output obj function is <x_new,f_new>+c
    #     % with A{end}, b{end} be standardized
    A=deepcopy(A1)
    b=deepcopy(b1)
    f=deepcopy(f1)
        i=length(b)
        z=nv(A[i])
        Index_pos=find(b[i].>z);
        Index_neg=find(b[i].<-z)
       
        A[i]=A[i][setdiff(1:end, [Index_neg;Index_pos]),:];
        b[i]=b[i][setdiff(1:end, [Index_neg;Index_pos]),:];
        v=0*f;
        v[Index_neg].=-1;
        v[Index_pos].=1;
        c=dot(f,v);

        f=f[setdiff(1:end, [Index_neg;Index_pos])];
 
    return A,b,f,c
end



function two_largest_elements_with_indices(vec)
    sorted_indices = sortperm(vec, rev=true)
    largest_index_1 = sorted_indices[1]
    largest_index_2 = sorted_indices[2]

    largest_element_1 = vec[largest_index_1]
    largest_element_2 = vec[largest_index_2]

    return   largest_index_1,largest_index_2
end
 
function Linf_Adv_attack_prob_generate(A,b,x0,delta)
    # For give BNN (A,b), and sample x0, L_inf norm delta  generate the problem std BNNV problem
    # the input data was bounded by -1~1
    x0=vec(x0)
    L=x0.-delta
    L=max.(L,-1)
    U=x0.+delta
    U=min.(U,1)
    n=[length(x0)]
    y,y0,o=runBNN(A,b,x0)
    largest_index_1,largest_index_2=two_largest_elements_with_indices(vec(y0))
    for i=1:(length(A)-1)
        n=[n;size(A[i],1)]
    end
    f0=A[end][largest_index_1,:]-A[end][largest_index_2,:]
    Rob_bound=b[end][largest_index_1]-b[end][largest_index_2]
    A=A[1:end-1];
    b=b[1:end-1];
    A1,b1=BNNV2stdBNNV_inf(A,b,L,U)
    # print size of a

    # t=x0.-((L.+U)./2).*(2)./(U.-L)
    
    A3,b3,f,c=preprocessingBNNForLastLayer(A1,b1,f0)
    Apop=A3
    bpop=b3
    println("Size of BNNs:")
    for i=1:length(A3)
        println(size(A3[i],2))
    end
    println(size(A3[end],1))
    println("The above list is the numbers of nodes of each of leayers A1")
    println("Rbound:")
    println(Rob_bound+c)

    opt,sol,data=POP4BNNstd(Apop, bpop, f) 
    opt=opt+Rob_bound+c
    return opt,sol,data,A3,b3,f,c,Rob_bound,f0,L,U
end

