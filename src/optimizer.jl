module Optimizer
# A module with optimization algorithms

using LinearAlgebra
using Random
using ..SymCubatures, ..OrthoPoly
#using IncompleteLU, SparseArrays, IterativeSolvers, Preconditioners

export pso, levenberg_marquardt, rosenbrock, rastrigin

"""
### SummationByParts.pso

Particle Swarm Optimization (PSO) for unconstrained and constrained problems

**Inputs**
* `fun`: function to be optimized
* `n`: number of dimensions (number of parameters to be found)
* `xinit`: initial guess (for random initial guess)
* `xL`: lower bound on the parameters 
* `xR`: upper bound on the parameters 
* `maxiter`: maximum number of iterations (default is 1000)
* `tol`: tolerance to stop iteration (default is 1e-14)

**Outputs**
* `fmin`: the optimized function value
* `xmin`: the minimizer (solution)
* `f_all`: all function evaluations
"""

function pso(fun::Function, ne; np::Int64=10, xinit=[], 
    xL::Float64=0.0, xR::Float64=1.0, maxiter::Int64=1000, tol::Float64=1e-4, save_iter=false, verbose=0)

    # set parameters
    ne = ne                 # number of elements in each particles
    np = np                 # number of particles
    vmax = 0.2*(xR-xL)      # maximum velocity allowed (20% of space range) 
    c1 = 1.5                # confidence in individual position
    c2 = 1.5                # confidence in group position
    w = 0.5                 # inertia weight

    # initialize the particles
    x = (xR-xL).* rand(np, ne) .+ xL
    if xinit!=[]
        x[1,:] = xinit
    end

    # initialize vectors
    v = zeros(np,ne)        # velocity vector
    fp = zeros(np,1)        # personal function values
    xpb = zeros(np,ne)      # personal best solution parameters

    for i = 1:np
        fp[i,1] = fun(x[i,:])
        xpb[i,:] = x[i,:]
    end
    fpb = copy(fp)          # personal best function values
    fgb = minimum(fpb)      # global best function value
    indx = argmin(fpb)[1]   # index of global best function value
    xgb = x[indx, :]        # global best solution parameters
    xp = zeros(np,ne)
    fgb1 = copy(fgb)
    cntr_perturb = 0

    # main PSO routine
    f_all=[]
    x_all=[]
    if save_iter
        push!(f_all,fun(xgb))
        push!(x_all,xgb)
    end
    iter = 1
    while (iter <= maxiter && fgb>=tol)
    # while (iter <= maxiter)
        for i = 1:np       
            # calculate velocity
            v[i,:] = w.*v[i,:] + c1.*rand(ne).*(xpb[i,:] - x[i,:]) + c2.*rand(ne).*(xgb - x[i,:])

            # limit the velocity
            for j = 1:ne
                if abs(v[i,j]) > vmax
                    v[i,j] = v[i,j].*vmax/abs(v[i,j])
                end
            end
             
            # update position
            xp[i,:] = x[i,:] + v[i,:]
            
            # enforce left and right bounds
            for j = 1:ne
                if (xp[i,j]>xR)
                    xp[i,j] = (xR-xL)*rand() + xL
                elseif (xp[i,j]<xL)
                    xp[i,j] = (xR-xL)*rand() + xL
                end
            end 
            
            # evaluate objective function for the new particle position
            fp[i] = fun(xp[i,:])
            
            # check personal best function and position
            if fp[i] < fpb[i]
                fpb[i] = fp[i]
                xpb[i,:] = xp[i,:]
            end
            
            # check global best objective function evaluation and position
            if fp[i] < fgb 
                fgb = fp[i] 
                xgb = xp[i,:]
            end   
            x[i,:] = xp[i,:]
        end
        
        if save_iter
            f_new = fgb; 
            if f_new < f_all[iter]
                push!(f_all,f_new)
                push!(x_all,xgb)
            else
                push!(f_all,f_all[iter])
                push!(x_all,x_all[iter])
            end
        end   
        
        # check if there is no change in the global best, if so perturb the current solution
        if mod(iter,round(Int,ne*2/(log10(ne)+1)))==0
            maxval = abs((fgb - fgb1))
            m = 1.0/1.0
            if fgb > 1.0
                maxval /= abs(fgb)
            end
            if (fgb>1e-6 && maxval<=1e-6)
                if abs(fgb) > 1.0
                    d = ((fgb^2-1.0)/(fgb^2+1.0))*(xR-xL)/m
                else
                    d =fgb*(xR-xL)/m
                end
                x += d .* rand(np, ne)
                for i = 1:np
                    fp[i,1] = fun(x[i,:])
                    xpb[i,:] = x[i,:]
                end
                fpb = copy(fp)
                fgb = minimum(fpb)
            elseif (fgb<=1e-6 && maxval<=1e-10)
                x += (fgb*(xR-xL)/m) .*rand(np,ne)
                for i = 1:np
                    fp[i,1] = fun(x[i,:])
                    xpb[i,:] = x[i,:]
                end
                fpb = copy(fp)
                fgb = minimum(fpb)
            end
            fgb1 = copy(fgb)
            cntr_perturb += 1
        end
        
        iter_show = 1000
        if verbose==1
            if mod(iter,iter_show)==0
                println(fgb)
            end
        elseif verbose==2
            if mod(iter,iter_show)==0
                println(fgb)
                println(xgb)
            end
        end

        iter = iter +1      
    end

    fmin = fgb 
    xmin = xgb
    return fmin, xmin, f_all, x_all, iter-1, cntr_perturb

end

"""
### SummationByParts.pso

Particle Swarm Optimization (PSO) for unconstrained and constrained problems

**Inputs**
* `fun`: function to be optimized
* `ne`: number of parameters to be found
* `cub`: cubature data
* `q`: the degree of the quadrature rule
* `mask`: a vector of the index of the parameters to be determined
* `np`: number of particles 
* `xinit`: initial parameter guess
* `xL`: lower bound on the parameters 
* `xR`: upper bound on the parameters 
* `delta1`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin <= 0.1
* `delta2`: purturbation constant that determines how much the parameters are perturbed
            in the case of stagnation with fmin > 0.1
* `maxiter`: maximum number of iterations (default is 1000)
* `tol`: tolerance to stop iteration (default is 1e-14)
* `save_iter`: boolean to save results of each iteration
* `verbose`: boolean to print results of iteration periodically

**Outputs**
* `fmin`: the optimized function value
* `xmin`: the minimizer (solution)
* `f_all`: all function evaluations
"""
function pso(fun::Function, ne::Int, cub::SymCub{T}, q::Int, mask::AbstractArray{Int64,1}; np::Int64=10, xinit::Array{T}=[], 
    xL::T=convert(T, 0.0), xR::T=convert(T, 2.0), delta1::Float64=1e-2, delta2::Float64=1e-2,
    maxiter::Int64=1000, tol=10*eps(typeof(real(one(T)))), save_iter=false, verbose=0) where {T}
    
    @assert(delta1>0.0 && delta1<1.0 && isa(delta1, Float64), "delta1 must be a positive Float64 in (0,1).")
    @assert(delta2>0.0 && delta2<1.0 && isa(delta2, Float64), "delta2 must be a positive Float64 in (0,1).")

    delta1_in = copy(delta1)
    delta2_in = copy(delta2)

    # set parameters
    ne = ne                 # number of elements in each particles
    np = np                 # number of particles
    vmax = 0.2*(xR-xL)      # maximum velocity allowed (20% of space range) 
    c1 = 1.5                # confidence in individual position
    c2 = 1.5                # confidence in group position
    w = 0.6                 # inertia weight
    compute_grad = false

    # initialize the particles
    x = zeros(T,np,ne)
    xrand = rand(T,np, ne)
    for i = 1:np
        x[i,1:cub.numparams] = cub.params
        x[i,cub.numparams+1:end] = cub.weights
        x[i,mask] = (xR-xL).* xrand[i,mask] .+ xL
    end

    if xinit!=[]
        for i = 1:np
            x[i,:] = xinit
        end
    end

    # initialize vectors
    v = zeros(T,np,ne)        # velocity vector
    fp = zeros(T,np,1)        # personal function values
    xpb = zeros(T,np,ne)      # personal best solution parameters

    for i = 1:np
        SymCubatures.setparams!(cub, x[i,:][1:cub.numparams])
        SymCubatures.setweights!(cub, x[i,:][cub.numparams+1:end])
        fp[i,1] = norm(fun(cub, q, compute_grad=compute_grad)[1])
        xpb[i,:] = x[i,:]
    end
    fpb = copy(fp)          # personal best function values
    fgb = minimum(fpb)      # global best function value
    indx = argmin(fpb)[1]   # index of global best function value
    xgb = x[indx, :]        # global best solution parameters
    xp = copy(x) #zeros(np,ne)
    fgb1 = copy(fgb)
    cntr_perturb = 0
    nperturb = 0
    
    # main PSO routine
    f_all=[]
    x_all=[]
    if save_iter
        SymCubatures.setparams!(cub, xgb[1:cub.numparams])
        SymCubatures.setweights!(cub, xgb[cub.numparams+1:end])
        fxgb = norm(fun(cub, q, compute_grad=compute_grad)[1])
        push!(f_all,fxgb)
        push!(x_all,xgb)
    end
    iter = 1
    while (iter <= maxiter && fgb>=tol)
        for i = 1:np     
            # calculate velocity
            v[i,:] = w.*v[i,:] + c1.*rand(ne).*(xpb[i,:] - x[i,:]) + c2.*rand(ne).*(xgb - x[i,:])

            # limit the velocity
            for j = 1:ne
                if abs(v[i,j]) > vmax
                    v[i,j] = v[i,j].*vmax/abs(v[i,j])
                end
            end
            if norm(v)<=1e-14
                v .= v .+1e-3
            end 
            # update position
            xp[i,mask] = x[i,mask] + v[i,mask]
            
            # enforce left and right bounds
            eps = 1e-6
            for j = 1:ne
                if (xp[i,j]>xR)
                    xp[i,j] = xR-(100*eps*rand()+eps) #(xR-xL)*rand() + xL
                elseif (xp[i,j]<xL)
                    xp[i,j] = 100*eps*rand()+eps #(xR-xL)*rand() + xL
                end
            end 

            # evaluate objective function for the new particle position
            SymCubatures.setparams!(cub, xp[i,:][1:cub.numparams])
            SymCubatures.setweights!(cub, xp[i,:][cub.numparams+1:end])
            fxp=norm(fun(cub, q, compute_grad=compute_grad)[1])
            fp[i] = fxp
            
            # check personal best function and position
            if fp[i] < fpb[i]
                fpb[i] = fp[i]
                xpb[i,:] = xp[i,:]
            end
            
            # check global best objective function evaluation and position
            if fp[i] < fgb 
                fgb = fp[i] 
                xgb = xp[i,:]
            end   
            x[i,:] = xp[i,:]
        end
        
        if save_iter
            f_new = fgb; 
            if f_new < f_all[iter]
                push!(f_all,f_new)
                push!(x_all,xgb)
            else
                push!(f_all,f_all[iter])
                push!(x_all,x_all[iter])
            end
        end   
       
        # check if there is no change in the global best, if so perturb the current solution
        if mod(iter,round(Int,4*ne/(log10(ne)+1)))==0
            maxval = abs((fgb - fgb1))
            
            delta1 = delta1_in
            delta2 = delta2_in 

            # if no improvement is shown after several purterbations, then use large perturbations
            if mod(cntr_perturb,5) == 0
                if fgb>1e-1
                    delta2 = 0.5
                else
                    delta1 = 0.2
                end
                cntr_perturb = 0
            end

            if fgb > 1.0
                maxval /= abs(fgb)
            end
            if (fgb>1e-1 && maxval<=1e-4)
                x[:,mask] = (1-delta2) .*x[:,mask] + delta2 .*rand(np, ne)[:,mask]
                for i = 1:np
                    SymCubatures.setparams!(cub, x[i,:][1:cub.numparams])
                    SymCubatures.setweights!(cub, x[i,:][cub.numparams+1:end])
                    fx=norm(fun(cub, q, compute_grad=compute_grad)[1])
                    fp[i,1] = fx
                    xpb[i,:] = x[i,:]
                end
                fpb = copy(fp)
                fgb = minimum(fpb)
                cntr_perturb += 1
                nperturb += 1
            elseif (fgb<=1e-1 && maxval<=1e-8)
                # x[:, mask] += (fgb*(xR-xL)/delta1) .*rand(np,ne)[:,mask]
                x[:,mask] = (1-delta1) .*x[:,mask] + delta1 .*rand(np, ne)[:,mask]
                for i = 1:np
                    SymCubatures.setparams!(cub, x[i,:][1:cub.numparams])
                    SymCubatures.setweights!(cub, x[i,:][cub.numparams+1:end])
                    fx=norm(fun(cub, q, compute_grad=compute_grad)[1])
                    fp[i,1] = fx
                    xpb[i,:] = x[i,:]
                end
                fpb = copy(fp)
                fgb = minimum(fpb)
                cntr_perturb += 1
                nperturb += 1
            end
            
            if fgb1<5e-15 && fgb <5e-15
                break
            end
            fgb1 = copy(fgb)
        end
        
        iter_show = 100
        if (verbose==1 || verbose==true)
            if mod(iter,iter_show)==0
                println(fgb)
            end
        elseif (verbose==2 || verbose==true)
            if mod(iter,iter_show)==0
                println(fgb)
                println(xgb)
            end
        end

        iter = iter +1      
    end
    fmin = fgb 
    xmin = xgb
    return fmin, xmin, f_all, x_all, iter-1, nperturb

end

"""
### SummationByParts.pso

Levenberg-Marquardt Algorithm (LMA) optimizer

**Inputs**
* `fun`: function to be optimized
* `cub`: cubature data
* `q`: the degree of the quadrature rule
* `mask`: a vector of the index of the parameters to be determined
* `np`: number of particles 
* `xinit`: initial parameter guess
* `xL`: lower bound on the parameters 
* `xR`: upper bound on the parameters 
* `nu`: parameter controling exploration (spliting between Newton's and steepest decent methods)
* `maxiter`: maximum number of iterations 
* `tol`: tolerance to stop iteration 
* `verbose`: boolean to print results of iteration periodically

**Outputs**
* `fmin`: the optimized function value
* `v`: the minimizer (solution)
* `iter`: number of itrations
"""
function levenberg_marquardt(fun::Function, cub::SymCub{T}, q::Int64, mask::AbstractArray{Int64,1}; xinit::Array{T}=[], 
    xL::T=convert(T, 0.0), xR::T=convert(T, 2.0), nu=1000.0, maxiter::Int64=1000, tol=10*eps(typeof(real(one(T)))), verbose=0) where{T}

    Jac = SymCubatures.calcjacobian(cub)

    # compute accuracy for initial guess 
    F, dF = fun(cub, q, compute_grad=true)
    res = norm(F)
    res_old = res
    verbose==1 ? print("solvecubature!:\n") : nothing
    verbose==1 ? print("res norm = ",res,"\n") : nothing
    # verbose==1 ? print("\titer ",0,": res norm = ",res,"\n") : nothing
    if (res < tol)
        v = zeros(T, (cub.numparams + cub.numweights))
        v[1:cub.numparams] = cub.params
        v[cub.numparams+1:end] = cub.weights
        return res, v, 0
    end

    # initialize the parameter vector
    v = zeros(T, (cub.numparams + cub.numweights))
    v[1:cub.numparams] = cub.params
    v[cub.numparams+1:end] = cub.weights
    if xinit!=[]
        v = xinit
    end

    alpha = 1.0
    res_lst = []
    iter_break = 1000
    dv = zeros(size(v))
    # P = Optimizer.preconditioner(Matrix(dF'))
    iter = 0
    for k = 1:maxiter
        J = dF*Jac
        JtJ = J'*J
        H = JtJ + nu*diagm(diag(JtJ))
        g = -J'*F

        # solve only for those parameters and weights that are in mask
        fill!(dv, zero(T))
        Hred = H[mask,mask]
        # dv[mask] = Hred\g[mask]
        dv[mask] = pinv(Hred,1e-14)*g[mask]
        # U,S,Vt=svd(Hred,full=true,alg=LinearAlgebra.QRIteration())
        # invS = zeros(length(S),length(S))
        # for i=1:length(axes(S,1))
        #     if abs(S[i])>1e-10
        #         invS[i,i] = 1.0/S[i]
        #     end
        # end 
        # pinvHred = Vt*invS*U'
        # dv[mask] = pinvHred*g[mask]

        # update cubature definition and check for convergence
        v += dv

        SymCubatures.setparams!(cub, v[1:cub.numparams])
        SymCubatures.setweights!(cub, v[cub.numparams+1:end])
        F, dF = fun(cub, q, compute_grad=true)
        res = norm(F)
        # print(res)
        # verbose==1 ? print("\titer ",k,": res norm = ",res,"\n") : nothing

        iter_show = 100
        if verbose==1
            if mod(k,iter_show)==0
                print("\titer ",k,": res norm = ",res,"\n")
            end
        end

        if res < tol
            return res, v, k
        end

        # trust-region like update
        if res > res_old
            v -= dv
            SymCubatures.setparams!(cub, v[1:cub.numparams])
            SymCubatures.setweights!(cub, v[cub.numparams+1:end])
            F, dF = fun(cub, q, compute_grad=true)
            nu *= 2.0
        else
            nu /= 2.0
            res_old = res
        end

        eps = 1e-4 
        for i=1:length(axes(v,1))
            if v[i]<xL
                v -= alpha*dv
                alpha = (eps - v[i])/(dv[i])
                v += alpha*dv
            elseif v[i]>xR
                v -= alpha*dv
                alpha = ((xR-eps) - v[i])/(dv[i])
                v += alpha*dv
            end
        end
        alpha=1.0
        iter+=1

        if mod(k,iter_show)==0
            push!(res_lst,res)
        end
        if mod(k,iter_break)==0
            res_break = res_lst[end-2:end] .- res
            if maximum(abs.(res_break)) < 1e-12
                break
            end
            iter_break += 1000
        end

    end
    # error("solvecubature failed to find solution in ",maxiter," iterations")
    return res, v, iter
end

function rosenbrock(x::Array{Float64}; xL::Float64=-Inf, xR::Float64=Inf)
    n = length(x)
    for i=1:n
        if (x[i]<xL || x[i]>xR)
            ErrorException("Some elements of vector x are not within the interval [xL, xR]")
        end
    end

    f = 0;
    for i=1:n-1
       f = f + 100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2 
    end

    gradf = zeros(n,1)
    gradf[1,1]= -400*x[1]*(x[2] - x[1]^2) - 2*(1-x[1])
    for i = 2:n-1
        gradf[i,1] = 200*(x[i]-x[i-1]^2)-400*x[i]*(x[i+1]-x[i]^2)-2*(1-x[i]);
    end    
    gradf[n,1] = 200*(x[n]-x[n-1]^2); 

    return f
end

function rastrigin(x::Array{Float64}; xL::Float64=-Inf, xR::Float64=Inf)    
    n = length(x)
    for i=1:n
        if (x[i]<xL || x[i]>xR)
            ErrorException("Some elements of vector x are not within the interval [xL, xR]")
        end
    end
     
    f = 0; 
    for i = 1:n
        f = f + (x[i]^2 - 10*cos(2*π*x[i]))
    end
    f += 10*n
          
    #gradient of the objective function 
    gradf = zeros(n,1)
    for i = 1:n
        gradf[i,1] = 2*x[i]+20*π*sin(2*π*x[i])
    end    
              
    return f
end

function preconditioner(V::Array{Float64}; ax::Int=1)
    @assert(ax>=1 && ax <=2)
    P = zeros(size(V)[ax],size(V)[ax])
    for i=1:length(axes(V,ax))
        if ax==1
            P[i,i] = 1.0/maximum(abs.(V[i,:]))
        elseif ax==2
            P[i,i] = 1.0/maximum(abs.(V[:,i]))
        end
    end
    return P
end

function preconditioner_saraswat(V::Array{Float64})
    m = size(V,1)
    A = zeros(m,1)
    B = zeros(m,1)
    for i = 1:m
        A[i] = norm(V[i,:])
        B[i] = norm(pinv(V,1e-2)'[i,:])
    end
    P = zeros(m,m)
    for i = 1:m
        P[i,i] = sqrt(B[i]/A[i])
    end

    return P
end

function preconditioner_shannon(V::Array{Float64})
    F = svd(V)
    P = F.V*inv(diagm(F.S))
    return P
end

# function preconditioner(P::Array{Float64})
#     for k=1:length(axes(P,1))
#         indx = argmax(abs.(P))
#         i = indx[1]
#         j = indx[2]

#         n = size(P, 1)
#         for r=1:n
#             if r != i
#                 ratio = P[r, j] / P[i, j]
#                 P[r, :] -= ratio * P[i, :]
#             end
#         end
#         P[i,:] ./= P[i,j]
#     end
#     return P
# end

end