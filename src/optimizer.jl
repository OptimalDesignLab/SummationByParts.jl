module Optimizer
# A module with optimization algorithms

using LinearAlgebra
using Random
using ..SymCubatures, ..OrthoPoly
#using IncompleteLU, SparseArrays, IterativeSolvers, Preconditioners

export pso, levenberg_marquardt #, rosenbrock, rastrigin

"""
### SummationByParts.pso

Particle Swarm Optimization (PSO) Algorithm

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
    xL::T=zero(T), xR::T=T(2), delta1::T=T(1e-2), delta2::T=T(1e-2),
    maxiter::Int64=1000, tol=10*eps(typeof(real(one(T)))), save_iter=false, verbose=0) where {T}
    
    @assert(delta1>zero(T) && delta1<one(T), "delta1 must be in (0,1).")
    @assert(delta2>zero(T) && delta2<one(T), "delta2 must be in (0,1).")

    delta1_in = copy(delta1)
    delta2_in = copy(delta2)

    # set parameters
    ne = ne                 # number of elements in each particles
    np = np                 # number of particles
    vmax = T(0.2)*(xR-xL)      # maximum velocity allowed (20% of space range) 
    c1 = T(1.5)                # confidence in individual position
    c2 = T(1.5)                # confidence in group position
    w = T(0.6)                 # inertia weight
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
            eps = T(1e-6)
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
                    delta1 = T(0.2)
                end
                cntr_perturb = 0
            end

            if fgb > one(T)
                maxval /= abs(fgb)
            end
            if (fgb>T(1e-1) && maxval<=T(1e-4))
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
### SummationByParts.levenberg_marquardt

Levenberg-Marquardt Algorithm (LMA) 

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
    xL::T=zero(T), xR::T=T(3), nu::T=T(1000), maxiter::Int64=1000, tol=10*eps(typeof(real(one(T)))), verbose=0) where{T}

    Jac = SymCubatures.calcjacobian(cub)

    # compute accuracy for initial guess 
    F, dF = fun(cub, q, compute_grad=true)
    res = norm(F)
    res_lst = T[]
    res_old = res
    verbose==1 ? print("solvecubature!:\n") : nothing
    verbose==1 ? print("res norm = ",res,"\n") : nothing
    if (res < tol)
        v = zeros(T, (cub.numparams + cub.numweights))
        v[1:cub.numparams] = cub.params
        v[cub.numparams+1:end] = cub.weights
        return res, v, 0, res_lst 
    end

    # initialize the parameter vector
    v = zeros(T, (cub.numparams + cub.numweights))
    v[1:cub.numparams] = cub.params
    v[cub.numparams+1:end] = cub.weights
    if xinit!=[]
        v = xinit
    end

    v_old = copy(v)
    alpha = one(T)
    iter_break = 1000
    dv = zeros(size(v))
    iter = 0
    for k = 1:maxiter
        J = dF*Jac
        JtJ = J'*J
        # nu = 10*(1/norm(diag(JtJ)))*(norm(F)^2)
        H = JtJ + nu*diagm(diag(JtJ))
        g = -J'*F

        # solve only for those parameters and weights that are in mask
        fill!(dv, zero(T))
        Hred = H[mask,mask]
        # dv[mask]=pinv(Hred,1e-14)*g[mask]
        xx = dv[mask]
        bb = g[mask]
        if !(any(isnan, Hred) || any(isinf, Hred) || any(isnan, bb) || any(isinf, bb))
            try
                xx = Hred\bb
            catch
                lambda = sqrt(tol)
                xx = (Hred + lambda * I)\bb
            end
        end
        dv[mask] = xx

        # update cubature definition and check for convergence
        v += dv

        eps = T(1e-10)
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

        if typeof(cub)==TriSymCub{T}
            vg = copy(v[1:cub.numparams])
            ng = 0
            vS21 = vg[1:cub.numS21]
            for i = 1:length(vS21)                
                if (v[ng+i] >= one(T))
                    alphaS21 = copy(alpha)
                    while (v[ng+i] >= one(T) && alphaS21>0)
                        v -= alphaS21*dv
                        alphaS21 = 0.5*alphaS21
                        v += alphaS21*dv
                    end
                end
            end
            ng+=cub.numS21 

            vS111 = vg[ng+1:ng+cub.numS111]
            for i = 1:2:length(vS111)
                if ((v[ng+i]+v[ng+i+1])/2 >= one(T))
                    alphaS111 = copy(alpha)
                    while ((v[ng+i]+v[ng+i+1])/2 >= one(T) && alphaS111>0)
                        v -= alphaS111*dv
                        alphaS111 = 0.5*alphaS111 
                        v += alphaS111*dv
                    end
                end
            end
            ng+=cub.numS111  
        end

        if typeof(cub)==TetSymCub{T}
            vg = copy(v[1:cub.numparams])
            ng = 0
            vS31 = vg[1:cub.numS31]
            for i = 1:length(vS31)                
                if (v[ng+i] >= one(T))
                    alphaS31 = copy(alpha)
                    while (v[ng+i] >= one(T) && alphaS31>0)
                        v -= alphaS31*dv
                        alphaS31 = 0.5*alphaS31
                        v += alphaS31*dv
                    end
                end
            end
            ng+=cub.numS31 

            vS22 = vg[ng+1:ng+cub.numS22]
            for i = 1:length(vS22)
                if (v[ng+i] >= one(T))
                    alphaS22 = copy(alpha)
                    while (v[ng+i] >= one(T) && alphaS22>0)
                        v -= alphaS22*dv
                        alphaS22 = 0.5*alphaS22 
                        v += alphaS22*dv
                    end
                end
            end
            ng+=cub.numS22 
            ng+=cub.numfaceS21 
            ng+=cub.numedge

            vS211 = vg[ng+1:ng+2*cub.numS211]
            for i = 1:2:length(vS211)
                if ((v[ng+i]+v[ng+i+1]/2)>= one(T))
                    alphaS211 = copy(alpha)
                    while ((v[ng+i]+v[ng+i+1]/2)>= one(T) && alphaS211>0)
                        v -= alphaS211*dv
                        alphaS211 = 0.5*alphaS211
                        v += alphaS211*dv
                    end
                end
            end
            ng+=2*cub.numS211 
            ng+=2*cub.numfaceS111

            vS1111 = vg[ng+1:ng+3*cub.numS1111] 
            for i = 1:3:length(vS1111)
                if ((v[ng+i]+v[ng+i+1]+v[ng+i+2])/2 >= one(T))
                    alphaS1111 = copy(alpha)
                    while ((v[ng+i]+v[ng+i+1]+v[ng+i+2])/2 >= one(T) && alphaS1111>0)
                        v -= alphaS1111*dv
                        alphaS1111 = 0.5*alphaS1111
                        v += alphaS1111*dv
                    end
                end
            end     
        end

        SymCubatures.setparams!(cub, v[1:cub.numparams])
        SymCubatures.setweights!(cub, v[cub.numparams+1:end])
        F, dF = fun(cub, q, compute_grad=true)
        res = norm(F)
        
        iter_show = 1
        if mod(k,iter_show)==0
            push!(res_lst,res)
        end

        if verbose==1
            if mod(k,iter_show)==0 
                print("\titer ",k,": res norm = ",res,"\n")
            end
        end

        if res < tol 
            return res, v, k, res_lst
        end

        if (res > T(1e3) || isnan(res))
            return res_old, v_old, 0, res_lst
        end

        # trust-region like update
        if res > res_old
            v -= dv
            SymCubatures.setparams!(cub, v[1:cub.numparams])
            SymCubatures.setweights!(cub, v[cub.numparams+1:end])
            F, dF = fun(cub, q, compute_grad=true)
            nu *= 5.0
        else
            nu /= 5.0
            res_old = res
        end
        alpha=one(T)
        iter+=1

        if mod(k,iter_break)==0
            res_break = res_lst[end-2:end] .- res
            if maximum(abs.(res_break)) < tol && norm(res) < tol
                return res, v, iter, res_lst
            end
            iter_break += 50
        end
        if k==maxiter && norm(res) > tol
            println("No solution found after a $k LMA iterations.")
        end
    end
    return res, v, iter
end

end
