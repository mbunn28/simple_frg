function magnus1!(RHS, V, dispersion, span=(10.0, 1E-6))
    #=
    # A Magnus based Integrator.
    # To first order the solution is y(h) = exp(h*A(y0,t0)) y0
    # Since we don't have A and don't want to calculate matrix exponentials we use
    # that A(y0,t0) * y0 = rhs(y0,t0). Generalizing this we define the function
    # exprhs(y0, yr) = A(y0,t0) yr
    =#

    upper_step = 0.08
    lower_step = 1E-5
    nk = size(V)[1]
    dV1 = StructArray(zeros(ComplexF64,size(V)))
    dV2 = StructArray(zeros(ComplexF64,size(V)))
    Lpp = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    Lph = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    lam = span[1]
    maxV = 1.0
    while lam > span[2]
        @show lam
        @show maximum(abs.(V))

        dlam = lam / maxV

        if dlam > upper_step * lam
            dlam = upper_step * lam
        elseif dlam < lower_step * lam
            dlam = lower_step * lam
        end
        # first calculate the Euler correction to the solution.
        RHS(dV1, V, dispersion, lam, Lpp, Lph)
        #next we apply the exprhs function once to the such obtained update vector from the rhs vector
        expRHS(dV2, V, dV1, dispersion, lam, Lpp, Lph)
        V .-= dlam .* dV1
        V .+= 0.5*dlam*dlam .* dV2

        # step size adaption, for now constant = f(l)
        lam -= dlam
        
        maxV = maximum(abs.(V))
        
        if maxV > 50
            break
        end
    end
end


function magnus2!(RHS, V, dispersion, span=(10.0, 1E-6))
    #=
    # A second order Magnus based Integrator. This one requires more data...
    =#
    
    upper_step = 0.1
    lower_step = 1E-5
    maxV = 1.0
    nk = size(V)[1]
    dV11 = StructArray(zeros(ComplexF64,size(V)))
    dV12 = StructArray(zeros(ComplexF64,size(V)))
    dV13 = StructArray(zeros(ComplexF64,size(V)))
    Vm = StructArray(zeros(ComplexF64,size(V)))
    dV21 = StructArray(zeros(ComplexF64,size(V)))
    dV22 = StructArray(zeros(ComplexF64,size(V)))
    dV23 = StructArray(zeros(ComplexF64,size(V)))    
    Lpp = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    Lph = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    lam = span[1]
    while lam > span[2]
        @show lam
        @show maximum(abs.(V))

        dlam = lam / maxV

        if dlam > upper_step * lam
            dlam = upper_step * lam
        elseif dlam < lower_step * lam
            dlam = lower_step * lam
        end
        # first calculate the Euler correction to the solution.
        RHS(dV11, V, dispersion, lam, Lpp, Lph)
        expRHS(dV12, V, dV11, dispersion, lam, Lpp, Lph)
        expRHS(dV13, V, dV12, dispersion, lam, Lpp, Lph)
        dlambda = dlam
        # approximate intermediate exponential
        Vm = V - dlambda/2 * dV11 + 0.5*(dlambda/2)*(dlambda/2)*dV12 - 1/6*(dlambda/2)*(dlambda/2)*(dlambda/2)*dV13
        expRHS(dV21, Vm, V, dispersion, lam + dlambda/2, Lpp, Lph)
        expRHS(dV22, Vm, dV21, dispersion, lam + dlambda/2, Lpp, Lph)
        expRHS(dV23, Vm, dV22, dispersion, lam + dlambda/2, Lpp, Lph)

        V .-= dlambda .* dV21
        V .+= 0.5*dlambda*dlambda .* dV22
        V .-= 1/6*dlambda*dlambda*dlambda .* dV23

        # step size adaption, for now constant = f(l)
        lam -= dlam
        maxV = maximum(abs.(V))

        if maxV > 50
            break
        end
    end
end


function fixed_step_integrator!(RHS, V, dispersion, span=(10.0, 4.0))
    #=
    # A fixed step-width integrator to check against the reference FRG results
    =#
    nk = size(V)[1]
    dV = StructArray(zeros(ComplexF64,size(V)))
    Lpp = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    Lph = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    lam = span[1]
    while lam > 1.1
        RHS(dV, V, dispersion, lam,Lpp,Lph)
        V .-= 0.1 .* dV
        @show lam
        lam -= 0.1
        @show maximum(abs.(V))
    end
end


function adaptive_euler!(RHS, V, dispersion, span=(50.0, 1E-6))
    #=
    # The basic adaptive euler integrator used in FRG calculations
    =#
    upper_step = 0.03
    lower_step = 1E-5
    nk = size(V)[1]
    dV = StructArray(zeros(ComplexF64, size(V)))
    Lpp = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    Lph = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    lam = span[1]
    maxV = 1.0
    while lam > span[2]
        #@show lam
        #@show maximum(abs.(V))

        dlam = lam / maxV

        if dlam > upper_step * lam
            dlam = upper_step * lam
        elseif dlam < lower_step * lam
            dlam = lower_step * lam
        end
        RHS(dV, V, dispersion, lam, Lpp, Lph)

        V .-= dlam .* dV
        lam -= dlam

        maxV = maximum(abs.(V))

        if maxV > 50
            break
        end
    end
    print(lam, "\n")
return lam

Base.@propagate_inbounds function Base.setindex!(s::StructArray{ComplexF64, <:Any, <:Any, Int}, vals::Float64, I::Int)
    @boundscheck checkbounds(s, I)
    StructArrays.foreachfield((col, val) -> (@inbounds col[I] = val), s, ComplexF64(vals))
    s
end

function library_integrator!(RHS, V, dispersion, integrator, span=(50.0, 1E-6))
    #=
    # Interface for DifferentialEquations.jl
    =#
    nk = size(V)[1]
    Lpp = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    Lph = StructArray(zeros(ComplexF64, nk, nk, nk, nk))
    divergence=50

    function condition(u,t,integrator)
        @show t
        return maximum(abs.(u)) > divergence
    end

    cb = DiscreteCallback(condition, terminate!)
    prob = ODEProblem((du, u,p,t)->RHS(du, u, dispersion, t, Lpp, Lph), V, span)

    sol = solve(prob, integrator, callback=cb, save_everystep=false, dtmin=1E-5)

    V .= sol.u[end]
end
