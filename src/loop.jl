function calc_Ldot_sharp(e0::Float64, e1::Float64, lam::Float64, type=-1.)
    return -type/(2.0*pi) *
        ( 1/(( type*im*lam - e0) * ( im*lam-e1))
        + 1/((-type*im*lam - e0) * (-im*lam-e1)) )
end

function calc_Ldot(e0::Float64, e1::Float64, lam::Float64, type =-1.) :: ComplexF64
    if abs(e0-type*e1)<1e-12
        return - (3*e0^2 - 4*abs(e0)*lam - lam^2) / (4* (abs(e0) + lam)^4)
    else
        return 1/(4*(e0 -type * e1))*
        (
                    e0 * (3*abs(e0) + lam) / (abs(e0) + lam)^3
            -type * e1 * (3*abs(e1) + lam) / (abs(e1) + lam)^3
        )
    end
end

function flip_fine(i, nk_fine, type)
    #=
    # This function should reverse the index (respective the midpoint) depending
    # on the channel (+i for particle-hole loop, -i for particle-particle loop)
    =#
    return type * i + (1-type)*ceil(nk_fine/2.)
end

function calc_L(L, lam, dispersion, type=-1., cutoff='o')::Nothing

    nk = size(dispersion)[1]
    nk_fine = size(dispersion)[3]

    L .= 0
    @sync begin
    for llx in 1:nk
    for lly in 1:nk
        for llpx in 1:nk
        for llpy in 1:nk
            Threads.@spawn begin
            for llx_fine in 1:nk_fine
            for lly_fine in 1:nk_fine

                lpx_fine::Int64 = flip_fine(llx_fine, nk_fine, type)
                lpy_fine::Int64 = flip_fine(lly_fine, nk_fine, type)
                L[llx, lly, llpx, llpy] +=
                    calc_Ldot(dispersion[llx,  lly,  llx_fine, lly_fine],
                                    dispersion[llpx, llpy, lpx_fine, lpy_fine],
                                    lam, type)
            end end
            end
        end end
    end end
    end

    L .*= (1.0 / nk_fine)^2

    return nothing

end
