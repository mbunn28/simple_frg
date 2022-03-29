struct hopping_parameters
    t::Float64
    tp::Float64
    mu::Float64
end

rhsevals = 0

function parse_parameters()
#=
# The main can now take command line arguments
# such as tp, mu. All other arguments will be ignored for now and are solely
# added for future use cases (currently we want to limit the different
# iterations of the code).
=#
    p = ArgParseSettings()
    @add_arg_table p begin
    "--nk"
        help = "Number of kpts"
        arg_type = Int
        default = 16
    "--nkf"
        help = "Number of fine kpts"
        arg_type = Int
        default = 9
    "--t"
        help = "Nearest neighbor hopping"
        arg_type = Float64
        default = 1.0
    "--tp"
        help = "Next-Nearest neighbor hopping"
        arg_type = Float64
        default = 0.0
    "--mu"
        help = "Chemical potential"
        arg_type = Float64
        default = 0.0
    "--U"
        help = "Interaction strength"
        arg_type = Float64
        default = 3.0
    end
    params = parse_args(p)
    return params["nk"], params["nkf"], params["U"],
        hopping_parameters(params["t"], params["tp"], params["mu"])
end

function main()
    nk, nkf, U, params = parse_parameters()
    lambdaspan=(20.0, 1E-4)
    divergence = 50.0
    filename = "adaeuler"

    @show nk
    @show nkf
    @show U
    @show params

    @assert(iseven(nk), "nk should be even to include HSP")
    @assert(isodd(nkf), "fine nk should be odd to respect symmetry")

    V = StructArray(initialize_interaction(nk, U))
    dispersion = initialize_dispersion(nk, nkf, params)

    # fixed_step_integrator!(RHS, V, dispersion)
    # magnus1!(RHS, V, dispersion)
    @time adaptive_euler!(RHS, V, dispersion)
    # library_integrator!(RHS,V,dispersion,Tsit5())
    @show maximum(abs.(V))

    println(analyze_vertex(V))

    #save out vertex at end of flow
    file = h5open(filename * ".h5", "w")
    file["hopping/t"] = params.t
    file["hopping/mu"] = params.mu
    file["hopping/tp"] = params.tp
    file["U"] = U
    file["nk"] = nk
    file["nkf"] = nkf
    file["reV"] = V.re
    file["imV"] = V.im
    file["global_results/rhsevals"] = rhsevals
    println("evals of rhs: ", rhsevals)
    close(file)
end


function e(kx::Int64, ky::Int64, nk::Int64,
            kxf::Int64=1, kyf::Int64=1, nkf::Int64=1,
            p::hopping_parameters = hopping_parameters(1.0, -0.1, 0.5)) :: Float64
    #=
    # The t-t' tight binding square lattice dispersion
    =#
    # This binds the index to be symmetric around the midpoint kx/ky.
    # Note that the refinement does not have periodicity
    fine_kx = kxf - ceil(nkf/2.)
    fine_ky = kyf - ceil(nkf/2.)

    kkx = (2*pi/nk) * (kx + fine_kx/nkf)
    kky = (2*pi/nk) * (ky + fine_ky/nkf)

    return -2*p.t*(cos(kkx) + cos(kky)) - 4*p.tp*cos(kkx)*cos(kky) - p.mu
end


function initialize_interaction(nk::Int64, U::Float64=3.) :: Array{ComplexF64,6}
    #=
    # Initialize the interaction: U at all points in the SU2 Hubbard Model
    =#
    return fill(U, (nk,nk,nk,nk,nk,nk))
end


function initialize_dispersion(nk::Int64, nk_fine::Int64, p::hopping_parameters) :: Array{Float64, 4}
    # Cache the dispersion to prevent reevaluations
    dispersion = 
    reshape([e(ix, iy, nk, ifx, ify, nk_fine, p)
                          for ify in 1:nk_fine for ifx in 1:nk_fine
                          for iy in 1:nk for ix in 1:nk],
                         nk, nk, nk_fine, nk_fine)
    return dispersion
end


function RHS(dV::StructArray{ComplexF64, 6}, V::StructArray{ComplexF64, 6}, dispersion::Array{Float64,4},
             lam::Float64, Lpp, Lph) :: Nothing
    #=
    # Evaluation of the RHS diagrams of the frg flow at scale lambda.
    # new dV returned in dV
    # Dispersion does not change over the duration of the flow!
    =#

    dV .= 0
    nk = size(V)[1]

    calc_L(Lpp, lam, dispersion, -1.0)
    calc_L(Lph, lam, dispersion, 1.0)

    @sync begin
        for k2x in 1:nk
        for k2y in 1:nk
            for k1x in 1:nk
            for k1y in 1:nk
                Threads.@spawn begin
                    @turbo for k0x in 1:nk
                    for k0y in 1:nk
                        for lx in 1:nk
                        for ly in 1:nk
                            #hack, mod1 currently not working with loopvectorization
                            lppx = mod(-lx+k0x+k1x-1, nk)+1
                            lppy = mod(-ly+k0y+k1y-1, nk)+1

                            lcphx = mod(lx+k1x-k2x-1, nk)+1
                            lcphy = mod(ly+k1y-k2y-1, nk)+1

                            lphx = mod(lx-k0x+k2x-1, nk)+1
                            lphy = mod(ly-k0y+k2y-1, nk)+1

                            dV.re[k0x,k0y,k1x,k1y,k2x,k2y] +=
                                # P
                                V.re[k0x,k0y,k1x,k1y,lx,ly] *
                                    Lpp.re[lx,ly,lppx, lppy] *
                                    V.re[lx, ly, lppx, lppy, k2x,k2y] +

                                - V.re[k0x,k0y,k1x,k1y,lx,ly] *
                                    Lpp.im[lx,ly,lppx, lppy] *
                                    V.im[lx, ly, lppx, lppy, k2x,k2y] +

                                - V.im[k0x,k0y,k1x,k1y,lx,ly] *
                                    Lpp.re[lx,ly,lppx, lppy] *
                                    V.im[lx, ly, lppx, lppy, k2x,k2y] +

                                - V.im[k0x,k0y,k1x,k1y,lx,ly] *
                                    Lpp.im[lx,ly,lppx, lppy] *
                                    V.re[lx, ly, lppx, lppy, k2x,k2y] +
                                # C
                                - V.re[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    V.re[lx, ly, k1x, k1y, k2x,k2y] +

                                + V.re[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.im[lx,ly,lcphx, lcphy] *
                                    V.im[lx, ly, k1x, k1y, k2x,k2y] +
                                + V.im[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    V.im[lx, ly, k1x, k1y, k2x,k2y] +
                                + V.im[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.im[lx,ly,lcphx, lcphy] *
                                    V.re[lx, ly, k1x, k1y, k2x,k2y] +
                                # D 1
                                +2.0 * V.re[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.re[lx, ly, k1x, k1y, lphx,lphy] +

                                -2.0 * V.re[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.im[lx, ly, k1x, k1y, lphx,lphy] +

                                -2.0 * V.im[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.im[lx, ly, k1x, k1y, lphx,lphy] +

                                -2.0 * V.im[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.re[lx, ly, k1x, k1y, lphx,lphy] +
                                # D 2
                                - V.re[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.re[k1x, k1y, lx, ly, lphx,lphy] +

                                + V.re[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.im[k1x, k1y, lx, ly, lphx,lphy] +

                                + V.im[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.im[k1x, k1y, lx, ly, lphx,lphy] +

                                + V.im[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.re[k1x, k1y, lx, ly, lphx,lphy] +
                                # D 3
                                - V.re[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.re[lx, ly, k1x, k1y, lphx,lphy] +

                                + V.re[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.im[lx, ly, k1x, k1y, lphx,lphy] +

                                + V.im[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.im[lx, ly, k1x, k1y, lphx,lphy] +

                                + V.im[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.re[lx, ly, k1x, k1y, lphx,lphy]



                            dV.im[k0x,k0y,k1x,k1y,k2x,k2y] +=
                                # P
                                V.im[k0x,k0y,k1x,k1y,lx,ly] *
                                    Lpp.re[lx,ly,lppx, lppy] *
                                    V.re[lx, ly, lppx, lppy, k2x,k2y] +

                                + V.re[k0x,k0y,k1x,k1y,lx,ly] *
                                    Lpp.im[lx,ly,lppx, lppy] *
                                    V.re[lx, ly, lppx, lppy, k2x,k2y] +

                                + V.re[k0x,k0y,k1x,k1y,lx,ly] *
                                    Lpp.re[lx,ly,lppx, lppy] *
                                    V.im[lx, ly, lppx, lppy, k2x,k2y] -

                                + V.im[k0x,k0y,k1x,k1y,lx,ly] *
                                    Lpp.im[lx,ly,lppx, lppy] *
                                    V.im[lx, ly, lppx, lppy, k2x,k2y] +
                                # C
                                - V.im[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    V.re[lx, ly, k1x, k1y, k2x,k2y] +

                                - V.re[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    V.im[lx, ly, k1x, k1y, k2x,k2y] +

                                - V.re[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    V.im[lx, ly, k1x, k1y, k2x,k2y] +

                                + V.im[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.im[lx,ly,lcphx, lcphy] *
                                    V.im[lx, ly, k1x, k1y, k2x,k2y] +
                                # D 1
                                +2.0 * V.im[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.re[lx, ly, k1x, k1y, lphx,lphy] +

                                +2.0 * V.re[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.re[lx, ly, k1x, k1y, lphx,lphy] +

                                +2.0 * V.re[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.im[lx, ly, k1x, k1y, lphx,lphy] +

                                -2.0 * V.im[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.im[lx, ly, k1x, k1y, lphx,lphy] +
                                # D 2
                                - V.im[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.re[k1x, k1y, lx, ly, lphx,lphy] +

                                - V.re[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.re[k1x, k1y, lx, ly, lphx,lphy] +
                                
                                - V.re[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.im[k1x, k1y, lx, ly, lphx,lphy] +
                                
                                + V.im[k0x,k0y,lphx,lphy,k2x,k2y] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.im[k1x, k1y, lx, ly, lphx,lphy] +
                                # D 3
                                - V.im[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.re[lx, ly, k1x, k1y, lphx,lphy] +

                                - V.re[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.re[lx, ly, k1x, k1y, lphx,lphy] +

                                - V.re[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    V.im[lx, ly, k1x, k1y, lphx,lphy] +

                                + V.im[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    V.im[lx, ly, k1x, k1y, lphx,lphy]


                                    
                        end end
                    end end
                end
            end end
        end end
    end

    # Weighting of the integral
    dV .*= (1.0 / nk)^2
    
    global rhsevals += 1
    return nothing
end


# TODO we should be able to combine these into one function calling the other?
function expRHS(dV::StructArray{ComplexF64, 6}, V::StructArray{ComplexF64, 6},
            VR::StructArray{ComplexF64, 6}, dispersion::Array{Float64,4},
             lam::Float64, Lpp, Lph) :: Nothing
    #=
    # new dV returned in dV
    # expRHS(V,V) = RHS(V)
    =#

    dV .= 0
    nk = size(V)[1]

    calc_L(Lpp,lam, dispersion, -1.0)
    calc_L(Lph,lam, dispersion, 1.0)

    @sync begin
        for k2x in 1:nk
        for k2y in 1:nk
            for k1x in 1:nk
            for k1y in 1:nk
                Threads.@spawn begin
                    @turbo for k0x in 1:nk
                    for k0y in 1:nk

                        for lx in 1:nk
                        for ly in 1:nk
                            lppx = mod(-lx+k0x+k1x-1, nk)+1
                            lppy = mod(-ly+k0y+k1y-1, nk)+1

                            lcphx = mod(lx+k1x-k2x-1, nk)+1
                            lcphy = mod(ly+k1y-k2y-1, nk)+1

                            lphx = mod(lx-k0x+k2x-1, nk)+1
                            lphy = mod(ly-k0y+k2y-1, nk)+1

                            dV.re[k0x,k0y,k1x,k1y,k2x,k2y] +=
                                # P
                                V.re[lx, ly, lppx, lppy, k2x,k2y] *
                                    Lpp.re[lx,ly,lppx, lppy] *
                                     VR.re[k0x,k0y,k1x,k1y,lx,ly] +

                                - V.re[lx, ly, lppx, lppy, k2x,k2y] *
                                    Lpp.im[lx,ly,lppx, lppy] *
                                     VR.im[k0x,k0y,k1x,k1y,lx,ly] +

                                - V.im[lx, ly, lppx, lppy, k2x,k2y] *
                                    Lpp.re[lx,ly,lppx, lppy] *
                                     VR.im[k0x,k0y,k1x,k1y,lx,ly] +

                                - V.im[lx, ly, lppx, lppy, k2x,k2y] *
                                    Lpp.im[lx,ly,lppx, lppy] *
                                     VR.re[k0x,k0y,k1x,k1y,lx,ly] +
                                # C
                                - V.re[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    VR.re[lx, ly, k1x, k1y, k2x,k2y] +

                                + V.re[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.im[lx,ly,lcphx, lcphy] *
                                    VR.im[lx, ly, k1x, k1y, k2x,k2y] +

                                + V.im[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    VR.im[lx, ly, k1x, k1y, k2x,k2y] +

                                + V.im[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.im[lx,ly,lcphx, lcphy] *
                                    VR.re[lx, ly, k1x, k1y, k2x,k2y] +
                                # D 1
                                +2.0 * V.re[lx, ly, k1x, k1y, lphx,lphy] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.re[k0x,k0y,lphx,lphy,k2x,k2y] +

                                -2.0 * V.re[lx, ly, k1x, k1y, lphx,lphy] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.im[k0x,k0y,lphx,lphy,k2x,k2y] +

                                -2.0 * V.im[lx, ly, k1x, k1y, lphx,lphy] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.im[k0x,k0y,lphx,lphy,k2x,k2y] +

                                -2.0 * V.im[lx, ly, k1x, k1y, lphx,lphy] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.re[k0x,k0y,lphx,lphy,k2x,k2y] +
                                # D 2
                                - V.re[k1x, k1y, lx, ly, lphx,lphy] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.re[k0x, k0y, lphx, lphy, k2x, k2y] +

                                + V.re[k1x, k1y, lx, ly, lphx,lphy] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.im[k0x, k0y, lphx, lphy, k2x, k2y] +

                                + V.im[k1x, k1y, lx, ly, lphx,lphy] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.im[k0x, k0y, lphx, lphy, k2x, k2y] +

                                + V.im[k1x, k1y, lx, ly, lphx,lphy] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.re[k0x, k0y, lphx, lphy, k2x, k2y] +
                                # D 3
                                - V.re[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.re[lx, ly, k1x, k1y, lphx,lphy] +

                                + V.re[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.im[lx, ly, k1x, k1y, lphx,lphy] +

                                + V.im[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.im[lx, ly, k1x, k1y, lphx,lphy] +

                                + V.im[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.re[lx, ly, k1x, k1y, lphx,lphy]

                                    dV.im[k0x,k0y,k1x,k1y,k2x,k2y] +=
                                # P
                                V.im[lx, ly, lppx, lppy, k2x,k2y] *
                                    Lpp.re[lx,ly,lppx, lppy] *
                                     VR.re[k0x,k0y,k1x,k1y,lx,ly] +

                                + V.re[lx, ly, lppx, lppy, k2x,k2y] *
                                    Lpp.im[lx,ly,lppx, lppy] *
                                     VR.re[k0x,k0y,k1x,k1y,lx,ly] +

                                + V.re[lx, ly, lppx, lppy, k2x,k2y] *
                                    Lpp.re[lx,ly,lppx, lppy] *
                                     VR.im[k0x,k0y,k1x,k1y,lx,ly] +

                                - V.im[lx, ly, lppx, lppy, k2x,k2y] *
                                    Lpp.im[lx,ly,lppx, lppy] *
                                     VR.im[k0x,k0y,k1x,k1y,lx,ly] +
                                # C
                                - V.im[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    VR.re[lx, ly, k1x, k1y, k2x,k2y] +

                                - V.re[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.im[lx,ly,lcphx, lcphy] *
                                    VR.re[lx, ly, k1x, k1y, k2x,k2y] +

                                - V.re[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.re[lx,ly,lcphx, lcphy] *
                                    VR.im[lx, ly, k1x, k1y, k2x,k2y] +

                                + V.im[k0x,k0y,lcphx,lcphy,lx,ly] *
                                    Lph.im[lx,ly,lcphx, lcphy] *
                                    VR.im[lx, ly, k1x, k1y, k2x,k2y] +
                                # D 1
                                +2.0 * V.im[lx, ly, k1x, k1y, lphx,lphy] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.re[k0x,k0y,lphx,lphy,k2x,k2y] +

                                +2.0 * V.re[lx, ly, k1x, k1y, lphx,lphy] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.re[k0x,k0y,lphx,lphy,k2x,k2y] +

                                +2.0 * V.re[lx, ly, k1x, k1y, lphx,lphy] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.im[k0x,k0y,lphx,lphy,k2x,k2y] +

                                -2.0 * V.im[lx, ly, k1x, k1y, lphx,lphy] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.im[k0x,k0y,lphx,lphy,k2x,k2y] +
                                # D 2
                                - V.im[k1x, k1y, lx, ly, lphx,lphy] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                     VR.re[k0x, k0y, lphx, lphy, k2x, k2y] +

                                - V.re[k1x, k1y, lx, ly, lphx,lphy] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                     VR.re[k0x, k0y, lphx, lphy, k2x, k2y] +

                                - V.re[k1x, k1y, lx, ly, lphx,lphy] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                     VR.im[k0x, k0y, lphx, lphy, k2x, k2y] +

                                + V.im[k1x, k1y, lx, ly, lphx,lphy] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                     VR.im[k0x, k0y, lphx, lphy, k2x, k2y] +
                                # D 3
                                - V.im[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.re[lx, ly, k1x, k1y, lphx,lphy] +

                                - V.re[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.re[lx, ly, k1x, k1y, lphx,lphy] +

                                - V.re[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.re[lx,ly,lphx, lphy] *
                                    VR.im[lx, ly, k1x, k1y, lphx,lphy] +

                                + V.im[k0x,k0y,lphx,lphy,lx,ly] *
                                    Lph.im[lx,ly,lphx, lphy] *
                                    VR.im[lx, ly, k1x, k1y, lphx,lphy]
                        end end
                    end end
                end
            end end
        end end
    end

    # Weighting of the integral
    dV .*= (1.0 / nk)^2
    global rhsevals += 1

    return nothing
end


function analyze_vertex(V::StructArray{ComplexF64, 6})
    #=
    # Function for a quick and dirty analysis of the FRG vertex to determine the
    # possible phases. We note that as there are only three possible phases in
    # the square lattice Hubbard model this is possible. Otherwise it would not
    # be the recommended approach.
    =#

    nk = size(V)[1]

    maximum_index = argmax(abs.(V))

    k0x = maximum_index[1]
    k0y = maximum_index[2]

    k1x = maximum_index[3]
    k1y = maximum_index[4]

    k2x = maximum_index[5]
    k2y = maximum_index[6]


    # TODO Is this correct?
    phase = ""
    # Check if Vmax is at q_P = k0 + k1 = [0,0]
    if (mod1(k0x + k1x, nk) == nk) && (mod1(k0y + k1y, nk) == nk)
        phase = "SC"
    # Check if Vmax is at q_C = k2 - k1 = [Pi, Pi]
    elseif (mod1(k2x - k1x, nk) == nk/2) && (mod1(k2y - k1y, nk) == nk/2)
        phase = "AFM"
    # Check if Vmax is at q_C = k2 - k1 = [0, 0]
    elseif (mod1(k2x - k1x, nk) == nk) && (mod1(k2y - k1y, nk) == nk)
        phase = "FM"
    else
        println("Could not identify the Phase. This might be an error")
    end

    return phase
end
