function distance(k1::Vector{Float64}, k2::Vector{Float64}) :: Float64
    return norm(k1-k2)
end

function closest(k::Vector{Float64}, mesh::Vector{Vector{Float64}}) :: Int64
    closest = pi
    closest_idx = 0

    for ii in eachindex(mesh)
        kk = mesh[ii]
        if distance(k,kk) < closest
            closest_idx = ii
            closest = distance(k,kk)
        end
    end

    return closest_idx
end


function get_hsp(pts::Vector{Vector{Float64}}) :: Tuple{Vector{Vector{Float64}}, Vector{Float64}, Vector{Float64}}
    
    N = 100

    path = [pts[1]]
    node_idx = [1]

    for ii in 1:length(pts)-1
        path_piece = collect(range(pts[ii], pts[ii+1]; length=N))
        path = vcat(path, path_piece[2:end])
        append!(node_idx, length(path))
    end

    dist_pts = cumsum(vcat(0.0,[distance(path[i-1], path[i]) for i in 2:length(path)]))

    dist_nodes = dist_pts[node_idx]

    return path, dist_pts, dist_nodes
end


function plot_along_hsp(dist_on_hsp::Vector{Float64}, dist_node::Vector{Float64}, labels::Vector{String}, energy::Vector{Float64}, title::String="",zero_line::Bool=false) :: Nothing
    fig, ax = subplots()
    ax.set_xlim(dist_node[1], dist_node[end-1])
    ax.set_xticks(dist_node)
    ax.set_xticklabels(labels)

    for n in 1:length(dist_node)
        ax.axvline(x=dist_node[n], linewidth=0.5, color="k")
    end
    ax.set_title(title)
    ax.set_xlabel("Path in k-space")
    ax.set_ylabel("Energy")

    ax.plot(dist_on_hsp, energy)
    if zero_line
        ax.plot(dist_on_hsp, np.zeros_like(data), "k--")
    end

    grid()
    fig.tight_layout()
    show()

    return nothing
end


function get_dispersion(disp_pts::Vector{Vector{Float64}}, mesh::Vector{Vector{Float64}}, energy::Vector{Float64}) :: Vector{Float64}
    disp = zeros(length(disp_pts))

    for ii in eachindex(disp_pts)
        kk = disp_pts[ii]
        nearest = closest(kk, mesh)
        disp[ii] = energy[nearest]
    end

    return disp
end


function HSP_dispersion(mesh :: Vector{Vector{Float64}}, disp::Vector{Float64}) ::Nothing 
    G = [0.0, 0.0]
    X = [pi, 0.0]
    M = [pi, pi]

    pts = [G, X, M, G]
    labels = ["\$\\Gamma\$", "\$X\$", "\$M\$", "\$\\Gamma\$"]

    path, dist_on_hsp, dist_nodes = get_hsp(pts)
    HSP_disp = get_dispersion(path, mesh, disp)

    plot_along_hsp(dist_on_hsp, dist_nodes, labels, HSP_disp)

    return nothing
end
