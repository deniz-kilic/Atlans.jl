function pore_pressure!(gw::GW where {GW<:GroundwaterColumn})
    @. gw.p = gw.γ_water * (gw.ϕ - (gw.z + 0.5 * gw.Δz))
    @. gw.dry = gw.p < 0.0
    gw.p[gw.dry] .= 0.0
    return
end

function plot(gw::GW where {GW<:GroundwaterColumn})
    plot(gw.ϕ, gw.z, color = :blue, label = "ϕ")
    plot!(gw.p / gw.γ_water, gw.z, color = :red, label = "p")
    ybot = gw.z .- 0.5 .* gw.Δz
    ytop = gw.z[end] + 0.5 * gw.Δz[end]
    hline!(ybot, color = :black, label = "")  # bottom
    hline!([ytop], color = :black, label = "") # top
end
