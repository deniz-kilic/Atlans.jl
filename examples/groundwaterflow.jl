# # Examples of groundwater flow
#
# This example demonstrates two concepts:
#
# * Linear head gradients between boundaries
# * One-dimensional vertical groundwater flow between boundaries

using Atlans

const Float = Float64

sg = SimpleGroundwater(
    collect(0.0:1.0:5.0),  # z
    fill(1.0, 5),  # Δz
    Int[1, 5],  # boundary
    Float[3.0, 4.0],  # boundary_ϕ
    fill(false, 5),  # dry
    fill(0.0, 5),  # ϕ
)
solve!(sg)
pore_pressure!(sg)
