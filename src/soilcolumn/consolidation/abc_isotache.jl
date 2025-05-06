
# compute intrinsic time (τ)
function set_τ0_ocr(abc::ABC, ocr::Float) where {ABC <: AbstractAbcIsotache}
    if abc.c < 1.0e-4
        τ = 1.0e-9
    else
        τ = τ_ref * ocr^((abc.b - abc.a) / abc.c)
    end
    return @set abc.τ = τ
end


function set_τ0_pop(abc::ABC, pop::Float) where {ABC <: AbstractAbcIsotache}
    ocr = (abc.σ′ + pop * 1.0e3) / abc.σ′
    return set_τ0_ocr(abc, ocr)
end


struct AbcIsotache <: AbstractAbcIsotache
    Δz::Float
    Δz_0::Float
    t::Float
    σ′::Float  # effective stress
    γ_wet::Float  # wet specific mass
    γ_dry::Float  # dry specific mass
    # Isotache parameters
    a::Float
    b::Float
    c::Float
    τ::Float
end


function AbcIsotache(Δz, γ_wet, γ_dry, a, b, c)
    return AbcIsotache(Δz, Δz, 0.0, NaN, γ_wet, γ_dry, a, b, c, NaN)
end


function consolidate!(abc::AbcIsotache, σ′::Float, Δt::Float)
    t = abc.t + Δt
    # Effective stress changes
    Δσ′ = σ′ - abc.σ′
    # τ changes
    τ = τ⃰ + Δt
    # consolidation
    strain = abc.c * log(abc.τ / τ⃰) + log(σ′ / (σ′ - loadstep))
    # Thickness should not go below 0
    consolidation = min(abc.Δz, strain * abc.Δz)
    # return new state
    return AbcIsotache(
        abc.Δz,
        abc.Δz_0,
        t,  # new
        σ′, # new
        abc.γ_wet,  # new
        abc.γ_dry,  # new
        abc.a,
        abc.b,
        abc.c,
        τ,  # new
    )
end


# 2.5x faster power method
"Faster method for exponentiation"
pow(x, y) = exp(y * log(x))

"""
Derivative of Qcreep with respect to head.
"""
function Qcreep_derivative(abc::AbcIsotache, σ′::Float, Δt::Float)
    exponent = (abc.b - abc.a) / abc.c
    numerator = γ_water * (abc.b - abc.a)
    denominator = -σ′ * (1.0 + abc.τ / Δt * pow((abc.σ′ / σ′), exponent))
    return numerator / denominator
end


function Qcreep(abc::AbcIsotache, σ′::Float, Δt)
    exponent = (abc.b - abc.a) / abc.c
    τ⃰ = abc.τ * pow(abc.σ / σ′, exponent)
    return abc.c * log((τ⃰ + Δt) / τ⃰)
end


#function Qcreep(abc::AbcIsotache, σ′::Float, Δt::Float)
#    exponent = (a - b) / c
#    return abc.c * log(pow(1.0 + Δt / abc.τ * abc.σ′ / σ′, exponent))
#end


"""
First term of Taylor expansion

f(x) ~= f(a) + f′(a) * (x - a)
"""
function formulate(abc::AbcIsotache, ϕ::Float, z::Float, σ::Float, Δt::Float)
    σ′ = σ - γ_water * (ϕ - z)
    deriv = Qcreep_derivative(abc, σ′, Δt)
    hcofterm = deriv
    rhsterm = Qcreep(abc, σ′, Δt) - deriv * ϕ
    return hcofterm, rhsterm
end

