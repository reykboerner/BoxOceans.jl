struct StommelLorenz84 < BoxModel end

function StommelLorenz84(u,p,t)
    x, y, z, T, S = u
    a, b, μ, ϵa, ϵf, F0, F1, G0, G1, θ0, θ1, σ0, σ1, x_mean, Δ_mean = p

    Δ = y^2 + z^2
    T_surf = θ0 + θ1 * (x - x_mean)/(sqrt(ϵf))
    S_surf = σ0 + σ1 * (Δ - Δ_mean)/(sqrt(ϵf))
    dx = 1/ϵf*(-Δ - a*(x - F0 - F1*T))
    dy = 1/ϵf*(x*y - b*x*z - (y - G0 + G1*T))
    dz = 1/ϵf*(b*x*y + x*z - z)
    dT = -1/ϵa*(T - T_surf) - T - μ*sign(S-T)*(S-T)*T
    dS = S_surf - S - μ*sign(S-T)*(S-T)*S

    SA[dx, dy, dz, dT, dS]
end

StommelLorenz84.default_params = [
    0.25,       # a
    4.0,        # b
    7.5,        # mu
    0.34,       # ϵa
    3e-4,       # ϵf
    8.0,        # F0
    0.1,        # F1
    1.0,        # G0
    0.0,        # G1
    1.0,        # θ0
    0.0195,     # θ1
    0.9,        # σ0
    0.00934,    # σ1
    1.0147,     # x_mean
    1.7463      # Δ_mean
]