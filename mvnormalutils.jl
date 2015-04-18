using Distributions

+(x::MvNormal, y::MvNormal) = MvNormal(x.μ + y.μ, x.Σ + y.Σ)
+(x::MvNormal, c::Vector) = MvNormal(x.μ + c, x.Σ)
+(c::Vector, x::MvNormal) = MvNormal(x.μ + c, x.Σ)

*(A::Matrix, x::MvNormal) = MvNormal(A * x.μ, PDMats.X_A_Xt(x.Σ, A))

getindex(x::MvNormal, i) = MvNormal(x.μ[i], cov(x)[i,i])

halfplanetail(sqmhd::Float64) = .5*erfc(sqrt(sqmhd/2))
halfplanecquantile(cp::Float64) = 2*erfcinv(2*cp)^2
halfplanetail(x::MvNormal, v::Vector) = ccdf(Normal(0, sqrt(PDMats.quad(x.Σ, v))), dot(v, v - x.μ))  # v normal to half-plane; PDMats.quad(x.Σ, v)) != sqmhd above

ellipsoidtail(dim::Int, sqmhd::Float64) = ccdf(Chisq(dim), sqmhd)
ellipsoidtail(x::MvNormal, v::Vector) = ellipsoidtail(length(x), sqmahal(x, v - x.μ))

function shiftnormal(x::MvNormal, v::Vector, sigma::Real = 1.0)
    x + v * sigma / sqrt(PDMats.invquad(x.Σ, v))
end

function prunenormal(x::MvNormal, hpvs::Vector{Vector{Float64}})
    μ_delta = zeros(x.μ)
    Σ_delta = zeros(full(x.Σ))
    for v in hpvs
        alpha = (dot(v,v) - dot(v,x.μ)) / sqrt(PDMats.quad(x.Σ, v))
        lambda = pdf(Normal(0, 1), alpha) / cdf(Normal(0,1), alpha)
        mu = dot(v, x.μ) + lambda * sqrt(PDMats.quad(x.Σ, v))
        sig2 = PDMats.quad(x.Σ, v) * (1 - lambda^2 + alpha*lambda)
        μ_delta += ((x.Σ * v) / PDMats.quad(x.Σ, v)) * (dot(v, x.μ) - mu)
        Σ_delta += ((x.Σ * v) / PDMats.quad(x.Σ, v)) * (PDMats.quad(x.Σ, v) - sig2) * ((x.Σ * v)' / PDMats.quad(x.Σ, v))
    end
    MvNormal(x.μ + μ_delta, x.Σ + Σ_delta)
end

# function kicknormal(x::MvNormal, v::Vector, sigma::Real = 1.0,
#                     P::Matrix = eye(dim(x.Σ)))
#     x + P*v
# end

# function kicknormal(x::MvNormal, v::Vector, sigma::Real = 1.0,
#                     P::Matrix = eye(dim(x.Σ)))
#     MvNormal(x.μ, sigma*x.Σ)
# end

function scalevariance(x::MvNormal, sf::Real = 1.0)
    MvNormal(x.μ, sf*x.Σ)
end