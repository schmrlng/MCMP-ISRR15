import Base.length
using StatsBase
using Distributions
using HypothesisTests
using ArrayViews
using Convex
using Devectorize
# redirect_stderr() # to shut CVX "inaccurate solution" up
include("mvnormalutils.jl")

### LQG typedef and utilities
type DiscreteLQG
    # Notation follows http://en.wikipedia.org/wiki/Linear-quadratic-Gaussian_control
    # Dynamics
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}      # output matrix for tracking controller
    Cws::Matrix{Float64}    # state space -> workspace
    
    # Cost function
    Q::Matrix{Float64}
    R::Matrix{Float64}
    F::Matrix{Float64}
    
    # Noise characteristics
    V::Matrix{Float64}
    W::Matrix{Float64}
    P0::Matrix{Float64}

    # To be computed
    P::Vector{Matrix{Float64}}
    S::Vector{Matrix{Float64}}
    K::Vector{Matrix{Float64}}
    L::Vector{Matrix{Float64}}
    dim::Int
    Ncomb::Vector{MvNormal}
    Σinv::Vector{PDMats.PDMat}
    Acomb::Matrix{Matrix{Float64}}
    T::Int
    dt::Float64
end

function DiscreteLQG(A,B,C,Cws,Q,R,F,V,W,P0,Tmax,dt)
    dim = size(A,1)
    P = Array(typeof(P0), Tmax+1)
    P[1] = P0
    S = Array(typeof(F), Tmax+1)
    S[1] = F                            # indexed backwards
    K = Array(typeof(A), Tmax)
    L = Array(typeof(A), Tmax)          # indexed backwards
    Ncomb = Array(MvNormal, Tmax+1)
    Σinv = Array(PDMats.PDMat, Tmax+1)
    Ncomb[1] = MvNormal(zeros(2*dim), [P0 zeros(A); zeros(A) sqrt(eps())*eye(A)]) # MvNormal requires PosDef covariance
    Σinv[1] = PDMats.inv(Ncomb[1][1:dim].Σ)
    for i in 1:Tmax
        P[i+1] = A*(P[i] - P[i]*C'*((C*P[i]*C' + W)\C)*P[i])*A' + V
        S[i+1] = A'*(S[i] - S[i]*B*((B'*S[i]*B + R)\B')*S[i])*A + Q
        K[i] = A*P[i]*C'/(C*P[i]*C' + W)
        L[i] = (B'*S[i]*B + R)\B'*S[i]*A
        Ncomb[i+1] = MvNormal(zeros(2*dim), [V zeros(A); zeros(A) K[i]*W*K[i]'])
        Σinv[i+1] = PDMats.inv(Ncomb[i+1].Σ)
    end
    Acomb = Array(typeof(A), 0, 0)
    DiscreteLQG(A,B,C,Cws,Q,R,F,V,W,P0,P,S,K,L,dim,Ncomb,Σinv,Acomb,0,dt)
end

function DiscreteLQG(SS::LinearQuadraticStateSpace, Wc, Vc, P0; dt = .025, nsf = 1.0, max_time = 20.)
    function GhettoIntegrate(f, t0, t1, N)   # Trapezoid rule would be a better name
        tt = linspace(t0, t1, N)
        ds = diff(tt)
        vv = map(f, tt)
        .5*(sum(map(*, vv[1:end-1], ds)) + sum(map(*, vv[2:end], ds)))
    end
    A = expm(SS.A*dt)
    B = GhettoIntegrate(s -> expm(SS.A*s), 0, dt, 20)*SS.B
    C = eye(SS.dim)
    Cws = SS.C
    Q = 5*eye(SS.dim)           # state regulator
    R = 1*eye(size(SS.B,2))     # control effort
    F = Q                       # final state penalty
    V = nsf*GhettoIntegrate(s -> expm(SS.A*s)*Vc*expm(SS.A*s)', 0, dt, 20)   # process noise
    W = nsf/dt*Wc                                                            # measurement noise
    P0 = nsf*P0                                                              # initial uncertainty
    DiscreteLQG(A,B,C,Cws,Q,R,F,V,W,P0,iceil(max_time/dt),dt)
end

function SingleIntegrator(dim = 2; dt = .025, nsf = 1.0, max_time = 2*dim)
    A = eye(dim)
    B = dt*eye(dim)    # assuming unit speed
    C = eye(dim)
    Cws = eye(dim)
    Q = 5*eye(dim)     # state regulator
    R = 1*eye(dim)     # control effort
    F = Q              # final state penalty
    if dim == 2
        V = 2*nsf*.2*dt*.01^2*[10 8; 8 10]         # process noise
        W = nsf*.2/(4*dt)*.01^2*[10 -8; -8 10]       # measurement noise
        P0 = nsf*.2*.25*.01^2*[10 -8; -8 10]       # initial uncertainty
    elseif dim == 3
        CM = full(blkdiag(sparse(eig([0 1; 1 0])[2]), sparse([1])))
        V = nsf*.2*dt*.01^2*CM*diagm([4,9,9])*CM'
        W = nsf*.2/(4*dt)*.01^2*CM*diagm([4,4,1])*CM'
        P0 = V/9
    end
    DiscreteLQG(A,B,C,Cws,Q,R,F,V,W,P0,iceil(max_time/dt),dt)
end

function Acomb(D::DiscreteLQG, t::Int, T::Int)
    [D.A          (-D.B*D.L[T-t+1]);
     D.K[t]*D.C   D.A-D.B*D.L[T-t+1]-D.K[t]*D.C]
end

function sethorizon(D::DiscreteLQG, T::Int, w=0)
    T > length(D.L) && error("DiscreteLQG horizon set beyond max horizon")
    D.Acomb = Array(typeof(D.A), T, w+1)    # Acomb[t,dt] = Ac_{t+dt-1} * ... * Ac_t
    for t in 1:T
        D.Acomb[t,1] = Acomb(D,t,T)
    end
    for s in 1:w, t in 1:T-s
        D.Acomb[t,s+1] = D.Acomb[t+1,s]*D.Acomb[t,1]
    end
    D.T = T
end

### Path and associated metadata for MC estimation
type COP                    # close obstacle point
    k::Int                  # waypoint index
    d2::Float64             # square mahalanobis distance
    v::Vector{Float64}       # obstacle point
    hpv::Vector{Float64}     # vector defining obstacle halfplane
    hpbp::Float64           # halfplane breach probability
end
type LQGPath
    path::Path
    D::DiscreteLQG
    pwu::Vector{MvNormal}  # uncertainty in workspace only
    CC::CollisionChecker
    cops::Vector{COP}
    pthresh::Float64

    function LQGPath(path::Path, D::DiscreteLQG)
        T = length(path) - 1
        T != D.T && sethorizon(D, T)
        combined_unc = Array(MvNormal, T+1)
        combined_unc[1] = D.Ncomb[1]
        for t in 1:T
            combined_unc[t+1] = D.Acomb[t,1] * combined_unc[t] + D.Ncomb[t+1]
        end
        pwu = [c[1:D.dim] for c in combined_unc]
        new(path, D, pwu)
    end
end
length(LP::LQGPath) = length(LP.path)

function noisify(LP::LQGPath; seed::Int = 0)
    path = LP.path
    D = LP.D
    T = length(path) - 1
    T != D.T && sethorizon(D, T)
    seed > 0 && srand(seed)
    xcomb = Array(Vector{Float64}, T+1)
    xcomb[1] = [rand(D.Ncomb[1][1:D.dim]), zeros(D.dim)]
    for t in 1:T
        xcomb[t+1] = D.Acomb[t,1] * xcomb[t] + rand(D.Ncomb[t+1])
    end
    path + Vector{Float64}[xc[1:D.dim] for xc in xcomb]
end

## LQG closest obstacle points
function nonoccluded_cops(p::AbstractVector, CC::CollisionChecker, W::AbstractMatrix, pthresh = 0.)
    cands = close(p, CC, W, halfplanecquantile(pthresh))    # sorted by Mahalanobis distance
    selector = trues(length(cands))
    for (i,dw) in enumerate(cands)
        if selector[i]
            w = dw[2]
            g = W*(w-p)
            d = dot(g, (w-p))
            for j in i+1:length(cands)
                selector[j] = selector[j] && (d > dot(g, (cands[j][2]-p)))
            end
        end
    end
    cands[selector]
end

function computecops(LP::LQGPath, CC::CollisionChecker, pthresh = 0.)
    LP.CC = CC
    LP.cops = vcat([[(W = sqrtm(full(LP.pwu[i].Σ));                     # likely the most ridiculous comprehension I've ever written
                      p = LP.path[i];
                      vf = p + W*pinv(LP.D.Cws*W)*(v-LP.D.Cws*p);       # d2f (hell yeah she is) is no different
                      COP(i, d2, vf, cop_to_hpv(LP.path[i], vf, LP.pwu[i]), halfplanetail(d2)))
                     for (d2,v) in nonoccluded_cops(LP.D.Cws*LP.path[i], CC, full(inv((LP.D.Cws*LP.pwu[i]).Σ)), pthresh)]
                    for i in 1:length(LP)]...)
    LP.pthresh = pthresh
end

## Patil et al. http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6224727&tag=1
function cop_to_hpv(path_pt::Vector, obs_pt::Vector, uncertainty::MvNormal) # hpv with path_pt at origin
    grad = inv(uncertainty.Σ)*(obs_pt-(path_pt+uncertainty.μ))
    vec(((obs_pt-path_pt)'*grad/(grad'*grad))[1]*grad)
end
function pointwise_pruned_uncertainty_CP_estimate(LP::LQGPath, CC::CollisionChecker)
    path = LP.path
    D = LP.D
    T = length(path) - 1
    T != D.T && sethorizon(D, T)
    combined_unc = Array(MvNormal, T+1)
    pruned_unc = Array(MvNormal, T+1)
    combined_unc[1] = D.Ncomb[1]
    CCP_estimate = 1.
    for t in 1:T+1
        copsws = nonoccluded_cops(D.Cws*(path[t]+combined_unc[t].μ[1:D.dim]), CC, full(inv((D.Cws*combined_unc[t][1:D.dim]).Σ)))
        W = sqrtm(full(combined_unc[t][1:D.dim].Σ))
        cops = [(path[t] + W*pinv(D.Cws*W)*(v-D.Cws*path[t])) for (d2,v) in copsws]
        hpvs = Vector{Float64}[[cop_to_hpv(path[t], cop, combined_unc[t][1:D.dim]), zeros(D.dim)] for cop in cops]
        CCP_estimate *= 1 - sum([halfplanetail(combined_unc[t], hpv) for hpv in hpvs])
        t > T && break
        pruned_unc[t] = prunenormal(combined_unc[t], hpvs)
        combined_unc[t+1] = D.Acomb[t,1] * pruned_unc[t] + D.Ncomb[t+1]
    end
    1 - CCP_estimate
end

## Pointwise probabilities
function half_plane_breach_probabilities(LP::LQGPath, CC::CollisionChecker, pthresh = 0.)
    (!isdefined(LP, :cops) || LP.pthresh != pthresh) && computecops(LP, CC, pthresh)
    Float64[cop.hpbp for cop in LP.cops]
end

function ellipsoid_breach_probabilities(LP::LQGPath, CC::CollisionChecker)
    ccop = [closest(LP.D.Cws*LP.path[i], CC, full(inv((LP.D.Cws*LP.pwu[i]).Σ))) for i in 1:length(LP)]
    [ellipsoidtail(length(LP.D.Cws*LP.pwu[i]), ccop[i][1]) for i in 1:length(LP)]
end

function half_plane_breach_count(LP::LQGPath, observed::Path)   # TODO: also try ellipsoid_breaches
    sum([dot(cop.hpv, observed[cop.k] - LP.path[cop.k] - cop.hpv) > 0 for cop in LP.cops])
end

type ISDistributionCache
    ds::Matrix{MvNormal}        # ds[i,t]
    wsizes::Vector{Int}
end

function ISDistributionCache(LP::LQGPath, cw = 0)
    path = LP.path
    T = length(path) - 1
    D = LP.D
    if T != D.T || cw != size(D.Acomb,2) - 1
        sethorizon(D, T, cw)    # Profile -> this could definitely take less time
    end
    ds = Array(MvNormal, length(LP.cops), T+1)
    wsizes = fill(cw, length(LP.cops))
    for (i, cop) in enumerate(LP.cops)
        k = cop.k
        for t in 1:T+1
            ds[i,t] = D.Ncomb[t]
        end
        bigA = hcat([D.Acomb[t,k-t] for t in max(k-cw,1):k-1]..., eye(2*D.dim))
        bigW = full(blkdiag([sparse(sqrtm(full(D.Ncomb[t].Σ))) for t in max(k-cw,1):k]...))
        bigV = bigW*pinv(bigA[1:D.dim,:]*bigW)*(cop.v - path[k])
        bigV = reshape(bigV, 2*D.dim, div(length(bigV), 2*D.dim))
        for s in max(k-cw,1):k
            ds[i,s] = ds[i,s] + bigV[:,s+1-max(k-cw,1)]
        end
        wsizes[i] = k - max(k-cw,1)
    end
    ISDistributionCache(ds, wsizes)
end

function noisify_with_kick(LP::LQGPath, alpha;
                           i = -1,              # we're sampling to make collision at i more likely
                           seed = 0,
                           cw = 0,              # the collision window i-cw:i over which we're adjusting the noise distribution
                           ISDC = ISDistributionCache(LP, cw))
    path = LP.path
    T = length(path) - 1
    D = LP.D
    if T != D.T || cw != size(D.Acomb,2) - 1
        sethorizon(D, T, cw)
    end
    seed > 0 && srand(seed)
    i == -1 && (i = rand(Categorical(alpha)))
    xcomb = Array(Vector{Float64}, T+1)

    ### pre-store noise for corresponding non-IS particle in xcomb
    xcomb[1] = rand(D.Ncomb[1][1:D.dim])
    for t in 1:T
        xcomb[t+1] = rand(D.Ncomb[t+1])
    end
    ###

    k = LP.cops[i].k
    for t in k-ISDC.wsizes[i]:k
        if t == 1
            xcomb[1] = rand(ISDC.ds[i,1])[1:D.dim]
        else
            xcomb[t] = rand(ISDC.ds[i,t])
        end
    end

    log_inv_likelihood_ratio = log(alpha)
    nominal_logpdfs = [PDMats.quad(D.Σinv[1], xcomb[1]), [PDMats.quad(D.Σinv[t], xcomb[t]) for t in 2:T+1]]
    diff1 = zeros(Float64, D.dim)
    diff2 = zeros(Float64, 2*D.dim)
    for j in 1:length(LP.cops)
        k = LP.cops[j].k
        for t in k-ISDC.wsizes[j]:k
            if t == 1
                for ii in 1:D.dim
                    diff1[ii] = xcomb[t][ii] - ISDC.ds[j,t].μ[ii]
                end
                log_inv_likelihood_ratio[j] += -0.5*(PDMats.quad(D.Σinv[t], diff1) - nominal_logpdfs[t])
            else
                for ii in 1:2*D.dim
                    diff2[ii] = xcomb[t][ii] - ISDC.ds[j,t].μ[ii]
                end
                log_inv_likelihood_ratio[j] += -0.5*(PDMats.quad(D.Σinv[t], diff2) - nominal_logpdfs[t])
            end
        end
    end

    xcomb[1] = [xcomb[1], zeros(xcomb[1])]
    for t in 1:T
        xcomb[t+1] = D.Acomb[t,1] * xcomb[t] + xcomb[t+1]
    end

    path + Vector{Float64}[xc[1:D.dim] for xc in xcomb], (1 / sum(exp(log_inv_likelihood_ratio)))
end