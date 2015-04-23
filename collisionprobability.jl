function shrink_state_space(SS::StateSpace, eps::Float64)   # ultimately unused - I can enforce boundaries with explicit obstacles if desired
    SSshrunk = copy(SS)
    SSshrunk.lo += eps
    SSshrunk.hi -= eps
    SSshrunk
end

function collision_probability_stats(P0::MPProblem, eps::Float64, LQG::DiscreteLQG, Nparticles = 1000; 
                                     cw::Int = -1,
                                     seed0::Int = abs(rand(Int)),
                                     defIS::Float64 = 0.,
                                     progressmod::Int = 1000,
                                     vis::Bool = {},
                                     VR::Bool = true,
                                     alphafilter = 1e-5)
    length(P0.V) < 2 && error("Near neighbor data needs to be prepopulated (e.g. from deterministic solution run)")
    CC0 = P0.CC
    CCI = inflate(P0.CC, eps)
    P0.CC = CCI
    tic()
    fmtstar!(P0, length(P0.V), connections = :R, r = P0.solution.metadata["r"])  # actually r doesn't matter - near neighbors sets are precomputed
    println("Planning Time $(toq())")
    local path::LQGPath
    try
        path = LQGPath(discretize_path(P0, LQG.dt), LQG)
    finally
        P0.CC = CC0
    end
    cw == -1 && (cw = length(path.path))
    local alpha::Vector{Float64}
    while true
        alpha = half_plane_breach_probabilities(path, CC0, alphafilter)
        isempty(alpha) && (alphafilter /= 10; continue)
        /(extrema(alpha)...) * Nparticles < .1 && (alphafilter *= 10; continue)
        break
    end
    theta = sum(alpha)
    alpha_normalized = [alpha, theta*defIS/(1-defIS)] / (theta + theta*defIS/(1-defIS))

    tic()
    particle_paths = [((mod(i, progressmod) == 0 && println(i));
                       noisify(path, seed=seed0+i))
                      for i in 1:Nparticles]
    f = Float64[~is_free_path(p, CC0) for p in particle_paths]
    h = Float64[half_plane_breach_count(path, p) for p in particle_paths]
    println("Naive MC: $(toq())")

    if VR
        tic()
        ISDC = ISDistributionCache(path, cw)
        IS_paths = [((mod(i, progressmod) == 0 && println(i));
                      noisify_with_kick(path, alpha_normalized,
                                        seed=seed0+i, cw=cw, ISDC=ISDC))
                    for i in 1:Nparticles]
        ISf = Float64[lr*(~is_free_path(p, CC0)) for (p, lr) in IS_paths]
        ISh = Float64[lr*half_plane_breach_count(path, p) for (p, lr) in IS_paths]
        println("Variance Reduction: $(toq())")
    else
        ISf = zeros(Nparticles)
        IS_paths = []
    end

    if vis
        plot(P0, sol=false, meta=false)
        plot_path([LQG.Cws*p for p in path.path], color="blue")
        eps > 0 && plot(CCI, P0.SS.lo, P0.SS.hi, alpha=.2)
        for t in 1:length(path.pwu)
            plot_ellipse(LQG.Cws*path.path[t], quantile(Chisq(2), .95)*cov(LQG.Cws*path.pwu[t]), color="purple", alpha=.2)
        end
        plot_line_segments(path.path[[c.k for c in path.cops]],
                           [c.v for c in path.cops], linewidth=.5, linestyle="-", zorder=1, color="black", alpha=0.8)
        plot_line_segments(path.path[[c.k for c in path.cops]],
                           path.path[[c.k for c in path.cops]] + [c.hpv for c in path.cops], linewidth=.5, linestyle="-", zorder=1, color="green", alpha=1.)

        # pwpu = pointwise_pruned_uncertainty(dpath, LQG, P0.obs)
        # for t in 1:1:length(pwpu)
        #     plot_ellipse(dpath[t], quantile(Chisq(2), .95)*cov(pwpu[t]), color="purple", alpha=.2)
        # end
    end

    prunedCPestimate = pointwise_pruned_uncertainty_CP_estimate(path, CC0)

    {
        "P0" => P0,
        "path" => path,
        "f" => f,
        "particle_paths" => particle_paths,
        "h" => h,
        "theta" => theta,
        "ISf" => ISf,
        "ISh" => ISh,
        "IS_paths" => IS_paths,
        "prunedCPestimate" => prunedCPestimate,
        "alpha" => alpha,
        "CP_VR" => mean(ISf) - (cov(ISf,ISh) / var(ISh)) * (mean(ISh) - theta) 
    }
end

function collision_probability(P0::MPProblem, eps::Float64, LQG::DiscreteLQG, Nparticles::Int=1000;
                               verbose::Bool = false,
                               progressmod::Int = 1000,
                               seed0::Int = abs(rand(Int)),
                               method::Symbol = :VR,
                               cw::Int = -1,
                               defIS::Float64 = 0.,
                               alphafilter::Float64 = 1e-5,
                               targeted::Bool = false,
                               CPgoal::Float64 = .01,
                               zhalt = 4.0)
    length(P0.V) < 2 && error("Near neighbor data needs to be prepopulated (e.g. from deterministic solution run)")
    CC0 = P0.CC
    CCI = inflate(P0.CC, eps)
    P0.CC = CCI
    tic()
    fmtstar!(P0, length(P0.V), connections = :R, r = P0.solution.metadata["r"])  # actually r doesn't matter - near neighbors sets are precomputed
    plan_time = toq()
    verbose && println("Planning Time $plan_time")
    local path::LQGPath
    try
        path = LQGPath(discretize_path(P0, LQG.dt), LQG)
    finally
        P0.CC = CC0
    end

    tic()
    i = 0
    if method == :NAIVE
        CP = 0.
        CPvarN = 0.
        CPstd = 0.
        for i in 1:Nparticles
            verbose && (mod(i, progressmod) == 0) && println(i)
            f = ~is_free_path(noisify(path, seed=seed0+i), CC0)
            CPvarN = CPvarN + (i - 1) * (f - CP)^2 / i
            CP = CP + (f - CP) / i
            CPstd = sqrt(CPvarN) / i
            targeted && i > 100 && abs(CP - CPgoal) > zhalt*CPstd && break
        end
    elseif method == :VR
        local alpha::Vector{Float64}
        while true
            alpha = half_plane_breach_probabilities(path, CC0, alphafilter)
            isempty(alpha) && (alphafilter /= 10; continue)
            /(extrema(alpha)...) * Nparticles < .1 && (alphafilter *= 10; continue)
            break
        end
        theta = sum(alpha)
        alpha_normalized = [alpha, theta*defIS/(1-defIS)] / (theta + theta*defIS/(1-defIS))
        cw == -1 && (cw = length(path.path))
        ISDC = ISDistributionCache(path, cw)
        CP = 0.
        CPvarN = 0.
        CPstd = 0.
        meanF = 0.
        meanH = 0.
        varFN = 0.
        varHN = 0.
        covFHN = 0.
        for i in 1:Nparticles
            verbose && (mod(i, progressmod) == 0) && println(i)
            p, lr = noisify_with_kick(path, alpha_normalized, seed=seed0+i, cw=cw, ISDC=ISDC)
            f = lr*(~is_free_path(p, CC0))
            h = lr*half_plane_breach_count(path, p)
            covFHN = covFHN + (i - 1) * (f - meanF) * (h - meanH) / i
            varFN = varFN + (i - 1) * (f - meanF)^2 / i
            varHN = varHN + (i - 1) * (h - meanH)^2 / i
            meanF = meanF + (f - meanF) / i
            meanH = meanH + (h - meanH) / i
            beta = covFHN / varHN
            CP = meanF - beta * (meanH - theta)
            CPvarN = max(varFN - 2*beta*covFHN + beta^2*varHN, 0.)   # numerical instability... crap
            CPstd = sqrt(CPvarN) / i
            targeted && i > 100 && abs(CP - CPgoal) > zhalt*CPstd && break
        end
    else
        error("Unsupported CP estimation method!")
    end
    {
        "P0" => P0,
        "CP" => CP,
        "CPstd" => CPstd,
        "path" => path,
        "Nparticles" => i,
        "plan_time" => plan_time,
        "cost" => P0.solution.cost,
        "MC_time" => toq()
    }
end

function binary_search_CP(P0::MPProblem, CPgoal::Float64, LQG::DiscreteLQG, Nparticles::Int=1000;
                          itermax = 20, lo::Float64 = 0., hi::Float64 = .04, reltol = .1, method::Symbol = :VR, verbose = false)
    tic()
    iter = 0
    local CPlo, CPhi, CPmid
    mid = 0.
    plan_time = 0.
    MC_time = 0.
    alphafilter = min(1e-5, CPgoal / Nparticles)    # not at all the right quantity, but it's something
    try
        CPlo = collision_probability(P0, lo, LQG, Nparticles, method = method, targeted = true, CPgoal = CPgoal, alphafilter = alphafilter)
        plan_time += CPlo["plan_time"]
        MC_time += CPlo["MC_time"]
    catch
        rethrow()
        error("No nominal solution at lo inflation $(lo)")
    end
    for i in 1:5
        try
            CPhi = collision_probability(P0, hi, LQG, Nparticles, method = method, targeted = true, CPgoal = CPgoal, alphafilter = alphafilter)
            plan_time += CPhi["plan_time"]
            MC_time += CPhi["MC_time"]
            CPhi["CP"] < CPgoal && break
            hi = lo + (hi-lo)*1.2   # super duper ad hoc
        catch
            hi = (lo+hi)/2
            if i > 4
                error("No nominal solution at hi inflation $(hi) (after $i halvings)")
            end
        end
    end
    if !(CPhi["CP"] < CPgoal < CPlo["CP"])    # eps is inversely related to CP, so CPhi < CPlo
        error("Initial bisection interval doesn't contain CPgoal: ", (CPhi["CP"], CPgoal, CPlo["CP"]))
    end
    verbose && @printf("Iteration %d: eps interval (%4f, %4f) CP interval (%4f, %4f, %4f) elapsed time %3fs\n", 
                       iter, lo, hi, CPhi["CP"], CPgoal, CPlo["CP"], toq())

    for iter in 1:itermax
        tic()
        iter += 1
        mid = (lo+hi)/2
        CPmid = collision_probability(P0, mid, LQG, Nparticles, method = method, targeted = true, CPgoal = CPgoal, alphafilter = alphafilter)
        plan_time += CPmid["plan_time"]
        MC_time += CPmid["MC_time"]
        verbose && @printf("Iteration %d: eps interval (%4f, %4f) CP interval (%4f, %4f, %4f) elapsed time %3fs\n", 
                           iter, lo, hi, CPhi["CP"], CPmid["CP"], CPlo["CP"], toq())
        abs(CPmid["CP"] - CPgoal) < reltol*CPgoal && break
        if (hi - lo) < 1e-4    #  also ad hoc; essentially all ad hoc decisions (including hi0)
            mid = hi           #  should be determined by noise characteristics
            CPmid = collision_probability(P0, mid, LQG, Nparticles, method = method, targeted = false)
            plan_time += CPmid["plan_time"]
            MC_time += CPmid["MC_time"]
            break
        end
        if CPmid["CP"] > CPgoal
            lo = mid
            CPlo = CPmid
        else
            hi = mid
            CPhi = CPmid
        end
    end

    alpha_ellipse = ellipsoid_breach_probabilities(CPmid["path"], P0.CC)
    {
        "ep" => mid,
        "CPdict" => CPmid,
        "plan_time" => plan_time,
        "MC_time" => MC_time,
        "cond_mult" => pointwise_pruned_uncertainty_CP_estimate(CPmid["path"], P0.CC),
        # "add" => sum(CPmid["alpha"]),
        # "mult" => (1-prod(1-CPmid["alpha"])),
        "adde" => sum(alpha_ellipse),
        "multe" => (1-prod(1-alpha_ellipse)),
        "CP" => CPmid["CP"],
        "CPstd" => CPmid["CPstd"],
        "cost" => CPmid["cost"]
    }
end

## Estimators

cummean(x) = cumsum(x) ./ [1:length(x)]
cumvar(x, m = cummean(x)) = cumsum((1 - 1 ./ [1:length(x)]).*(x - [0, m[1:end-1]]).^2) ./ [1:length(x)]
cumcov(x, y, mx = cummean(x), my = cummean(y)) = cumsum((1 - 1 ./ [1:length(x)]).*(x - [0, mx[1:end-1]]).*(y - [0, my[1:end-1]])) ./ [1:length(x)]

function CP_estimates(CPS)
    f = CPS["f"]
    h = CPS["h"]
    θ_true = CPS["theta"]
    ISf = CPS["ISf"]
    ISh = CPS["ISh"]

    Nparticles = length(f)
    p = cummean(f)
    θ = cummean(h)
    β = cumcov(f,h) ./ cumvar(h)
    p_β = p - β .* (θ - θ_true)
    stdp = sqrt(cumvar(f) ./ [1:Nparticles])
    stdp_β = sqrt(abs((cumvar(f) + β.^2 .* cumvar(h) - 2 * β .* cumcov(f,h)) ./ [1:Nparticles]))

    p_q = cummean(ISf)
    θ_q = cummean(ISh)
    β_q = cumcov(ISf,ISh) ./ cumvar(ISh)
    p_β_q = p_q - β_q .* (θ_q - θ_true)
    stdp_q = sqrt(cumvar(ISf) ./ [1:Nparticles])
    stdp_β_q = sqrt(abs((cumvar(ISf) + β_q.^2 .* cumvar(ISh) - 2 * β_q .* cumcov(ISf,ISh)) ./ [1:Nparticles]))

    {
        "Naive MC" => (p, stdp),
        "CV only" => (p_β, stdp_β),
        "IS only" => (p_q, stdp_q),
        "CV+IS" => (p_β_q, stdp_β_q)
    }
end

function plot_series_with_error(x, e, color, label)
    plt.plot([1:length(x)], x, color=color, label=label)
    fill_between([1:length(x)], x-e, x+e, alpha = .3, edgecolor="none", facecolor=color)
end

function plot_CP_stats(CPS; bw = .25)
    Nparticles = length(CPS["f"])
    CP = CP_estimates(CPS)
    plot_series_with_error(CP["Naive MC"]..., "red", "Naive MC")
    plot_series_with_error(CP["CV only"]..., "blue", "MC with CV")
    plot_series_with_error(CP["IS only"]..., "green", "MC with IS")
    plot_series_with_error(CP["CV+IS"]..., "black", "MC with CV+IS")
    # plt.plot([1:Nparticles], CP["Naive MC"][1], color="red", label="Naive MC")
    # plt.plot([1:Nparticles], CP["CV only"][1], color="blue", label="MC with CV")
    # plt.plot([1:Nparticles], CP["IS only"][1], color="green", label="MC with IS")
    # plt.plot([1:Nparticles], CP["CV+IS"][1], color="black", label="MC with CV+IS")
    plt.plot([1:Nparticles], CPS["theta"]*ones(Nparticles), color="orange", label="Additive")
    plt.plot([1:Nparticles], (1-prod(1-CPS["alpha"]))*ones(Nparticles), color="brown", label="Multiplicative")
    plt.plot([1:Nparticles], CPS["prunedCPestimate"]*ones(Nparticles), color=(.85,.85,0.), label="Cond Mult")
    axis([0, Nparticles, CP["CV+IS"][1][end] - .025*bw, CP["CV+IS"][1][end] + .025*bw])
    legend(loc="center left", bbox_to_anchor=(1, 0.5))
end