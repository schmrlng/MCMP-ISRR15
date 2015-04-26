using PyPlot
using MotionPlanning
using ImmutableArrays
using JSON
!isdefined(:ISRR_POLY) && include(Pkg.dir("MotionPlanning")*"/test/obstaclesets/2D.jl")
!isdefined(:BOXES3D) && include(Pkg.dir("MotionPlanning")*"/test/obstaclesets/ND.jl")
include("discreteLQG.jl")
include("collisionprobability.jl")
include("printhelpers.jl")

function double_integrator_noise(d)
    srand(1)
    V0 = rand(d,d) |> x-> x'*x
    V1 = rand(d,d) |> x-> x'*x
    W0 = rand(d,d) |> x-> x'*x
    W1 = rand(d,d) |> x-> x'*x
    srand("/dev/random")   # unix only in julia 0.3
    pos_v = .05^2
    vel_v = .05^2
    Vc = [pos_v*V0 zeros(d,d); zeros(d,d) vel_v*V1]
    # V = diagm([pos_v^2*ones(d), vel_v^2*ones(d)])
    pos_w = .01^2
    vel_w = .005^2
    # W = diagm([pos_w^2*ones(d), vel_w^2*ones(d)])
    Wc = [pos_w*W0 zeros(d,d); zeros(d,d) vel_w*W1]
    P0 = diagm([pos_v^2*ones(d), vel_v^2*ones(d)])

    Vc, Wc, P0
end

!isdefined(:DOUBLE_INTEGRATOR_2D) && (const DOUBLE_INTEGRATOR_2D = DoubleIntegrator(2, vmax = 0.5))
!isdefined(:DOUBLE_INTEGRATOR_3D) && (const DOUBLE_INTEGRATOR_3D = DoubleIntegrator(3, vmax = 0.5))

function ISRR_test_problems(setup = "SI2GOOD", SIdt = 0.015, DIdt = .05)
    if setup == "SI2GOOD"
        P = MPProblem(UnitHypercube(2),
                      Vector2(.1,.1),
                      PointGoal([.9, .9]),
                      PointRobot2D(ISRR_POLY))
        DLQG = SingleIntegrator(2, nsf=0.4, dt = SIdt)
        lo, hi = 0.001, 0.04
    elseif setup == "SI2BAD"
        P = MPProblem(UnitHypercube(2),
                      Vector2(.1,.1),
                      PointGoal([.9, .9]),
                      PointRobot2D(ISRR_POLY_WITH_SPIKE))
        DLQG = SingleIntegrator(2, nsf=0.4, dt = SIdt)
        lo, hi = 0.001, 0.04
    elseif setup == "DI2"
        P = MPProblem(DOUBLE_INTEGRATOR_2D,
                      [.1, .1, 0., 0.],
                      PointGoal([.9, .9, 0., 0.]),
                      PointRobot2D(ISRR_POLY))
        P.SS.dist.cmax = 1.
        DLQG = DiscreteLQG(P.SS, double_integrator_noise(2)..., nsf=0.6, dt = DIdt)
        lo, hi = 0.001, 0.04
    elseif setup == "SI3"
        P = MPProblem(UnitHypercube(3),
                      [.1,.1,.1],
                      PointGoal([.9,.5,.1]),
                      PointRobotNDBoxes(BOXES3D))
        DLQG = SingleIntegrator(3, nsf=0.5, dt = SIdt)
        lo, hi = 0.001, 0.04
    elseif setup == "DI3"
        P = MPProblem(DOUBLE_INTEGRATOR_3D,
                      [.1,.1,.1,0.,0.,0.],
                      PointGoal([.9,.5,.1,0.,0.,0.]),
                      PointRobotNDBoxes(BOXES3D))
        P.SS.dist.cmax = 1.5        # TODO: better way of integrating with FMT (really, solver should set this)
        DLQG = DiscreteLQG(P.SS, double_integrator_noise(3)..., nsf=0.6, dt = .1)
        lo, hi = 0.001, 0.03
    end
    P, DLQG, lo, hi
end

function run_tests(setup, CPgoal, M=10, N=20; verbose=true, writefile=false)
    P, DLQG, lo, hi = ISRR_test_problems(setup)
    plan_cache_times = Float64[]
    results = {}

    for m in 1:M
        println("======== $(setup) Sample Set $(m) ========")
        P.V = defaultNN(P.SS, P.init)
        tic()
        isa(P.SS, RealVectorMetricSpace) && fmtstar!(P, 5000, connections = :R, rm = 1.5)
        setup == "DI2" && fmtstar!(P, 2500, connections = :R, r = 1.)
        setup == "DI3" && fmtstar!(P, 3500, connections = :R, r = 1.5)
        push!(plan_cache_times, toq())
        println("Planner Cache Time: $(plan_cache_times[end])")

        for n in 1:N
            println("=> RUN $n")
            push!(results, binary_search_CP(P, CPgoal, DLQG, 500, lo = lo, hi = hi, verbose = verbose)[1])
        end
    end
    P.V = defaultNN(P.SS, P.init)   # might help with a memory leak?

    ret = plan_cache_times, results
    if writefile
        timestr = strftime("%Y-%m-%d_%H:%M:%S", time())
        s = open("DATA/"*setup*"_$(timestr)", "a")
        write(s, "$(JSON.json(ret))")
        flush(s)
    end
    ret
end

if !isinteractive()
    run_tests(ARGS[1], parsefloat(ARGS[2]), 1, 1, verbose = false, writefile = false)
    run_tests(ARGS[1], parsefloat(ARGS[2]), parseint(ARGS[3]), parseint(ARGS[4]), writefile = true)
end