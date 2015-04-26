using Interact

function visualize_CP_evolution(cpds)
    flush(STDOUT)
    f = figure()
    @manipulate for i in slider(1:length(cpds), value=1)
        withfig(f) do
            plot_path_uncertainty_visualization(cpds[i]["P0"], cpds[i]["LQG"], cpds[i]["path"], cpds[i]["eps"])
            cp = cpds[i]["CP"]
            cpstd = cpds[i]["CPstd"]
            eps = cpds[i]["eps"]
            title(latexstring("CP = $(round(100*cp,2))\% \$\\pm\$ $(round(100*cpstd,2))\%, eps = $(round(eps,5))"))
        end
    end
end

function test_run(setup, CPgoal)
    P, DLQG, lo, hi = ISRR_test_problems(setup)
    isa(P.SS, RealVectorMetricSpace) && fmtstar!(P, 5000, connections = :R, rm = 1.5)
    setup == "DI2" && fmtstar!(P, 2500, connections = :R, r = 1.)
    setup == "DI3" && fmtstar!(P, 3500, connections = :R, r = 1.5)
    visualize_CP_evolution(binary_search_CP(P, CPgoal, DLQG, 500, lo = lo, hi = hi, verbose = true, vis = true)[2])
end