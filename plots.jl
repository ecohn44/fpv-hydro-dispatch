using LaTeXStrings
using Plots

function gen_plots(path, T, ps, PS, a, ph, PH, PF)
    
    ## Plot 1: Solar generation 
    plot1 = plot(1:T, ps, label="FPV Generation", lw=2, legend = :outertopright)
    # Add maximum solar threshold 
    hline!(plot1, [PS], color=:red, linestyle=:dash, label="Nameplate")
    # Add solar capacity 
    plot!(plot1, 1:T, a*PS, label = "Solar Cap.")
    xlabel!(plot1, "Hour")
    ylabel!(plot1, "Ps (MWh)")
    title!(plot1, "Hourly FPV Generation")
    savefig(plot1, path * "/FPV_gen.png");

    ## Plot 2: Hydropower generation 
    plot2 = plot(1:T, ph, xlabel="Hour", ylabel="Ph (MW)", title="Hydropower Generation and Available Capacity via Hydraulic Head", label="Hydropower Generation", lw=2,legend = :topright)
    # Add maximum threshold from hydraulic head
    plot!(plot2, 1:T, PH, color=:red, linestyle=:dash, label="Hydraulic Head Capacity")
    xlabel!(plot2, "Hour")
    ylabel!(plot2, "MWh")
    title!(plot2, "Hourly Hydropower Generation")
    savefig(plot2,  path * "/HP_gen.png");

    ## Plot 3: Hydropower generation (limit removed)
    plot3 = plot(1:T, ph, xlabel="Hour", ylabel="Ph (MW)", title="Hydropower Generation", label="Hydro", lw=2)
    xlabel!(plot3, "Hour")
    ylabel!(plot3, "Ph (MWh)")
    title!(plot3, "Hourly Hydropower Generation")
    savefig(plot3,  path * "/HP_gen_no_th.png");

    ## Plot 4: Joint generation 
    plot4 = plot(1:T, ph + ps, xlabel="Hour", ylabel="P (MW)", title="Total Generation", label="Generation", lw=2)
    hline!(plot4, [PF], color=:red, linestyle=:dash, label="Feeder Cap.")
    xlabel!(plot4, "Hour")
    ylabel!(plot4, "P (MWh)")
    title!(plot4, "Hourly System Generation")
    ylims!(plot4, 0*PF, 1.1*PF)
    savefig(plot4,  path * "/total_gen.png");
    
    plot(plot1, plot3, plot4, layout=(1, 3))

end

function release_plots(path, T, u, min, max)
    plot1 = plot(1:T, u, xlabel="Hour", ylabel="Release (m3)", title="Hourly Water Release", label="Water", lw=2, legend = :topright)
    # Add maximum th 
    hline!(plot1, [max], color=:red, linestyle=:dash, label="Max Release")
    # Add min th
    hline!(plot1, [min], color=:red, linestyle=:dash, label="Min Release")
    xlabel!(plot1, "Hour")
    ylabel!(plot1, "Water Release (m3)")
    title!(plot1, "Hourly Water Release")
    savefig(plot1,  path * "/release.png");
end

function ramp_rate(path, T, u)
    plot(u[2:end]-u[1:end-1], label = "Simulated Ramp Rate", legend = :outertopright)   
    hline!([RR_up], color=:red, linestyle=:dash, label="Max Ramp Rate")
    hline!([RR_dn], color=:red, linestyle=:dash, label="Min Ramp Rate")
    title!("Relaxed Ramp Rate for Simulated Policy")
    savefig(path * "/ramp_rate.png");
end 

function release_overlay_plots(path, T, u_b, u, min, max)
    # plot baseline
    plot1 = plot(1:T, u_b, xlabel="Hour", ylabel="Release (m3)", title="Hourly Water Release of Baseline and Relaxed Policy", label="Baseline Policy", lw=2, legend = :topright)
    # plot relaxed
    plot!(1:T, u, label="Relaxed Policy", color =:purple)
    # Add maximum th 
    hline!(plot1, [max], color=:red, linestyle=:dash, label="Max Release")
    # Add min th
    hline!(plot1, [min], color=:red, linestyle=:dash, label="Min Release")
    xlabel!(plot1, "Hour")
    ylabel!(plot1, "Water Release (m3)")
    title!(plot1, "Hourly Water Release")
    savefig(plot1,  path * "/release_overlay.png");
end

function generation_overlay_plots(path, T, p_b, p, type, max)
    # plot baseline
    plot1 = plot(1:T, p_b, label="Baseline Policy", lw=2, legend = :outertopright)
    # plot relaxed
    plot!(1:T, p, label="Relaxed Policy", color =:purple)
    # Add maximum th 
    hline!(plot1, [max], color=:red, linestyle=:dash, label="Capacity")
    xlabel!(plot1, "Hour")
    ylabel!(plot1, type*" Generation (MWh)")
    title!(plot1, "Comparison of " *type* " Generation Policies")
    savefig(plot1,  path * "/"*type*"_generation_overlay.png");
end

function release_plots_LMP(path, T, u, min, max, L)
    plot(1:T, u, xlabel="Hour", ylabel="Release (m3)", label="Water Release", lw=2,legend = :outertopright,title="Hourly Water Release vs Local Marginal Price",titlefont=font(12))
    # Add maximum th 
    hline!([max], color=:red, linestyle=:dash, label="Max Release Rate")
    # Add min th
    hline!([min], color=:red, linestyle=:dash, label="Min Release Rate")
    # Add LMP with different yaxis
    plot!(twinx(), 1:T, L, ylabel="\$/MWh", label="Hourly Marginal Price", color=:green,legend = :outerbottomright)
    savefig(path * "/LMPrelease.png");
end

function gen_plots_LMP(path, T, p, max, L)
    plot(1:T, p, xlabel="Hour", ylabel="Total Generation (MWh)", label="Power Generation", lw=2,legend = :outertopright,title="Total System Generation vs Local Marginal Price",titlefont=font(12))
    # Add maximum th 
    hline!([max], color=:red, linestyle=:dash, label="Feeder Capacity")
    # Add LMP with different yaxis
    plot!(twinx(), 1:T, L, ylabel="\$/MWh", label="Hourly Marginal Price", color=:green,legend = :outerbottomright)
    savefig(path * "/LMPgeneration.png");
end
    
function flow_plot(path, T, Q)
    plot1 = plot(1:T, Q, lw=2, legend = :outertopright, label = "Inflow")
    xlabel!(plot1, "Day")
    ylabel!(plot1, "Inflow (m3)")
    title!(plot1, "Monthly System Inflow")
    savefig(plot1,  path * "/inflow.png");
end

function LMP_plot(path, T, L)
    plot1 = plot(1:T, L, lw=2, legend = :outertopright, label = "LMP")
    xlabel!(plot1, "Hour")
    ylabel!(plot1, "Marginal Price (\$/MWh)")
    title!(plot1, "Hourly LMP")
    savefig(plot1,  path * "/LMP.png");
end

function policy_plot(path, T, u_pi)
    # Plot water release from dual-value driven policy
    plot(1:T, u_pi, xlabel="Hour", ylabel="Policy Water Release (m3/s)", label="Lagrangian Policy", lw=2,legend = :topright,title="Lagrangian Water Release Simulations",titlefont=font(12))
    savefig(path * "/policy_waterrelease.png");
end

function overlay_policy_plot(path, T, u_pi, u_star)
    # Plot water release from dual-value driven policy
    plot(1:T, u_pi, xlabel="Hour", ylabel="Policy Water Release (m3/s)", label="Lagrangian Policy", lw=2,legend = :outertopright,title="Optimal and Lagrangian Water Release Simulations",titlefont=font(12))
    # Add optimal water release with different yaxis
    plot!(twinx(), 1:T, u_star, ylabel="Optimal Water Release (m3/s)", label="Optimal Policy", color=:green, legend = :outerbottomright)
    savefig(path * "/policy_waterrelease_overlay.png");
end

function overlay_policy_plot_solar(path, T, u_pi, ps_star)
    # Plot water release from dual-value driven policy
    plot(1:T, u_pi, xlabel="Hour", ylabel="Policy Water Release (m3/s)", label="Lagrangian Policy", lw=2,legend = :outertopright,title="Lagrangian Water Release & FPV Generation",titlefont=font(12))
    # Add optimal water release with different yaxis
    plot!(twinx(), 1:T, ps_star, ylabel="Optimal FPV Generation (MWh)", label="FPV Generation", color=:green,legend = :outerbottomright)
    savefig(path * "/policy_waterrelease_solar.png");
end

function duals_plot(path, T, x, dv_name, eq_name)
    dv_print = replace(dv_name, r"$|\\|[{}]" => "");
    println(dv_print)
    plot1 = plot(1:T, x, lw=2, legend = :topright, label = dv_name, size=(850,800))
    xlabel!(plot1, "Hour")
    ylabel!(plot1, "Dual Value " * dv_name)
    title!(plot1, eq_name * " Dual Value over Historical Simulation")
    savefig(plot1,  path * "/" * dv_print * "duals.png");
end

function iters_plot(path, T, x, dv_name)
    dv_print = replace(dv_name, r"$|\\|[{}]" => "");
    println(dv_print)
    plot1 = plot(1:T, x, lw=2, legend = :topright, label = dv_name * ": Water Contract", size=(850,800))
    xlabel!(plot1, "Iteration (k)")
    ylabel!(plot1, "Dual Value " * dv_name)
    title!(plot1, "Convergence of the Dual Value over BST Algorithm")
    savefig(plot1,  path * "/" * dv_print * "duals.png");
end

function hhead_plots(path, T, hh, hh_d)
    
    ## Plot 1: Hydraulic Head
    plot1 = plot(1:T, hh, label="Hydraulic Head", lw=2, legend = :topright)
    xlabel!(plot1, "Hour")
    ylabel!(plot1, "Hydraulic Head [m]")
    title!(plot1, "Available Hydraulic Head [m]")
    savefig(plot1, path * "/hh.png");

    ## Plot 2: Hydraulic Head Derivative
    plot2 = plot(1:T, hh_d, label="Derivative of Hydraulic Head", lw=2, legend = :topright)
    xlabel!(plot2, "Hour")
    ylabel!(plot2, "Derivative of Hydraulic Head")
    title!(plot2, "Derivative of Hydraulic Head")
    savefig(plot2, path * "/hh_d.png");

end

function netflow_plot()
    mode1 = 837620904/1e6
    mode2 = 1889976499/1e6
    delta1 = mode1 - 0.05*mode1
    delta2 = mode1 + 0.05*mode1

    net = water_contract - total_inflow
    plot(net/1e6, label = "Netflow")
    hline!([mode1], color=:red, linestyle=:dash, label="Max Release")
    hline!([delta1], color=:red, linestyle=:dash, label="+5% Max Release")
    hline!([delta2], color=:red, linestyle=:dash, label="-5% Max Release")
    title!("System Net Flow")
    xaxis!("Month")
    yaxis!("Volume (Million m3)")
end

function inflow_plot()
    #plot(total_inflow/1e6, label = "Inflow")
    plot(water_contract/1e6, label = "Water Contract (Outflow)")
    hline!([max_release/1e6], color=:red, linestyle=:dash, label="Max Release")
    title!("Total Monthly Water Contract Vs Inflow")
    xaxis!("Month")
    yaxis!("Volume (Million m3)")
end