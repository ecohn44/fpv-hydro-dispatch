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
