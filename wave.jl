include("preamble.jl")
import Distributions 

T = [0.0 0.0; 1.0 -1.0]
c = [1.0; 0.0]
model = BoundedFluidQueue(T,c)

nodes = 0.0:1.0:10.0

x_lims = (-0.1,10.1)

norm_pdf(x) = Distributions.pdf(Distributions.Normal(2.5,0.75),x)
norm_cdf(x) = Distributions.cdf(Distributions.Normal(2.5,0.75),x)
## First, lets examin closing vectors for QBDRAP 
orders = 1:2:21
CDF0s = [
    (z,i)->Float64(z>=0.0)*(i==1)*(z>=0.0)
    (z,i)->Float64(z>0.5)*(i==1)*(z>=0.0);
    (z,i)->((z<0.5)*2*z+(z>=0.5))*(i==1)*(z>=0.0);
    # (z,i)->(z>0.5)*(z^2/2-z/2+1/8)*8;
    (z,i)->(z*(z<1.0) + Float64(z>=1.0))*(i==1)*(z>=0.0);
    (z,i)->(z^3/10^3)*(i==1)*(z>=0.0);
    # (z,i)->-2*z^3+3*z^2;
    (z,i)->((norm_cdf(z)-norm_cdf(0.0))/(norm_cdf(10.0)-norm_cdf(0.0)))*(i==1)*(z>=0.0);
    (z,i)->(2*(z/2 - sin(2*pi*z)/(4*pi))*(z<1.0) + Float64(z>=1.0))*(i==1)*(z>=0.0);
    (z,i)->((2*(z/2 - sin(2*pi*z)/(4*pi))+z)/10.0)*(i==1)*(z>=0.0);
    (z,i)->cdf(ConcentratedMatrixExponential(5,mean=2.5))(z)*(i==1)*(z>=0.0)
    ]
PDF0s = [
    (z,i)->0.0
    (z,i)->NaN;
    (z,i)->2*Float64(z<0.5)*(i==1)*(z>=0.0);
    # (z,i)->(z>0.5)*8*(z-0.5);
    (z,i)->1.0*(z<1.0)*(i==1)*(z>=0.0);
    (z,i)->z^2*3/10^3*(i==1)*(z>=0.0);
    # (z,i)->-6*z^2+6*z;
    (z,i)->(norm_pdf(z)/(norm_cdf(10.0)-norm_cdf(0.0)))*(i==1)*(z>=0.0);
    (z,i)->((cos(2*π*(z+0.5))+1)*(z<1.0))*(i==1)*(z>=0.0);
    (z,i)->(cos(2*π*(z+0.5))+1)/10*(i==1)*(z>=0.0);
    (z,i)->pdf(ConcentratedMatrixExponential(5,mean=2.5))(z)*(i==1)*(z>=0.0)
    ]
CDFs = [(z,j)->CDF0s[i](z-4.0,j) for i in 1:length(CDF0s)]
PDFs = [(z,j)->PDF0s[i](z-4.0,j) for i in 1:length(PDF0s)]

init_distns = [x->left_point_mass(2,x);
         x->interior_point_mass(0.5,1,x);
    [ x->SFMDistribution(PDF0s[i],x,TrapezoidRule;fun_evals=2001) for i in 3:length(PDFs)]]
    # x->SFMDistribution(PDFs[2],x,TrapezoidRule;fun_evals=2001);
    # x->SFMDistribution(PDFs[3],x,TrapezoidRule;fun_evals=2001);
    # [ x->SFMDistribution(PDFs[i],x) for i in 4:length(PDFs)]]
# x->SFMDistribution((z,i)->1.0,x),
# x->SFMDistribution((z,i)->2.0*z,x),
# x->SFMDistribution((z,i)->-6*z^2+6*z,x),
# x->SFMDistribution((z,i)->3*exp(-3*z)/(1-exp(-3)),x),
# x->SFMDistribution((z,i)->cos(2*π*(z+0.5))+1,x),
# x->SFMDistribution((z,i)->cos(4*π*(z+0.5))+1,x),
# x->SFMDistribution((z,i)->pdf(ConcentratedMatrixExponential(5,mean=0.5))(z),x),
# ]
m_loop_vec = (DGMesh,DGMesh,FRAPMesh)
for (func_count,init_dist_fun) in enumerate(init_distns)
    pth = mkpath((@__FILE__)[1:end-3]*"/func_count_"*string(func_count)*"/mesh_comp")
    mkpath(pth*"/data")
    mkpath(pth*"/figs")
    
    L1_cdf_errors = DataFrame(DG = Float64[], DGₗ = Float64[], QBDRAP = Float64[])
    ks_errors = DataFrame(DG = Float64[], DGₗ = Float64[], QBDRAP = Float64[])
    L1_pdf_errors = DataFrame(DG = Float64[], DGₗ = Float64[], QBDRAP = Float64[])
    L2_pdf_errors = DataFrame(DG = Float64[], DGₗ = Float64[], QBDRAP = Float64[])
    for plt_type in (cdf,pdf)
        p = plot(layout=(length(m_loop_vec),length(orders)),legend=false)
        linetypes = [:solid,:solid,:solid]
        c_order = 0
        for order in orders
            c_order += 1
            c_mesh = 0

            L1_cdf_error_row = zeros(1,length(m_loop_vec))
            ks_error_row = zeros(1,length(m_loop_vec))
            L1_pdf_error_row = zeros(1,length(m_loop_vec))
            L2_pdf_error_row = zeros(1,length(m_loop_vec))
            for mtype in m_loop_vec
                c_mesh += 1
                mesh = mtype(nodes,order)

                dq = DiscretisedFluidQueue(model,mesh)
                B = build_full_generator(dq)

                d0 = init_dist_fun(dq)

                dt = ((mtype==DGMesh) && (c_mesh==1)) ? (limit_str=" NoLimit ";integrate_time(d0,B,4.0,StableRK4(0.01); limiter=NoLimiter)) : (limit_str=" ";integrate_time(d0,B,4.0,StableRK4(0.01)))

                rec = (mtype==FRAPMesh) ? plt_type(dt, eval(Symbol(:normalised_closing_operator_,plt_type))) : plt_type(dt)

                ## errors 
                (plt_type==cdf) && (L1_cdf_error_row[c_mesh] = log10(DiscretisedFluidQueues.Lp(x->rec(x,1),x->CDFs[func_count](x,1),range(-eps(),10+eps(),length=2001))))
                (plt_type==cdf) && (ks_error_row[c_mesh] = log10(DiscretisedFluidQueues.kolmogorov_smirnov(x->rec(x,1),x->CDFs[func_count](x,1),range(-eps(),10+eps(),length=2001))))
                (plt_type==pdf) && (L1_pdf_error_row[c_mesh] = log10(DiscretisedFluidQueues.Lp(x->rec(x,1),x->PDFs[func_count](x,1),range(-eps(),10+eps(),length=2001))))
                (plt_type==pdf) && (L2_pdf_error_row[c_mesh] = log10(DiscretisedFluidQueues.Lp(x->rec(x,1),x->PDFs[func_count](x,1),range(-eps(),10+eps(),length=2001),2)))
                
                ## plotting...
                y_lim_vals = (plt_type==pdf) ? (-2.0,5.0) : (-0.2,1.4)
                y_tick_vals = ((plt_type==pdf) ? (y_lim_vals[1]:3:y_lim_vals[2]-1) : (0:0.5:1))
                std_plot(args...; kwargs...) = plot!(
                    p.layout.grid[c_mesh,c_order],
                    x->rec(x,1),
                    nodes[1]+eps(), nodes[end]+eps();
                    ylim=y_lim_vals,
                    linestyle=linetypes[c_mesh],
                    xticks=false,
                    yticks=false,
                    grid=false,
                    legend=false,
                    kwargs...,
                )
                if c_order == 1
                    std_plot(;label=false,ylabel=string(mtype)[1:findfirst("M",string(mtype))[1]-1],legend=false)
                    yticks!(p.layout.grid[c_mesh,c_order], y_tick_vals)
                elseif c_order == length(orders)
                    std_plot(;legend=false)#:outertopright)
                else
                    std_plot(;label=false,legend=false)
                end
                if c_mesh == 1 
                    plot!(p.layout.grid[c_mesh,c_order]; title=string(order),legend=false)
                elseif c_mesh == length(m_loop_vec)
                    xticks!(p.layout.grid[c_mesh,c_order], 0:3:11)
                    plot!(p.layout.grid[c_mesh,c_order]; xlabel="x",legend=false)
                end
                std_plot2(args... ; kwargs...) = if plt_type==cdf
                    plot!(
                        p.layout.grid[c_mesh,c_order],
                        x->CDFs[func_count](x,1),
                        nodes[1], nodes[end],
                        args...;
                        linestyle=:dot,
                        legend=false,
                        kwargs...,
                    )
                elseif plt_type==pdf
                    plot!(
                        p.layout.grid[c_mesh,c_order],
                        x->PDFs[func_count](x,1),
                        nodes[1], nodes[end],
                        args...;
                        subplot=c_order,
                        linestyle=:dot,
                        legend=false,
                        kwargs...
                    )
                end
                if c_order==length(orders)
                    std_plot2(;legend=false)#:outertopright)
                else
                    std_plot2(;label=false)
                end
            end
            (plt_type==cdf) && push!(L1_cdf_errors,L1_cdf_error_row)
            (plt_type==cdf) && push!(ks_errors,ks_error_row)
            (plt_type==pdf) && push!(L1_pdf_errors,L1_pdf_error_row)
            (plt_type==pdf) && push!(L2_pdf_errors,L2_pdf_error_row)
        end
        #display(p)
        # display(pth)
        savefig(p,pth*"/figs/mesh_"*string(plt_type)*"_func_count_"*string(func_count)*".svg")
    end
    # #display(L1_cdf_errors)
    file = pth*"/data/meshs_l1_cdf_"*"func_count_"*string(func_count)
    CSV.write(file*".csv",L1_cdf_errors)
    q = plot()
    for reconstruction in names(L1_cdf_errors)
        plot!(orders,L1_cdf_errors[:,reconstruction],label=reconstruction)
    end
    plot!(xlabel="Order", ylabel="log₁₀ Error", 
        title="L¹ error between the true CDF and approximations",
        legend=:right)
    #display(q)
    # #display(pth)
    file = pth*"/figs/meshs_l1_cdf_"*"func_count_"*string(func_count)
    savefig(q,file*".svg")

    # #display(ks_errors)
    file = pth*"/data/meshs_ks_"*"func_count_"*string(func_count)
    CSV.write(file*".csv",ks_errors)
    q = plot()
    for reconstruction in names(ks_errors)
        plot!(orders,ks_errors[:,reconstruction],label=reconstruction)
    end
    plot!(xlabel="Order", ylabel="log₁₀ Error", 
        title="KS error between the true CDF and approximations",
        legend=:right)
    #display(q)
    file = pth*"/figs/meshs_ks_"*"func_count_"*string(func_count)
    savefig(q,file*".svg")

    # #display(L1_pdf_errors)
    file = pth*"/data/meshs_l1_pdf_"*"func_count_"*string(func_count)
    CSV.write(file*".csv",L1_pdf_errors)
    q = plot()
    for reconstruction in names(L1_pdf_errors)
        plot!(orders,L1_pdf_errors[:,reconstruction],label=reconstruction)
    end
    plot!(xlabel="Order", ylabel="log₁₀ Error", 
        title="L¹ error between the true PDF and approximations",
        legend=:right)
    #display(q)
    file = pth*"/figs/meshs_l1_pdf_"*"func_count_"*string(func_count)
    savefig(q,file*".svg")

    # #display(L2_pdf_errors)
    file = pth*"/data/meshs_l2_pdf_"*"func_count_"*string(func_count)
    CSV.write(file*".csv",L2_pdf_errors)
    q = plot()
    for reconstruction in names(L2_pdf_errors)
        plot!(orders,L2_pdf_errors[:,reconstruction],label=reconstruction)
    end
    plot!(xlabel="Order", ylabel="log₁₀ Error", 
        title="L² error between the true PDF and approximations",
        legend=:right)
    #display(q)
    file = pth*"/figs/meshs_l2_pdf_"*"func_count_"*string(func_count)
    savefig(q,file*".svg")

end