include("preamble.jl")

T = fill(0.0,1,1)
c = [1.0]
model = BoundedFluidQueue(T,c)

nodes = collect(0.0:1.0:10.0);

errors = DataFrame()
x_lims = (-0.1,3.0)

## First, lets examin closing vectors for QBDRAP 

for plt_type in (cdf,pdf)
    p = plot(layout=(3,4))
    linetypes = [:solid,:solid,:solid]
    c_order = 0
    for order in 1:2:7
        c_order += 1
        c_mesh = 0
        rec_types = (eval(Symbol(:unnormalised_closing_operator_,plt_type)), 
            eval(Symbol(:naive_normalised_closing_operator_,plt_type)),
            eval(Symbol(:normalised_closing_operator_,plt_type)))
        for reconstruction in rec_types
            c_mesh += 1
            mesh = FRAPMesh(nodes,order)

            dq = DiscretisedFluidQueue(model,mesh)

            d0 = interior_point_mass(1.0+eps(),1,dq)
            V = DiscretisedFluidQueues.vandermonde(order)

            rec = plt_type(d0, reconstruction)
            y_lim_vals = (plt_type==pdf) ? (-3.0,9.0) : (-0.2,1.4)
            y_tick_vals = ((plt_type==pdf) ? (y_lim_vals[1]:3:y_lim_vals[2]) : (y_lim_vals[1]:0.4:y_lim_vals[2]))
            std_plot(args...; kwargs...) = plot!(
                p.layout.grid[c_mesh,c_order],
                x->rec(x,1),
                x_lims[1], x_lims[2];
                ylim=y_lim_vals,
                linestyle=linetypes[c_mesh],
                xticks=false,
                yticks=false,
                grid=false,
                kwargs...,
            )
            if c_order == 1
                std_plot(;label=false,ylabel=uppercase(string(plt_type)),legend=false)
                yticks!(p.layout.grid[c_mesh,c_order], y_tick_vals)
            elseif c_order == 4
                std_plot(;label=string(reconstruction)[1:findfirst("_",string(reconstruction))[1]-1],legend=:outertopright)
            else
                std_plot(;label=false,legend=false)
            end
            if c_mesh == 1 
                plot!(p.layout.grid[c_mesh,c_order]; title="Order "*string(order))
            elseif c_mesh == 3
                xticks!(p.layout.grid[c_mesh,c_order], nodes)
                plot!(p.layout.grid[c_mesh,c_order]; xlabel="x")
            end
            std_plot2(args... ; kwargs...) = if plt_type==cdf
                plot!(
                    p.layout.grid[c_mesh,c_order],
                    [x_lims[1];1+eps();1+eps();x_lims[2]],
                    [0.0;0.0;1.0;1.0],
                    args...;
                    linestyle=:dot,
                    kwargs...,
                )
            elseif plt_type==pdf
                plot!(
                    p.layout.grid[c_mesh,c_order],
                    [x_lims[1];1+eps();1+2*eps();x_lims[2]],
                    [0.0;100;0.0;0.0],
                    args...;
                    subplot=c_order,
                    linestyle=:dot,
                    kwargs...
                )
            end
            if c_order==4
                std_plot2(;label="truth",legend=:outertopright)
            else
                std_plot2(;label=false)
            end
        end
    end
    display(p)
end

# for plt_type in (cdf,pdf)
#     p = plot(layout=(3,4))
#     linetypes = [:solid,:solid,:solid]
#     c_order = 0
#     for order in 1:2:7
#         c_order += 1
#         c_mesh = 0
#         for mtype in (DGMesh, FVMesh, FRAPMesh)
#             c_mesh += 1
#             mesh = mtype(nodes,order)

#             dq = DiscretisedFluidQueue(model,mesh)

#             d0 = interior_point_mass(1.0+eps(),1,dq)
#             V = DiscretisedFluidQueues.vandermonde(order)
#             # (mtype<:DGMesh)&&display(d0.coeffs)

#             rec = plt_type(d0)
#             y_lim_vals = (plt_type==pdf) ? (-3.0,9.0) : (-0.2,1.4)
#             y_tick_vals = ((plt_type==pdf) ? (y_lim_vals[1]:3:y_lim_vals[2]) : (y_lim_vals[1]:0.4:y_lim_vals[2]))
#             std_plot(args...; kwargs...) = plot!(
#                 p.layout.grid[c_mesh,c_order],
#                 x->rec(x,1),
#                 x_lims[1], x_lims[2];
#                 ylim=y_lim_vals,
#                 linestyle=linetypes[c_mesh],
#                 xticks=false,
#                 yticks=false,
#                 grid=false,
#                 kwargs...,
#             )
#             if c_order == 1
#                 std_plot(;label=false,ylabel=uppercase(string(plt_type)),legend=false)
#                 yticks!(p.layout.grid[c_mesh,c_order], y_tick_vals)
#             elseif c_order == 4
#                 std_plot(;label=string(mtype),legend=:outertopright)
#             else
#                 std_plot(;label=false,legend=false)
#             end
#             if c_mesh == 1 
#                 plot!(p.layout.grid[c_mesh,c_order]; title="Order "*string(order))
#             elseif c_mesh == 3
#                 xticks!(p.layout.grid[c_mesh,c_order], nodes)
#                 plot!(p.layout.grid[c_mesh,c_order]; xlabel="x")
#             end
#             std_plot2(args... ; kwargs...) = if plt_type==cdf
#                 plot!(
#                     p.layout.grid[c_mesh,c_order],
#                     [x_lims[1];1+eps();1+eps();x_lims[2]],
#                     [0.0;0.0;1.0;1.0],
#                     args...;
#                     linestyle=:dot,
#                     kwargs...,
#                 )
#             elseif plt_type==pdf
#                 plot!(
#                     p.layout.grid[c_mesh,c_order],
#                     [x_lims[1];1+eps();1+2*eps();x_lims[2]],
#                     [0.0;100;0.0;0.0],
#                     args...;
#                     subplot=c_order,
#                     linestyle=:dot,
#                     kwargs...
#                 )
#             end
#             if c_order==4
#                 std_plot2(;label="truth",legend=:outertopright)
#             else
#                 std_plot2(;label=false)
#             end
#         end
#     end
#     display(p)
# end