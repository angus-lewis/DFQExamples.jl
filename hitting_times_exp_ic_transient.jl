include("preamble.jl")

include("hitting_times_model/reflecting_model/model_def.jl")
#model.P_lwr[1], model.P_lwr[2] = [0.0 1.0]
#model.P_upr[1], model.P_upr[2] = [1.0 0.0]
function exp_initial_condition(x,i)
    return (i∈(1:2))*exp(-x)/(1-exp(-model.b))/2
end

t = 2.0
function run_dg(t)
    n_bases_per_cell = 21
    nodes = 0.0:10.0
    mesh = DGMesh(nodes, n_bases_per_cell)

    dq = DiscretisedFluidQueue(model,mesh)
    B = build_full_generator(dq)

    d0 = SFMDistribution(exp_initial_condition, dq)

    t_step = 0.005
    t_integration_scheme = StableRK4(t_step)

    dt = integrate_time( d0, B, t, t_integration_scheme; limiter=NoLimiter)
    return dt
end

n_bases_per_cell = 15
nodes = 0.0:10.0
mesh = DGMesh(nodes, n_bases_per_cell)

dq = DiscretisedFluidQueue(model,mesh)
B = build_full_generator(dq)

d0 = SFMDistribution(exp_initial_condition, dq)# TrapezoidRule, fun_evals=10_001)

dt = run_dg(t)
#dt.coeffs[:] = d0.coeffs'*exp(B*2)
dt_cdf = cdf(dt)

x = range(eps(),10-10*eps(),length=10001)
dt_cdf_eval = dt_cdf.( x, (1:2)')

sum(abs.(dt_cdf_eval[:,1].-sim[1]./sim[1][end]*dt_cdf_eval[end,1]))*10/length(x) + sum(abs.(dt_cdf_eval[:,2].-sim[2]./sim[2][end]*dt_cdf_eval[end,2]))*10/length(x)

## 

sim = JSON.parsefile(
    "/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_model/"*
    "reflecting_model/transient_distribution/exp/data/sims/cdf_evaluated.json"
)

plot!( x, sim[2]-dt_cdf_eval[:,2])
plot!( x, sim[1]-dt_cdf_eval[:,1])

tau = fixed_time(t)

rng = StableRNGs.StableRNG(16092021)#StableRNGs.StableRNG(1)
function initial_condition_generator(n, rng)
    i = rand(rng,1:2,n) 
    X = Vector{Float64}(undef,n)
    for k in 1:n
        x = Inf
        while x>=model.b
            x = -log(rand(rng))
        end
        X[k]=x
    end
    # X=5.0*ones(n)#10.0*rand(rng,n)
    return (X=X,φ=i)
end
n_sims = 5_000_000
sim_initial_condition = initial_condition_generator(n_sims, rng)

sims = simulate(model, tau, sim_initial_condition)

const idx1 = sims.φ.==1
const idx2 = sims.φ.==2

sim_cdf1(x) = sum(sims.X[idx1].<x)./n_sims
sim_cdf2(x) = sum(sims.X[idx2].<x)./n_sims

cdf1_eval = sim_cdf1.(x)
cdf2_eval = sim_cdf2.(x)

const idx1 = sim_initial_condition.φ.==1
const idx2 = sim_initial_condition.φ.==2

n1 = sum(idx1)
n2 = sum(idx2)

for i in 1:n_sims
    if n1>n2
        if idx1[i] 
            idx1[i]=0
            n1 -= 1
        end
    elseif n2>n1 
        if idx2[i]
            idx2[i]=0
            n2 -= 1
        end
    else
        break
    end
end

n_new = n1+n2

sim_cdf1(x) = sum(sim_initial_condition.X[idx1].<x)./n_new
sim_cdf2(x) = sum(sim_initial_condition.X[idx2].<x)./n_new

cdf1_eval = sim_cdf1.(x)
cdf2_eval = sim_cdf2.(x)

plot( x, cdf2_eval+cdf1_eval-dt_cdf_eval[:,2]-dt_cdf_eval[:,1],label="2")
plot!( x, cdf1_eval-dt_cdf_eval[:,1],label="1")

sum(abs.(dt_cdf_eval[:,1].-cdf1_eval))*10/length(x) + sum(abs.(dt_cdf_eval[:,2].-cdf2_eval))*10/length(x)

x2 = (x[1:end-1]+x[2:end])./2
plot!( x2, diff(cdf2_eval)*length(x)/10-pdf(d0).(x2,2),label="2")
plot!( x2, diff(cdf1_eval)*length(x)/10-pdf(d0).(x2,1),label="1")

plot( x2, diff(sim_cdf2.(x))*length(x)/10-pdf(dt).(x2,2),label="2")
plot!( x2, diff(sim_cdf1.(x))*length(x)/10-pdf(dt).(x2,1),label="1")
#plot( x, log10.(abs.(cdf(d0).(x,1)-exp_ic_cdf.(x, 1))))
#plot!( x, log10.(abs.(cdf(d0).(x,2)-exp_ic_cdf.(x, 2))))

exp_ic_cdf(x,i) = (1-exp(-x))/(1-exp(-model.b))/2

2*sum(abs.(cdf(d0).(x,1)-exp_ic_cdf.(x, 1)))*10.0/10_001
2*sum(abs.(pdf(d0).(x,1)-exp_initial_condition.(x, 1)))*10.0/10_001

sum(abs.(dt_cdf_eval[:,1]-cdf1_eval))*10.0/10_001+sum(abs.(dt_cdf_eval[:,2]-cdf2_eval))*10.0/10_001
sum(abs.(dt_cdf_eval[:,1]-exp_ic_cdf.(x,1)))*10.0/10_001+
    sum(abs.(dt_cdf_eval[:,2]-exp_ic_cdf.(x,2)))*10.0/10_001

plot(x,abs.(dt_cdf_eval[:,1]-exp_ic_cdf.(x,1)))
plot!(x,abs.(dt_cdf_eval[:,2]-exp_ic_cdf.(x,2)))

plot(x,dt_cdf_eval[:,1])
plot!(x,dt_cdf_eval[:,2])

plot((x[1:end-1]+x[2:end])./2,diff(cdf1_eval)*10_001/10)
plot!(x,pdf(dt).(x,1))
plot!(x,pdf(dt).(x,2))


## qbd
function run1(t,n)
    n_bases_per_cell = n
    nodes = range(0,10,length=11)
    mesh = FRAPMesh(nodes, n_bases_per_cell)
    @show n_bases_per_phase(mesh)

    dq = DiscretisedFluidQueue(model,mesh)
    B = build_full_generator(dq)
    #display(Matrix(B))

    d0 = SFMDistribution(exp_initial_condition, dq, TrapezoidRule; fun_evals = 20_001)
    t_step = diff(zglj(n,0,0))[1]/2/2 # smallest diff between nodes on [-1,1], divided by 2 to get to [0,1], then divided by 2 to be safe
    t_integration_scheme = StableRK4(t_step)
    dt = integrate_time( d0, B, t, t_integration_scheme)
    return dt
end
function run2(t,n)
    mesh2 = DGMesh(range(0,10,length=11), n)
    @show n_bases_per_phase(mesh2)

    dq2 = DiscretisedFluidQueue(model,mesh2)
    B2 = build_full_generator(dq2)

    #display(Matrix(B2))
    d02 = SFMDistribution(exp_initial_condition, dq2)
    t_step = diff(zglj(n,0,0))[1]/2/2 # smallest diff between nodes on [-1,1], divided by 2 to get to [0,1], then divided by 2 to be safe
    t_integration_scheme = StableRK4(t_step)
    dt2 = integrate_time( d02, B2, t, t_integration_scheme; limiter = NoLimiter)
    return dt2
end
function run_all_qbd()
    plot(layout=(2,1))
    dt_pdf_eval_qbd = zeros(10_001,2)
    dt_pdf_eval_dg = zeros(10_001,2)
    for n in 49
        tqbd = 2.0
        
        dt_dg, dt_qbd = run1(tqbd,n), run2(tqbd,n)

        dt_cdf_qbd = cdf(dt_qbd)

        dt_cdf_dg = cdf(dt_dg)
        dt_pdf_qbd = pdf(dt_qbd)

        dt_pdf_dg = pdf(dt_dg)

        x = range(eps(),10-10*eps(),length=10_001)
        dt_pdf_eval_qbd .= dt_pdf_qbd.( x, (1:2)')
        dt_pdf_eval_dg .= dt_pdf_dg.( x, (1:2)')

        println(n)
        println(sum(abs.(dt_pdf_eval_qbd[:,1].-dt_pdf_eval_dg[:,1]))*10/length(x) + 
            sum(abs.(dt_pdf_eval_qbd[:,2].-dt_pdf_eval_dg[:,2]))*10/length(x))
        # display(plot!((abs.(dt_pdf_eval_qbd[:,1].-dt_pdf_eval_dg[:,1])) + 
        #     (abs.(dt_pdf_eval_qbd[:,2].-dt_pdf_eval_dg[:,2]))))
        display(plot!(dt_pdf_eval_dg, subplot=1))
        display(plot!(dt_pdf_eval_qbd, subplot=2))
    end
    return dt_pdf_eval_qbd, dt_pdf_eval_dg
end

qbd_d, dg_d = run_all_qbd()

n=25

plot(x,dt_cdf_eval_qbd.-[sim[1][:] sim[2][:]])
plot!(x,dt_cdf_eval.-[sim[1][:] sim[2][:]])


exp_err(p::Float64, n) = sqrt(2/pi)*sqrt(p*(1-p)/n) #sqrt(1/(2*pi*(1-p)*p*n)) 

n = 5_000_000

sims_pm_transient = JSON.parsefile(
    "/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_model/reflecting_model/"*
    "transient_distribution/point_mass/data/sims/cdf_evaluated.json"
)

p_pm_transient = sims_pm_transient[1][end]

exp_err(p_pm_transient,n)

sims_exp_transient = JSON.parsefile(
    "/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_model/"*
    "reflecting_model/transient_distribution/exp/data/sims/cdf_evaluated.json"
)

p_exp_transient = sims_exp_transient[1][end]

exp_err(p_exp_transient,n)

sims_exp_hitting = CSV.read(
    "/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_model/hitting_times/data/exp/sims.csv",
    DataFrame
)

p_exp_hitting = sum(sims_exp_hitting.φ.==1)/nrow(sims_exp_hitting)

exp_err(p_exp_hitting,n)

sims_pm_hitting = CSV.read(
    "/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/hitting_times_model/hitting_times/data/point_mass/sims.csv",
    DataFrame
)

p_pm_hitting = sum(sims_pm_hitting.φ.==1)/nrow(sims_pm_hitting)

exp_err(p_pm_hitting,n)

sims_first_return_discont = CSV.read(
    "/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/fluidfluid_discontinuous/data/sims_cdf_evaluated.csv",
    DataFrame
)

p_first_return_discont = sims_first_return_discont[end,end]

exp_err(p_first_return_discont,n)

sims_first_return_cont = CSV.read(
    "/Users/anguslewis/Documents/SFFMProject/DFQExamples.jl/fluidfluid_continuous/data/sims_cdf_evaluated.csv",
    DataFrame
)

p_first_return_cont = sims_first_return_cont[end,end]

exp_err(p_first_return_cont,n)



function integrateD(evals,params)
    # evals is an integer specifying how many points to eval the function at
    # params is a CMEParams dictionary entry, i.e. CMEParams[3]
    N = 2*params["n"]+1 # ME order

    α = zeros(N)
    α[1] = params["c"]
    a = params["a"]
    b = params["b"]
    ω =  params["omega"]
    for k in 1:params["n"]
        kω = k*ω
        α[2*k] = (1/2)*( a[k]*(1+kω) - b[k]*(1-kω) )/(1+kω^2)
        α[2*k+1] = (1/2)*( a[k]*(1-kω) + b[k]*(1+kω) )/(1+kω^2)
    end

    period = 2*π/ω # the orbit repeats after this time
    edges = range(0,period,length=evals+1) # points at which to evaluate the fn
    h = period/(evals)

    orbit_LHS = α
    orbit_RHS = zeros(N)
    v_RHS = zeros(N)
    v_RHS[1] = 1
    v_LHS = ones(N)
    D = zeros(N,N)
    for t in edges[2:end]
        orbit_RHS[1] = α[1]
        for k in 1:params["n"]
            kωt = k*ω*t
            idx = 2*k
            idx2 = idx+1
            temp_cos = cos(kωt)
            temp_sin = sin(kωt)
            orbit_RHS[idx] = α[idx]*temp_cos + α[idx2]*temp_sin
            orbit_RHS[idx2] = -α[idx]*temp_sin + α[idx2]*temp_cos
            v_RHS[idx] = temp_cos - temp_sin
            v_RHS[idx2] = temp_sin + temp_cos
        end
        orbit_RHS = orbit_RHS./sum(orbit_RHS)
        orbit = (orbit_LHS+orbit_RHS)./2

        v = exp(-(t-h))*(v_LHS - exp(-h)*v_RHS)

        Dᵢ = v*orbit'
        D += Dᵢ

        orbit_LHS = copy(orbit_RHS)
        v_LHS = copy(v_RHS)
    end
    D = (1/(1-exp(-period)))*D
    return D
end

# numerical approximation of D
k = 100_000_000
T=[]
N=[]
n=184
evals = Int(ceil(k/(n/3)))# number of function evals
CMEParams[n]["D"], t = @timed integrateD(evals,CMEParams[n])
CMEParams[n]["intDevals"] = evals
push!(N,n)
push!(T,t)