const t=10.0

function initial_condition_generator(rng,n,model)
    i = rand(rng,1:2,n) 
    X = Vector{Float64}(undef,n)
    for k in 1:n
        x = Inf
        while x>=1
            x = -log(rand(rng))
        end
        X[k]=x
    end
    return (X=X,φ=i)
end

const _first_exit = first_exit_x(0.5*eps(),1.0)
function first_exit_or_t10(model,SFM,SFM0)
    t1 = _first_exit(model,SFM,SFM0)
    i2 = SFM.t>10.0
    return (Ind=(i2||t1.Ind), SFM=t1.SFM)
end

const tau = first_exit_or_t10

const x = range(0,10,length=6_001)