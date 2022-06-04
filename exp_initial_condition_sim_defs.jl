const t=2.0

function initial_condition_generator(rng,n,model)
    i = rand(rng,1:2,n) 
    X = Vector{Float64}(undef,n)
    for k in 1:n
        x = Inf
        while x>=model.b
            x = -log(rand(rng))
        end
        X[k]=x
    end
    return (X=X,φ=i)
end

const tau = fixed_time(t)

const x = range(0,10,length=10_001)