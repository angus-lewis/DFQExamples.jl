const t=10.0

@inline initial_condition_generator(rng,n,model) = (X=[eps()],φ=[1])

const _first_exit = first_exit_x(0.5*eps(),1.0)
function first_exit_or_t10(model,SFM,SFM0)
    t1 = _first_exit(model,SFM,SFM0)
    i2 = SFM.t>10.0
    i2 && SFM = (t=10.0, X=SFM.X, φ=SFM.φ, n=SFM.m)
    return (Ind=(i2||t1.Ind), SFM=t1.SFM)
end

const tau = first_exit_or_t10

const x = range(0,10,length=6_001)