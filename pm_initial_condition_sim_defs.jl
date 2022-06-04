const t=2.0

@inline initial_condition_generator(rng,n,model) = (X=[0.0],φ=[1])

const tau = fixed_time(t)

const x = range(0,10,length=10_001)