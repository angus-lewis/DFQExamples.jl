@inline initial_condition_generator(rng,n,model) = (φ=[4], X=[2.0], Y=[0.0])

const tau = first_exit_y(0.0,model.dq.model.b+1.0)

const x = range(0,1.6,length=10_001)