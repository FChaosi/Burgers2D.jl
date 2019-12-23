module Burgers2D

include("grids.jl")

include("RHS.jl")


#Structures
struct Domain
     a1
     a2
     b1
     b2
     nx
     ny
 end

 struct Inputs
     d :: Domain
     f0 :: Function
     g0 :: Function
     tstep
     tfin
 end


#Burgers 2D Algorithm
function B2D(input :: Inputs)

#Grid
g = Burgers2D.grid(input.d.a1, input.d.a2, input.d.b1, input.d.b2, input.d.nx, input.d.ny)

#Initial Conditions f and g
n1 = input.d.nx
n2 = input.d.ny

ft0 = zeros(n1, n2)
gt0 = zeros(n1, n2)

for i in 1:n1, j in 1:n2
    ft0[i,j] = input.f0(g.X[i,j], g.Y[i,j])
    gt0[i,j] = input.g0(g.X[i,j], g.Y[i,j])
end

uvt0 = Data(ft0, gt0)

#Further Initial Conditions
uvt1 = euler(uvt0, g.K, g.L, input.tstep)
uvt2 = euler(uvt1, g.K, g.L, input.tstep)

#Iteration
nt = round( (input.tfin) / (input.tstep) )
update = Triple(uvt0, uvt1, uvt2)

for i in 1:(nt-2)
    update = AB3(update, g.K, g.L, input.tstep)
end

return update.hnn
end




end # module
