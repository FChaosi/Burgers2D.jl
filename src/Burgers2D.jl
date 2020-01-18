module Burgers2D

using FFTW

include("grids.jl")

include("RHS.jl")

include("AB3.jl")

#Structures
struct Domain
     a1 :: Float64
     b1 :: Float64
     a2 :: Float64
     b2 :: Float64
     nx :: Int64
     ny :: Int64
 end

#Burgers 2D Algorithm
function B2D(d::Domain, u0, v0, tstep, tfin, method::String)

#Grid
gr = Grid(d.a1, d.b1, d.a2, d.b2, d.nx, d.ny)

#Initial Conditions f and g
f = @. u0(gr.X, gr.Y)
g = @. v0(gr.X, gr.Y)

uv_t0 = (f, g)

nt = round(tfin / tstep)

if method=="AB3"
#Further Initial Conditions
uv_t1 = Burgers2D.euler(uv_t0, gr.K, gr.L, tstep)
uv_t2 = Burgers2D.euler(uv_t1, gr.K, gr.L, tstep)

#Iteration
AB3data = (uv_t0, uv_t1, uv_t2)

for i in 1:(nt - 2)
    AB3data = AB3(AB3data[1], AB3data[2], AB3data[3], gr.K, gr.L, tstep)
end

(u, v) = AB3data[3]

elseif method=="RK4"

RK4data = uv_t0

for i in 1:nt
    RK4data = RK4(RK4data, gr.K, gr.L, tstep)
end

(u_fin, v_fin) = RK4data

end

end




end # module
