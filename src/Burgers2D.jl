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

dom = Burgers2D.Domain(-pi, pi, -2pi, 2pi, 16, 64)

gr = Grid.(dom.a1, dom.b1, dom.a2, dom.b2, dom.nx, dom.ny)

cx = 0.01
sx = 0.1

cy = 0.3
sy = 0.2

tstep1 = 1e-5
tstep2 = 1e-4
tfin = 1e-2

gaussx(x) = cx * exp(-x^2/(2*cx^2))

gaussy(y) = cy * exp(-y^2/(2*cy^2))

gauss(x, y) = gaussx(x) * gaussy(y)

gaussv = gauss.(gr.X, gr.Y)

function anx(x, t)
    st = 2*t + sx^2
    return cx * sx/sqrt(st) * exp(-x^2/(2*st))
end

function any(y, t)
    st = 2*t + sy^2
    return cy * sy/sqrt(st) * exp(-y^2/(2*st))
end

an(x, y, t) = anx(x, t) * any(y, t)

num_val1 = B2D(dom, gauss, gauss, tstep1, tfin, "AB3")[1]

num_val2 = B2D(dom, gauss, gauss, tstep2, tfin, "AB3")[1]

an_val = an.(gr.X, gr.Y, tfin)

diff1 = maximum(abs.(num_val1 - an_val))

diff2 = maximum(abs.(num_val2 - an_val))

println("Error with tstep=1e-5 is:")
println(diff1)

println("Error with tstep=1e-4 is:")
println(diff2)


end # module
