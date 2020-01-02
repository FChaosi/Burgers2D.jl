module Burgers2D

using FFTW

include("grids.jl")

include("RHS.jl")

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
function B2D(d::Domain, u0, v0, tstep, tfin)

#Grid
gr = Burgers2D.grid(d.a1, d.b1, d.a2, d.b2, d.nx, d.ny)

#Initial Conditions f and g
f = @. u0(gr.X, gr.Y)
g = @. v0(gr.X, gr.Y)

#Convert f and g to Fourier space
fh = rfft(f)
gh = rfft(g)

uvh_t0 = (fh, gh)

#Further Initial Conditions
uvh_t1 = euler(uvh_t0, gr.K, gr.L, tstep)
uvh_t2 = euler(uvh_t1, gr.K, gr.L, tstep)

#Iteration
nt = round(tfin / tstep)

AB3data = (uvh_t0, uvh_t1, uvh_t2)

for i in 1:(nt - 2)
    AB3data = AB3(AB3data[1], AB3data[2], AB3data[3], gr.K, gr.L, tstep)
end

(uh, vh) = AB3data[3]

(lengthx, lengthy) = size(uh)
nx = 2*lengthx - 2

u = irfft(uh, nx)
v = irfft(vh, nx)

return (u, v)
end




end # module
