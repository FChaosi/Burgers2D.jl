module Burgers2D

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
f = zeros(d.nx, d.ny)
g = zeros(d.nx, d.ny)

for i in 1:(d.nx), j in 1:(d.ny)
    f[i, j] = u0(gr.X[i, j], gr.Y[i, j])
    g[i, j] = v0(gr.X[i, j], gr.Y[i, j])
end

uvt0 = (f, g)

#Further Initial Conditions
uvt1 = euler(uvt0, gr.K, gr.L, tstep)
uvt2 = euler(uvt1, gr.K, gr.L, tstep)

#Iteration
nt = round(tfin / tstep)

AB3data = (uvt0, uvt1, uvt2)

for i in 1:(nt-2)
    AB3data = AB3(AB3data[1], AB3data[2], AB3data[3], gr.K, gr.L, tstep)
end

return AB3data[3]
end


function tester(tstep, tfin)

#Amplitude coefficients for u
c1, d1 = 0.3, 1

#Amplitude coefficients for v
c2, d2 = 1, 0.5

#Standard deviations for u
s1, e1 = 1, 0.2

#Standard deviations for v
s2, e2 = 1, 0.6

#Diffusivity coefficient
k=1

#Sigma functions for u
sig11t(t) = sqrt(2*k*t + (s1)^2)
sig12t(t) = sqrt(2*k*t + (e1)^2)

#Sigma functions for v
sig21t(t) = sqrt(2*k*t + (s2)^2)
sig22t(t) = sqrt(2*k*t + (e2)^2)

#Analytic solution for u
gaussianuX(x, t) = c1 * (s1 / sig11t(t)) * exp(-x^2/(2 * ( sig11t(t) )^2 ))
gaussianuY(y, t) = d1 * (e1 / sig12t(t)) * exp(-y^2/(2 * ( sig12t(t) )^2 ))

u_analytic(x, y, t) = gaussianuX(x, t) * gaussianuY(y, t)

#Analytic solution for v
gaussianvX(x, t) = c2 * (s2 / sig21t(t)) * exp(-x^2/(2 * ( sig21t(t) )^2 ))
gaussianvY(y, t) = d2 * (e2 / sig22t(t)) * exp(-y^2/(2 * ( sig22t(t) )^2 ))

v_analytic(x, y, t) = gaussianvX(x, t) * gaussianvY(y, t)

#Initial Conditions
u0(x, y) = u_analytic(x, y, 0)
v0(x, y) = v_analytic(x, y, 0)

#Set up grid
d = Burgers2D.Domain(-2pi, 2pi, -4pi, 4pi, 128, 256)
g = Burgers2D.grid(-2pi, 2pi, -4pi, 4pi, 128, 256)

(u_numerical, v_numerical) = B2D(d, u0, v0, tstep, tfin)

nt = round(tfin / tstep)

usol_analytic = zeros(128, 256)
vsol_analytic = zeros(128, 256)

for i in 1:128, j in 1:256
    usol_analytic[i, j] = u_analytic(g.X[i, j], g.Y[i, j], tfin)
    vsol_analytic[i, j] = v_analytic(g.X[i, j], g.Y[i, j], tfin)
end

uerror = maximum(abs.(usol_analytic - u_numerical))
verror = maximum(abs.(vsol_analytic - v_numerical))

return (uerror, verror)

end




end # module
