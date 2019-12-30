using Burgers2D, Test, FFTW

#Test for x, y grid-spacing
function test_gridxdy()
  Lx, Ly = 200, 4pi
  nx, ny = 100, 200
  g = Burgers2D.grid(-Lx/2, Lx/2, -Ly/2, Ly/2, nx, ny)
  dx, dy = Lx/nx, Ly/ny

  return isapprox(g.X[3, 1] - g.X[2, 1], dx, rtol = 1e-12) & isapprox(g.Y[1, 3]-g.Y[1, 2], dy, rtol=1e-12)
end

#Test for x-, y-wavenumber spacing
function test_wavenumbers()
  Lx, Ly = 2*pi, 3.5*pi
  nx, ny = 100, 300
  g = Burgers2D.grid(-Lx/2, Lx/2, -Ly/2, Ly/2, nx, ny)
  dk, dl = 2*pi/Lx, 2*pi/Ly

 return (isapprox(g.K[2,3] - g.K[1, 3], dk, rtol=1e-12) & isapprox(g.L[3,2] - g.L[3, 1], dl, rtol=1e-12) )
end

#Test for spectral derivative function
function test_spectralderivative()
  g = Burgers2D.grid(-2pi, 2pi, -4pi, 4pi, 100, 300)
  u = @. sin(g.X)*cos(g.Y)
  uh = rfft(u)
  uxx_analytic = @. -sin(g.X)*cos(g.Y)
  uxx_numerical = Burgers2D.specd(uh, g.K, 2)
  return isapprox(uxx_analytic, uxx_numerical, rtol=1e-12)
end






function test_diffusion_AB3()
err = zeros(2, 2)
tstep_array = [1e-4, 1e-5]
tfin_array = [1e-7, 1e-7]

#Amplitude coefficients for u
c1, d1 = 20, 16

#Amplitude coefficients for v
c2, d2 = 17, 13

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

for i in 1:2
tstep = tstep_array[i]
tfin = tfin_array[i]

(u_numerical, v_numerical) = Burgers2D.B2D(d, u0, v0, tstep, tfin)

usol_analytic = @. u_analytic(g.X, g.Y, tfin)
vsol_analytic = @. v_analytic(g.X, g.Y, tfin)

(uerror, verror) = (maximum(abs.(usol_analytic - u_numerical)), maximum(abs.(vsol_analytic - v_numerical)))

err[i, 1] = uerror
err[i, 2] = verror
end

return (isapprox(err[1,1]/err[2, 1], 100, rtol=1e-12) & isapprox(err[1,2]/err[2, 2], 100, rtol=1e-12) )
end



# begin tests
@testset "Grid Tests" begin
  @test test_gridxdy()
  @test test_wavenumbers()
end

@testset "Derivative Tests" begin
  @test test_spectralderivative()
end

@testset "Diffusion Tests" begin
 @test test_diffusion_AB3()
end
