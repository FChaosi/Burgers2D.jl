using Burgers2D, Test, FFTW

#Test for x, y grid-spacing
function test_gridxdy()
  Lx, Ly = 20.2, 4pi
  nx, ny = 200, 30
  g = Burgers2D.Grid(-Lx/2, Lx/2, -Ly/2, Ly/2, nx, ny)
  dx, dy = Lx/nx, Ly/ny

  return isapprox(g.X[3, 1] - g.X[2, 1], dx, rtol = 1e-12) & isapprox(g.Y[1, 3]-g.Y[1, 2], dy, rtol=1e-12)
end

#Test for x-, y-wavenumber spacing
function test_wavenumbers()
  Lx, Ly = 2*pi, 3.5*pi
  nx, ny = 20, 30
  g = Burgers2D.Grid(-Lx/2, Lx/2, -Ly/2, Ly/2, nx, ny)
  dk, dl = 2*pi/Lx, 2*pi/Ly

 return (isapprox(g.K[2,3] - g.K[1, 3], dk, rtol=1e-12) & isapprox(g.L[3,2] - g.L[3, 1], dl, rtol=1e-12) )
end

#Test for spectral derivative function
function test_spectralderivative()
  Lx, Ly = 4pi, 7pi
  g = Burgers2D.Grid(-2pi, 2pi, -4pi, 4pi, 40, 60)

  x, y = g.X, g.Y
  k0, l0 = 2pi/Lx, 2pi/Ly

  u = @. sin(2k0*x)*cos(l0*y)
  uxx_analytic  = @. -(2k0)^2*u
  uxx_numerical = Burgers2D.specd(u, g.K, 2)

  return isapprox(uxx_analytic, uxx_numerical, rtol=1e-12)
end

function test_laplacian()
  Lx, Ly = 3pi, 2pi
  g = Burgers2D.Grid(-Lx/2, Lx/2, -Ly/2, Ly/2, 40, 60)
  x, y = g.X, g.Y
  k0, l0 = 2pi/Lx, 2pi/Ly

  f = @. cos(4k0*x)*sin(3l0*y)

  fxx =  @. -(4*k0)^2*f
  fyy =  @. -(3*l0)^2*f

  laplacian_analytic = @. fxx + fyy
  laplacian_numeric = Burgers2D.laplacian(f, g.K, g.L)

  isapprox(laplacian_analytic, laplacian_numeric, rtol=1e-12)
end

function test_advection()
  Lx, Ly = 3pi, 2pi
  g = Burgers2D.Grid(-Lx/2, Lx/2, -Ly/2, Ly/2, 40, 60)
  x, y = g.X, g.Y
  k0, l0 = 2pi/Lx, 2pi/Ly

  u = @. sin(2k0*x)*cos(l0*y)
  v = @. sin(3k0*x)^2*cos(3l0*y)
  f = @. cos(4k0*x)*sin(4l0*y)

  fx =  @. -4*k0*sin(4k0*x)*sin(4l0*y)
  fy =  @.  4*l0*cos(4k0*x)*cos(4l0*y)

  advection_analytic = @. u*fx + v*fy
  advection_numeric = Burgers2D.advection((u, v), f, g.K, g.L)

  isapprox(advection_analytic, advection_numeric, rtol=1e-12)
end

function test_diffusion_sin(timestepper, errorscale)
  tstep_arr = [1e-4 1e-5]
  tfin = 1e-3
  
  u_err_arr = zeros(2)
  v_err_arr = zeros(2)

  Lx, Ly = 2pi, 4pi
  g = Burgers2D.Grid(-Lx/2, Lx/2, -Ly/2, Ly/2, 64, 128)
  d = Burgers2D.Domain(-Lx/2, Lx/2, -Ly/2, Ly/2, 64, 128)
  
  nu = 1.0
  
  x, y = g.X, g.Y

  u_analytic(x, y, t) = sin(x)*sin(y)*exp(-2t)
  v_analytic(x, y, t) = cos(x)*cos(y)*exp(-2t)

  u0(x, y) = u_analytic(x, y, 0)
  v0(x, y) = v_analytic(x, y, 0)

  u_analytic_val = u_analytic.(g.X, g.Y, tfin)
  v_analytic_val = v_analytic.(g.X, g.Y, tfin)

  for i in 1:2
    (u_numeric, v_numeric) = Burgers2D.B2D(d, u0, v0, tstep_arr[i], tfin, timestepper, nu)

    u_err_arr[i] = maximum( abs.(u_numeric - u_analytic_val) )
    v_err_arr[i] = maximum( abs.(v_numeric - v_analytic_val) )
  end

  #Test for quadratic scaling of the error
  u_err_scale = u_err_arr[1]/u_err_arr[2]
  v_err_scale = v_err_arr[1]/v_err_arr[2]

  isapprox(u_err_scale, errorscale, rtol=1e-1) & isapprox(v_err_scale, errorscale, rtol=1e-1)
end

#New test
function test_diffusion_gauss(timestepper, errorscale)
  d = Burgers2D.Domain(-pi, pi, -3.0, 3.0, 128, 130)
  g = Burgers2D.Grid(d)
  nu = 0.01

  tstep_arr = [1e-2 1e-3]
  tfin = 2e-1
  u_err_arr = zeros(2)
  v_err_arr = zeros(2)
  
  gaussian(x, σ, ampl) = ampl*exp(-x^2/(2*σ^2))  
  
  function analytic_solution(x, σ, ampl, t, nu)
    σt = sqrt(2*nu*t + σ^2)
    return gaussian(x, σt, ampl*σ/σt) 
  end
  
  σx_u, σy_u, ampl_u = 0.20, 0.19, 0.01
  σx_v, σy_v, ampl_v = 0.21, 0.22, 0.008
  
  u_analytic(x, y, t) = ampl_u*analytic_solution(x, σx_u, 1, t, nu)*analytic_solution(y, σy_u, 1, t, nu)
  v_analytic(x, y, t) = ampl_v*analytic_solution(x, σx_v, 1, t, nu)*analytic_solution(y, σy_v, 1, t, nu)

  u0(x, y) = u_analytic(x, y, 0)
  v0(x, y) = v_analytic(x, y, 0)

  u_analytic_val = u_analytic.(g.X, g.Y, tfin)
  v_analytic_val = v_analytic.(g.X, g.Y, tfin)

  for i in 1:2
    (u_numeric, v_numeric) = Burgers2D.B2D(d, u0, v0, tstep_arr[i], tfin, timestepper, nu)

    u_err_arr[i] = maximum( abs.(u_numeric - u_analytic_val) )
    v_err_arr[i] = maximum( abs.(v_numeric - v_analytic_val) )
  end

  #Test for quadratic scaling of the error
  u_err_scale = u_err_arr[1]/u_err_arr[2]
  v_err_scale = v_err_arr[1]/v_err_arr[2]

  println("Gaussian")
  println(u_err_scale)
  println(v_err_scale)

  println("Error with tstep=1e-5 is:")
  println(u_err_arr[1])

  println("Error with tstep=1e-4 is:")
  println(u_err_arr[2])

  isapprox(u_err_scale, errorscale, rtol=1e-1) & isapprox(v_err_scale, errorscale, rtol=1e-1)
end

# begin tests
@testset "Grid Tests" begin
  @test test_gridxdy()
  @test test_wavenumbers()
end

@testset "Derivative Tests" begin
  @test test_spectralderivative()
  @test test_laplacian()
end

@testset "Advection Test" begin
  @test test_advection()
end

@testset "Diffusion Test" begin
  @test test_diffusion_sin("ForwardEuler", 10)
  @test test_diffusion_sin("AB3", 100)
  @test test_diffusion_gauss("ForwardEuler", 10)
  @test test_diffusion_gauss("AB3", 100)
  @test test_diffusion_gauss("RK4", 10000)
end
