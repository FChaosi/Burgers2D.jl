using Burgers2D, Test

#Test for x, y grid-spacing
function test_gridxdy()
  Lx, Ly = 20.2, 4pi
  nx, ny = 200, 30
  g = Burgers2D.grid(-Lx/2, Lx/2, -Ly/2, Ly/2, nx, ny)
  dx, dy = Lx/nx, Ly/ny

  return isapprox(g.X[3, 1] - g.X[2, 1], dx, rtol = 1e-12) & isapprox(g.Y[1, 3]-g.Y[1, 2], dy, rtol=1e-12)
end

#Test for x-, y-wavenumber spacing
function test_wavenumbers()
  Lx, Ly = 2*pi, 3.5*pi
  nx, ny = 20, 30
  g = Burgers2D.grid(-Lx/2, Lx/2, -Ly/2, Ly/2, nx, ny)
  dk, dl = 2*pi/Lx, 2*pi/Ly

 return (isapprox(g.K[2,3] - g.K[1, 3], dk, rtol=1e-12) & isapprox(g.L[3,2] - g.L[3, 1], dl, rtol=1e-12) )
end

#Test for spectral derivative function
function test_spectralderivative()
  Lx, Ly = 4pi, 7pi
  g = Burgers2D.grid(-2pi, 2pi, -4pi, 4pi, 40, 60)
  
  x, y = g.X, g.Y
  k0, l0 = 2pi/Lx, 2pi/Ly
  
  u = @. sin(2k0*x)*cos(l0*y)
  uxx_analytic  = @. -(2k0)^2*u
  uxx_numerical = Burgers2D.specd(u, g.K, 2)
  
  return isapprox(uxx_analytic, uxx_numerical, rtol=1e-12)
end

function test_advection()
  Lx, Ly = 3pi, 2pi
  g = Burgers2D.grid(-Lx/2, Lx/2, -Ly/2, Ly/2, 40, 60)
  x, y = g.X, g.Y
  k0, l0 = 2pi/Lx, 2pi/Ly
  
  f = @. cos(4k0*x)*sin(3l0*y)
  
  fxx =  @. -(4*k0)^2*f
  fyy =  @. -(3*l0)^2*f
  
  laplacian1 = @. fxx + fyy
  laplacian2 = Burgers2D.laplacian(f, (g.K, g.L))
  
  isapprox(laplacian1, laplacian2, rtol=1e-12)
end

function test_advection()
  Lx, Ly = 3pi, 2pi
  g = Burgers2D.grid(-Lx/2, Lx/2, -Ly/2, Ly/2, 40, 60)
  x, y = g.X, g.Y
  k0, l0 = 2pi/Lx, 2pi/Ly
  
  u = @. sin(2k0*x)*cos(l0*y)
  v = @. sin(3k0*x)^2*cos(3l0*y)
  f = @. cos(4k0*x)*sin(4l0*y)
  
  fx =  @. -4*k0*sin(4k0*x)*sin(4l0*y)
  fy =  @.  4*l0*cos(4k0*x)*cos(4l0*y)
  
  advection1 = @. u*fx + v*fy
  advection2 = Burgers2D.advection((u, v), f, (g.K, g.L))
  
  isapprox(advection1, advection2, rtol=1e-12)
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

