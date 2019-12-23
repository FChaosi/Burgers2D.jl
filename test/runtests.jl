using Burgers2D, Test

#Test for x spacing
function test_gridx()
  g = Burgers2D.grid(-100, 100, -2pi, 2pi, 100, 200)
  dx = 200/100

  return isapprox(g.X[3, 1] - g.X[2, 1], dx, rtol = 1e-12)
end

#Test for y spacing
function test_gridy()
  g = Burgers2D.grid(-pi, pi, 0, 5, 10, 20)
  dy = 5/20

  return isapprox(g.Y[1, 3]-g.Y[1, 2], dy, rtol=1e-12)
end

#Test for x-, y-wavenumber spacing
function test_wavenumbers()
  Lx, Ly = 2*pi, 3.5*pi
  g = Burgers2D.grid(-Lx/2, Lx/2, -Ly/2, Ly/2, 100, 300)

  dk, dl = 2*pi/Lx, 2*pi/Ly

 return (isapprox(g.K[2,3] - g.K[1, 3], dk, rtol=1e-12) & isapprox(g.L[3,2] - g.L[3, 1], dl, rtol=1e-12) )
end

#Test for spectral derivative function
function test_spectralderivative()
  g = Burgers2D.grid(-2pi, 2pi, -4pi, 4pi, 100, 300)
  u = @. sin(g.X)*cos(g.Y)
  uxx_analytic = @. -sin(g.X)*cos(g.Y)
  uxx_numerical = Burgers2D.specd(u, g.K, 2)
  return isapprox(uxx_analytic, uxx_numerical, rtol=1e-12)
end

# begin tests
@testset "Grid Tests" begin
  @test test_gridx()
  @test test_gridy()
  @test test_wavenumbers()
end

@testset "Derivative Tests" begin
  @test test_spectralderivative()
end
