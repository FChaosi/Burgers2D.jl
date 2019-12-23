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
function test_wavex()
  Lx, Ly = 2*pi, 3.5*pi
  g = Burgers2D.grid(-Lx/2, Lx/2, -Ly/2, Ly/2, 100, 300)

  dk, dl = 2*pi/Lx, 2*pi/Ly

 return (isapprox(g.K[2,3] - g.K[1, 3], dk, rtol=1e-12) & isapprox(g.L[3,2] - g.L[3, 1], dl, rtol=1e-12) )
end

#Test for spectral derivative function
function testspec()
  g = Burgers2D.grid(-2pi, 2pi, -4pi, 4pi, 100, 300)

  u = zeros(100, 300)
  actual = zeros(100, 300)

  for i in 1:100, j in 1:300
      u[i,j] = sin(g.X[i,j])*cos(g.Y[i,j])
      actual[i,j] = -sin(g.X[i,j])*cos(g.Y[i,j])
  end

  spec = Burgers2D.specd(u, g.K, 2)

  err = maximum( abs.(actual - spec) )

  return isapprox(err, 0, atol=1e-10)
end








@test test_gridx()
@test test_gridy()
@test test_wavex()
@test testspec()
