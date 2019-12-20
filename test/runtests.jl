using Burgers2D, Test

function test_grid()
  return isapprox(1, 1.0000001, rtol=1e-4)
end

@test test_grid()
