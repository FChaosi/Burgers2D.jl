struct Grid
    X
    Y
    K
    L
end

function grid(a1, b1, a2, b2, nx, ny)

  Lx = b1 - a1
  Ly = b2 - a2

  dx = Lx/nx
  dy = Ly/ny

  x = reshape(Array(range(a1, stop = b1 - dx, length = nx)), (nx, 1))
  y = reshape(Array(range(a2, stop = b2 - dy, length = ny)), (1, ny))

  k = (2pi/Lx)*[0:nx/2; -(nx/2 - 1):-1]
  l = (2pi/Ly)*[0:ny/2; -(ny/2 - 1):-1]

  X, Y = x .+ 0*y, 0*x .+ y

  K, L = zeros(nx, ny), zeros(nx, ny)

  for i in 1:nx, j in 1:ny
    K[i, j] = k[i]
    L[i, j] = l[j]
  end

  return Grid(X, Y, K, L)

end
