function grid(a1,b1,a2,b2,nx,ny)

  L1 = b1 - a1
  L2 = b2 - a2

  dx = L1/nx
  dy = L2/ny

  x = range(a1, stop = b1 - dx, length = nx)
  y = range(a2, stop = b2 - dy, length = ny)

  k = [0:nx/2; -(nx/2 - 1):-1]
  l = [0:ny/2; -(ny/2 - 1):-1]

  X, Y = zeros(nx, ny), zeros(nx, ny)

  K, L = zeros(nx, ny), zeros(nx, ny)

  for i in 1:nx, j in 1:ny
    X[i,j] = x[i]
    Y[i,j] = y[j]

    K[i,j] = k[i]
    L[i,j] = l[j]

  end

return (X,Y,K,L)

end
