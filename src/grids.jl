struct Grid
    X :: Array{Float64, 2}
    Y :: Array{Float64, 2}
    K :: Array{Float64, 2}
    L :: Array{Float64, 2}
end

function grid(a1, b1, a2, b2, nx, ny)
  Lx, Ly = b1 - a1, b2 - a2
  dx, dy = Lx/nx, Ly/ny

  x = reshape(Array(range(a1, stop = b1 - dx, length = nx)), (nx, 1))
  y = reshape(Array(range(a2, stop = b2 - dy, length = ny)), (1, ny))

  k = reshape( Array( (2pi/Lx)*cat(0:nx/2, dims=1) ), (Int(nx/2 + 1), 1) )
  l = reshape( Array( (2pi/Ly)*cat(0:ny/2, -(ny/2 - 1):-1, dims=1) ), (1, ny) )

  X, Y = x .+ 0*y, 0*x .+ y
  K, L = k .+ 0*l, 0*k .+ l

  return Grid(X, Y, K, L)
end
