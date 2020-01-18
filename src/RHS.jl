#Spectral Method for Derivative
function specd(phys_values, wave, n)
  (nx, ny)  = size(phys_values)
  wavepower = (im * wave).^n
  spectral_deriv = wavepower .* rfft(phys_values)
  return irfft(spectral_deriv, nx)
end

"""
    advection((u, v), f)

Returns the advection of field f with flow (u, v):
u*∂f/∂x + v*∂f/∂y
"""
function advection((u, v), f, wave_x, wave_y)

  fx = specd(f, wave_x, 1)
  fy = specd(f, wave_y, 1)
  advection = @. u*fx + v*fy
end

"""`    `
    laplacian(f, (kx, ky))

Returns the Laplacian of field f: ∂²f/∂x² + ∂²f/∂y²
"""
function laplacian(f, wave_x, wave_y)
 return specd(f, wave_x, 2) + specd(f, wave_y, 2)
end

#Calculate RHS
function RHS((un, vn), wave_x, wave_y)
  #Calculate Laplacian of un and vn
  laplun = laplacian(un, wave_x, wave_y)
  laplvn = laplacian(vn, wave_x, wave_y)

  #Calculate advection terms for un and vn
  advectionu = advection((un, vn), un, wave_x, wave_y)
  advectionv = advection((un, vn), vn, wave_x, wave_y)

  un_t = laplun - (0*advectionu)
  vn_t = laplvn - (0*advectionv)

  return (un_t, vn_t)
end

#Forward Euler Time Step
function euler((un, vn), wave_x, wave_y, tstep)
    (un_t, vn_t) = RHS((un, vn), wave_x, wave_y)

    return (un, vn) .+ ( tstep .* (un_t, vn_t) )
end

#AB3 Method
function AB3((un, vn), (unn, vnn), (unnn, vnnn), wave_x, wave_y, tstep)

    inc1 = 23 .*  RHS((unnn, vnnn), wave_x, wave_y)
    inc2 = -16 .* RHS((unn, vnn), wave_x, wave_y)
    inc3 = 5 .* RHS((un, vn), wave_x, wave_y)

    increment = (tstep/12) .* (inc1 .+ inc2 .+ inc3)

    return ((unn, vnn), (unnn, vnnn), (unnn, vnnn) .+ increment )
end

function RK4((un,vn), wave_x, wave_y, tstep)

    k1 = tstep .* RHS((un,vn), wave_x, wave_y)
    k2 = (2 * tstep) .* RHS((un,vn) .+ (0.5 .* k1), wave_x, wave_y)
    k3 = (2 * tstep) .* RHS((un,vn) .+ (0.5 .* k2), wave_x, wave_y)
    k4 = tstep .* RHS((un,vn) .+ k3, wave_x, wave_y)

    increment = 1/6 .* (k1 .+ k2 .+ k3 .+ k4)
    (unn,vnn) = (un,vn) .+ increment
    return (unn, vnn)
end
