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

"""`    `
    RHS((un, vn), wave_x, wave_y)

Returns the a tuple with the RHS for 2D Burgers equation.
"""
function RHS((un, vn), wave_x, wave_y, nu)
  #Calculate Laplacian of un and vn
  diffu = nu*laplacian(un, wave_x, wave_y)
  diffv = nu*laplacian(vn, wave_x, wave_y)

  #Calculate advection terms for un and vn
  advectionu = advection((un, vn), un, wave_x, wave_y)
  advectionv = advection((un, vn), vn, wave_x, wave_y)

  un_t = diffu - (0*advectionu)
  vn_t = diffv - (0*advectionv)

  return (un_t, vn_t)
end
