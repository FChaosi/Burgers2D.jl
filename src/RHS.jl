#Spectral Method for Derivative
function specd(values, wave, n)
  (nx, ny)  = size(values)
  wavepower = (im * wave).^n
  spectral = wavepower .* rfft(values)
  return irfft(spectral, nx)
end

"""
    advection((u, v), f)

Returns the advection of field f with flow (u, v):
u*∂f/∂x + v*∂f/∂y
"""
function advection((u, v), f, (kx, ky))
  fx = specd(f, kx, 1)
  fy = specd(f, ky, 1)
  advection = @. u*fx + v*fy
end

"""
    laplacian(f, (kx, ky))

Returns the Laplacian of field f: ∂²f/∂x² + ∂²f/∂y²
"""
function laplacian(f, (kx, ky))
  laplacianf = specd(f, kx, 2) + specd(f, ky, 2)
end

#Calculate RHS
function RHS((un, vn), wavex, wavey)
  #Calculate Laplacian of un and vn
  laplun = laplacian(un, (kx, ky))
  laplvn = laplacian(vn, (kx, ky))

  #Calculate advection terms for un and vn
  advectionu = advection((un, vn), un, (wavex, wavey))
  advectionv = advection((un, vn), vn, (wavex, wavey))

  forcingu = -(laplun - advectionu)
  forcingv = -(laplvn - advectionv)

  ut = laplun - advectionu + forcingu
  vt = laplvn - advectionv + forcingv

  return (ut, vt)
end

#Forward Euler Time Step in Fourier Space
function euler((unh, vnh), wavex, wavey, tstep)
    (ut_h, vt_h) = RHS((unh, vnh), wavex, wavey)

    (unn_h, vnn_h) = (unh, vnh) .+ (tstep .* (ut_h, vt_h))

    return (unn_h, vnn_h)
end

#AB3 Method in Fourier Space
function AB3((unh, vnh), (unn_h, vnn_h), (unnn_h, vnnn_h), wavex, wavey, tstep)

    inc1 = RHS((unh, vnh), wavex, wavey)
    inc2 = RHS((unn_h, vnn_h), wavex, wavey)
    inc3 = RHS((unnn_h, vnnn_h), wavex, wavey)

    increment = (tstep/12) .* (23. .* inc1 .+ (-16) .* inc2 .+ 5 .* inc3)

    (unh, vnh) = (unn_h, vnn_h)
    (unn_h, vnn_h) = (unnn_h, vnnn_h)
    (unnn_h, vnnn_h) = (unn_h, vnn_h) .+ increment

    return ((unh, vnh), (unn_h, vnn_h), (unnn_h, vnnn_h))
end
