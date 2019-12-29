using FFTW



#Spectral Method for Derivative
function specd(values, wave, n)
(nx, ny)  = size(wave)
wavepower = (im * wave).^n
index = Int(nx/2 + 1)
realwave = wavepower[1:index, :]
  spectral = realwave .* rfft(values)
  return irfft(spectral, nx)
end

#Calculate RHS
function RHS((un, vn), wavex, wavey)
    #Calculate Laplacian of un and vn
    laplun = specd(un, wavex, 2) + specd(un, wavey, 2)
    laplvn = specd(vn, wavex, 2) + specd(vn, wavey, 2)

    unx = specd(un, wavex, 1)
    uny = specd(un, wavey, 1)

    vnx = specd(vn, wavex, 1)
    vny = specd(vn, wavey, 1)

    advectionu = (un).*unx + (vn).*uny
    advectionv = (un).*vnx + (vn).*vny

    forcingu = -(laplun - advectionu)
    forcingv = -(laplvn - advectionv)

    ut = laplun - advectionu + forcingu
    vt = laplvn - advectionv + forcingv

    return (ut, vt)
end

#Forward Euler Time Step
function euler((un, vn), wavex, wavey, tstep)
    (ut, vt) = RHS((un, vn), wavex, wavey)

    (unn, vnn) = (un, vn) .+ (tstep .* (ut, vt))

    return (unn, vnn)
end

#AB3 Method
function AB3((un, vn), (unn, vnn), (unnn, vnnn), wavex, wavey, tstep)

    inc1 = RHS((un, vn), wavex, wavey)
    inc2 = RHS((unn, vnn), wavex, wavey)
    inc3 = RHS((unnn, vnnn), wavex, wavey)

    increment= (tstep/12) .* (23. .* inc1 .+ (-16) .* inc2 .+ 5 .* inc3)

    (un, vn) = (unn, vnn)
    (unn, vnn) = (unnn, vnnn)
    (unnn, vnnn) = (unn, vnn) .+ increment

    return ((un, vn), (unn, vnn), (unnn, vnnn))
end
