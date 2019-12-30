using FFTW



#Spectral Method for Derivative
function specd(specvalues, wave, n)
  (lengthx, lengthy) = size(specvalues)
  nx = 2*lengthx - 2
  scaling = (im * wave).^n
  spectral = scaling .* specvalues
  return irfft(spectral, nx)
end

#Laplacian
function lapln(specvalues, wavex, wavey)
valuesxx = specd(specvalues, wavex, 2)
valuesyy = specd(specvalues, wavey, 2)
  return valuesxx + valuesyy
end

#Calculate RHS
function RHS((unh, vnh), wavex, wavey)

    #Laplacian
    laplun = lapln(unh, wavex, wavey)
    laplvn = lapln(vnh, wavex, wavey)

    #First Derivatives
    unx = specd(unh, wavex, 1)
    uny = specd(unh, wavey, 1)

    vnx = specd(vnh, wavex, 1)
    vny = specd(vnh, wavey, 1)

    #Conversion to physical space
    (lengthx, lengthy) = size(unh)
    nx = 2*lengthx - 2

    un = irfft(unh, nx)
    vn = irfft(vnh, nx)

    #Advection
    advectionu = (un).*unx + (vn).*uny
    advectionv = (un).*vnx + (vn).*vny

    ut = laplun - (0 .* advectionu)
    vt = laplvn - (0 .* advectionv)

    ut_h = rfft(ut)
    vt_h = rfft(vt)

    return (ut_h, vt_h)
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
