using FFTW

struct Data
  un
  vn
end

struct Triple
  h
  hn
  hnn
end

#Spectral Method for Derivative
function specd(values, wave, n)
  spectral = ((im*wave).^n).*fft(values)
  return real(ifft(spectral))
end


#Calculate RHS
function RHS(data::Data, wavex, wavey)
    #Calculate Laplacian of un and vn
    laplun = specd(data.un, wavex, 2) + specd(data.un, wavey, 2)
    laplvn = specd(data.vn, wavex, 2) + specd(data.vn, wavey, 2)

    unx = specd(data.un, wavex, 1)
    uny = specd(data.un, wavey, 1)

    vnx = specd(data.vn, wavex, 1)
    vny = specd(data.vn, wavey, 1)

    advectionu = (data.un).*unx + (data.vn).*uny
    advectionv = (data.un).*vnx + (data.vn).*vny

    ut = laplun - 0*advectionu
    vt = laplvn - 0*advectionv

    return Data(ut, vt)
end

#Forward Euler Time Step
function euler(data::Data, wavex, wavey, tstep)
    dd = RHS(data, wavex, wavey)

    unn = data.un + tstep*(dd.un)
    vnn = data.vn + tstep*(dd.vn)

    return Data(unn, vnn)
end

#AB3 Method
function AB3(triple :: Triple, wavex, wavey, tstep)

    data1 = triple.h
    data2 = triple.hn
    data3 = triple.hnn

    inc1 = RHS(data1, wavex, wavey)
    inc2 = RHS(data2, wavex, wavey)
    inc3 = RHS(data3, wavex, wavey)

    incrementu = (tstep/12)*( 23*(inc1.un) -16*(inc2.un) + 5*(inc3.un) )
    incrementv = (tstep/12)*( 23*(inc1.vn) -16*(inc2.vn) + 5*(inc3.vn) )

    data1new = Data(data2.un, data2.vn)
    data2new = Data(data3.un, data3.vn)
    data3new = Data(data3.un + incrementu, data3.vn + incrementv)

    return Triple(data1new, data2new, data3new)
end
