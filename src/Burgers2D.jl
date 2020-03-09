module Burgers2D

using FFTW

include("grids.jl")
include("RHS.jl")
include("timesteppers.jl")


#Burgers 2D Algorithm
function B2D(d::Domain, u0, v0, tstep, tfin, method::String, nu)

  #Grid
  gr = Grid(d)

  #Initial Conditions f and g
  f = @. u0(gr.X, gr.Y)
  g = @. v0(gr.X, gr.Y)

  uv_t0 = (f, g)

  nt = Int(round(tfin / tstep))

  if method=="ForwardEuler" 
    uv_t = uv_t0
    #Iteration
    for i in 1:nt
      uv_t = Burgers2D.euler(uv_t, gr.K, gr.L, tstep, nu)
    end
    
    (u, v) = (uv_t[1], uv_t[2])
    
  elseif method=="AB3"
    
    #Further Initial Conditions
    uv_t1 = Burgers2D.euler(uv_t0, gr.K, gr.L, tstep, nu)
    uv_t2 = Burgers2D.euler(uv_t1, gr.K, gr.L, tstep, nu)

    #Iteration
    AB3data = (uv_t0, uv_t1, uv_t2)

    for i in 1:(nt - 2)
      AB3data = AB3(AB3data[1], AB3data[2], AB3data[3], gr.K, gr.L, tstep, nu)
    end

    (u, v) = AB3data[3]
    
  elseif method=="RK4"

    RK4data = uv_t0

    for i in 1:nt
      RK4data = RK4(RK4data, gr.K, gr.L, tstep, nu)
    end

    (u_fin, v_fin) = RK4data
  end

end


end # module
