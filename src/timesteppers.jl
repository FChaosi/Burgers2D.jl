#Forward Euler Time Step
function euler((un, vn), wave_x, wave_y, tstep, nu)
  (un_t, vn_t) = RHS((un, vn), wave_x, wave_y, nu)
  
  return (un, vn) .+  tstep .* (un_t, vn_t) 
end

#AB3 Method
function AB3((un, vn), (unn, vnn), (unnn, vnnn), wave_x, wave_y, tstep, nu)
  inc1 =  23 .* RHS((unnn, vnnn), wave_x, wave_y, nu)
  inc2 = -16 .* RHS((unn, vnn), wave_x, wave_y, nu)
  inc3 =   5 .* RHS((un, vn), wave_x, wave_y, nu)

  increment = (tstep/12) .* (inc1 .+ inc2 .+ inc3)

  return ((unn, vnn), (unnn, vnnn), (unnn, vnnn) .+ increment )
end

function RK4((un,vn), wave_x, wave_y, tstep, nu)
  k1 = tstep .* RHS((un, vn), wave_x, wave_y, nu)
  k2 = tstep .* RHS((un, vn) .+ (0.5 .* k1), wave_x, wave_y, nu)
  k3 = tstep .* RHS((un, vn) .+ (0.5 .* k2), wave_x, wave_y, nu)
  k4 = tstep .* RHS((un, vn) .+ k3, wave_x, wave_y, nu)

  increment = 1/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
  
  return (un, vn) .+ increment
end
