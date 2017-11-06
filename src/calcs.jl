function constructSimModelMesh(model::SimModel3D{T},grid_size) where T
  origin = model.cell_origins[1]
  a_vec = model.cell_origins[end,1,1]-origin+Point3D{T}(model.prim_vec[1,:])
  b_vec = model.cell_origins[1,end,1]-origin+Point3D{T}(model.prim_vec[2,:])
  c_vec = model.cell_origins[1,1,end]-origin+Point3D{T}(model.prim_vec[3,:])
  points = Array{T,3}(grid_size,grid_size,grid_size)
  for (ia,a) in enumerate(linspace(0,1,grid_size)),(ib,b) in enumerate(linspace(0,1,grid_size)),(ic,c) in enumerate(linspace(0,1,grid_size))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    w = zero(Complex{T})
    for wfc in model.wfcs
      w+=wave_equation(point,wfc)
    end
    points[ia,ib,ic] = real(w)
  end
  return points
end

function calculateDipoleMesh(model::SimModel3D{T},cs,grid_size) where T
  origin = model.cell_origins[1]
  a_vec = model.cell_origins[end,1,1]-origin+Point3D{T}(model.prim_vec[1,:])
  b_vec = model.cell_origins[1,end,1]-origin+Point3D{T}(model.prim_vec[2,:])
  c_vec = model.cell_origins[1,1,end]-origin+Point3D{T}(model.prim_vec[3,:])
  dip_mesh = Array{Tuple{Point3D{T},Point3D{T}},3}(grid_size,grid_size,grid_size)
  for (ia,a) in enumerate(linspace(0,1,grid_size)),(ib,b) in enumerate(linspace(0,1,grid_size)),(ic,c) in enumerate(linspace(0,1,grid_size))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    dip = Point3D{T}(0.0)
    for (wfc1,c1) in zip(model.wfcs,cs)
      for (wfc2,c2) in zip(model.wfcs,cs)
        dip += conj(c1*wave_equation(point,wfc1))*c2*wave_equation(point,wfc2)*point
      end
    end
    dip_mesh[ia,ib,ic] = (point,dip)
  end
  return dip_mesh
end

function calculateTotDipole(model::SimModel3D{T},cs,grid_size) where T
  origin = model.cell_origins[1]
  a_vec = model.cell_origins[end,1,1]-origin+Point3D{T}(model.prim_vec[1,:])
  b_vec = model.cell_origins[1,end,1]-origin+Point3D{T}(model.prim_vec[2,:])
  c_vec = model.cell_origins[1,1,end]-origin+Point3D{T}(model.prim_vec[3,:])
  dip_x = zero(typeof(cs[1]))
  dip_y = zero(typeof(cs[1]))
  dip_z = zero(typeof(cs[1]))
  n1 = zero(T)
  n2 = zero(T)
  for (ia,a) in enumerate(linspace(0,1,grid_size)),(ib,b) in enumerate(linspace(0,1,grid_size)),(ic,c) in enumerate(linspace(0,1,grid_size))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    for i1=1:length(model.wfcs)
      wfc1,c1 = model.wfcs[i1],cs[i1]
      w1 = wave_equation(point,wfc1)
      n1 += norm(w1)^2
      for i2=1:length(model.wfcs)
        wfc2,c2 = model.wfcs[i2],cs[i2]
        w2 = wave_equation(point,wfc2)
        n2+=norm(w2)^2
        fac = conj(c1*w1)*c2*w2
        dip_x += fac*point.x
        dip_y += fac*point.y
        dip_z += fac*point.z
      end
    end
  end
  n = sqrt(n1*n2)
  return Point3D{T}(real(dip_x)/n,real(dip_y)/n,real(dip_z)/n)
end