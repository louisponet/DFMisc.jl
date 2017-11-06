function write_SIM_model(filename::String,model::SimModel3D{T},grid_size) where T
  origin = model.cell_origins[1]
  a_vec = model.cell_origins[end,1,1]-origin+Point3D{T}(model.prim_vec[1,:])
  b_vec = model.cell_origins[1,end,1]-origin+Point3D{T}(model.prim_vec[2,:])
  c_vec = model.cell_origins[1,1,end]-origin+Point3D{T}(model.prim_vec[3,:])
  points = Array{WfPoint3D{T},3}(grid_size,grid_size,grid_size)
  for (ia,a) in enumerate(linspace(0,1,grid_size)),(ib,b) in enumerate(linspace(0,1,grid_size)),(ic,c) in enumerate(linspace(0,1,grid_size))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    w = zero(Complex{T})
    for wfc in model.wfcs
      w+=wave_equation(point,wfc)
    end
    points[ia,ib,ic] = WfPoint3D{T}(w,point)
  end
  tmp_wfc = Wfc3D(points,Point3D{T}[],PhysAtom())
  write_xsf_file(filename,tmp_wfc)
end