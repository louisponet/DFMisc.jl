struct SimWfc3D{T<:AbstractFloat}
  character::Tuple{T,Int64,Symbol,Symbol}
  unit_cell::Tuple{Int64,Int64,Int64}
  center::Point3D{T}
end

function SimWfc3D(n::Int64,character::Symbol,sub_character::Symbol,unit_cell,center::Point3D{T},Z) where T
  if character==:s
    return SimWfc3D{T}((Z,n,character,sub_character),(unit_cell...),center)
  elseif character==:p
    return SimWfc3D{T}((Z,n,character,sub_character),(unit_cell...),center)
  end
end

mutable struct SimModel3D{T<:AbstractFloat}
  cell_origins::Array{Point3D{T},3}
  prim_vec::Matrix{T}
  wfcs::Array{SimWfc3D{T},1}
end

function SimModel3D(supercell,characters::Array{Tuple{Int64,Symbol,Symbol}},centers,Z_s,unit_cells,prim_vec::Matrix{T}) where T
  #we create the supercell + a bounding box so that we don't get boudary effects in the calculations
  tot_cells = supercell+2
  model_origin = -(div(tot_cells-1,2))*Point3D{T}(prim_vec[1,:]+prim_vec[2,:]+prim_vec[3,:])
  #we store the origin of each of the unit cells
  #right now indexing runs from 1->tot_cell, would it be nicer to have -div(tot_cells,2):div(tot_cells,2)?
  cell_origins = [model_origin+Point3D{T}(a*prim_vec[1,:]+b*prim_vec[2,:]+c*prim_vec[3,:]) for a=0:tot_cells-1,b=0:tot_cells-1,c=0:tot_cells-1]
  #convert relative unit_cells to absolute ones
  abs_unit_cells = [cell + div(tot_cells-1,2) for cell in unit_cells]
  abs_centers = [model_origin + center + Point3D{T}(unit_cell[1]*prim_vec[1,:]+unit_cell[2]*prim_vec[2,:]+unit_cell[3]*prim_vec[3,:]) for (center,unit_cell) in zip(centers,abs_unit_cells)]
  wfcs = [SimWfc3D(character[1],character[2],character[3],unit_cell,center,Z) for (character,unit_cell,center,Z) in zip(characters,unit_cells,abs_centers,Z_s)]
  return SimModel3D(cell_origins,prim_vec,wfcs)
end
