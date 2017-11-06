
function construct1S(Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    R = 2*Z^1.5 * exp(-Z*norm(point-center))
    Y = sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

function construct2S(Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    rho = 2*Z*norm(point-center)/2
    R = 1/(2*sqrt(2)) *(2 - rho) * Z^1.5 *exp(-rho/2)
    Y = sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

function construct3S(Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    rho = 2*Z*norm(point-center)/3
    R = 1/(9*sqrt(3)) *(6-6rho+rho^2) * Z^1.5 *exp(-rho/2)
    Y = sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

function construct4S(Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    rho = 2*Z*norm(point-center)/4
    R = 1/(96) *(24-36rho+12rho^2-rho^3) * Z^1.5 *exp(-rho/2)
    Y = sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

function construct5S(Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    rho = 2*Z*norm(point-center)/5
    R = 1/(300sqrt(5)) *(120-240rho+120rho^2-20rho^3+rho^4) * Z^1.5 *exp(-rho/2)
    Y = sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

function construct2P(direction,Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    rho = 2*Z*norm(point-center)/2
    R = 1/2sqrt(6) * rho * Z^1.5 *exp(-rho/2)
    Y = sqrt(3)*getfield(point-center,direction)/norm(point-center)*sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

function construct3P(direction,Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    rho = 2*Z*norm(point-center)/3
    R = 1/9sqrt(6) * rho*(4-rho) * Z^1.5 *exp(-rho/2)
    Y = sqrt(3)*getfield(point-center,direction)/norm(point-center)*sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

function construct4P(direction,Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    rho = 2*Z*norm(point-center)/4
    R = 1/(32*sqrt(15)) * rho*(20-10rho+rho^2) * Z^1.5 *exp(-rho/2)
    Y = sqrt(3)*getfield(point-center,direction)/norm(point-center)*sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

function construct5P(direction,Z,center::Point3D{T},cell::Array{Point3D{T}},dims,supercells) where T
  n_a,n_b,n_c = dims
  origin = -(div(supercells-1,2)+1)*sum(cell)
  a_vec,b_vec,c_vec = cell.*supercells
  points = Array{WfPoint3D{T},3}(n_a,n_b,n_c)
  R,Y = zero(T),zero(T)
  for (ia,a) in enumerate(linspace(0,1,n_a)),(ib,b) in enumerate(linspace(0,1,n_b)),(ic,c) in enumerate(linspace(0,1,n_c))
    point = origin+(a*a_vec+b*b_vec+c*c_vec)
    rho = 2*Z*norm(point-center)/5
    R = 1/(150*sqrt(30)) * rho*(120-910rho+18rho^2-rho^3) * Z^1.5 *exp(-rho/2)
    Y = sqrt(3)*getfield(point-center,direction)/norm(point-center)*sqrt(1/4pi)
    points[ia,ib,ic] = WfPoint3D{T}(R*Y+0.0im,origin+point)
  end
  return Wfc3D{T}(points,cell,PhysAtom(T))
end

@inline function S(point::Point3D{T},n,center::Point3D{T},Z)::T where T
  rho = 2*Z*norm(point-center)/n
  Y = sqrt(1/4pi)
  if n==1
    R = 2*Z^1.5 * exp(-Z*norm(point-center))
  elseif n==2
    R = 1/2sqrt(2) *(2 - rho) * Z^1.5 *exp(-rho/2)
  elseif n==3
    R = 1/9sqrt(3) *(6-6*rho+rho^2) * Z^1.5 *exp(-rho/2)
  elseif n==4
    R = 1/96 * (24-36rho+12rho^2-rho^3) * Z^1.5 *exp(-rho/2)
  elseif n==5
    R = 1/300sqrt(5) * (120-240rho+120rho^2-20rho^3+rho^4) * Z^1.5 *exp(-rho/2)
  end
  return R*Y
end

@inline function P(point::Point3D{T},n,kind,center::Point3D{T},Z)::T where T
  rho = 2*Z*norm(point-center)/n
  if n==2
    R = 1/2sqrt(6) * rho * Z^1.5 *exp(-rho/2)
    Y = sqrt(3)*getfield(point-center,kind)/norm(point-center)*sqrt(1/4pi)
  elseif n==3
    R = 1/9sqrt(6) * rho * (4-rho) * Z^1.5 *exp(-rho/2)
    Y = sqrt(3)*getfield(point-center,kind)/norm(point-center)*sqrt(1/4pi)
  elseif n==4
    R = 1/32*sqrt(15) * rho * (20-10rho+rho^2) * Z^1.5 *exp(-rho/2)
    Y = sqrt(3)*getfield(point-center,kind)/norm(point-center)*sqrt(1/4pi)
  elseif n==5
    R = 1/150sqrt(30) * rho * (120-910rho+18rho^2-rho^3) * Z^1.5 *exp(-rho/2)
    Y = sqrt(3)*getfield(point-center,kind)/norm(point-center)*sqrt(1/4pi)
  end
  return Y*R
end

@inline function wave_equation(point::Point3D{T},wfc::SimWfc3D{T})::T where T
  if wfc.character[3] == :s
    return S(point,wfc.character[2],wfc.center,wfc.character[1])
  elseif wfc.character[3] == :p
    return P(point,wfc.character[2],wfc.character[4],wfc.center,wfc.character[1])
  end
end