using DFWannier
using Optim

function change_lsoc(model,indices,lambda)
  for wfc in model.wfcs[indices[1] : indices[2]]
    wfc.atom = PhysAtom(wfc.atom.center,lambda)
  end
end

function optimize_l()
  T = Float64
  
  cell_param = 4.3392*T[ 0.8424444 -0.4863855 -3.7e-8;-1.5e-8 0.972771 6.0e-9;-5.7e-8 6.0e-9 1.5226926]
  atom_pos   = [cell_param*[1/3 ,2/3, 0.3194],cell_param*[2/3,1/3,0.7363],cell_param*[0.0,0.0,0.0]]
  atoms = [[PhysAtom(Point3D(atom_pos[1]),0.0) for i=1:4];[PhysAtom(Point3D(atom_pos[2]),-0.8) for i=1:4];[PhysAtom(Point3D(atom_pos[3]),-1.0) for i=1:4]]
  x = WannierModel{T}("/home/ponet/Documents/PhD/BiTeI/NSOC/","/home/ponet/Documents/PhD/BiTeI/NSOC/2BiTeI_bands.out",atoms);
  
  dfbands = read_qe_bands_file("/home/ponet/Documents/PhD/BiTeI/SOC/2BiTeI_bands.out",T)[21:44]

  lambda1_s = -0.462138
  lambda2_s = 0.355636
  lambda3_s = -0.768408
  function f(lambdas)
    change_lsoc(x,[1,4],lambdas[1])
    change_lsoc(x,[5,8],lambdas[2])
    change_lsoc(x,[9,12],lambdas[3])
    tbbands = calculate_eig_soc_bloch(x)
    out_diff = 0.0
    for i=1:length(tbbands)
      for (calc,exact) in zip(tbbands[i].eigvals[80:120],dfbands[i].eigvals[80:120])
        out_diff += abs(calc-exact)
      end
    end
    return out_diff
  end
  return optimize(f,[lambda1_s,lambda2_s,lambda3_s],show_trace = true,show_every = 1,extended_trace=true)
end
results = optimize_l()
println(results)
