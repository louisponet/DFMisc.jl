using DFTools
using Optim

function optimize_l()
  x = create_TB_model("/home/ponet/GeTe_2/nonrel/","/home/ponet/GeTe_2/rel/GeTe_bands.out",[[PhysAtom(0.0,0.0,-0.0239129,-0.018) for i=1:4]...,[PhysAtom(0.0,0.0,5.5540692,0.2995182189) for i=1:4]...],Float64);
  # x = create_TB_model("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/","/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out",[[PhysAtom(0.0,0.0,-0.0239129,-0.018) for i=1:4]...,[PhysAtom(0.0,0.0,5.5540692,0.2995182189) for i=1:4]...],Float64);
  dfbands = read_bands_file("/home/ponet/GeTe_2/rel/GeTe_bands.out",Float64)[21:30]
  lambda1_s = -0.018
  lambda2_s = 0.297346
  function f(lambdas)
    changeLambdaSOC!(x,[1,4],lambdas[1])
    changeLambdaSOC!(x,[5,8],lambdas[2])
    tbbands,tmp = calculate_eig_SOC(x)
    out_diff = 0.0
    for i=1:10
      for (calc,exact) in zip(tbbands[i].eigvals,dfbands[i].eigvals)
        out_diff += abs(calc-exact)
      end
    end
    return out_diff
  end
  return optimize(f,[lambda1_s,lambda2_s],show_trace = true,show_every = 1,extended_trace=true,time_limit=3600)
end
results = optimize_l()
println(results)
