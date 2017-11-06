using DFTools

filename = ARGS[1]
output = ARGS[2]
curdir = pwd()

if filename[1]!="/"
  filename = curdir*"/"*filename
end

if output[1]!="/"
  output = curdir*"/"*output
end

zmat,yticks = read_KPDOS(filename)
heatmap(zmat_ge,legend=false,color=:BuPu,yticks=(yticks[1][5:5:end]+10,["","","",""]),xticks=([0,101,201],["A","Z","U"])
savefig(output)
