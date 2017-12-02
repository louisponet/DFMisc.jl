using Plots
using LaTeXStrings
using DFWannier
pyplot()
plot_font=font(15,"DejaVu Sans")

pyplot(lab="",yguidefont=plot_font,ytickfont=plot_font,xtickfont=plot_font,legendfont=plot_font)
T=Float64
# x = create_TB_model("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/","/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out",[[PhysAtom(0.0,0.0,-0.0239129,-0.155854) for i=1:4]...,[PhysAtom(0.0,0.0,5.5540692,0.318205) for i=1:4]...],Float64);

# x = WannierModel{T}("/home/ponet/Documents/PhD/GeTe/NSOC/paperxsf/","/home/ponet/Documents/PhD/GeTe/SOC/GeTe_bands.out",[[PhysAtom(T[0.0,0.0,-0.0239129,-0.155854]...) for i=1:4]...,[PhysAtom(T[0.0,0.0,5.5540692,0.318205]...) for i=1:4]...]);
x = WannierModel{T}("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/test1","/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out",[[PhysAtom(T[0.0,0.0,-0.0239129,-0.155854]...) for i=1:4]...,[PhysAtom(T[0.0,0.0,5.5540692,0.318205]...) for i=1:4]...]);
calculate_angmom(x.wfcs[6],x.wfcs[8])
points = [WfcPoint3D(-p.w,p.p) for p in x.wfcs[8].points]
x.wfcs[8] = Wfc3D(points,x.wfcs[8].cell,x.wfcs[8].atom)
x.wfcs[1].cell
kz = 0.6049734
ky_max = 1.0
kx_max = 1.0
kx_max = 1.184212/10
ky_max = 1.0255575/10

recip = [1.184212  0.000000  0.403316;-0.592106  1.025558  0.403316;-0.592106 -1.025558  0.403316]
kp_cart = [[kx,ky,kz] for kx=linspace(-kx_max,kx_max,30),ky=linspace(-ky_max,ky_max,30)]
kp_rec = [inv(recip')*kc for kc in kp_cart]
kp_rec = reshape(kp_rec,(900,1))
k_input = reshape(kp_rec,900,1)
k_plot = reshape(kp_cart,(length(kp_cart)))
tbbands=calculate_eig_cm_angmom(x,k_input);
tbbandssoc=calculate_eig_cm_angmom_soc(x,k_input);

kxs = [k[1] for k in k_plot]
kys = [k[2] for k in k_plot]
# kzs = [k[3] for k in k_input]
u = Float64[]
v = Float64[]
c= Float64[]
# w = Float[]
for (i,spinx) in enumerate(tbbandssoc[10].spins)
  lx,ly,lz = Array(spinx[2])
  push!(u,lx/50)
  push!(v,ly/50)
  push!(c,lz)
end
println(maximum(c))
quiver(kxs,kys,quiver=(u,v),color=c)
plot(tbbandssoc[6],:angmom2_x)
test_points = [[ka,kb,kc] for (ka,kb,kc) in zip(linspace(0.5,0.25,201), linspace(1.0,1.0,201),linspace(0.0,0.25,201))]
dfbands = read_qe_bands_file("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/rest/bands.out",T)
# dfbands = read_qe_bands_file("/home/ponet/Documents/PhD/GeTe/NSOC/xsftest/GeTe_bands_notrot.out",T)
dfbandssoc = read_qe_bands_file("/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out",T);
# dfbandsnsoc = read_qe_bands_file("/home/ponet/Documents/PhD/GeTe/NSOC/bands.out",T)
# potential=read_potential_file("/home/ponet/Documents/PhD/GeTe/NSOC/p_1.0/density.out")
tbbandssocf=calculate_eig_soc(x)
tbbandsf = calculate_eig_angmom(x)
tbbands=calculate_eig_cm_angmom(x,90:0.2:110);
tbbands=calculate_eig_cm_angmom(x,k_input);
tbbandssoc = calculate_eig_cm_angmom_soc(x,90:0.2:110);
tbbandssoc = calculate_eig_cm_angmom_soc(x);
tbbandssoc = calculate_eig_cm_angmom_soc(x,k_input);
plot(plot(tbbandssoc[9],:angmom2_x),plot(tbbands[5],:angmom2_x),plot(tbbandssoc[9],:angmom2_y),plot(tbbands[5],:angmom2_y))
#----- paper plots
# two small and 2 big
plot(plot([tbbandssoc[9:10]...,tbbands[5]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[9:10]...,tbbands[5]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)"), plot(plot([tbbandssoc[9:10]...,tbbands[5]],:angmom2_x,label = ["SOC" "SOC" "No SOC"],leg=false,title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[9:10]...,tbbands[5]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title=""),layout=(2,1)), plot(plot([tbbandssoc[9:10]...,tbbands[5]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[9:10]...,tbbands[5]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),layout=(2,1)),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"))
plot(plot([tbbandssoc[5:6]...,tbbands[3]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[5:6]...,tbbands[3]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)"), plot(plot([tbbandssoc[5:6]...,tbbands[3]],:angmom2_x,label = ["SOC" "SOC" "No SOC"],leg=false,title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[5:6]...,tbbands[3]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title=""),layout=(2,1)), plot(plot([tbbandssoc[5:6]...,tbbands[3]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[5:6]...,tbbands[3]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),layout=(2,1)),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"))
#all same size
plot(plot([tbbandssoc[9:10]...,tbbands[5]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(2.5:0.2:4)),plot([tbbandssoc[9:10]...,tbbands[5]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)"), plot([tbbandssoc[9:10]...,tbbands[5]],:angmom2_x,label = ["SOC" "SOC" "No SOC"],leg=false,title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(-1:0.2:0.2)),plot([tbbandssoc[9:10]...,tbbands[5]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[9:10]...,tbbands[5]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title="",yticks=collect(-0.2:0.2:0.6)),plot([tbbandssoc[9:10]...,tbbands[5]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"),layout=(3,2))
savefig("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/rest/oamvseigvalsv1.pdf")
plot(plot([tbbandssoc[5:6]...,tbbands[3]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(3.4:0.2:5)),plot([tbbandssoc[5:6]...,tbbands[3]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)",yticks=collect(-3:0.2:-1)), plot([tbbandssoc[5:6]...,tbbands[3]],:angmom2_x,leg=false,label = ["SOC" "SOC" "No SOC"],title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(-0.5:0.2:1)),plot([tbbandssoc[5:6]...,tbbands[3]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[5:6]...,tbbands[3]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title=""),plot([tbbandssoc[5:6]...,tbbands[3]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"),layout=(3,2))
savefig("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/rest/oamvseigvalsv3.pdf")
plot(plot([tbbandssoc[7:8]...,tbbands[4]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[7:8]...,tbbands[4]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)"), plot([tbbandssoc[7:8]...,tbbands[4]],:angmom2_x,leg=false,label = ["SOC" "SOC" "No SOC"],title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[7:8]...,tbbands[4]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[7:8]...,tbbands[4]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title=""),plot([tbbandssoc[7:8]...,tbbands[4]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"),layout=(3,2))
savefig("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/rest/oamvseigvalsv2.pdf")

#-----------plot 3 bandstructures-------------9.2879
plot(plot(dfbands[11:19],tbbandsf,title="No SOC",label="",fermi=9.2879,xticks=([0,101,201],["A","Z","U"]),subplot_index=1,ylim=[-13,9]),plot(dfbands[21:37],tbbandssoc,title="SOC",yl="",fermi=9.2879,xticks=([0,101,201],["A","Z","U"]),subplot_index=2,ylim=[-13,9],leg=false),plot(dfbands[21:37],tbbandssoc,title="SOC",yl="",fermi=9.2879,xticks=([0,101,201],["A","Z","U"]),subplot_index=3,ylim=[-5,5]),titlefont = font(30,"Courier"),size=(1800,600),layout=(1,3))
plot(plot(dfbands[11:19],tbbandsf,title="No SOC",label="",fermi=9.2879,subplot_index=1,ylim=[-13,9]),plot(dfbandssoc[21:37],tbbandssoc,title="SOC",yl="",fermi=9.2879,subplot_index=2,ylim=[-13,9],leg=false),plot(dfbandssoc[21:37],tbbandssoc,title="SOC",yl="",fermi=9.2879,subplot_index=3,ylim=[-5,5]),titlefont = font(30,"Courier"),size=(1800,600),layout=(1,3))
#-----------plot heatmaps-----------------
begin
  T=Float64
  dfbandssoc = read_bands_file("/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out",T);
  zmat_te,yticks=read_KPDOS("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/rest/pwo.pdos_atm#1(Te)_wfc#2(p)",fermi=9.2879);
  zmat_ge,yticks=read_KPDOS("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/rest/pwo.pdos_atm#2(Ge)_wfc#2(p)",fermi=9.2879);
  # zmat_ge,yticks=read_KPDOS("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/pwo.pdos_atm#2(Ge)_wfc#2(p)",fermi=9.2879)
  k_lim = [-norm(dfbandssoc[1].k_points[101]-dfbandssoc[1].k_points[1]),0.0,norm(dfbandssoc[1].k_points[101]-dfbandssoc[1].k_points[end])]
end
plot(plot(dfbandssoc,fermi=9.2879,color=:red,xticks=(k_lim,["A","Z","U"]),ylims=(-5-9.2879,18.01-9.2879),width=1,label="",title=""),heatmap(zmat_ge,legend=false,color=:BuPu,yticks=(yticks[1][5:5:end]+10,["","","",""]),xticks=([0,101,201],["A","Z","U"])),heatmap(zmat_te,color=:BuPu,xticks=([0,101,201],["A","Z","U"]),yticks=(yticks[1][5:5:end]+10,["","","",""])),layout=(1,3),size=(1200,500),right_margin=5mm)
savefig("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/paperxsf/rest/heatmap_bs.png")
#-----
recip = [1.184212  0.000000  0.403316;
        -0.592106  1.025558  0.403316;
        -0.592106 -1.025558  0.403316]
k_cart = []
for k in x.k_points
  push!(k_cart,recip'*k)
end
x.k_points[end]
k_cart[end]
#------Calculation of H_eff----------------------------------------------------------------#
function Hd(c,Ez,band)
  recip = [1.184212  0.000000  0.403316;
          -0.592106  1.025558  0.403316;
          -0.592106 -1.025558  0.403316]
  out = []
  for i in eachindex(band.k_points)
    push!(out,c*Ez*(recip'*band.k_points[i])[2]*(band.spins[i][2].x+band.angmoms[i][2].x))
  end
  return out
end

function Hso(l_so,band)
  out = []
  for i in eachindex(band.k_points)
    push!(out,l_so*(band.angmoms[i][2].x*band.spins[i][2].x))
  end
  return out
end

H0 = tbbands[5].eigvals

H10=H0+Hd(0.1,-0.5,tbbandssoc[10])+Hso(0.29,tbbandssoc[10])
H9=H0+Hd(0.1,-0.5,tbbandssoc[9])+Hso(0.29,tbbandssoc[9])
H5=tbbands[3].eigvals+Hd(9.0,-0.5,tbbandssoc[5])+Hso(0.29,tbbandssoc[5])
H6=tbbands[3].eigvals+Hd(9.0,-0.5,tbbandssoc[6])+Hso(0.29,tbbandssoc[6])
plot([H10,H9,tbbandssoc[6].eigvals,tbbandssoc[5].eigvals],line=[:blue :red])
plot(H0)
#---------------Proof of J rather than L------------#
\lambda
tbbandssoc[10].angmoms[50][2]
tbbandssoc[9].angmoms[50][2]
tbbandssoc[6].angmoms[50][2]
tbbandssoc[5].angmoms[50][2]
tbbandssoc[10].cms[49]-tbbandssoc[10].cms[50]
tbbandssoc[9].cms[49]-tbbandssoc[9].cms[50]
tbbandssoc[6].cms[49]-tbbandssoc[6].cms[50]
tbbandssoc[5].cms[49]-tbbandssoc[5].cms[50]
tbbandssoc[10].angmoms[50][2]+tbbandssoc[10].spins[50][2]
tbbandssoc[9].angmoms[50][2]+tbbandssoc[9].spins[50][2]
tbbandssoc[6].angmoms[50][2]+tbbandssoc[6].spins[50][2]
tbbandssoc[5].angmoms[50][2]+tbbandssoc[5].spins[50][2]

 #-------------CM stuff-------------------
Ge_s_b = construct_bloch_sum(x.wfcs[1],[0.5,0.5,0.5])
Ge_pz_b = construct_bloch_sum(x.wfcs[2],[0.5,0.5,0.5])
Ge_px_b = construct_bloch_sum(x.wfcs[3],[0.5,0.5,0.5])
Ge_py_b = construct_bloch_sum(x.wfcs[4],[0.5,0.5,0.5])
Te_s_b = construct_bloch_sum(x.wfcs[5],[0.5,0.5,0.5])
Te_pz_b = construct_bloch_sum(x.wfcs[6],[0.5,0.5,0.5])
Te_px_b = construct_bloch_sum(x.wfcs[7],[0.5,0.5,0.5])
Te_py_b = construct_bloch_sum(x.wfcs[8],[0.5,0.5,0.5])
Ge_s_b = construct_bloch_sum(x.wfcs[1],[0.,0.,0.])
Ge_pz_b = construct_bloch_sum(x.wfcs[2],[0.4,0.55,0.5])
Ge_px_b = construct_bloch_sum(x.wfcs[3],[0.4,0.55,0.5])
Ge_py_b = construct_bloch_sum(x.wfcs[4],[0.4,0.55,0.5])
Te_s_b = construct_bloch_sum(x.wfcs[5],[0.4,0.55,0.5])
Te_pz_b = construct_bloch_sum(x.wfcs[6],[0.4,0.55,0.5])
Te_px_b = construct_bloch_sum(x.wfcs[7],[0.4,0.55,0.5])
Te_py_b = construct_bloch_sum(x.wfcs[8],[0.4,0.55,0.5])
write_xsf_file("/Users/ponet/Downloads/test.xsf",Ge_s_b)
Ge_s = x.wfcs[1]
Ge_pz = x.wfcs[2]
Ge_px = x.wfcs[3]
Ge_py = x.wfcs[4]
Te_s = x.wfcs[5]
Te_pz = x.wfcs[6]
Te_px = x.wfcs[7]
Te_py = x.wfcs[8]
calculate_cm(Te_pz,Te_pz)

#---------------Eq 11-------------
term1,term2,wfc = eq11(x,9,98)
write_dipole_mesh("/Users/ponet/Downloads/term2test.xsf",term2,:z)
write_xsf_file("/Users/ponet/Downloads/wfctest.xsf",wfc)
#--------------different polarizations------

x = create_TB_model("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.5/","/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out",[[PhysAtom(0.0,0.0,-0.0239129,-0.155854) for i=1:4]...,[PhysAtom(0.0,0.0,5.5540692,0.318205) for i=1:4]...],Float64);
tbbandssocf,tmp=calculate_eig_SOC(x);
tbbands,tmp=calculate_eig_cm(x);
tbbands,dip1=calculate_eig_cm_angmom(x,90:0.2:110);
tbbandssoc,dip= calculate_eig_cm_angmom_SOC(x,90:0.2:110);
dfbands = read_qe_bands_file("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.0/bands.out",Float64);
plotlyjs()
dfbandssoc = read_qe_bands_file("/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out");
plot(plot([tbbandssoc[9:10]...,tbbands[5]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(0:0.2:10)),plot([tbbandssoc[9:10]...,tbbands[5]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)"), plot([tbbandssoc[9:10]...,tbbands[5]],:angmom2_x,label = ["SOC" "SOC" "No SOC"],leg=false,title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(-1:0.2:1)),plot([tbbandssoc[9:10]...,tbbands[5]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[9:10]...,tbbands[5]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title="",yticks=collect(-1:0.2:1)),plot([tbbandssoc[9:10]...,tbbands[5]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"),layout=(3,2))
savefig("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.5/oamvseigvalsv1.pdf")
plot(plot([tbbandssoc[5:6]...,tbbands[3]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(0:0.2:10)),plot([tbbandssoc[5:6]...,tbbands[3]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)",yticks=collect(-3:0.2:-1)), plot([tbbandssoc[5:6]...,tbbands[3]],:angmom2_x,leg=false,label = ["SOC" "SOC" "No SOC"],title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(-1:0.2:1)),plot([tbbandssoc[5:6]...,tbbands[3]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[5:6]...,tbbands[3]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title=""),plot([tbbandssoc[5:6]...,tbbands[3]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"),layout=(3,2))
savefig("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.5/oamvseigvalsv3.pdf")
plot(dfbandssoc[21:37],tbbandssocf,title="SOC",yl="",fermi=9.2879,subplot_index=2,ylim=[-13,9],leg=false)
dfbands = read_bands_file("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.0/bands.out")
plot(tbbandssocf,legend=false,ylim=[2,10])
pyplot()
begin
  T=Float64
  dfbandssoc = read_bands_file("/Users/ponet/Documents/Fysica/PhD/GeTe/fullrel/GeTe_bands.out",T);
  zmat_te,yticks=read_KPDOS("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.0/pwo.pdos_atm#1(Te)_wfc#2(p)",fermi=9.2879);
  zmat_ge,yticks=read_KPDOS("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.0/pwo.pdos_atm#2(Ge)_wfc#2(p)",fermi=9.2879);
  # zmat_ge,yticks=read_KPDOS("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/pwo.pdos_atm#2(Ge)_wfc#2(p)",fermi=9.2879)
  k_lim = [-norm(dfbandssoc[1].k_points[101]-dfbandssoc[1].k_points[1]),0.0,norm(dfbandssoc[1].k_points[101]-dfbandssoc[1].k_points[end])]
end
zmat_ges,yticks=read_KPDOS("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.0/pwo.pdos_atm#2(Ge)_wfc#1(s)",fermi=9.2879);
zmat_tes, = read_KPDOS("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.0/pwo.pdos_atm#1(Te)_wfc#1(s)",fermi=9.2879);
plot(heatmap(zmat_ge),heatmap(zmat_te),heatmap(zmat_ges),heatmap(zmat_tes))
savefig("/Users/ponet/Documents/Fysica/PhD/GeTe/colin/p_0.0/heatmap_p.png")

#----------------total angmom-----------------#
#notes: Weird things happen with the soc calculation. Nonetheless the non-soc calculation seems to work out fine, and results
#       in a much larger Lx than is seen around the atoms. The same dynamics seem to still apply though. The issue with the sign
#       of Ly still remains
x = create_TB_model("/home/ponet/Documents/PhD/GeTe/NSOC/paperxsf/","/home/ponet/Documents/PhD/GeTe/SOC/GeTe_bands.out",[[PhysAtom(0.0,0.0,-0.0239129,-0.155854) for i=1:4]...,[PhysAtom(0.0,0.0,5.5540692,0.318205) for i=1:4]...],Float64);

tbbands=calculate_eig_cm_angmom(x,90:0.2:110);
tbbandssoc= calculate_eig_cm_angmom_SOC(x,90:0.2:110);

total_angmom = []
total_angmom_soc = []
for i = 1:length(tbbands[5].eigvec)
  vec = tbbands[5].eigvec[i]
  vec_soc = tbbandssoc[10].eigvec[i]
  wfc = sum([v*w for (v,w) in zip(vec,x.wfcs)])
  wfc_soc = sum([v*w for (v,w) in zip(vec_soc,[x.wfcs;x.wfcs])])
  center = calculate_cm(wfc,wfc)
  center_soc = calculate_cm(wfc_soc,wfc_soc)
  push!(total_angmom,calculate_angular_momentum(wfc,wfc,center))
  push!(total_angmom_soc,calculate_angular_momentum(wfc_soc,wfc_soc,center_soc))
end
total_angmom[1][4]

plot(map(x->real(x[1]),total_angmom))


using Images, GeometryTypes, GLVisualize, Reactive, GLWindow, Colors
using GLAbstraction, GLFW
import GLAbstraction: singlepressed, imagespace
import GLVisualize: moving_average

if !isdefined(:runtests)
    window = glscreen()
end

description = """
Drawing lines and then adding symmetrie to them.
"""

window.color = RGBA(0f0, 0f0, 0f0, 1f0)

const linebuffer = Signal(fill(Point2f0(NaN), 4))
const middle = Point2f0(widths(window)) / 2
const symmetries = 18
const unsmoothed = copy(value(linebuffer))


function angle_between(a, b)
    a = normalize(a)
    b = normalize(b)
    atan2(cross(a, b), dot(a, b))
end


symmetry_lines = foldp(Point2f0[], linebuffer) do v0, lines
    resize!(v0, symmetries * length(lines))
    i = 1;
    for fi = linspace(0, 2pi, symmetries)
        for point in lines
            diff = middle - point
            radius = norm(diff)
            _fi = (fi + angle_between(diff, Point2f0(0, 1))) % 2pi
            v0[i] = radius * Point2f0(
                sin(_fi + pi),
                cos(_fi + pi)
            ) + middle
            i += 1
        end
    end
    v0
end

const lineobj = visualize(
    symmetry_lines, :lines,
    color = RGBA(0.8f0, 0.8f0, 0.8f0, 1f0),
    thickness = 1.5f0,
    boundingbox = AABB{Float32}(value(window.area)) # boundingbox for center!
)
_view(lineobj, window, camera = :orthographic_pixel)

@materialize mouseposition, mouse_buttons_pressed, mouseinside, buttons_pressed = window.inputs

camera = window.cameras[:orthographic_pixel]
const history = Point2f0[]
s = map(mouseposition, mouse_buttons_pressed, init = nothing) do mp, mbp
    l0 = value(linebuffer)
    if singlepressed(mbp, GLFW.MOUSE_BUTTON_LEFT) && value(mouseinside) && isempty(value(buttons_pressed))
        p = imagespace(mp, camera)
        keep, p = moving_average(p, 1.5f0, history)
        if keep
            push!(linebuffer, push!(l0, p))
        end
    else
        if !isnan(last(l0)) # only push one NaN to seperate
            empty!(history) # reset
            push!(linebuffer, push!(l0, Point2f0(NaN)))
        end
    end
    nothing
end
# preserve signals, so that it doesn't get garbage collected.
preserve(s)
# we need to define init, because otherwise it will be initialised by calling
# the function one time, which

if !isdefined(:runtests)
    renderloop(window)
end

using GLVisualize, GeometryTypes, GLAbstraction
using Colors, Reactive, FileIO
if !isdefined(:runtests)
    window = glscreen()
    timesignal = bounce(linspace(0f0, 1f0, 360))
end

description = """
Showing off the flexibility of the particle system by animating
all kind of atributes for an arbitrary mesh as the particle.
"""

cat    = GLNormalMesh(loadasset("cat.obj"))
sphere = GLNormalMesh(Sphere{Float32}(Vec3f0(0), 1f0), 12)

function scale_gen(v0, nv)
    l = length(v0)
    @inbounds for i=1:l
        v0[i] = Vec3f0(1,1,sin((nv*l)/i))/2
    end
    v0
end
function color_gen(v0, t)
    l = length(v0)
    @inbounds for x=1:l
        v0[x] = RGBA{N0f8}(x/l,(cos(t)+1)/2,(sin(x/l/3)+1)/2.,1.)
    end
    v0
end

t            = const_lift(x->x+0.1, timesignal)
ps           = sphere.vertices
scale_start  = Vec3f0[Vec3f0(1,1,rand()) for i=1:length(ps)]
scale        = foldp(scale_gen, scale_start, t)
colorstart   = color_gen(zeros(RGBA{N0f8}, length(ps)), value(t))
color        = foldp(color_gen, colorstart, t)
rotation     = sphere.normals
cats = visualize((cat, ps), scale = scale, color = color, rotation = rotation)

_view(cats, window)

if !isdefined(:runtests)
    renderloop(window)
end