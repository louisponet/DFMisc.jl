using LaTeXStrings

using Plots
using DFWannier

T=Float64
cell_param = 4.3392*T[ 0.8424444 -0.4863855 -3.7e-8;-1.5e-8 0.972771 6.0e-9;-5.7e-8 6.0e-9 1.5226926]
atom_pos   = [cell_param*[1/3 ,2/3, 0.3194],cell_param*[2/3,1/3,0.7363],cell_param*[0.0,0.0,0.0]]
# atoms = [[PhysAtom(Point3D(atom_pos[1]),0.210752) for i=1:4];[PhysAtom(Point3D(atom_pos[2]),0.296619) for i=1:4];[PhysAtom(Point3D(atom_pos[3]), -0.770879) for i=1:4]]
atoms = [[PhysAtom(Point3D(atom_pos[1]),-0.0462138) for i=1:4];[PhysAtom(Point3D(atom_pos[2]),0.355636) for i=1:4];[PhysAtom(Point3D(atom_pos[3]), -0.768408) for i=1:4]]
BiTeI = WannierModel{T}("/home/ponet/Documents/PhD/BiTeI/NSOC/","/home/ponet/Documents/PhD/BiTeI/NSOC/2BiTeI_bands.out",atoms);
# x = WannierModel{T}("/Users/ponet/Documents/Fysica/PhD/BiTeI/NSOC/","/Users/ponet/Documents/Fysica/PhD/BiTeI/NSOC/2BiTeI_bands.out",atoms);

calculate_angmom(BiTeI.wfcs[6],BiTeI.wfcs[7])
calculate_angmom(BiTeI.wfcs[6],BiTeI.wfcs[8])

tbbands_soc_tot = calculate_eig_soc(x)
df_bands_soc = read_qe_bands_file("/home/ponet/Documents/PhD/BiTeI/SOC/2BiTeI_bands.out",T)
plot(df_bands_soc,tbbands_soc_tot,ylim=[5,12])
plot(df_bands_soc,ylim=[-5,5],fermi=8.1176)
minimum(df_bands_soc[39].eigvals)-maximum(df_bands_soc[38].eigvals)
for (i,wfc) in enumerate(x.wfcs)
  wfc.atom = atoms[i]
end

tbbands_tot = calculate_eig_cm(x);
df_bands = read_qe_bands_file("/home/ponet/Documents/PhD/BiTeI/NSOC/2BiTeI_bands.out",T)
tbbands    = calculate_eig_cm_angmom(x,90:0.2:110);
tbbandssoc = calculate_eig_cm_angmom_soc(x,test_poitns);
tbbandssoc = calculate_eig_cm_angmom_soc(x,90:0.2:110);

plot(tbbandssoc)
)plot(df_bands_soc[21:44]
#----------all the plots--------
plot(plot([tbbandssoc[9:10]...,tbbands[5]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(2.5:0.2:4)),plot([tbbandssoc[9:10]...,tbbands[5]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)"), plot([tbbandssoc[9:10]...,tbbands[5]],:angmom1_x,label = ["SOC" "SOC" "No SOC"],leg=false,title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(-1:0.2:0.2)),plot([tbbandssoc[9:10]...,tbbands[5]],:spin1_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[9:10]...,tbbands[5]],:angmom1_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title="",yticks=collect(-0.2:0.2:0.6)),plot([tbbandssoc[9:10]...,tbbands[5]],:spin1_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"),layout=(3,2))
plot(plot([tbbandssoc[5:6]...,tbbands[3]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(3.4:0.2:5)),plot([tbbandssoc[5:6]...,tbbands[3]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)",yticks=collect(-3:0.2:-1)), plot([tbbandssoc[5:6]...,tbbands[3]],:angmom2_x,leg=false,label = ["SOC" "SOC" "No SOC"],title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""]),yticks=collect(-0.5:0.2:1)),plot([tbbandssoc[5:6]...,tbbands[3]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[5:6]...,tbbands[3]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title=""),plot([tbbandssoc[5:6]...,tbbands[3]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"),layout=(3,2))
plot(plot([tbbandssoc[7:8]...,tbbands[4]],:cm_z,label=["SOC" "SOC" "No SOC"],leg=false,xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[7:8]...,tbbands[4]],:eigvals,fermi=9.2879,label=["SOC" "SOC" "No SOC"],xticks=([-0.05,0.0,0.05],["","",""]),yguide=L"E-E_f"*" (eV)"), plot([tbbandssoc[7:8]...,tbbands[4]],:angmom2_x,leg=false,label = ["SOC" "SOC" "No SOC"],title="OAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[7:8]...,tbbands[4]],:spin2_x,fermi=9.2879,leg=false,title="SAM around Te",xticks=([-0.05,0.0,0.05],["","",""])),plot([tbbandssoc[7:8]...,tbbands[4]],:angmom2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),xguide=L"|$\mathbf{k}_r$|",title=""),plot([tbbandssoc[7:8]...,tbbands[4]],:spin2_y,fermi=9.2879,leg=false,xticks=([-0.05,0.0,0.05],[0.05,"Z",0.05]),title="",xguide=L"|$\mathbf{k}_r$|"),size=(1200,1200),guidefont=font(20,"DejaVu Sans"),titlefont=font(20,"DejaVu Sans"),layout=(3,2))
#-----------tbbands vs dfbands-----------------#

dfbands = read_bands_file("/Users/ponet/Documents/Fysica/PhD/BiTeI/NSOC/02BiTeI_bands.out",Float64);
plot(dfbands[8:20],tbbands,legend=false)
outputs = pull_outputs(BiTeI_soc)

BiTeI_nsoc = load_server_job("ponet@10.255.9.115","BiTeI/NSOC","/home/ponet/Documents/PhD/BiTeI/NSOC")
BiTeI_soc = load_server_job("ponet@10.255.9.115","BiTeI/SOC","/home/ponet/Documents/PhD/BiTeI/SOC")
 print_flow(BiTeI_nsoc)
add_calculation!(BiTeI_nsoc,default_inputs[:projwfc],4)
#  BiTeI_nsoc = load_job(phd_dir*"BiTeI/NSOC",server = default_server,server_dir = "BiTeI/NSOC")
 BiTeI_soc = load_server_job("ponet@10.255.9.115","BiTeI/SOC","/home/ponet/Documents/PhD/BiTeI/SOC")
atomic_positions = get_data(BiTeI_soc,"bands",:atomic_positions)
Te_pos = atomic_positions[:Te]
I_pos  = atomic_positions[]
atomic_positions[:Te] = I_pos
atomic_positions[:I] = Te_pos
atomic_positions[:I] = [Point3D{Float32}(1/3,2/3,0.3076)]
atomic_positions[:Te] = [Point3D{Float32}(2/3,1/3,0.7482)]
change_data!(BiTeI_soc, ["bands","scf"],:atomic_positions, atomic_positions)
cell_param = Float32[0.842444445  -0.486385500  -0.000000037;-0.000000015   0.972770975   0.000000006;-0.000000057   0.000000006   1.522692506]
change_flags!(BiTeI_nsoc,Dict(:dis_win_min => 0.f0,:dis_win_max => 10.f0, :dis_froz_min => 0.f0,:dis_froz_max => 10.f0))
outputs = pull_outputs(BiTeI_nsoc,extras=["*.xsf","*r.dat"])

bands = read_qe_bands_file(BiTeI_nsoc.local_dir*"2BiTeI_bands.out")
using Plots


change_flow!(BiTeI_soc,[("vc",false),("nscf",false),("bands",true)])
fermi = read_fermi_from_qe_file(BiTeI_nsoc.local_dir*"1BiTeI_scf.out")
set_flags!(job_nsoc,:projwfc,Dict(:kresolveddos=>true))
set_flags!(job_nsoc,:system,Dict(:occupations => "'smearing'",:smearing => "'cold'",:degauss => 0.01))
change_flow!(job_nsoc,[(5,true),(6,false),(7,false)])
print_flag(job_nsoc,:A)
print_flow(job
plot(bands,ylim=[-10,15])
submit_job(BiTeI_nsoc)

zmat_Bi,ticks = read_qe_kpdos(BiTeI_nsoc.local_dir*"pwo.pdos_atm#3(Bi)_wfc#2(p)")
zmat_Te,ticks = read_qe_kpdos(BiTeI_nsoc.local_dir*"pwo.pdos_atm#2(Te)_wfc#2(p)")
zmat_I,ticks = read_qe_kpdos(BiTeI_nsoc.local_dir*"pwo.pdos_atm#1(I)_wfc#2(p)")
plot(heatmap(zmat_Bi,title = "Bi_p"), heatmap(zmat_Te, title="Te_p"), heatmap(zmat_I,title = "I_p"),ticks=ticks)