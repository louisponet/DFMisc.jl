using DFControl
# using DFVisualize
using Plots
bialo3 = load_server_job("BiAlO3/SOC",phd_dir*"BiAlO3/SOC")
bialo3.server_dir = "BiAlO3/SOC/tetra/"
bialo3.local_dir  = "BiAlO3/SOC/tetra/"
cell = Float32[1.0000000   0.0000000   0.0000000;0.0000000   1.0000000   0.0000000;0.0000000   0.0000000   1.0315789]
change_flags!(bialo3,:A => 3.8f0)
change_cell_parameters!(bialo3,cell)
atoms = Dict{Symbol,Array{Point3D{Float32},1}}(:Bi => Point3D{Float32}[Float32[0,0,0]],:Al => Point3D{Float32}[Float32[0.5,0.5,0.4247]],:O => Point3D{Float32}[Float32[0.5,0.5,0.8899],Float32[0.5,0.,0.365],Float32[0.,0.5,0.365]])
change_atoms!(bialo3,atoms)
push_job(bialo3)
submit_job(bialo3)
plot(read_qe_output(outputs[2])[:bands],fermi=10.4557,ylims=[-8,8],legend=false)
test_outputs = pull_outputs(bialo3_nsoc)
gui()
bialo3= load_server_job("BiAlO3/SOC",phd_dir*"BiAlO3/SOC")
savefig("/Users/ponet/Downloads/tmp.png")
bialo3_nsoc = load_server_job("BiAlO3/NSOC",phd_dir*"BiAlO3/NSOC")
bialo3_nsoc.name ="BiAlO3_NSOC"
bialo3.name = "BiAlO3_SOC"
print_data(bialo3,:atomic_species)
print_data(bialo3,:k_points)
print_data(bialo3,:cell_parameters)
print_data(bialo3,:atomic_positions)
change_flags!(bialo3_nsoc,:ecutwfc=>45)
cell = get_data(bialo3,"bands",:cell_parameters)
add_flags!(bialo3_nsoc,:system,:occupations => "'smearing'",:degauss=>0.01f0)
cell[3,1:3] = Float32[0,0,3.92/3.80]
change_cell_parameters!(bialo3_nsoc,cell)
atoms = get_data(bialo3,"bands",:atomic_positions)
atoms[:Al] = [Point3D{Float32}(0.5,0.5,0.4247)]
atoms[:O]  = Point3D{Float32}[[0.5f0,0.5f0,0.8899f0],[0.5f0,0.0f0,0.3650f0]]
change_atoms!(bialo3_nsoc,atoms,pseudo_set=:sssp)
bialo3_test = deepcopy(bialo3_nsoc)
bialo3_test.server_dir = "BiAlO3/NSOC/test/"
bialo3_test.local_dir = phd_dir*"BiAlO3/NSOC/test/"
bialo3_test.name = "BiAlO3_nsoc_test"
change_flags!(bialo3,:A=>5.37546f0)
add_flags!(bialo3_test,:system,:smearing=>"'cold'")
print_flags(bialo3)
atoms = get_data(read_qe_input("/Users/ponet/Downloads/01BiAlO3_scf.in"),:atomic_positions)
cell = Float32[ 0.577350269189626   0.000000000000000   0.830523651309221;-0.288675134594813   0.500000000000000   0.830523651309221;-0.288675134594813  -0.500000000000000   0.830523651309221]
change_cell_parameters!(bialo3_test,cell)
atoms=Dict{Symbol,Array{Point3D{Float32},1}}( )

change_atoms!(bialo3_test,atoms)
change_flow!(bialo3,"scf"=>true,"bands"=>true)
print_data(bialo3_test,:atomic_species)
k_path = get_data(bialo3,"bands",:k_points)
k_path[1] = Float32[0.0,0.0,0.5,100]
k_path[2] = Float32[0.0,0.5,0.5,100]
k_path[3] = Float32[0.5,0.5,0.5,1]
change_data!(bialo3_nsoc,"bands",:k_points,k_path)
change_flow!(bialo3,"_scf"=>true,"_nscf"=>false,"wan"=>false,"bands"=>true)
print_flags(bialo3_nsoc)
submit_job(bialo3)
change_cell_parameters!(bialo3,cell)
change_atoms!(bialo3,atoms)
outputs=pull_outputs(bialo3)
using Plots

# BiAlO3_nsoc = load_server_job("ponet@10.255.9.115","BiAlO3/NSOC","/home/ponet/Documents/PhD/BiAlO3/NSOC")
BiAlO3_nsoc = load_job("/home/ponet/Documents/PhD/BiAlO3/NSOC",server="ponet@10.255.9.115",server_dir="BiAlO3/NSOC")
outputs = pull_outputs(BiAlO3_nsoc,extras=["*.xsf"])

zmat_Bi,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#2(Bi)_wfc#2(p)")
zmat_Al,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#1(Al)_wfc#2(p)")
zmat_O1,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#3(O)_wfc#2(p)")
zmat_O2,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#4(O)_wfc#2(p)")
zmat_O3,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#5(O)_wfc#2(p)")
zmat_Bi,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#2(Bi)_wfc#1(s)")
zmat_Al,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#1(Al)_wfc#1(s)")
zmat_O1,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#3(O)_wfc#1(s)")
zmat_O2,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#4(O)_wfc#1(s)")
zmat_O3,ticks = read_qe_kpdos(BiAlO3_nsoc.local_dir*"pwo.pdos_atm#5(O)_wfc#1(s)")
change_flags!(BiAlO3_nsoc,Dict(:dis_froz_min => 4.0f0,:dis_froz_max => 22.5f0))
set_flags!(BiAlO3_nsoc,Dict(:dis_win_min => -5.0f0))
using Plots
plot(heatmap(zmat_Bi,title="Bi_p"),heatmap(zmat_Al,title="Al_p"),heatmap(zmat_O1,title="O1_p"),heatmap(zmat_O2,title="O2_p"),heatmap(zmat_Bi,title="O3_p"),ticks=ticks)
