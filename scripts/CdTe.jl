using DFControl

CdTe = create_job("CdTe",phd_dir*"CdTe/NSOC",default_inputs[:scf],default_inputs[:bands],server_dir = "CdTe/NSOC")
atoms = Dict(:Te => [Point3D{Float32}(0.3333333,0.6666667,0.40000),Point3D{Float32}(0.6666667,.3333333,0.87500)],:Cd =>[Point3D{Float32}(0.3333333,0.6666667,0.0000000),Point3D{Float32}(0.6666667,0.3333333,0.5000000)])
change_atoms!(CdTe,atoms,pseudo_set=:pbesol,pseudo_fuzzy="paw")
change_data!(CdTe,"bands",:k_points,[[0.5f0,0.5f0,0.0f0,100],[0.0f0,0.0f0,0.0f0,100f0],[0.5f0,0.0f0,0.0f0,1]])
cell = Float32[0.8660254  -0.5000000   0.0000000;0.0000000   1.0000000   0.0000000;0.0000000   0.0000000   1.6367615]*Float32[4.5700000 0.0 0.0 ;0.0 4.5700000 0.0;0.0 0.0 7.4800000]
change_cell_parameters!(CdTe,cell)
remove_flags!(CdTe,:nbnd)
add_flags!(CdTe,:system,:A => 1.0f0)
change_data_option!(CdTe,:cell_parameters,:alat)
submit_job(CdTe)

outputs = pull_outputs(CdTe,extras=["*s_j0.5*"])
using  Plots
plot(read_qe_bands_file(outputs[2]),fermi=read_fermi_from_qe_file(outputs[1]),ylims=[-5,5],legend=false)
set_default_input(:scf,get_input(CdTe,"scf"))
set_default_input(:bands,get_input(CdTe,"bands"))

set_default_input(:projwfc,read_qe_input(DFControl.assets_dir*"inputs/qe/projwfc.in",run_command="mpirun -np 24 ~/bin/projwfc.x"))
add_calculation!(CdTe,default_inputs[:projwfc])
change_flow!(CdTe,[("scf",false),("bands",false)])
heatmap(read_qe_kpdos(outputs[3])[1])



CdTe_soc = deepcopy(CdTe)
atoms = Dict(:Te => [Point3D{Float32}(0.3333333,0.6666667,0.385000),Point3D{Float32}(0.6666667,.3333333,0.86500)],:Cd =>[Point3D{Float32}(0.3333333,0.6666667,0.0000000),Point3D{Float32}(0.6666667,0.3333333,0.5000000)])
change_atoms!(CdTe_soc,atoms,pseudo_set=:pbesolrel,pseudo_fuzzy="paw")
add_flags!(CdTe_soc,:system,:lspinorb=>true,:noncolin=>true)
CdTe_soc.server_dir = "CdTe/SOC/"
CdTe_soc.local_dir  = phd_dir*"CdTe/SOC/"
change_flow!(CdTe_soc,[("scf",true),("bands",true)])
submit_job(CdTe_soc)

outputs = pull_outputs(CdTe_soc)
plot(read_qe_bands_file(outputs[2]),fermi=read_fermi_from_qe_file(outputs[1]),ylims=[-5,5],legend=false)