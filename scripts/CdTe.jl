using DFControl

# CdTe = create_job("CdTe",phd_dir*"CdTe/NSOC",default_inputs[:scf],default_inputs[:bands],server_dir = "CdTe/NSOC")
cdte = load_server_job("CdTe/NSOC",phd_dir*"CdTe/NSOC")
add_calculation!(cdte,deepcopy(get_input(cdte,"scf")))
cdte.calculations[end].filename = "nscf2.in"

add_calculation!(cdte,deepcopy(get_input(cdte,"nscf2")),filename="nscf3.in")
print_info(cdte)
change_flags!(cdte,"nscf.in",:calculation => "'nscf'")
change_flags!(cdte,"nscf3.in",:gdir => 3)
change_data!(cdte,["nscf1.in","nscf2.in","nscf3.in"],:k_points, gen_k_grid(12,12,12,:nscf))
change_flow!(cdte,"scf.in"=>true,"bands.in"=>false,"projwfc"=>false)
change_data_option!(cdte,"nscf",:k_points,:crystal)
remove_flags!(cdte,:smearing)
change_flags!(cdte,:occupations=>"'fixed'",:degauss=>0.0f0)
atoms = Dict(:Te => [Point3D{Float32}(0.3333333,0.6666667,0.40000),Point3D{Float32}(0.6666667,.3333333,0.87500)],:Cd =>[Point3D{Float32}(0.3333333,0.6666667,0.0000000),Point3D{Float32}(0.6666667,0.3333333,0.5000000)])
change_atoms!(cdte,new_atoms,pseudo_set=:pbesol,pseudo_fuzzy="paw")
change_data!(CdTe,"bands",:k_points,[[0.5f0,0.5f0,0.0f0,100],[0.0f0,0.0f0,0.0f0,100f0],[0.5f0,0.0f0,0.0f0,1]])
cell = Float32[0.8660254  -0.5000000   0.0000000;0.0000000   1.0000000   0.0000000;0.0000000   0.0000000   1.6367615]*Float32[4.5700000 0.0 0.0 ;0.0 4.5700000 0.0;0.0 0.0 7.4800000]
change_cell_parameters!(CdTe,cell)
remove_flags!(CdTe,:nbnd)
add_flags!(CdTe,:system,:A => 1.0f0)
change_data_option!(CdTe,:cell_parameters,:alat)
replace_header_word!(cdte,"frontend","defpart")
submit_job(cdte)

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
orig_atoms = Dict(:Te => [Point3D{Float32}(0.3333333,0.6666667,0.375000),Point3D{Float32}(0.6666667,.3333333,0.87500)],:Cd =>[Point3D{Float32}(0.3333333,0.6666667,0.0000000),Point3D{Float32}(0.6666667,0.3333333,0.5000000)])
new_atoms = copy(orig_atoms)
for (key,val) in orig_atoms
    if key == :Te
        println(val[1])
        for i=1:length(val)

        new_atoms[key][i] = val[i] .+ Point3D{Float32}(0.03,0.03,0.03)
        end
    end
end
change_atoms!(CdTe_soc,atoms,pseudo_set=:pbesolrel,pseudo_fuzzy="paw")
add_flags!(CdTe_soc,:system,:lspinorb=>true,:noncolin=>true)
CdTe_soc.server_dir = "CdTe/SOC/"
CdTe_soc.local_dir  = phd_dir*"CdTe/SOC/"
change_flow!(CdTe_soc,[("scf",true),("bands",true)])
submit_job(CdTe_soc)

outputs = pull_outputs(CdTe_soc)
plot(read_qe_bands_file(outputs[2]),fermi=read_fermi_from_qe_file(ou7puts[1]),ylims=[-5,5],legend=false)