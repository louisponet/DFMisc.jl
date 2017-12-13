using DFControl

# CdTe = create_job("CdTe",phd_dir*"CdTe/NSOC",default_inputs[:scf],default_inputs[:bands],server_dir = "CdTe/NSOC")
cdte = load_job(phd_dir*"CdTe/NSOC")
pull_outputs(cdte)
add_calculation!(cdte,deepcopy(get_input(cdte,"scf")),filename="nscf1.in")
change_data!(cdte,["nscf1.in","nscf2.in","nscf3.in"],:k_points,[6,6,6,1,1,1])
add_flags!(cdte,["nscf1.in","nscf2.in","nscf3.in"],:control,:nppstr=>6)
change_data_option!(cdte,["nscf1.in","nscf2.in","nscf3.in"],:k_points,:automatic)
print_flags(cdte,"nscf3")
remove_flags!(cdte,Symbol("starting_magnetization(2)"))
cdte.calculations[end].filename = "nscf2.in"
add_flags!(cdte,"nscf1.in",:control,:lberry=>true,:gdir=>1,:nppstr=>12)
add_calculation!(cdte,deepcopy(get_input(cdte,"nscf1")),filename="nscf3.in")
add_flag
print_info(cdte)
change_flags!(cdte,:nppstr => 10)
change_flags!(cdte,"nscf3.in",:gdir => 1)
add_flags!(cdte,"nscf1.in",:control,:nppstr => 10)
change_data!(cdte,["nscf3.in"],:k_points, gen_k_grid(6,6,6,:nscf))
print_flow(cdte)
cdte.server_dir = "CdTe/NSOC/"
print_flow(cdte)
change_flow!(cdte,"nscf1.in"=>true,"nscf2.in"=>true,"scf"=>true)
change_data_option!(cdte,"nscf",:k_points,:crystal)
remove_flags!(cdte,"nscf3.in",:nppstr)
change_flags!(cdte,:occupations=>"'fixed'",:degauss=>0.0f0)
atoms = Dict(:Te => [Point3D{Float32}(0.3333333,0.6666667,0.40000),Point3D{Float32}(0.6666667,.3333333,0.87500)],:Cd =>[Point3D{Float32}(0.3333333,0.6666667,0.0000000),Point3D{Float32}(0.6666667,0.3333333,0.5000000)])
change_atoms!(cdte,atoms,pseudo_set=:pbesol,pseudo_fuzzy="paw")
change_data!(CdTe,"bands",:k_points,[[0.5f0,0.5f0,0.0f0,100],[0.0f0,0.0f0,0.0f0,100f0],[0.5f0,0.0f0,0.0f0,1]])
cell = Float32[0.8660254  -0.5000000   0.0000000;0.0000000   1.0000000   0.0000000;0.0000000   0.0000000   1.6367615]*Float32[4.5700000 0.0 0.0 ;0.0 4.5700000 0.0;0.0 0.0 7.4800000]
change_cell_parameters!(CdTe,cell)
remove_flags!(CdTe,:nbnd)
add_flags!(CdTe,:system,:A => 1.0f0)
change_data_option!(CdTe,:cell_parameters,:alat)
replace_header_word!(cdte,"frontend","defpart")
cdte.server_dir = "CdTe/NSOC/"
change_flow!(cdte,"nscf2.in"=>true)
change_flags!(cdte,:ecutwfc=>35)
submit_job(cdte)
cdte.server_dir

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
for x=0.1:0.1:0.9
  for (key,val) in orig_atoms
    if key == :Te
      for i=1:length(val)
        new_atoms[key][i] = val[i] .+ Point3D{Float32}(0.0,0.0,x*0.06)
      end
    end
  end
  change_atoms!(cdte,new_atoms,pseudo_set=:pbesol,pseudo_fuzzy="paw")
  cdte.server_dir = "CdTe/NSOC/berryz$x/"
  cdte.local_dir  = phd_dir*"CdTe/NSOC/berryz$x/"
  submit_job(cdte)
end
cdte.local_dir  = phd_dir*"CdTe/NSOC/"
cdte.server_dir = "CdTe/NSOC/"
outputs = pull_outputs(cdte)
polarizations = Array{Tuple{Point3D{Float64},Float64},1}()
for x=0.1:0.1:0.8
  cdte.local_dir  = phd_dir*"CdTe/NSOC/berryz$x/"
  cdte.server_dir = "CdTe/NSOC/berryz$x/"
  files = pull_outputs(cdte)[end-2:end]
  t     = Point3D{Float64}()
  mod   = zero(Float64)
  for f in files
    t_,mod =read_qe_polarization(f,Float64)
    t+=t_
  end
  push!(polarizations,(t,mod))
end
plot(getfield.(getindex.(polarizations,1),:x))
add_flags!(CdTe_soc,:system,:lspinorb=>true,:noncolin=>tru
CdTe_soc.server_dir = "CdTe/SOC/"
CdTe_soc.local_dir  = phd_dir*"CdTe/SOC/"
change_flow!(CdTe_soc,[("scf",true),("bands",true)])
submit_job(CdTe_soc)
cdte_soc = load_server_job("CdTe/SOC/",phd_dir*"CdTe/SOC/")
print_data(cdte_soc,"bands")
change_atoms!(cdte_soc,orig_atoms)
outputs = pull_outputs(cdte_soc)
plot(read_qe_bands_file(outputs[2]),fermi=read_fermi_from_qe_file(outputs[1]),ylims=[-0.2,0.2],legend=false)
print_flow(cdte_soc)
submit_job(cdte_soc)