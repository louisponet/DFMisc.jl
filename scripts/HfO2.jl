using DFControl
set_default_pseudo_dir(:pbesolrel,"/home/ponet/pseudos/pbesolrel/PSEUDOPOTENTIALS/")
set_default_pseudo_dir(:pbesol,"/home/ponet/pseudos/pbesol/PSEUDOPOTENTIALS/")
set_default_pseudo_dir(:sssp,"/home/ponet/pseudos/sssp/")
configure_default_pseudos()
set_default_input(:scf,read_qe_input(DFControl.assets_dir*"inputs/qe/scf.in",run_command="mpirun -np 24 ~/bin/pw.x"))
set_default_input(:bands,read_qe_input(DFControl.assets_dir*"inputs/qe/bands.in",run_command="mpirun -np 24 ~/bin/pw.x"))
set_default_input(:nscf,read_qe_input(DFControl.assets_dir*"inputs/qe/nscf.in",run_command="mpirun -np 24 ~/bin/pw.x"))
default_inputs[:scf]
atoms = get_data(HfO2,"scf",:atomic_positions)
change_atoms!(HfO2,atoms,pseudo_set=:sssp)
HfO2 = create_job("HfO2",phd_dir*"HfO2/NSOC",default_inputs[:scf],default_inputs[:bands],default_inputs[:projwfc],default_inputs[:nscf],server_dir = "HfO2/NSOC")
HfO2 = load_job(phd_dir*"HfO2/NSOC",server_dir="HfO2/NSOC")
pull_outputs(HfO2)
using Plots
bands = read_qe_bands_file(HfO2.local_dir*"/bands.out")
plot(read_qe_bands_file(HfO2.local_dir*"/bands.out"),fermi=read_fermi_from_qe_file(HfO2.local_dir*"scf.out"),ylims=[-2,2],legend=false)
remove_flags!(HfO2,:occupations,:smearing,:degauss)
submit_job(HfO2)
change_flow!(HfO2,[("nscf",false)])

using DFControl
hfo2    = load_server_job("HfO2/NSOC/",phd_dir*"HfO2/NSOC")
outputs = pull_outputs(hfo2)
using Plots
plot(read_qe_bands_file(outputs[2]),ylims=[-18,2])
print_flags(hfo2,"scf")
add_flags!(hfo2,:system,:nbnd=>40)
print_data(hfo2,"scf")
remove_flags!(hfo2,:occupations,:smearing,:degauss)
submit_job(hfo2)
hfo2.job_header