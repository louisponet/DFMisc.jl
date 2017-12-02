using DFControl
set_default_pseudo_dir(:pbesolrel,"pseudos/pbesolrel/PSEUDOPOTENTIALS/")
set_default_pseudo_dir(:pbesol,"pseudos/pbesol/PSEUDOPOTENTIALS/")
configure_default_pseudos()
HfO2 = create_job("HfO2",phd_dir*"HfO2/NSOC",server_dir = "HfO2/NSOC")

pull_file("GeTe_2/rel/GeTe_scf.in",phd_dir*"GeTe/fullrel/")
scf_calc = read_qe_input(phd_dir*"GeTe/fullrel/GeTe_scf.in")
scf_calc.name 
change_atoms!(scf_calc,Dict(:Hf => [Point3D{Float32}(0.9668,0.7337,0.1231)],:O=>[Point3D{Float32}(0.1350,0.0666,0.2672),Point3D{Float32}(0.7288,0.4633,0.8744)]),pseudo_set = :pbesolrel,pseudo_fuzzy="paw")
set_default_server("ponet@10.255.9.115")