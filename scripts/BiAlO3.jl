using DFControl

const phd_dir = "/home/ponet/Documents/PhD/"
job_soc = load_server_job("ponet@10.255.9.115","BiAlO3/SOC","/home/ponet/Documents/PhD/BiAlO3/SOC")
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
