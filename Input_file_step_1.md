# displacement simulation, lammps input for relaxation# 
variable read_file1 string    /home/mosab/resources/data.v_bcc26266   #ok
variable pot_file1 string    /home/mosab/resources/Fesdasdas_DD-BN.eam.fs  #ok
variable atom_symbol1 string mo   #ok
variable temp equal        30  #ok

units  metal
boundary  p p p
atom_style atomic
atom_modify map array
read_data   ${read_file1}

pair_style eam/fs
pair_coeff * * ${pot_file1} ${atom_symbol1}

group tungusten  type 1

velocity  all create  ${temp}  1 mom  yes rot yes  dist gaussian

timestep 0.001 #corresponding to 1 fs

thermo_style custom vol
thermo 100

dump 1 all custom 10 dump.w-0${temp}k id type x y z vx vy vz
 
fix 1 all npt temp ${temp} ${temp} 0.1 iso 0.0 0.0 1.0 
restart 50000 restart.w-${temp}K
run 100000
minimize 1.0e-16 1.0e-16  2000 6000
