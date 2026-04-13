""" This code will generate 
        1) fchk for gaussian chk files
        2) optimised xyz for these generated fchk files
        3) orca input files for all optimised xyz
        4) an array slurm job submission file 
        5) a configuration file for the array job
    Following this, the jobs will be automatically submitted.

    The user can adjust 
        1) the slurm_template for their HPC system
        2) the calculation parameters in the main function and the orca_template"""


import os, sys, subprocess

def convert_chk_to_fchk(project_path):
    #set up folders
    log_path  = f"{project_path}/gaussian/log/"
    chk_path  = f"{project_path}/gaussian/chk/"
    fchk_path = f"{project_path}/gaussian/fchk/"
    os.makedirs(fchk_path, exist_ok=True)

    #Run gaussian in terminal to convert .chk to .fchk
    molecule_list = [file[:-4] for file in os.listdir(log_path) if file.endswith('.log')]

    for index, molecule in enumerate(molecule_list):
        sys.stdout.write(f"\rReading {molecule}.chk and writing {molecule}.fchk ({index + 1}/{len(molecule_list)})                              ")
        sys.stdout.flush()

        subprocess.run(["bash", "-lc", f"module load gaussian; formchk {chk_path}{molecule}.chk {fchk_path}{molecule}.fchk"], check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print(f".fchk created in {fchk_path}")

    return molecule_list, fchk_path

def convert_fchk_to_xyz(project_path, fchk_path, molecule_list):
    #setup folder
    xyz_path  = f"{project_path}/xyz_optimised/"
    os.makedirs(xyz_path, exist_ok=True)

    #run openbabel in terminal to convert .fchk to .xyz
    for index, molecule in enumerate(molecule_list):
        sys.stdout.write(f"\rReading {molecule}.fchk and writing {molecule}.xyz ({index + 1}/{len(molecule_list)})                              ")
        sys.stdout.flush()

        subprocess.run(["bash", "-lc", f"module load obabel; obabel -ifchk {fchk_path}{molecule}.fchk -oxyz -O {xyz_path}{molecule}.xyz"], check=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print(f".xyz created in {xyz_path}")

    return xyz_path

def create_input_files(project_path, xyz_path, molecule_list, active_space_electrons, active_space_orbitals, multiplicity, roots, nproc, maxcore):
    #setup folder
    orca_path = os.path.join(f"{project_path}/orca/")
    os.makedirs(orca_path, exist_ok=True)
    input_path = os.path.join(f"{orca_path}/inp/")
    os.makedirs(input_path, exist_ok=True)
    
    for index,molecule in enumerate(molecule_list):
        sys.stdout.write(f"\rReading {molecule}.xyz and writing Orca {molecule}.inp ({index + 1}/{len(molecule_list)})                            ")
        sys.stdout.flush()
        
        #read .xyz file
        with open(f"{xyz_path}{molecule}.xyz", 'r') as f:
            lines = f.readlines()
        xyz_coords = "".join(lines[2:]).strip()
    
        #orca input files        
        #orca template
        orca_template = f"""! def2-TZVPP AutoAux LARGEPRINT

%casscf
   nel {active_space_electrons}        #active space number of electrons
   norb {active_space_orbitals}       #active space number of orbitals
   MULT   {multiplicity}   #Multiplicity (Triplet, Singlet)
   NRoots {roots}   #Number of excited states for each multiplicity

   PTMethod SC_NEVPT2
   trafostep ri  #RI approximation for CASSCF and NEVPT2

   PTSettings
      NThresh 1e-6
      D4Step  Fly
      D4Tpre  1e-10
      D3Tpre  1e-14
      EWIN  -3,1000
      TSMallDenom 1e-2
      CanonStep 1  
      QDType  QD_VANVLECK  
   end
end

* xyz 0 1
{xyz_coords}
*

%maxcore {maxcore}

%pal
   nprocs {nproc}
end
"""
        #Complete the template with the desired values.
        #content = orca_template.format(active_space_electrons, active_space_orbitals, multiplicity, roots, xyz_coords=xyz_coords, nproc=nproc)
        
        with open(f"{input_path}{molecule}.inp", "w") as inp_file:
            inp_file.write(orca_template)
       
    print(f"\nOrca .inp files created in {input_path}\n")
    
    return input_path, orca_path

def create_sh_file(project_path, molecule_list, input_path, orca_path, nproc, partition):
    #folder setup
    slurm_output_path = os.path.join(f"{project_path}/slurm_out/")
    os.makedirs(slurm_output_path, exist_ok=True)
    output_path = os.path.join(f"{orca_path}out/")
    os.makedirs(output_path, exist_ok=True)
    
    #configuration file for array job            
    with open (f"{project_path}/config_orca.txt", "w") as config_file:
        config_file.write("ArrayTaskID    Sample\n")
        for index, molecule in enumerate(molecule_list):
            config_file.write(f"{index}               {molecule}\n")
    
    #.sh file
    slurm_template = f"""#!/bin/bash
# Propagate environment variables to the compute node
#SBATCH --export=ALL
#SBATCH --partition={partition}
#SBATCH --account=your-account
#SBATCH --ntasks={nproc}
#SBATCH --time=48:00:00
#SBATCH --job-name=orca
#SBATCH --output={slurm_output_path}orca-%j.out

#SBATCH --array=0-{len(molecule_list)-1}
config={project_path}/config_orca.txt
molecule=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {{print $2}}' $config)
 
module load orca/6.1.0
 
#=========================================================
# Prologue script to record job details 
/opt/software/scripts/job_prologue.sh
#----------------------------------------------------------
 
$ORCABIN/orca {input_path}${{molecule}}.inp > {output_path}${{molecule}}.out
 
#=========================================================
# Epilogue script to record job endtime and runtime 
/opt/software/scripts/job_epilogue.sh
#----------------------------------------------------------
"""
    with open(f"{project_path}/orca.sh", "w") as f:
        f.write(slurm_template)
    print(f"SLURM script and config file created in {project_path}")

def main():
    #Define parameters here
    project_path = os.getcwd()

    active_space_electrons = "8"
    active_space_orbitals = "8"
    multiplicity = "3,1" #triplet, singlet
    roots = "3,3" #number of excited states for each multiplicity
    nproc = 40
    maxcore = 3750
    partition = "standard"

    molecule_list, fchk_path = convert_chk_to_fchk(project_path)
    xyz_path = convert_fchk_to_xyz(project_path, fchk_path, molecule_list)
    input_path, orca_path = create_input_files(project_path, xyz_path, molecule_list, active_space_electrons, active_space_orbitals, multiplicity, roots, nproc, maxcore)
    create_sh_file(project_path, molecule_list, input_path, orca_path, nproc, partition)
    
    os.system(f"sbatch {project_path}/orca.sh")

if __name__ == "__main__":
    main()

