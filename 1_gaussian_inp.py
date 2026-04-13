""" This code will generate 
        1) gaussian com files for all geometric files in a folder or specified .xyz file
        2) an array slurm job submission file 
        3) a configuration file for the array job
    Following this, the jobs will be automatically submitted.

    The user can adjust 
        1) the slurm_template for their HPC system
        2) the calculation parameters in the main function"""

import os, sys, glob, argparse

def get_xyz_list(input_path):
    if os.path.isfile(input_path):
        return [os.path.splitext(os.path.basename(input_path))[0]], os.path.dirname(input_path)
    elif os.path.isdir(input_path):
        xyz_files = [file for file in os.listdir(input_path) if file.endswith(".xyz")]
        xyz_list = [os.path.splitext(f)[0] for f in xyz_files]
        return xyz_list, input_path
    else:
        raise ValueError("Input must be a valid .xyz file or directory")

def create_input_files(project_path, xyz_path, xyz_list, nproc, memory, jobtype, functional, basis_set):
    #folder and xyz setup
    gaussian_path = os.path.join(f"{project_path}/gaussian/")
    os.makedirs(gaussian_path, exist_ok=True)
    com_path = os.path.join(f"{gaussian_path}/com/")
    os.makedirs(com_path, exist_ok=True)
    chk_path = os.path.join(f"{gaussian_path}/chk/")
    os.makedirs(chk_path, exist_ok=True)
    
    for index,xyz_file in enumerate(xyz_list):
        sys.stdout.write(f"\rReading {xyz_file}.xyz and writing Gaussian {xyz_file}.com ({index + 1}/{len(xyz_list)})                              ")
        sys.stdout.flush()
        
        #read .xyz file
        try:
            with open(f"{xyz_path}/{xyz_file}.xyz", "r") as f:
                lines = f.readlines()
            if len(lines) < 3:
                raise ValueError("The xyz file does not have the expected format (at least 3 lines).")
                #Read in from line 3 onwards, as the first two lines in a standard .xyz contain number of atoms and a comment
            xyz_coords = "".join(lines[2:]).strip()
        except Exception as e:
            print(f"Error reading {xyz_path}/{xyz_file}: {e}")
            return
    
        #gaussian input files        
        #gaussian template
        gaussian_template = f"""%nprocshared={nproc}
%mem={memory}GB
%chk={chk_path}{xyz_file}.chk
#p {jobtype}{functional}/{basis_set}

{xyz_file}

0, 1
{xyz_coords}

"""
        
        with open(f"{com_path}{xyz_file}.com", "w") as com_file:
            com_file.write(gaussian_template)
            
    
    print(f"\nGaussian .com files created in {com_path}\n")
    
    return com_path, gaussian_path

def create_sh_file(project_path, xyz_list, com_path, nproc, partition):
    #folder setup
    slurm_output_path = os.path.join(f"{project_path}/slurm_out/")
    os.makedirs(slurm_output_path, exist_ok=True)
    
    #configuration file for array job             
    with open (f"{project_path}/config_gaussian.txt", "w") as config_file:
        config_file.write("ArrayTaskID    Sample\n")
        for index, molecule in enumerate(xyz_list):
            config_file.write(f"{index}               {molecule}\n")
    
    #.sh file
    slurm_template = f"""#!/bin/bash
# Propagate environment variables to the compute node
#SBATCH --export=ALL
#SBATCH --partition={partition}
#SBATCH --account=your_account
#SBATCH --ntasks={nproc}
#SBATCH --time=48:00:00
#SBATCH --job-name=gaussian
#SBATCH --output={slurm_output_path}gaussian-%j.out

#SBATCH --array=0-{len(xyz_list)-1}
config={project_path}/config_gaussian.txt
molecule=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {{print $2}}' $config)
 
module purge
module load gaussian
 
#=========================================================
# Prologue script to record job details 
/opt/software/scripts/job_prologue.sh
#----------------------------------------------------------
 
g16 {com_path}/${{molecule}}.com
 
#=========================================================
# Epilogue script to record job endtime and runtime 
/opt/software/scripts/job_epilogue.sh
#----------------------------------------------------------
"""

        
    with open(f"{project_path}/gaussian.sh", "w") as f:
        f.write(slurm_template)
    print(f"SLURM script and config file created in {project_path}")

def main():
    #Define parameters here
    project_path = os.getcwd()
    
    nproc=16
    memory=72
    partition="standard"
    jobtype="opt freq " #optimisation = "opt " (note the space after opt), single point = ""
    functional="b3lyp"
    basis_set="6-31G*"

    #cli setup
    parser = argparse.ArgumentParser(description="Process XYZ files")
    parser.add_argument("input", help="Path to .xyz file or directory containing .xyz files")
    args = parser.parse_args()

    #create files
    xyz_list, xyz_path = get_xyz_list(args.input)
    com_path, gaussian_path=create_input_files(project_path, xyz_path, xyz_list, nproc, memory, jobtype, functional, basis_set)
    create_sh_file(project_path, xyz_list, com_path, nproc, partition)

    os.system(f"sbatch {project_path}/gaussian.sh")
        
if __name__ == "__main__":
    main()

