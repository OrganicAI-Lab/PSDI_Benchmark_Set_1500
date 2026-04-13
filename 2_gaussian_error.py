""" This code will generate 
        1) gaussian com files for log files with error l103, l502 and l9999
        2) an array slurm job submission file for these error com files 
        3) a configuration file for this error files array job
    Following this, the jobs will be automatically submitted.

    The user can adjust 
        1) the slurm_template for their HPC system
        2) the calculation parameters in the main function"""

import os, os.path, shutil, fileinput,sys

def gaussian_error_determination(project_path, nproc, memory, jobtype, functional, basis_set):  
    #folder and list setup
    com_path=f"{project_path}/gaussian/com/"
    log_path=f"{project_path}/gaussian/log/"
    os.makedirs(log_path, exist_ok=True)
    error_path=f"{project_path}/gaussian/error/"
    os.makedirs(error_path, exist_ok=True)
    
    molecule_list = [file[:-4] for file in os.listdir(f"{com_path}") if file.endswith(".log")]
    successful_molecule_list=[]; fixed_error_molecule_list=[]; error_molecule_list=[]; error_lists = {"103": [], "301": [], "301_H": [], "301_basis_set": [], "502": [], "9999": []}

    #define log files depending on errors in file or normal termination
    for molecule in range(len(molecule_list)):
        sys.stdout.write(f"\rSorting error files {molecule + 1}/{len(molecule_list)}        ")
        sys.stdout.flush()

        with open(f"{com_path}{molecule_list[molecule]}.log", "r") as log_file:
            log_file_content = log_file.read()

        lines = [line.strip() for line in log_file_content.splitlines() if line.strip()]
        last_line = lines[-1]

        if last_line.startswith("Normal termination of Gaussian"):
                shutil.move(f"{com_path}{molecule_list[molecule]}.log", f"{log_path}{molecule_list[molecule]}.log")
                successful_molecule_list.append(molecule_list[molecule])
                error="None"

        #l103 () error molecules (moving, fixing+listing) 
        elif ("Error termination via Lnk1e in /opt/software/gaussian/g16/l103.exe" in log_file_content or "Error termination via Lnk1e in /opt/software/gaussian/g16_avx/g16/l103.exe" in log_file_content):
            error="103"

        #l301 (missing H) or atom not in basis set error molecules, moving+listing !!ONLY!!) 
        elif ("Error termination via Lnk1e in /opt/software/gaussian/g16/l301.exe" in log_file_content or "Error termination via Lnk1e in /opt/software/gaussian/g16_avx/g16/l301.exe" in log_file_content):
            if "The combination of multiplicity" in log_file_content:
                error="301_H"
            elif "Atomic number out of range" in log_file_content:
                error="301_basis_set"
            else:
                error="301"

        #l502 () error molecules (moving, fixing+listing) 
        elif ("Error termination via Lnk1e in /opt/software/gaussian/g16/l502.exe" in log_file_content or "Error termination via Lnk1e in /opt/software/gaussian/g16_avx/g16/l502.exe" in log_file_content):
            error="502"

        #l9999 (need longer/more cycles, only give extra cycles once, if not fixed the issue might be something else) error molecules (moving, fixing+listing) 
        elif ("Error termination via Lnk1e in /opt/software/gaussian/g16/l9999.exe" in log_file_content or "Error termination via Lnk1e in /opt/software/gaussian/g16_avx/g16/l9999.exe" in log_file_content):
            error="9999"

        #other error molecules molecules and molecules still running (moving+listing !!ONLY!!)                                                                                 
        else:
            #shutil.move(f"{com_path}{molecule_list[molecule]}.log", f"{error_path}{molecule_list[molecule]}.log")
            error_molecule_list.append(molecule_list[molecule])
            error="other"

        #move and list log file based on error
        if error in ["103", "301", "301_H", "301_basis_set", "502", "9999"]:
            error_subfolder_path=f"{error_path}/l{error}/"
            os.makedirs(error_subfolder_path, exist_ok=True)
            shutil.move(f"{com_path}{molecule_list[molecule]}.log", f"{error_subfolder_path}{molecule_list[molecule]}.log")
            error_lists[error].append(molecule_list[molecule])
                    
            #edit current energy calculation to try fix the error and change chk to chk_2
            if error in ["103", "301_H", "502"]:
                with fileinput.input(f"{com_path}{molecule_list[molecule]}.com", inplace=True) as com_file:
                    for line in com_file:
                        if error == "103":
                            line = line.replace("opt", "opt=cartesian")
                        if error == "502":
                            line = line.replace(f" {functional}/{basis_set} ", f" {functional}/{basis_set} SCF=Fermi SCF=Noincfock SCF=QC ")
                        line = line.replace(f".chk", f"_2.chk")
                        print(line, end="")
            #error 9999 needs a new com file without the coordinates
            elif error == "9999":
                gaussian_template = f"""%nprocshared={nproc}
%mem={memory}GB
%oldchk={project_path}/gaussian/chk/{molecule_list[molecule]}.chk
%chk={project_path}/gaussian/chk/{molecule_list[molecule]}_2.chk
#p {jobtype}{functional}/{basis_set} guess=read geom=check

{molecule_list[molecule]}_2

0, 1

"""
                with open(f"{com_path}{molecule_list[molecule]}.com", "w") as com_file:
                    com_file.write(gaussian_template)

            #list the errors that were automatically attempted to be fixed and do not need manual inspection
            if error in ["103", "502", "9999"]:
                fixed_error_molecule_list.append(molecule_list[molecule])

    print(f"\nERROR l301 (missing H & other) - manually add H & check = {error_lists['301_H']} & {error_lists['301']}")
    print(f"ERROR l301 (atoms not in basis set) - adjust basis set = {len(error_lists['301_basis_set'])}")
    print(f"ERROR (other) - check log files in {com_path}= {error_molecule_list}")
    print(f"Number of l9999, l502 & l103 - already rerunning, check if fixed once job completed = {len(error_lists['9999'])}, {len(error_lists['502'])} & {len(error_lists['103'])}")
    print(f"Number of successful molecules = {len(successful_molecule_list)}")

    return fixed_error_molecule_list, com_path

def create_sh_file(project_path, molecule_list, com_path, nproc, error_name, partition):
    #folder setup
    slurm_output_path = os.path.join(f"{project_path}/slurm_out/")
    
    #configuration file for array job             
    with open (f"{project_path}/config_gaussian{error_name}.txt", "w") as config_file:
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
#SBATCH --job-name=gaussian{error_name}
#SBATCH --output={slurm_output_path}/gaussian{error_name}-%j.out

#SBATCH --array=0-{len(molecule_list)-1}
config={project_path}/config_gaussian{error_name}.txt
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

        
    with open(f"{project_path}/gaussian{error_name}.sh", "w") as f:
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

    fixed_error_molecule_list, com_path = gaussian_error_determination(project_path, nproc, memory, jobtype, functional, basis_set)
    if len(fixed_error_molecule_list) > 0:
        error_name = "_error" #for naming the jobs
        create_sh_file(project_path, fixed_error_molecule_list, com_path, nproc, error_name, partition)

        os.system(f"sbatch {project_path}/gaussian{error_name}.sh")
        
if __name__ == "__main__":
    main()
