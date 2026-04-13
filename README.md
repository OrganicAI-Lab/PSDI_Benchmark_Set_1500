# PSDI_Benchmark_Set_1500
This repository contains the code associated with the publication "Multireference Excited-State Screening Reveals Hidden Candidate Space in Organic Semiconductors" by Malin Zollner and Tahereh Nematiaram, Department of Pure and Applied Chemistry, University of Strathclyde, 295 Cathedral Street, Glasgow G1 1XL, UK. 
The code is released under the Creative Commons Attribution 4.0 International (CC BY 4.0) licence. Any use of the code or associated data must be accompanied by appropriate citation of the publication.

0. Required input: 
        1) Copy the 4 python codes (1_gaussian_inp.py, 2_gaussian_error.py, 3_orca_inp.py, 4_cclib.py) into your desired working directory.
        2) In your terminal, go to your desired working directory.
        3) Have a folder with xyz files or the path to a specific .xyz file ready.
1. Run "python 1_gaussian_inp.py molecule.xyz" or "python 1_gaussian_inp.py path/to/molecule/folder" in your terminal. This code will generate 
        1) gaussian com files for all geometric files in a folder or specified .xyz file.
        2) an array slurm job submission file .
        3) a configuration file for the array job.
    Following this, the jobs will be automatically submitted.

    The user can adjust 
        1) the slurm_template for their HPC system.
        2) the calculation parameters in the main function.


2. Once calculations are completed run "python 2_gaussian_error.py" in your terminal. This code will determine if any calculations produced the errors l103, l502 or l9999 and attempt to fix the error by generating
        1) gaussian com files for log files with error l103, l502 and l9999.
        2) an array slurm job submission file for these error com files .
        3) a configuration file for this error files array job.
    Following this, the jobs will be automatically submitted.

    The user can adjust 
        1) the slurm_template for their HPC system.
        2) the calculation parameters in the main function.
    
    If the calculations did not produce an error, they are ready for the next step instead.


3. Once all calculations are completed and successful, run "python 3_orca_inp.py" in your terminal. This code will generate 
        1) fchk for gaussian chk files.
        2) optimised xyz for these generated fchk files.
        3) orca input files for all optimised xyz.
        4) an array slurm job submission file.
        5) a configuration file for the array job.
    Following this, the jobs will be automatically submitted.

    The user can adjust 
        1) the slurm_template for their HPC system.
        2) the calculation parameters in the main function and the orca_template.


4. Run "python 4_cclib.py" in your terminal. This code will generate 
        1) Excited_states_energy.csv with ID, job cpu time, SMILES, InChI, molecular formula, number of atoms, and E(S1), E(S2), E(T1), E(T2), f(S1), f(S2) for each theory level.

