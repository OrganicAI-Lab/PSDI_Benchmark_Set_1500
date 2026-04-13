""" This code will generate 
        1) Excited_states_energy.csv with ID, job cpu time, SMILES, InChi, molecular formula, number of atoms, and E(S1),E(S2), E(T1), E(T2), f(S1), f(S2) for each theory level
        """

import cclib, sys, os, re
import pandas as pd

def process_out_files(project_path):
    out_path = f"{project_path}/orca/out/"
    if os.path.exists(f"{out_path}/.out"):
        os.remove(f"{out_path}/.out") #sometimes an empty file is created, this will remove it
    molecule_list = [file[:-4] for file in os.listdir(f"{out_path}") if file.endswith(".out")]
    
    all_data = []

    for index, molecule in enumerate(molecule_list):
        sys.stdout.write(f"\rExtracting results from {molecule}.out ({index + 1}/{len(molecule_list)})        ")
        sys.stdout.flush()
            
        with open(f"{out_path}{molecule}.out", "r") as f:
            full_text = f.read()
                
        #Extract CPU time, smiles, number of atoms and molecular formula
        total_cpu_time = run_time(full_text)
        smiles, number_of_atoms, formula, inchi=smiles_from_xyz(project_path, molecule)
        
        #Extract results
        casscf_energies = extract_energies(full_text, "SA-CASSCF TRANSITION ENERGIES", "DENSITY MATRIX")
        nevpt2_energies = extract_energies(full_text, "NEVPT2 TRANSITION ENERGIES", "NEVPT2 CORRECTION TO THE TRANSITION ENERGY")
        casscf_oscillator_strengths = extract_oscillator_strengths(full_text, "CASSCF UV, CD spectra")
        nevpt2_oscillator_strengths = extract_oscillator_strengths(full_text, "CASSCF (NEVPT2 diagonal energies) UV, CD spectra")

        #save data in dictionary for each molecule
        row = {
            "ID": molecule,
            "Job CPU Time (s)": total_cpu_time,
            "SMILES": smiles,
            "InChi": inchi,
            "formula": formula,
            "#atoms": number_of_atoms,
            "SA_CASSCF_E(S1)": casscf_energies["E(S1)"],
            "NEVPT2_E(S1)": nevpt2_energies["E(S1)"],
            "SA_CASSCF_E(S2)": casscf_energies["E(S2)"],
            "NEVPT2_E(S2)": nevpt2_energies["E(S2)"],
            "SA_CASSCF_E(T1)": casscf_energies["E(T1)"],
            "NEVPT2_E(T1)": nevpt2_energies["E(T1)"],
            "SA_CASSCF_E(T2)": casscf_energies["E(T2)"],
            "NEVPT2_E(T2)": nevpt2_energies["E(T2)"],
            "SA_CASSCF_f(S1)": casscf_oscillator_strengths["fosc_0-1A->1-1A"],
            "NEVPT2_f(S1)": nevpt2_oscillator_strengths["fosc_0-1A->1-1A"],
            "SA_CASSCF_f(S2)": casscf_oscillator_strengths["fosc_0-1A->2-1A"],
            "NEVPT2_f(S2)": nevpt2_oscillator_strengths["fosc_0-1A->2-1A"],
        }
        all_data.append(row)
    
    #save results for all molecules
    df = pd.DataFrame(all_data)
    df.to_csv(f"{project_path}/Excited_states_energy.csv", index=False)
    print(f"\nResults saved in {project_path}/Excited_states_energy.csv")

def run_time(full_text):
    total_run_match = re.search(
                r"TOTAL RUN TIME:\s*(\d+)\s*days\s*(\d+)\s*hours\s*(\d+)\s*minutes\s*(\d+)\s*seconds\s*([\d\.]+)\s*msec",
                full_text, re.IGNORECASE)
    if total_run_match:
        days    = int(total_run_match.group(1))
        hours   = int(total_run_match.group(2))
        minutes = int(total_run_match.group(3))
        secs    = int(total_run_match.group(4))
        msec    = float(total_run_match.group(5))
        total_cpu_time = days*86400 + hours*3600 + minutes*60 + secs + msec/1000
    else:
        total_cpu_time = None
    return total_cpu_time

def extract_energies(full_text, method_keyword, end_keyword=None):
    transition_pattern = re.compile(
        r"^\s*(\d+):\s+\d+\s+(\d+)\s+([\d\.\+\-Ee]+)\s+([\d\.\+\-Ee]+)\s+([\d\.\+\-Ee]+)",
        re.MULTILINE
    )

    #Searches for a block starting with method_keyword and, if provided, stops at end_keyword
    if method_keyword in full_text:
        section = full_text.split(method_keyword, 1)[1]
        if end_keyword and end_keyword in section:
            section = section.split(end_keyword, 1)[0]
        #Attempt to locate the header line of the table:
        header_match = re.search(r"STATE\s+ROOT\s+MULT", section)
        if header_match:
            block_text = section[header_match.start():]
        else:
            block_text = section
        matches = transition_pattern.findall(block_text)
        singlets = []
        triplets = []
        for match in matches:
            #match[0]: state number; match[1]: multiplicity; match[3]: energy (eV)
            state_num = int(match[0])
            multiplicity = int(match[1])
            try:
                energy_eV = float(match[3])
            except Exception:
                continue
            if multiplicity == 1:
                singlets.append((energy_eV, state_num))
            elif multiplicity == 3:
                triplets.append((energy_eV, state_num))
        singlets.sort(key=lambda x: x[0])
        triplets.sort(key=lambda x: x[0])
        energies = {"E(S1)": singlets[0][0], "E(S2)": singlets[1][0], "E(T1)": triplets[0][0], "E(T2)": triplets[1][0]}
        return energies
    else:
        return {"E(S1)": None, "E(S2)": None, "E(T1)": None, "E(T2)": None}

def extract_oscillator_strengths(full_text, method_header):
    result = {"fosc_0-1A->1-1A": None, "fosc_0-1A->2-1A": None}
    if method_header not in full_text:
        return result

    #Delimit the section corresponding to the method
    section = full_text.split(method_header, 1)[1]

    #Extract the absorption section until the CD block starts (or until the end of the text)
    absorption_match = re.search(
        r"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS(.*?)(?:CD SPECTRUM VIA|\Z)",
        section,
        re.DOTALL
    )
    if not absorption_match:
        return result

    block = absorption_match.group(1)
    #Pattern to extract the transition line and the fosc(D2) value
    fosc_pattern = re.compile(
        r"^\s*0-1A\s*->\s*(1-1A|2-1A)\s+[\d\.\+\-Ee]+\s+[\d\.\+\-Ee]+\s+[\d\.\+\-Ee]+\s+([\d\.\+\-Ee]+)",
        re.MULTILINE
)
    matches = fosc_pattern.findall(block)
    for trans, fosc in matches:
        result[f"fosc_0-1A->{trans}"] = float(fosc)
    return result

def smiles_from_xyz(project_path,molecule):
    from openbabel import openbabel

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "smi")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, f"{project_path}/xyz_optimised/{molecule}.xyz")

    mol.AddHydrogens()

    smiles=(obConversion.WriteString(mol))
    smiles=smiles.split()[0]
    number_of_atoms=mol.NumAtoms()
    formula=mol.GetFormula()
    inchi=Chem.MolToInchi(mol)
    return smiles, number_of_atoms, formula, inchi
    

def clean_directory(project_path):
    for filename in os.listdir(f"{project_path}/orca/inp/"):
        if not filename.endswith(".inp"):
            os.remove(f"{project_path}/orca/inp/{filename}")

#%%main
def main():
    #Define parameters here
    project_path = os.getcwd()

    process_out_files(project_path)
    clean_directory(project_path)

if __name__ == "__main__":
    main()

