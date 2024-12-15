import requests
from rdkit import Chem
from rdkit.Chem import AllChem

def get_pubchem_cid(chemical_name):
    """
    Fetch the PubChem CID for a given chemical name.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/cids/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if "IdentifierList" in data and "CID" in data["IdentifierList"]:
            return data["IdentifierList"]["CID"][0]
    return None

def get_smiles(cid):
    """
    Fetch the SMILES string for a given CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    return None

def smiles_to_mol2(smiles):
    """
    Converts a SMILES string to MOL2 format.
    """
    try:
        # Generate molecule object from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError("Invalid SMILES string.")
        
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        # Export to MOL2 format
        mol2_data = Chem.MolToMol2Block(mol)
        return mol2_data
    except Exception as e:
        print(f"Error converting SMILES to MOL2: {e}")
        return None

# Example Usage
if __name__ == "__main__":
    chemical_name = input("Enter the name of the chemical: ")
    
    # Step 1: Fetch CID
    cid = get_pubchem_cid(chemical_name)
    if not cid:
        print(f"Chemical '{chemical_name}' not found in PubChem.")
    else:
        print(f"PubChem CID: {cid}")
        
        # Step 2: Fetch SMILES
        smiles = get_smiles(cid)
        if not smiles:
            print(f"SMILES not found for CID {cid}.")
        else:
            print(f"SMILES: {smiles}")
            
            # Step 3: Convert SMILES to MOL2
            mol2_data = smiles_to_mol2(smiles)
            if mol2_data:
                with open(f"{chemical_name}.mol2", "w") as file:
                    file.write(mol2_data)
                print(f"MOL2 file created: {chemical_name}.mol2")
