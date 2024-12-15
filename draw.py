import requests
import streamlit as st

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

def get_mol2(cid):
    """
    Fetch the MOL2 file for a given CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/MOL2"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    return None

def get_chebi_mol2(chemical_name):
    """
    Fetch the MOL2 file for a given chemical name using ChEBI API.
    """
    # Step 1: Search for the compound by name
    search_url = f"https://www.ebi.ac.uk/chebi/webServices.do?operation=searchByName&name={chemical_name}"
    response = requests.get(search_url)
    if response.status_code == 200 and "<chebiId>" in response.text:
        chebi_id = response.text.split("<chebiId>")[1].split("</chebiId>")[0]
        # Step 2: Get MOL2 file by ChEBI ID
        structure_url = f"https://www.ebi.ac.uk/chebi/webServices.do?operation=getStructureById&chebiId={chebi_id}"
        mol2_response = requests.get(structure_url)
        if mol2_response.status_code == 200:
            return mol2_response.text
    return None

# Streamlit App
def main():
    st.title("Chemical Structure Finder (PubChem + ChEBI)")
    st.write("Enter the name of a chemical to fetch its SMILES and MOL2 file.")

    # Input: Chemical Name
    chemical_name = st.text_input("Enter Chemical Name", placeholder="e.g., Aspirin")

    if chemical_name:
        st.write(f"Searching for: **{chemical_name}**")
        
        # Step 1: Get CID from PubChem
        cid = get_pubchem_cid(chemical_name)
        mol2_data = None

        if cid:
            st.success(f"Chemical found! PubChem CID: **{cid}**")
            
            # Step 2: Get SMILES from PubChem
            smiles = get_smiles(cid)
            if smiles:
                st.write(f"**SMILES**: `{smiles}`")
                st.download_button("Download SMILES", smiles, file_name=f"{chemical_name}_smiles.txt", mime="text/plain")

            # Step 3: Try to get MOL2 from PubChem
            mol2_data = get_mol2(cid)
            if mol2_data:
                st.download_button("Download MOL2 File (PubChem)", mol2_data, file_name=f"{chemical_name}.mol2", mime="chemical/x-mol2")
            else:
                st.warning("MOL2 file not available on PubChem.")
        
        # Step 4: Fallback to ChEBI API if MOL2 is not found on PubChem
        if not mol2_data:
            st.info("Searching ChEBI database for MOL2 file...")
            mol2_data = get_chebi_mol2(chemical_name)
            if mol2_data:
                st.download_button("Download MOL2 File (ChEBI)", mol2_data, file_name=f"{chemical_name}_chebi.mol2", mime="chemical/x-mol2")
            else:
                st.error("MOL2 file not available on both PubChem and ChEBI.")

if __name__ == "__main__":
    main()
