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

# Streamlit App
def main():
    st.title("PubChem Chemical Structure Finder")
    st.write("Enter the name of a chemical to fetch its SMILES and MOL2 file from PubChem.")

    # Input: Chemical Name
    chemical_name = st.text_input("Enter Chemical Name", placeholder="e.g., Aspirin")

    if chemical_name:
        st.write(f"Searching for: **{chemical_name}**")
        
        # Step 1: Get CID
        cid = get_pubchem_cid(chemical_name)
        if cid:
            st.success(f"Chemical found! PubChem CID: **{cid}**")
            
            # Step 2: Get SMILES
            smiles = get_smiles(cid)
            if smiles:
                st.write(f"**SMILES**: `{smiles}`")
                st.download_button("Download SMILES", smiles, file_name=f"{chemical_name}_smiles.txt", mime="text/plain")

            # Step 3: Get MOL2
            mol2_data = get_mol2(cid)
            if mol2_data:
                st.download_button("Download MOL2 File", mol2_data, file_name=f"{chemical_name}.mol2", mime="chemical/x-mol2")
            else:
                st.warning("MOL2 format not available for this compound.")
        else:
            st.error("Chemical not found in PubChem database.")

if __name__ == "__main__":
    main()
