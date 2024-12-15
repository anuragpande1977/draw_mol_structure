import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os

def smiles_to_structure(smiles):
    """
    Converts SMILES to 3D structure and returns PNG image and MOL2 data.
    
    Parameters:
        smiles (str): Input SMILES string of the molecule.
        
    Returns:
        img (BytesIO): PNG image of the molecule.
        mol2_block (str): MOL2 format representation of the molecule.
    """
    try:
        # Create RDKit Mol object from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError("Invalid SMILES string. Could not parse molecule.")
        
        # Add hydrogens and generate 3D conformation
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        # Create PNG image
        img = Draw.MolToImage(mol, size=(300, 300))
        
        # Convert to MOL2 format
        mol2_block = Chem.MolToMol2Block(mol)
        return img, mol2_block

    except Exception as e:
        st.error(f"Error: {e}")
        return None, None


# Streamlit app
def main():
    st.title("SMILES to 3D Structure Converter")
    st.write("Input a SMILES string to generate a 3D structure and download its PNG and MOL2 files.")

    # User input for SMILES
    smiles = st.text_input("Enter SMILES string", placeholder="e.g., C1=CC=CC=C1")

    if smiles:
        # Generate outputs
        img, mol2_block = smiles_to_structure(smiles)
        
        if img and mol2_block:
            # Display the 3D structure image
            st.image(img, caption="3D Structure", use_column_width=False)

            # Offer downloads for PNG and MOL2 files
            st.download_button(
                label="Download PNG",
                data=img.tobytes(),
                file_name="structure.png",
                mime="image/png"
            )
            st.download_button(
                label="Download MOL2",
                data=mol2_block,
                file_name="structure.mol2",
                mime="text/plain"
            )

if __name__ == "__main__":
    main()
