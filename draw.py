import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from io import BytesIO

def smiles_to_structure(smiles):
    """
    Converts a SMILES string to 3D structure and generates PNG and MOL2 outputs.
    
    Parameters:
        smiles (str): Input SMILES string.

    Returns:
        mol2_str (str): MOL2 format string.
        img_bytes (BytesIO): PNG image in bytes.
    """
    try:
        # Create RDKit Mol object
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError("Invalid SMILES string.")
        
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        # Create PNG image
        img = Draw.MolToImage(mol, size=(300, 300))
        img_bytes = BytesIO()
        img.save(img_bytes, format="PNG")
        img_bytes.seek(0)
        
        # Convert to MOL2 format
        mol2_str = Chem.MolToMol2Block(mol)
        
        return mol2_str, img_bytes

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
        # Process SMILES and generate outputs
        mol2_str, img_bytes = smiles_to_structure(smiles)

        if mol2_str and img_bytes:
            # Display PNG image
            st.image(img_bytes, caption="3D Structure", use_column_width=False)

            # Provide download options
            st.download_button(
                label="Download PNG",
                data=img_bytes,
                file_name="structure.png",
                mime="image/png"
            )
            st.download_button(
                label="Download MOL2",
                data=mol2_str,
                file_name="structure.mol2",
                mime="text/plain"
            )

if __name__ == "__main__":
    main()
