import streamlit as st
from openbabel import pybel
from PIL import Image, ImageDraw
import os

def smiles_to_mol2_and_png(smiles):
    """
    Converts a SMILES string to MOL2 and PNG formats.
    
    Parameters:
        smiles (str): Input SMILES string of the molecule.

    Returns:
        mol2_data (str): MOL2 format string of the molecule.
        img (Image): PNG image of the molecule.
    """
    try:
        # Create molecule from SMILES
        mol = pybel.readstring("smi", smiles)
        mol.addh()  # Add hydrogens
        mol.make3D()  # Generate 3D structure

        # Convert to MOL2
        mol2_data = mol.write("mol2")

        # Generate PNG image
        mol.draw(show=False, filename="temp.png")
        img = Image.open("temp.png")
        os.remove("temp.png")  # Clean up temporary file
        return mol2_data, img
    except Exception as e:
        st.error(f"Error: {e}")
        return None, None

# Streamlit app
def main():
    st.title("SMILES to 3D Structure Converter (Open Babel)")
    st.write("Input a SMILES string to generate a 3D structure and download its PNG and MOL2 files.")

    # User input for SMILES
    smiles = st.text_input("Enter SMILES string", placeholder="e.g., C1=CC=CC=C1")

    if smiles:
        mol2_data, img = smiles_to_mol2_and_png(smiles)

        if mol2_data and img:
            # Display molecule image
            st.image(img, caption="3D Structure", use_column_width=False)

            # Download buttons for MOL2 and PNG
            st.download_button(
                label="Download MOL2",
                data=mol2_data,
                file_name="structure.mol2",
                mime="text/plain"
            )
            st.download_button(
                label="Download PNG",
                data=img.tobytes(),
                file_name="structure.png",
                mime="image/png"
            )

if __name__ == "__main__":
    main()
