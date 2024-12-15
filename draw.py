import streamlit as st
from PIL import Image, ImageDraw

def generate_placeholder_image(smiles):
    """
    Generates a placeholder image for the given SMILES string.
    """
    img = Image.new("RGB", (400, 300), color="white")
    draw = ImageDraw.Draw(img)
    draw.text((10, 150), f"SMILES: {smiles}", fill="black")
    return img

# Streamlit app
def main():
    st.title("SMILES Structure Viewer")
    st.write("Enter a SMILES string to generate a placeholder PNG image.")

    # Input for SMILES
    smiles = st.text_input("Enter SMILES string", placeholder="e.g., C1=CC=CC=C1 (benzene)")

    if smiles:
        # Generate placeholder image
        img = generate_placeholder_image(smiles)

        # Display image
        st.image(img, caption="Generated Structure", use_column_width=False)

        # Provide download button
        img_bytes = BytesIO()
        img.save(img_bytes, format="PNG")
        img_bytes.seek(0)

        st.download_button(
            label="Download PNG",
            data=img_bytes,
            file_name="molecule.png",
            mime="image/png"
        )

if __name__ == "__main__":
    main()
