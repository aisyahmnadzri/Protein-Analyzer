import streamlit as st
import requests
import matplotlib.pyplot as plt
import networkx as nx
import py3Dmol

# -------------------------------
# Function to retrieve protein data from UniProt
# -------------------------------
def get_protein_data(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = requests.get(url)

    if response.status_code != 200:
        st.error("Invalid UniProt ID or unable to fetch data.")
        return {}

    protein_data = {"ID": uniprot_id}
    lines = response.text.splitlines()

    for i, line in enumerate(lines):
        if line.startswith("ID"):
            fields = line.split()
            protein_data["Name"] = " ".join(fields[1:]) if len(fields) >= 2 else "Not available"

        elif line.startswith("SQ") and i + 1 < len(lines):
            weight_line = lines[i + 1]
            weight = weight_line.split()[-1]
            protein_data["Weight"] = weight

        elif line.startswith("DE   RecName: Full="):
            function = line.split("Full=")[1].split(";")[0]
            protein_data["Function"] = function

        elif line.startswith("DR   SUPFAM"):
            structure = line.split(";")[1]
            protein_data["Structure"] = structure

        elif line.startswith("FT   MOD_RES"):
            parts = line.split()
            if len(parts) >= 3:
                protein_data["Length"] = int(parts[2])
            ptms = line.split(";")[1:]
            if ptms:
                protein_data["PTMs"] = [ptm.strip() for ptm in ptms]

        elif line.startswith("CC   -!- SUBCELLULAR LOCATION:"):
            subcellular_location = line.split("CC   -!- SUBCELLULAR LOCATION:")[1].strip()
            protein_data["Subcellular Location"] = subcellular_location

        elif line.startswith("DR   Reactome"):
            pathway = line.split(";")[1]
            protein_data["Pathway"] = pathway

        elif line.startswith("DR   MIM"):
            disease = line.split(";")[1]
            protein_data["Disease"] = disease

        elif line.startswith("DR   PDB"):
            pdb_id = line.split(";")[1].strip()
            protein_data["PDB_ID"] = pdb_id

        elif line.startswith("DR   AlphaFoldDB"):
            af_id = line.split(";")[1].strip()
            protein_data["AlphaFold_ID"] = af_id

    return protein_data

# -------------------------------
# Function to display protein characteristics
# -------------------------------
def display_protein_characteristics(protein_data):
    st.subheader("Protein Characteristics")
    st.write(f"**UniProt ID:** {protein_data.get('ID', 'N/A')}")
    st.write(f"**Name:** {protein_data.get('Name', 'N/A')}")
    st.write(f"**Length:** {protein_data.get('Length', 'Not available')}")
    st.write(f"**Molecular Weight:** {protein_data.get('Weight', 'N/A')}")
    st.write(f"**Function:** {protein_data.get('Function', 'N/A')}")

    if "Structure" in protein_data:
        st.write(f"**Structure:** {protein_data['Structure']}")
    if "Subcellular Location" in protein_data:
        st.write(f"**Subcellular Location:** {protein_data['Subcellular Location']}")
    if "PTMs" in protein_data:
        st.write("**Post-Translational Modifications (PTMs):**")
        for ptm in protein_data["PTMs"]:
            st.write(f"- {ptm}")
    if "Pathway" in protein_data:
        st.write(f"**Pathway Involvement:** {protein_data['Pathway']}")
    if "Disease" in protein_data:
        st.write(f"**Disease Association:** {protein_data['Disease']}")

# -------------------------------
# Real protein-protein interaction network from STRING DB
# -------------------------------
def get_ppi_network(uniprot_id, species=9606):
    url = f"https://string-db.org/api/json/network?identifiers={uniprot_id}&species={species}"
    response = requests.get(url)

    if response.status_code != 200:
        st.error("Failed to fetch PPI data from STRING DB.")
        return []

    data = response.json()
    interactions = []
    for entry in data:
        protein1 = entry.get("preferredName_A")
        protein2 = entry.get("preferredName_B")
        if protein1 and protein2:
            interactions.append((protein1, protein2))
    return interactions

def display_ppi_network(interactions):
    st.subheader("Protein-Protein Interaction Network")

    if not interactions:
        st.warning("No interactions found for this protein in STRING DB.")
        return

    # Draw the graph
    fig, ax = plt.subplots()
    G = nx.Graph()
    G.add_edges_from(interactions)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color="skyblue", node_size=1500,
            font_size=10, font_color="black", edge_color="black", linewidths=1)
    plt.title("Protein-Protein Interaction Network")
    st.pyplot(fig)

    # Dynamic explanation under the graph
    st.markdown(
        f"""
        **Explanation:**  
        This network contains **{len(G.nodes)} proteins** and **{len(G.edges)} interactions** retrieved from STRING DB.  
        - Each node represents a protein.  
        - Each edge indicates a known or predicted interaction, based on STRING’s evidence sources (experiments, curated databases, text mining, co-expression).  
        - The spring layout clusters proteins with more connections closer together.  

        Exploring these interactions helps researchers uncover functional relationships, biological pathways, 
        and potential therapeutic targets for the queried protein.
        """
    )

# -------------------------------
# Function to display 3D structure (PDB or AlphaFold)
# -------------------------------
def display_structure(protein_data):
    st.subheader("3D Structure Viewer")

    pdb_id = protein_data.get("PDB_ID")
    af_id = protein_data.get("AlphaFold_ID")

    if pdb_id:
        st.write(f"Showing PDB structure: **{pdb_id}**")
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        pdb_data = requests.get(pdb_url).text
    elif af_id:
        st.write(f"Showing AlphaFold predicted structure for UniProt ID: **{af_id}**")
        pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{af_id}-F1-model_v4.pdb"
        pdb_data = requests.get(pdb_url).text
    else:
        st.warning("No PDB or AlphaFold structure available for this protein.")
        return

    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(pdb_data, "pdb")
    viewer.setStyle({"cartoon": {"color": "spectrum"}})
    viewer.zoomTo()
    st.components.v1.html(viewer._make_html(), height=500)

# -------------------------------
# Main Streamlit App
# -------------------------------
def main():
    st.title("🔬 Protein Data Retrieval, Interaction Network & 3D Structure")

    uniprot_id = st.text_input("Enter UniProt ID (e.g., P69905 for Hemoglobin alpha):")

    if st.button("Retrieve Data"):
        if uniprot_id:
            protein_data = get_protein_data(uniprot_id)
            if protein_data:
                display_protein_characteristics(protein_data)
                ppi_network = get_ppi_network(uniprot_id)  # Real STRING DB data
                display_ppi_network(ppi_network)
                display_structure(protein_data)

if __name__ == "__main__":
    main()
