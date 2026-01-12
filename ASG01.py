import streamlit as st
import requests
import matplotlib.pyplot as plt
import networkx as nx

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
# Dummy protein-protein interaction network
# -------------------------------
def get_ppi_network(uniprot_id):
    proteins = ["ProteinA", "ProteinB", "ProteinC", "ProteinD", "ProteinE",
                "ProteinF", "ProteinG", "ProteinH", "ProteinI", "ProteinJ"]
    interactions = [(proteins[i], proteins[j]) for i in range(len(proteins)) for j in range(i + 1, len(proteins))]
    return interactions

def display_ppi_network(interactions):
    st.subheader("Protein-Protein Interaction Network")
    fig, ax = plt.subplots()
    G = nx.Graph()
    G.add_edges_from(interactions)
    pos = nx.spring_layout(G)
    nx.draw(
        G, pos, with_labels=True,
        node_color="skyblue", node_size=1500,
        font_size=10, font_color="black",
        edge_color="black", linewidths=1
    )
    plt.title("Protein-Protein Interaction Network")
    st.pyplot(fig)

# -------------------------------
# Main Streamlit App
# -------------------------------
def main():
    st.title("🔬 Protein Data Retrieval and Analysis")

    uniprot_id = st.text_input("Enter UniProt ID (e.g., P69905 for Hemoglobin alpha):")

    if st.button("Retrieve Data"):
        if uniprot_id:
            protein_data = get_protein_data(uniprot_id)
            if protein_data:
                display_protein_characteristics(protein_data)
                ppi_network = get_ppi_network(uniprot_id)  # Dummy implementation
                display_ppi_network(ppi_network)

if __name__ == "__main__":
    main()
