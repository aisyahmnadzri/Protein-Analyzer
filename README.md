# 🔬 Protein Data Retrieval, Interaction Network & 3D Structure

Play around with the app [here](https://protein-analyzer-aisyahmnadzri.streamlit.app/)

This Streamlit app lets you:
- Retrieve protein characteristics from **UniProt**.
- Visualize real **protein-protein interaction networks** from **STRING DB**.
- Explore **3D structures** from **PDB** or **AlphaFold DB** interactively.

---

## 🚀 Features
- Enter a UniProt ID (e.g., `P69905` for Hemoglobin alpha).
- See protein metadata (function, length, weight, pathways, etc.).
- Generate a network graph of protein-protein interactions.
- View 3D structures in an interactive viewer.

---

## 📂 Project Structure
```
protein-app/
│
├── app.py            # Main Streamlit app
├── requirements.txt  # Dependencies
└── README.md         # Documentation
```

---

## ⚙️ Installation
Clone the repository and install dependencies:

```
git clone https://github.com/aisyahmnadzri/Protein-Analyzer.git
cd protein-app
pip install -r requirements.txt
```

## ▶️ Run Locally
Start the app with:

```
streamlit run app.py
```

## ☁️ Deployment
You can deploy easily on Streamlit Cloud:
- Push this repo to GitHub.
- Go to Streamlit Cloud.
- Create a new app → select your repo → choose app.py as the entry point.

## 📋 Requirements
Dependencies are listed in requirements.txt:
```
streamlit
requests
matplotlib
networkx
py3Dmol
```
## 🧪 Example Usage
- Input: P69905 (Hemoglobin subunit alpha).
Output:
- Protein characteristics (function, weight, etc.).
- STRING DB network of interactions.
- 3D structure from PDB/AlphaFold.
