import streamlit as st
from Bio import SeqIO
from io import StringIO
import requests
from Bio.Blast import NCBIXML
import csv
from Bio.SeqUtils import ProtParam

# Function to load proteome data from uploaded FASTA files
def load_proteome_data(proteome_file):
    try:
        proteome_content = proteome_file.read().decode("utf-8")
        proteome_stream = StringIO(proteome_content)
        proteome_sequences = list(SeqIO.parse(proteome_stream, "fasta"))
        return proteome_sequences
    except Exception as e:
        st.error(f"Error loading proteome data from {proteome_file.name}: {e}")
        return None

# Function to run BLASTp analysis using NCBI BLAST API
def run_blastp_online(query_proteins, subject_proteins):
    try:
        query_sequences = "\n".join([f">{seq.id}\n{seq.seq}" for seq in query_proteins])
        subject_sequences = "\n".join([f">{seq.id}\n{seq.seq}" for seq in subject_proteins])

        # NCBI BLAST API URL
        blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

        # Send query sequences to NCBI BLAST API
        response_query = requests.post(blast_url, data={
            'CMD': 'Put',
            'PROGRAM': 'blastp',
            'DATABASE': 'nr',
            'QUERY': query_sequences
        })

        if response_query.status_code != 200:
            st.error(f"Error submitting query sequences to BLAST API: {response_query.status_code}")
            return None, None

        rid_query = None
        for line in response_query.text.split("\n"):
            if line.startswith("RID"):
                rid_query = line.split("=")[1].strip()
                break

        if not rid_query:
            st.error("Error retrieving RID for query sequences")
            return None, None

        # Send subject sequences to NCBI BLAST API
        response_subject = requests.post(blast_url, data={
            'CMD': 'Put',
            'PROGRAM': 'blastp',
            'DATABASE': 'nr',
            'QUERY': subject_sequences
        })

        if response_subject.status_code != 200:
            st.error(f"Error submitting subject sequences to BLAST API: {response_subject.status_code}")
            return None, None

        rid_subject = None
        for line in response_subject.text.split("\n"):
            if line.startswith("RID"):
                rid_subject = line.split("=")[1].strip()
                break

        if not rid_subject:
            st.error("Error retrieving RID for subject sequences")
            return None, None

        # Wait for the BLAST results
        st.info("Waiting for BLAST results. This may take a few minutes.")

        while True:
            result_query = requests.get(blast_url, params={'CMD': 'Get', 'RID': rid_query, 'FORMAT_OBJECT': 'SearchInfo'})
            if 'Status=READY' in result_query.text:
                break
            st.info("Waiting for query BLAST results...")

        while True:
            result_subject = requests.get(blast_url, params={'CMD': 'Get', 'RID': rid_subject, 'FORMAT_OBJECT': 'SearchInfo'})
            if 'Status=READY' in result_subject.text:
                break
            st.info("Waiting for subject BLAST results...")

        # Fetch the BLAST results
        result_query = requests.get(blast_url, params={'CMD': 'Get', 'RID': rid_query, 'FORMAT_TYPE': 'XML'})
        result_subject = requests.get(blast_url, params={'CMD': 'Get', 'RID': rid_subject, 'FORMAT_TYPE': 'XML'})

        return result_query.text, result_subject.text

    except Exception as e:
        st.error(f"Error running BLASTp: {e}")
        return None, None

# Function to prioritize essential proteins based on BLASTp results
def prioritize_essential_proteins(blastp_results_query, blastp_results_subject, pathogen_proteins):
    try:
        query_records = NCBIXML.read(StringIO(blastp_results_query))
        subject_records = NCBIXML.read(StringIO(blastp_results_subject))

        non_homologous_proteins = set()
        for query_record in query_records:
            if not query_record.alignments:
                non_homologous_proteins.add(query_record.query.split()[0])

        prioritized_proteins = [protein for protein in pathogen_proteins if protein.id in non_homologous_proteins]

        st.success(f"Prioritized {len(prioritized_proteins)} essential proteins based on BLASTp results.")
        return prioritized_proteins

    except Exception as e:
        st.error(f"Error prioritizing essential proteins: {e}")
        return []

# Function to calculate protein properties
def calculate_properties(sequence):
    valid_residues = set("ACDEFGHIKLMNPQRSTVWY")
    if not set(sequence).issubset(valid_residues):
        raise ValueError(f"Invalid sequence: {sequence}. Contains non-standard characters.")

    protein_analysis = ProtParam.ProteinAnalysis(sequence)

    mw = protein_analysis.molecular_weight()
    theoretical_pI = protein_analysis.isoelectric_point()
    gravy_score = protein_analysis.gravy()
    instability_index = protein_analysis.instability_index()

    return mw, theoretical_pI, gravy_score, instability_index

# Function to get KEGG ID from UniProt ID
def get_kegg_id(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = requests.get(url)
    if response.status_code == 200:
        for line in response.text.splitlines():
            if line.startswith("DR   KEGG;"):
                return line.split(";")[1].strip().split(":")[1]
    return None

# Function to fetch KEGG pathways for a given KEGG ID
def fetch_kegg_pathways(kegg_id):
    url = f"http://rest.kegg.jp/link/pathway/{kegg_id}"
    response = requests.get(url)
    pathways = []
    if response.status_code == 200 and response.text:
        for line in response.text.splitlines():
            pathways.append(line.split('\t')[1].strip())
    return pathways

# Streamlit app setup
def main():
    st.title("Proteomics Analysis App")
    st.write("---")

    st.header("Select Pathogen KEGG ID")
    popular_organisms = {
        "Escherichia coli (k12)": "eco",
        "Bacillus subtilis": "bsu",
        "Homo sapiens": "hsa",
        "Saccharomyces cerevisiae": "sce",
        "Arabidopsis thaliana": "ath",
        "Enterococcus faecium": "efa",
        "Staphylococcus aureus": "sau",
        "Klebsiella pneumoniae": "kpn",
        "Acinetobacter baumannii": "aba",
        "Pseudomonas aeruginosa": "pae"
    }
    custom_option = "Custom (Enter KEGG ID)"
    selected_organism = st.selectbox("Select Organism", options=list(popular_organisms.keys()) + [custom_option])

    if selected_organism == custom_option:
        kegg_id = st.text_input("Enter Custom KEGG ID:")
    else:
        kegg_id = popular_organisms.get(selected_organism)

    if not kegg_id:
        st.warning("Please select or enter a valid KEGG ID.")
        st.stop()

    st.write(f"Selected KEGG ID: {kegg_id}")
    st.write("---")

    st.header("Upload Proteome Data")
    pathogen_proteome_file = st.file_uploader("Upload Pathogen Proteome FASTA", type=["fasta"], key="pathogen")
    host_proteome_file = st.file_uploader("Upload Host Proteome FASTA", type=["fasta"], key="host")

    if st.button("Run Complete Analysis"):
        if pathogen_proteome_file and host_proteome_file:
            pathogen_proteins = load_proteome_data(pathogen_proteome_file)
            host_proteins = load_proteome_data(host_proteome_file)

            if pathogen_proteins and host_proteins:
                st.header("Running BLASTp Analysis")
                blastp_results_query, blastp_results_subject = run_blastp_online(pathogen_proteins, host_proteins)
                if blastp_results_query and blastp_results_subject:
                    st.success("BLASTp analysis completed successfully.")
                else:
                    st.error("BLASTp analysis failed.")
                    st.stop()

                st.header("Prioritizing Essential Proteins")
                prioritized_proteins = prioritize_essential_proteins(blastp_results_query, blastp_results_subject, pathogen_proteins)

                st.header("Calculating Protein Properties")
                protein_properties = {}
                for protein in prioritized_proteins:
                    sequence = str(protein.seq)
                    organism = protein.description.split("OS=")[1].split(" ")[0]
                    try:
                        mw, theoretical_pI, gravy, instability_index = calculate_properties(sequence)
                        kegg_id = get_kegg_id(protein.id)
                        kegg_pathways = fetch_kegg_pathways(kegg_id) if kegg_id else []
                        protein_properties[protein.id] = {
                            "Organism": organism,
                            "Molecular Weight": mw,
                            "Theoretical pI": theoretical_pI,
                            "GRAVY Score": gravy,
                            "Instability Index": instability_index,
                            "KEGG Pathways": ", ".join(kegg_pathways) if kegg_pathways else "None"
                        }
                    except Exception as e:
                        st.warning(f"Skipping protein {protein.id} due to error: {e}")

                st.header("Protein Properties")
                st.write(protein_properties)

                st.header("Download Results")
                csv_buffer = StringIO()
                csv_writer = csv.writer(csv_buffer)
                csv_writer.writerow(["Protein ID", "Organism", "Molecular Weight", "Theoretical pI", "GRAVY Score", "Instability Index", "KEGG Pathways"])
                for protein_id, properties in protein_properties.items():
                    csv_writer.writerow([
                        protein_id,
                        properties["Organism"],
                        properties["Molecular Weight"],
                        properties["Theoretical pI"],
                        properties["GRAVY Score"],
                        properties["Instability Index"],
                        properties["KEGG Pathways"]
                    ])

                csv_buffer.seek(0)
                st.download_button(
                    label="Download Protein Properties CSV",
                    data=csv_buffer,
                    file_name="protein_properties.csv",
                    mime="text/csv"
                )

                st.success("Analysis completed successfully.")
            else:
                st.error("Failed to load proteome data.")
        else:
            st.warning("Please upload both pathogen and host proteome FASTA files.")

if __name__ == "__main__":
    main()
