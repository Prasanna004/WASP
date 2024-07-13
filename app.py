import streamlit as st
from Bio import SeqIO
from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import os
import tempfile
import csv
from Bio.SeqUtils import ProtParam
import requests
import zipfile
import shutil

# Function to load proteome data from uploaded FASTA files
def load_proteome_data(proteome_file):
    try:
        # Read file content
        proteome_content = proteome_file.read().decode("utf-8")
        
        # Use StringIO to handle FASTA content
        proteome_stream = StringIO(proteome_content)
        
        # Parse sequences
        proteome_sequences = list(SeqIO.parse(proteome_stream, "fasta"))
        
        return proteome_sequences
    except Exception as e:
        st.error(f"Error loading proteome data from {proteome_file.name}: {e}")
        return None

# Function to run BLASTp analysis
def run_blastp(query_proteins, subject_proteins, output_file):
    try:
        # Write query and subject proteins to temporary FASTA files
        temp_dir = tempfile.TemporaryDirectory()
        query_file = os.path.join(temp_dir.name, "query_proteins.fasta")
        subject_file = os.path.join(temp_dir.name, "subject_proteins.fasta")

        SeqIO.write(query_proteins, query_file, "fasta")
        SeqIO.write(subject_proteins, subject_file, "fasta")

        # Run BLASTp
        blastp_cline = NcbiblastpCommandline(query=query_file, subject=subject_file, outfmt=5, out=output_file)
        stdout, stderr = blastp_cline()

        # Remove temporary directory
        temp_dir.cleanup()

        if stderr:
            st.error(stderr)
            return False
        else:
            st.success(f"BLASTp results saved to {output_file}")
            return True
    except Exception as e:
        st.error(f"Error running BLASTp: {e}")
        return False

# Function to prioritize essential proteins based on BLASTp results
def prioritize_essential_proteins(blastp_results_file, pathogen_proteins):
    try:
        # Parse non-homologous proteins from BLASTp results
        with open(blastp_results_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            non_homologous_proteins = set()
            for blast_record in blast_records:
                if not blast_record.alignments:
                    non_homologous_proteins.add(blast_record.query.split()[0])

        # Filter and prioritize essential proteins
        prioritized_proteins = [protein for protein in pathogen_proteins if protein.id in non_homologous_proteins]

        st.success(f"Prioritized {len(prioritized_proteins)} essential proteins based on BLASTp results.")
        return prioritized_proteins

    except Exception as e:
        st.error(f"Error prioritizing essential proteins: {e}")
        return []

# Function to calculate protein properties
def calculate_properties(sequence):
    # Check if the sequence contains only valid amino acid residues
    valid_residues = set("ACDEFGHIKLMNPQRSTVWY")
    if not set(sequence).issubset(valid_residues):
        raise ValueError(f"Invalid sequence: {sequence}. Contains non-standard characters.")

    # Create a ProteinAnalysis object
    protein_analysis = ProtParam.ProteinAnalysis(sequence)

    # Calculate properties
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

    # Step 0: Select KEGG ID
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

    # Step 1: Upload proteome data
    st.header("Upload Proteome Data")
    pathogen_proteome_file = st.file_uploader("Upload Pathogen Proteome FASTA", type=["fasta"], key="pathogen")
    host_proteome_file = st.file_uploader("Upload Host Proteome FASTA", type=["fasta"], key="host")

    # Run the complete analysis pipeline
    if st.button("Run Complete Analysis"):
        if pathogen_proteome_file and host_proteome_file:
            pathogen_proteins = load_proteome_data(pathogen_proteome_file)
            host_proteins = load_proteome_data(host_proteome_file)

            if pathogen_proteins and host_proteins:
                # Step 2: Run BLASTp analysis
                st.header("Running BLASTp Analysis")
                output_result_file = "blastp_results.xml"
                if run_blastp(pathogen_proteins, host_proteins, output_result_file):
                    st.success("BLASTp analysis completed successfully.")
                else:
                    st.error("BLASTp analysis failed.")
                    st.stop()

                # Step 3: Prioritize essential proteins based on BLASTp results
                st.header("Prioritizing Essential Proteins")
                prioritized_proteins = prioritize_essential_proteins(output_result_file, pathogen_proteins)

                # Step 4: Calculate protein properties and fetch KEGG pathways
                st.header("Calculating Protein Properties and Fetching KEGG Pathways")
                protein_properties = {}
                for protein in prioritized_proteins:
                    sequence = str(protein.seq)
                    organism = protein.description.split("OS=")[1].split(" ")[0]  # Extract organism name
                    try:
                        mw, theoretical_pI, gravy, instability_index = calculate_properties(sequence)
                        kegg_id = get_kegg_id(protein.id)
                        pathways = fetch_kegg_pathways(kegg_id) if kegg_id else ["No KEGG Pathway found"]
                        protein_properties[protein.id] = {
                            "Organism": organism,
                            "Molecular Weight": mw,
                            "Theoretical pI": theoretical_pI,
                            "GRAVY Score": gravy,
                            "Instability Index": instability_index,
                            "KEGG Pathways": "; ".join(pathways)
                        }
                    except ValueError as e:
                        st.error(f"Error calculating properties for {protein.id}: {str(e)}")

                # Step 5: Save results to a folder
                st.header("Saving Results")
                result_folder = "proteomics_results"
                os.makedirs(result_folder, exist_ok=True)

                # Save BLASTp results
                shutil.copy(output_result_file, os.path.join(result_folder, "blastp_results.xml"))

                # Save prioritized essential proteins
                prioritized_file = os.path.join(result_folder, "prioritized_essential_proteins.fasta")
                SeqIO.write(prioritized_proteins, prioritized_file, "fasta")

                # Save protein properties with KEGG pathways
                properties_file = os.path.join(result_folder, "protein_properties.csv")
                with open(properties_file, mode="w", newline="") as csvfile:
                    fieldnames = ["Protein ID", "Organism", "Molecular Weight", "Theoretical pI", "GRAVY Score", "Instability Index", "KEGG Pathways"]
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    for protein_id, properties in protein_properties.items():
                        row = {"Protein ID": protein_id}
                        row.update(properties)
                        writer.writerow(row)

                # Create a zip file of the result folder
                result_zip = "proteomics_results.zip"
                with zipfile.ZipFile(result_zip, 'w') as zipf:
                    for root, _, files in os.walk(result_folder):
                        for file in files:
                            zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), result_folder))

                # Provide download link
                with open(result_zip, "rb") as zip_file:
                    st.download_button(label="Download Results", data=zip_file, file_name=result_zip, mime="application/zip")

                st.success("Results saved and ready for download.")
            else:
                st.error("Error loading proteome data.")
        else:
            st.warning("Please upload both pathogen and host proteome FASTA files.")

if __name__ == "__main__":
    main()
