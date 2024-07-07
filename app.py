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

# Function definitions (paste your functions here)

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

# Home page function
def home():
    st.title("WASP: Web-based Approach for Subtractive Proteomics")
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
    if st.button("Run Analysis"):
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

                # Step 4: Calculate protein properties
                st.header("Calculating Protein Properties")
                protein_properties = {}
                for protein in prioritized_proteins:
                    sequence = str(protein.seq)
                    organism = protein.description.split("OS=")[1].split(" ")[0]  # Extract organism name
                    try:
                        mw, theoretical_pI, gravy, instability_index = calculate_properties(sequence)
                        protein_properties[protein.id] = {
                            "Organism": organism,
                            "Molecular Weight": mw,
                            "Theoretical pI": theoretical_pI,
                            "GRAVY Score": gravy,
                            "Instability Index": instability_index
                        }
                    except ValueError as e:
                        st.error(f"Error calculating properties for {protein.id}: {str(e)}")

                # Step 5: Fetch KEGG pathways
                st.header("Fetching KEGG Pathways")
                pathways = {}
                for protein_id, _ in protein_properties.items():
                    kegg_id = get_kegg_id(protein_id)
                    if kegg_id:
                        pathways[protein_id] = fetch_kegg_pathways(kegg_id)
                    else:
                        pathways[protein_id] = ["Unclassified"]

                # Step 6: Save results to a folder
                st.header("Saving Results")
                result_folder = "proteomics_results"
                os.makedirs(result_folder, exist_ok=True)

                # Save BLASTp results
                blastp_results_path = os.path.join(result_folder, output_result_file)
                shutil.move(output_result_file, blastp_results_path)
                st.write(f"BLASTp results saved to {blastp_results_path}")

                # Save prioritized essential proteins to FASTA
                prioritized_proteins_file = os.path.join(result_folder, "prioritized_essential_proteins.fasta")
                SeqIO.write(prioritized_proteins, prioritized_proteins_file, "fasta")
                st.write(f"Prioritized essential proteins saved to {prioritized_proteins_file}")

                # Save protein properties to CSV
                protein_properties_file = os.path.join(result_folder, "protein_properties.csv")
                with open(protein_properties_file, "w", newline="") as csvfile:
                    fieldnames = ["Protein ID", "Organism", "Molecular Weight", "Theoretical pI", "GRAVY Score", "Instability Index"]
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    for protein_id, properties in protein_properties.items():
                        row = {"Protein ID": protein_id}
                        row.update(properties)
                        writer.writerow(row)
                st.write(f"Protein properties saved to {protein_properties_file}")

                # Save pathways to CSV
                pathways_file = os.path.join(result_folder, "protein_pathways.csv")
                with open(pathways_file, "w", newline="") as csvfile:
                    fieldnames = ["Protein ID", "KEGG Pathways"]
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    for protein_id, kegg_pathways in pathways.items():
                        row = {"Protein ID": protein_id, "KEGG Pathways": "; ".join(kegg_pathways)}
                        writer.writerow(row)
                st.write(f"KEGG pathways saved to {pathways_file}")

                # Create a ZIP archive of the results
                zip_filename = "analysis_results.zip"
                with zipfile.ZipFile(zip_filename, "w") as zipf:
                    for foldername, subfolders, filenames in os.walk(result_folder):
                        for filename in filenames:
                            file_path = os.path.join(foldername, filename)
                            zipf.write(file_path, os.path.relpath(file_path, result_folder))

                st.success(f"Results saved to {zip_filename}")

                # Step 7: Provide download link for the results
                st.header("Download Results")
                with open(zip_filename, "rb") as f:
                    st.download_button("Download ZIP", f, file_name=zip_filename)
        else:
            st.warning("Please upload both pathogen and host proteome FASTA files.")

# About page function
def about():
    st.title("About WASP")
    st.write(
        """
        # About the WASP Code

        WASP (Web-based Approach for Subtractive Proteomics) is a comprehensive bioinformatics tool designed to aid researchers in identifying potential drug targets by comparing the proteome of a pathogen with its host. The application is built using Python and the Streamlit framework, providing an interactive and user-friendly interface for conducting subtractive proteomics analyses.

        ## Code Structure and Features

        ### Main Components

        1. **Home Page**: The main interface where users can:
           - Select a KEGG ID for the pathogen.
           - Upload proteome data for both the pathogen and the host.
           - Run the complete analysis pipeline, including BLASTp analysis, prioritization of essential proteins, calculation of protein properties, and fetching KEGG pathways.

        2. **About Page**: Provides an overview of the WASP tool, its features, and functionality.

        3. **Help Page**: Offers guidance on how to use the WASP tool, including step-by-step instructions and answers to frequently asked questions.

        ### Functions

        1. **load_proteome_data**:
           - Loads and parses proteome data from uploaded FASTA files.
           - Handles file content and parses sequences using BioPython.

        2. **run_blastp**:
           - Executes BLASTp analysis by comparing pathogen proteins against host proteins.
           - Writes query and subject proteins to temporary FASTA files and runs BLASTp using the `NcbiblastpCommandline`.

        3. **prioritize_essential_proteins**:
           - Prioritizes essential proteins based on BLASTp results.
           - Identifies non-homologous proteins and filters essential proteins accordingly.

        4. **calculate_properties**:
           - Calculates various properties of proteins, including molecular weight, theoretical pI, GRAVY score, and instability index.
           - Uses BioPython's `ProtParam` module to perform the calculations.

        5. **get_kegg_id**:
           - Retrieves the KEGG ID for a given UniProt ID by querying the UniProt database.

        6. **fetch_kegg_pathways**:
           - Fetches KEGG pathways associated with a given KEGG ID by querying the KEGG database.

        ### Analysis Pipeline

        1. **Select Pathogen KEGG ID**:
           - Users can select a KEGG ID from a predefined list of popular organisms or enter a custom KEGG ID.

        2. **Upload Proteome Data**:
           - Users upload FASTA files for both the pathogen and host proteomes.

        3. **Run BLASTp Analysis**:
           - Compares pathogen proteins against host proteins to identify non-homologous proteins.

        4. **Prioritize Essential Proteins**:
           - Identifies and prioritizes essential proteins based on BLASTp results.

        5. **Calculate Protein Properties**:
           - Calculates various physicochemical properties of the prioritized essential proteins.

        6. **Fetch KEGG Pathways**:
           - Retrieves KEGG pathway information for the prioritized proteins.

        7. **Save and Download Results**:
           - Saves analysis results to a folder and provides an option to download them as a ZIP file.

        ### Usage

        1. **Select Pathogen KEGG ID**: Choose from a dropdown menu or enter a custom KEGG ID.
        2. **Upload Proteome Data**: Upload FASTA files for the pathogen and host proteomes.
        3. **Run Analysis**: Click the 'Run Analysis' button to start the analysis pipeline.
        4. **Download Results**: After the analysis completes, download the results as a ZIP file.

        For more information, please visit the [GitHub repository](https://github.com/username/WASP).
        """
    )

# Help page function
def help_page():
    st.title("Help")
    st.write("Here you can provide help content and FAQs for your users.")

# Main app
def main():
    st.sidebar.title("Navigation")
    pages = {
        "Home": home,
        "About": about,
        "Help": help_page
    }
    selection = st.sidebar.radio("Go to", list(pages.keys()))
    page = pages[selection]
    page()

if __name__ == "__main__":
    main()
