About WASP
About the WASP Code
WASP (Web-based Approach for Subtractive Proteomics) is a comprehensive bioinformatics tool designed to aid researchers in identifying potential drug targets by comparing the proteome of a pathogen with its host. The application is built using Python and the Streamlit framework, providing an interactive and user-friendly interface for conducting subtractive proteomics analyses.

Code Structure and Features
Main Components
Home Page: The main interface where users can:

Select a KEGG ID for the pathogen.
Upload proteome data for both the pathogen and the host.
Run the complete analysis pipeline, including BLASTp analysis, prioritization of essential proteins, calculation of protein properties, and fetching KEGG pathways.
About Page: Provides an overview of the WASP tool, its features, and functionality.

Help Page: Offers guidance on how to use the WASP tool, including step-by-step instructions and answers to frequently asked questions.

Functions
load_proteome_data:

Loads and parses proteome data from uploaded FASTA files.
Handles file content and parses sequences using BioPython.
run_blastp:

Executes BLASTp analysis by comparing pathogen proteins against host proteins.
Writes query and subject proteins to temporary FASTA files and runs BLASTp using the NcbiblastpCommandline.
prioritize_essential_proteins:

Prioritizes essential proteins based on BLASTp results.
Identifies non-homologous proteins and filters essential proteins accordingly.
calculate_properties:

Calculates various properties of proteins, including molecular weight, theoretical pI, GRAVY score, and instability index.
Uses BioPython's ProtParam module to perform the calculations.
get_kegg_id:

Retrieves the KEGG ID for a given UniProt ID by querying the UniProt database.
fetch_kegg_pathways:

Fetches KEGG pathways associated with a given KEGG ID by querying the KEGG database.
Analysis Pipeline
Select Pathogen KEGG ID:

Users can select a KEGG ID from a predefined list of popular organisms or enter a custom KEGG ID.
Upload Proteome Data:

Users upload FASTA files for both the pathogen and host proteomes.
Run BLASTp Analysis:

Compares pathogen proteins against host proteins to identify non-homologous proteins.
Prioritize Essential Proteins:

Identifies and prioritizes essential proteins based on BLASTp results.
Calculate Protein Properties:

Calculates various physicochemical properties of the prioritized essential proteins.
Fetch KEGG Pathways:

Retrieves KEGG pathway information for the prioritized proteins.
Save and Download Results:

Saves analysis results to a folder and provides an option to download them as a ZIP file.
Usage
Select Pathogen KEGG ID: Choose from a dropdown menu or enter a custom KEGG ID.
Upload Proteome Data: Upload FASTA files for the pathogen and host proteomes.
Run Analysis: Click the 'Run Analysis' button to start the analysis pipeline.
Download Results: After the analysis completes, download the results as a ZIP file.
