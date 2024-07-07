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
