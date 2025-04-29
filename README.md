# stringdb_api

# Overview
STRING-DB Physical Protein Interaction Tool

This tool queries the STRING database for physical protein-protein interactions, with
special emphasis on experimental and database evidence. It produces multiple output files
including detailed interaction data, a network visualization, and a list of unique interactors.

# INSTALLATION:
```
#Clone this repository
git clone https://github.com/feltus/stringdb_api.git

#Install dependencies
#Note: You may want to create a conda or virtual environment (venv) for this code.
pip install requests pandas matplotlib networkx numpy
```

# USAGE:
```
python string_ppi.py --input_list GENE_LIST_FILE [OPTIONS]

REQUIRED ARGUMENTS:
  --input_list FILE       Path to comma-delimited file containing gene symbols

OPTIONS:
  --species INT           NCBI taxonomy ID (default: 9606 for human)
  --score_threshold FLOAT Minimum confidence score (0-1) (default: 0.7)
  --evidence_contribution_threshold FLOAT
                          Minimum contribution of experimental/database evidence to total
                          score (0-1) (default: 0.2)
  --required_evidence_types STRING
                          Comma-separated list of required evidence types
                          (default: experimental,database)
  --include_second_shell BOOL
                          Whether to include second shell interactions (True/False)
                          (default: True)
  --additional_nodes INT  Number of additional proteins to include in second shell
                          (default: 1000)
  --output_prefix STRING  Prefix for output files (default: string_output)
  --node_spacing FLOAT    Spacing factor for nodes in network visualization (default: 2.0)
                          Higher values produce more spread-out networks
  -h, --help              Show this help message and exit

INPUT FILE FORMAT:
  A comma-separated text file containing gene symbols, for example:
  TP53, MDM2, CDKN1A, BAX, BBC3, BCL2

OUTPUT FILES:
  [prefix]_interactions.csv       All interactions with detailed evidence scores
  [prefix]_unique_interactors.csv List of all unique proteins with ENSEMBL IDs and confidence scores
  [prefix]_network.png            Network visualization of protein interactions
  [prefix]_no_interactions.csv    Genes with no identified interactions

EXAMPLES:
  # Basic usage with default parameters
  python string_ppi.py --input_list genes.txt

  # Custom configuration with higher confidence and spacing
  python string_ppi.py --input_list genes.txt --score_threshold 0.9 --node_spacing 3.0

  # Running with more second shell proteins and lower evidence threshold
  python string_ppi.py --input_list genes.txt --additional_nodes 2000 --evidence_contribution_threshold 0.1
```
