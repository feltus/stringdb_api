import requests
import pandas as pd
import json
from typing import List, Dict, Any, Optional, Set, Tuple
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import argparse
import sys
import os

class STRINGPhysicalInteractionFetcher:
    """Class to retrieve direct physical protein interactions from STRING database"""
    
    BASE_URL = "https://string-db.org/api"
    
    def __init__(self, species: int = 9606):
        """
        Initialize with species ID
        
        Args:
            species: NCBI taxonomy ID (default: 9606 for human)
        """
        self.species = species
        
    def get_physical_interactions(self, 
                                gene_list: List[str],
                                score_threshold: float = 0.7,
                                evidence_contribution_threshold: float = 0.2,
                                required_evidence_types: List[str] = ["experimental", "database"],
                                include_second_shell: bool = True,
                                additional_nodes: int = 1000) -> tuple:
        """
        Get direct physical protein interactions for a list of genes
        
        Args:
            gene_list: List of HUGO gene symbols
            score_threshold: Minimum overall confidence score (0-1)
            evidence_contribution_threshold: Minimum contribution of evidence types to total score
            required_evidence_types: Evidence types to consider
            include_second_shell: Whether to include 2nd shell interactions
            additional_nodes: Number of additional proteins to include in 2nd shell
            
        Returns:
            Tuple of (DataFrame with filtered interactions, Set of genes with no interactions)
        """
        # STRING API endpoint for network
        url = f"{self.BASE_URL}/json/network"
        
        # Parameters for the API request
        params = {
            "identifiers": "%0d".join(gene_list),
            "species": self.species,
            "caller_identity": "python_string_app",
            "required_score": int(score_threshold * 1000),  # Convert to STRING format (0-1000)
            "network_type": "physical"  # Physical interactions only
        }
        
        # Add second shell parameter if requested
        if include_second_shell:
            params["add_nodes"] = additional_nodes
        
        try:
            # Make the API request
            response = requests.post(url, data=params)
            response.raise_for_status()
            
            # Parse the JSON response
            data = response.json()
            
            if not data:
                print(f"No interactions found for the provided genes with score ≥ {score_threshold}")
                # If no interactions found, all genes have no interactions
                return pd.DataFrame(), set(gene_list)
            
            # Convert to DataFrame
            network_df = pd.DataFrame(data)
            
            if network_df.empty:
                return pd.DataFrame(), set(gene_list)
            
            # Print available columns and a sample row to help diagnose issues
            print(f"\nAvailable columns in response: {', '.join(network_df.columns)}")
            print("\nSample interaction row:")
            print(network_df.iloc[0][['preferredName_A', 'preferredName_B', 'score'] + 
                                   [col for col in network_df.columns if col.endswith('_score')]])
                
            # Normalize scores to 0-1 scale
            for col in network_df.columns:
                if col.endswith('_score') or col == 'score':
                    network_df[f"{col}_normalized"] = network_df[col] / 1000
            
            # Keep interactions with overall score above threshold 
            filtered_df = network_df[network_df['score'] >= int(score_threshold * 1000)]
            print(f"\nAfter overall score filtering: {len(filtered_df)} interactions")
            
            # Calculate the contribution of experimental OR database evidence to total score
            if 'experimental_score' in filtered_df.columns and 'database_score' in filtered_df.columns:
                # Show distribution of evidence types
                print("\nEvidence score distribution (mean values):")
                for col in [c for c in filtered_df.columns if c.endswith('_score')]:
                    print(f"{col}: {filtered_df[col].mean():.2f}")
                
                # Calculate normalized scores (0-1 scale) for each evidence type
                for evidence_type in required_evidence_types:
                    score_col = f"{evidence_type}_score"
                    if score_col in filtered_df.columns:
                        filtered_df[f"{score_col}_normalized"] = filtered_df[score_col] / filtered_df['score']
                
                # CHANGED: Instead of summing contributions, check if ANY required evidence type meets threshold
                evidence_mask = False  # Start with False
                for evidence_type in required_evidence_types:
                    norm_score_col = f"{evidence_type}_score_normalized"
                    if norm_score_col in filtered_df.columns:
                        # Use logical OR to combine masks
                        evidence_mask = evidence_mask | (filtered_df[norm_score_col] >= evidence_contribution_threshold)
                
                # Show the evidence contribution distribution
                if len(filtered_df) > 0:
                    print(f"\nEvidence contribution stats:")
                    for evidence_type in required_evidence_types:
                        norm_score_col = f"{evidence_type}_score_normalized"
                        if norm_score_col in filtered_df.columns:
                            print(f"{evidence_type} min: {filtered_df[norm_score_col].min():.2f}")
                            print(f"{evidence_type} max: {filtered_df[norm_score_col].max():.2f}")
                            print(f"{evidence_type} mean: {filtered_df[norm_score_col].mean():.2f}")
                
                # Keep interactions where experimental OR database contribution meets the threshold
                before_filtering = len(filtered_df)
                filtered_df = filtered_df[evidence_mask]
                print(f"\nAfter evidence contribution filtering: {len(filtered_df)} interactions " +
                      f"(removed {before_filtering - len(filtered_df)})")
                print(f"Filter logic: Kept interactions with EITHER {' OR '.join(required_evidence_types)} " +
                      f"evidence contribution ≥ {evidence_contribution_threshold}")
            
            # Identify which proteins are from the original set vs second shell
            if 'preferredName_A' in filtered_df.columns and 'preferredName_B' in filtered_df.columns:
                original_proteins = set(gene_list)
                filtered_df['A_original'] = filtered_df['preferredName_A'].isin(original_proteins)
                filtered_df['B_original'] = filtered_df['preferredName_B'].isin(original_proteins)
                filtered_df['interaction_type'] = 'second_shell'
                
                # Mark interactions between original proteins
                both_original = filtered_df['A_original'] & filtered_df['B_original']
                filtered_df.loc[both_original, 'interaction_type'] = 'original'
                
                # Count different types of interactions
                original_count = sum(filtered_df['interaction_type'] == 'original')
                second_shell_count = sum(filtered_df['interaction_type'] == 'second_shell')
                
                # Get unique second shell proteins
                second_shell_proteins = set()
                for _, row in filtered_df[filtered_df['interaction_type'] == 'second_shell'].iterrows():
                    if not row['A_original']:
                        second_shell_proteins.add(row['preferredName_A'])
                    if not row['B_original']:
                        second_shell_proteins.add(row['preferredName_B'])
            
            # Find genes with no interactions
            genes_with_interactions = set()
            if 'preferredName_A' in filtered_df.columns and 'preferredName_B' in filtered_df.columns:
                for _, row in filtered_df.iterrows():
                    if row['A_original']:
                        genes_with_interactions.add(row['preferredName_A'])
                    if row['B_original']:
                        genes_with_interactions.add(row['preferredName_B'])
            
            genes_with_no_interactions = set(gene_list) - genes_with_interactions
            
            # Print summary
            print(f"\nFound {len(data)} total interactions")
            print(f"After filtering for physical interactions with experimental/database evidence: {len(filtered_df)} interactions")
            
            if 'interaction_type' in filtered_df.columns:
                print(f"Original set interactions: {original_count}")
                print(f"Second shell interactions: {second_shell_count}")
                print(f"Number of unique second shell proteins: {len(second_shell_proteins)}")
                
                if second_shell_proteins:
                    print(f"Second shell proteins (top 5): {list(second_shell_proteins)[:5]}")
                    if len(second_shell_proteins) > 5:
                        print(f"... and {len(second_shell_proteins)-5} more")
            
            print(f"\nGenes with no interactions: {len(genes_with_no_interactions)}")
            if genes_with_no_interactions:
                print(f"Examples: {list(genes_with_no_interactions)[:5]}")
                if len(genes_with_no_interactions) > 5:
                    print(f"... and {len(genes_with_no_interactions)-5} more")
            
            return filtered_df, genes_with_no_interactions
            
        except requests.exceptions.RequestException as e:
            print(f"Error querying STRING API: {e}")
            return pd.DataFrame(), set(gene_list)
    
    def visualize_network(self, network_df: pd.DataFrame, output_file: Optional[str] = None,
                         node_spacing: float = 2.0) -> None:
        """
        Create a network visualization of protein interactions with adjustable node spacing
        
        Args:
            network_df: DataFrame with filtered interactions
            output_file: Optional file path to save the visualization
            node_spacing: Factor to control spacing between nodes (higher = more spread out)
        """
        if network_df.empty:
            print("No interactions to visualize")
            return
            
        try:
            # Create graph
            G = nx.Graph()
            
            # Add nodes and edges
            for _, row in network_df.iterrows():
                source = row['preferredName_A']
                target = row['preferredName_B']
                
                # Calculate edge width based on overall score
                score = row.get('score_normalized', 0.5)
                width = score * 2  # Scale for visibility
                
                # Determine edge color based on interaction type
                if 'interaction_type' in row:
                    if row['interaction_type'] == 'original':
                        color = 'blue'  # Original set interactions
                    else:
                        color = 'green'  # Second shell interactions
                else:
                    # Determine edge color based on dominant evidence type as fallback
                    exp_score = row.get('experimental_score_normalized', 0)
                    db_score = row.get('database_score_normalized', 0)
                    color = 'purple' if exp_score > db_score else 'blue'
                
                # Add nodes with attributes
                is_original_a = row.get('A_original', False)
                is_original_b = row.get('B_original', False)
                
                if source not in G:
                    G.add_node(source, original=is_original_a)
                if target not in G:
                    G.add_node(target, original=is_original_b)
                
                # Add edge with attributes
                G.add_edge(source, target, weight=width, color=color)
            
            # Create layout - use spring layout for smaller networks, kamada_kawai for larger ones
            if len(G.nodes) > 50:
                # MODIFIED: Added k parameter to adjust node spacing in kamada_kawai_layout
                pos = nx.kamada_kawai_layout(G, scale=node_spacing)
            else:
                # MODIFIED: Added k parameter to adjust node spacing in spring layout
                pos = nx.spring_layout(G, k=node_spacing/len(G.nodes), seed=42)
            
            # Draw nodes with different colors for original vs second shell
            plt.figure(figsize=(14, 12))
            
            # Get node lists by type
            original_nodes = [n for n, attr in G.nodes(data=True) if attr.get('original', False)]
            second_shell_nodes = [n for n, attr in G.nodes(data=True) if not attr.get('original', False)]
            
            # Draw nodes by type
            nx.draw_networkx_nodes(G, pos, nodelist=original_nodes, node_color='red', 
                                  node_size=700, alpha=0.8, label='Original proteins')
            nx.draw_networkx_nodes(G, pos, nodelist=second_shell_nodes, node_color='lightgreen', 
                                  node_size=500, alpha=0.6, label='Second shell proteins')
            
            # Draw edges with varying colors and widths
            edges = G.edges()
            if edges:  # Check if there are any edges
                edge_colors = [G[u][v]['color'] for u, v in edges]
                edge_widths = [G[u][v]['weight'] for u, v in edges]
                
                nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color=edge_colors, alpha=0.7)
            
            # Draw labels with different fonts for original vs second shell
            nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', 
                                   font_color='black', labels={n: n for n in original_nodes})
            nx.draw_networkx_labels(G, pos, font_size=8, font_weight='normal',
                                   font_color='darkgreen', labels={n: n for n in second_shell_nodes})
            
            # Add legend
            plt.plot([0], [0], color='blue', linewidth=2, label='Original set interactions')
            plt.plot([0], [0], color='green', linewidth=2, label='Second shell interactions')
            plt.legend(loc='upper right')
            
            plt.title("Physical Protein-Protein Interactions")
            plt.axis('off')
            
            if output_file:
                plt.savefig(output_file, bbox_inches='tight', dpi=300)
                print(f"Network visualization saved to {output_file}")
            else:
                plt.show()
                
        except Exception as e:
            print(f"Error visualizing network: {e}")
            print(f"Error details: {str(e)}")
    
    def export_to_csv(self, network_df: pd.DataFrame, output_file: str) -> None:
        """Export the interaction data to CSV with all available scores and details"""
        if network_df.empty:
            print("No interactions to export")
            return
        
        # Define all possible columns we want to include, in order of importance
        all_possible_columns = [
            # Protein identifiers
            'preferredName_A', 'preferredName_B', 'stringId_A', 'stringId_B',
            # Basic interaction information
            'score', 'score_normalized', 'interaction_type', 'evidence_contribution',
            # Normalized evidence scores
            'experimental_score_normalized', 'database_score_normalized', 
            'textmining_score_normalized', 'coexpression_score_normalized',
            'neighborhood_score_normalized', 'fusion_score_normalized', 
            'cooccurrence_score_normalized',
            # Raw evidence scores
            'experimental_score', 'database_score', 'textmining_score', 
            'coexpression_score', 'neighborhood_score', 'fusion_score', 
            'cooccurrence_score',
            # Original protein flags
            'A_original', 'B_original',
            # Any other useful columns
            'nscore', 'ascore', 'escore', 'dscore', 'tscore'
        ]
        
        # Use only columns that exist in the dataframe
        available_columns = [col for col in all_possible_columns if col in network_df.columns]
        export_df = network_df[available_columns].copy()
        
        # Rename columns for better readability
        column_mapping = {
            # Protein identifiers
            'preferredName_A': 'Protein_A',
            'preferredName_B': 'Protein_B',
            'stringId_A': 'StringID_A',
            'stringId_B': 'StringID_B',
            # Scores and types
            'score': 'Combined_Score_Raw',
            'score_normalized': 'Combined_Score',
            'interaction_type': 'Interaction_Type',
            'evidence_contribution': 'Experimental_Database_Contribution',
            # Normalized evidence scores
            'experimental_score_normalized': 'Experimental_Evidence',
            'database_score_normalized': 'Database_Evidence',
            'textmining_score_normalized': 'TextMining_Evidence',
            'coexpression_score_normalized': 'Coexpression_Evidence',
            'neighborhood_score_normalized': 'Neighborhood_Evidence',
            'fusion_score_normalized': 'Fusion_Evidence',
            'cooccurrence_score_normalized': 'Cooccurrence_Evidence',
            # Raw scores
            'experimental_score': 'Experimental_Raw',
            'database_score': 'Database_Raw',
            'textmining_score': 'TextMining_Raw',
            'coexpression_score': 'Coexpression_Raw',
            'neighborhood_score': 'Neighborhood_Raw',
            'fusion_score': 'Fusion_Raw',
            'cooccurrence_score': 'Cooccurrence_Raw',
            # Flags
            'A_original': 'A_In_Input_List',
            'B_original': 'B_In_Input_List'
        }
        
        # Rename only columns that exist in our export dataframe
        rename_dict = {old: new for old, new in column_mapping.items() if old in available_columns}
        export_df.rename(columns=rename_dict, inplace=True)
        
        # Add metadata row at the top
        with open(output_file, 'w') as f:
            f.write(f"# STRING database physical protein-protein interactions\n")
            f.write(f"# Export date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Total interactions: {len(export_df)}\n")
            f.write(f"# Score threshold: {network_df['score'].min()/1000 if 'score' in network_df.columns else 'N/A'}\n")
            f.write("\n")
        
        # Export to CSV, appending to the file with the metadata
        export_df.to_csv(output_file, index=False, mode='a')
        print(f"Interaction data exported to {output_file} with {len(available_columns)} columns of details")
    
    def export_no_interactions(self, genes_with_no_interactions: Set[str], output_file: str) -> None:
        """Export the list of genes with no interactions to a file"""
        if not genes_with_no_interactions:
            print("No genes without interactions to export")
            return
        
        # Create a dataframe with the genes
        no_interactions_df = pd.DataFrame(list(genes_with_no_interactions), columns=['Gene'])
        
        # Export to CSV
        no_interactions_df.to_csv(output_file, index=False)
        print(f"Genes with no interactions exported to {output_file}")
    
    # NEW METHOD: Export unique interactors list with ENSEMBL IDs and scores
    def export_unique_interactors(self, network_df: pd.DataFrame, output_file: str) -> None:
        """
        Export a list of all unique interactors with their ENSEMBL IDs and confidence scores
        
        Args:
            network_df: DataFrame with filtered interactions
            output_file: Path to output file
        """
        if network_df.empty:
            print("No interactions to export unique interactors from")
            return
        
        # Create sets to hold unique proteins
        unique_proteins = set()
        
        # Collect all unique proteins and their metadata
        protein_data = {}
        
        # Check for necessary columns
        has_ensembl = 'stringId_A' in network_df.columns and 'stringId_B' in network_df.columns
        has_scores = 'score' in network_df.columns
        has_origin = 'A_original' in network_df.columns and 'B_original' in network_df.columns
        
        # Process each row to extract unique proteins
        for _, row in network_df.iterrows():
            # Process protein A
            protein_a = row['preferredName_A']
            unique_proteins.add(protein_a)
            
            # Initialize protein data if not already present
            if protein_a not in protein_data:
                protein_data[protein_a] = {
                    'ensembl_id': '',
                    'interaction_count': 0, 
                    'avg_score': 0.0,
                    'max_score': 0.0,
                    'is_original': False
                }
            
            # Update protein A data
            protein_data[protein_a]['interaction_count'] += 1
            if has_ensembl:
                # Extract ENSEMBL ID from STRING ID (format: species.ENSEMBL)
                string_id = row['stringId_A']
                if '.' in string_id:
                    protein_data[protein_a]['ensembl_id'] = string_id.split('.')[1]
            if has_scores:
                score = row['score'] / 1000  # Normalize to 0-1
                protein_data[protein_a]['avg_score'] += score
                protein_data[protein_a]['max_score'] = max(protein_data[protein_a]['max_score'], score)
            if has_origin:
                protein_data[protein_a]['is_original'] = protein_data[protein_a]['is_original'] or row['A_original']
                
            # Process protein B
            protein_b = row['preferredName_B']
            unique_proteins.add(protein_b)
            
            # Initialize protein data if not already present
            if protein_b not in protein_data:
                protein_data[protein_b] = {
                    'ensembl_id': '',
                    'interaction_count': 0, 
                    'avg_score': 0.0,
                    'max_score': 0.0,
                    'is_original': False
                }
            
            # Update protein B data
            protein_data[protein_b]['interaction_count'] += 1
            if has_ensembl:
                # Extract ENSEMBL ID from STRING ID (format: species.ENSEMBL)
                string_id = row['stringId_B']
                if '.' in string_id:
                    protein_data[protein_b]['ensembl_id'] = string_id.split('.')[1]
            if has_scores:
                score = row['score'] / 1000  # Normalize to 0-1
                protein_data[protein_b]['avg_score'] += score
                protein_data[protein_b]['max_score'] = max(protein_data[protein_b]['max_score'], score)
            if has_origin:
                protein_data[protein_b]['is_original'] = protein_data[protein_b]['is_original'] or row['B_original']
        
        # Calculate average scores
        for protein in protein_data:
            if protein_data[protein]['interaction_count'] > 0:
                protein_data[protein]['avg_score'] /= protein_data[protein]['interaction_count']
        
        # Convert to DataFrame
        interactors_df = pd.DataFrame([
            {
                'Gene': protein,
                'ENSEMBL_ID': protein_data[protein]['ensembl_id'],
                'Avg_Confidence_Score': round(protein_data[protein]['avg_score'], 3),
                'Max_Confidence_Score': round(protein_data[protein]['max_score'], 3),
                'Interaction_Count': protein_data[protein]['interaction_count'],
                'In_Input_List': protein_data[protein]['is_original']
            }
            for protein in unique_proteins
        ])
        
        # Sort by interaction count (descending), then by protein name
        interactors_df = interactors_df.sort_values(
            by=['Interaction_Count', 'Gene'], 
            ascending=[False, True]
        )
        
        # Add metadata to the file
        with open(output_file, 'w') as f:
            f.write(f"# STRING database unique protein interactors\n")
            f.write(f"# Export date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Total unique proteins: {len(unique_proteins)}\n")
            f.write("\n")
        
        # Export to CSV
        interactors_df.to_csv(output_file, index=False, mode='a')
        print(f"Unique interactors exported to {output_file}")

# REINSTATED ALTERNATIVE QUERY METHOD
def query_string_directly(gene_list: List[str], 
                         species: int = 9606,
                         score_threshold: float = 0.7,
                         include_second_shell: bool = True,
                         additional_nodes: int = 1000) -> tuple:
    """
    Alternative approach using different STRING API endpoints
    
    Returns:
        Tuple of (DataFrame with filtered interactions, Set of genes with no interactions)
    """
    # Map gene symbols to STRING IDs first
    url = "https://string-db.org/api/json/get_string_ids"
    params = {
        "identifiers": "\r".join(gene_list),
        "species": species,
        "limit": 1,
        "echo_query": 1
    }
    
    try:
        print("\nTrying alternative query method...")
        print("Step 1: Mapping gene symbols to STRING IDs...")
        
        response = requests.post(url, data=params)
        response.raise_for_status()
        mapping_data = response.json()
        
        # Extract STRING IDs and create mapping
        string_ids = []
        id_to_gene = {}
        gene_to_id = {}
        unmapped_genes = set(gene_list)
        
        for item in mapping_data:
            if "stringId" in item and "queryItem" in item:
                string_id = item["stringId"]
                gene = item["queryItem"]
                string_ids.append(string_id)
                id_to_gene[string_id] = gene
                gene_to_id[gene] = string_id
                if gene in unmapped_genes:
                    unmapped_genes.remove(gene)
        
        if not string_ids:
            print("Could not map gene symbols to STRING IDs")
            return pd.DataFrame(), set(gene_list)
        
        print(f"Successfully mapped {len(string_ids)} genes to STRING IDs")
        if unmapped_genes:
            print(f"Warning: {len(unmapped_genes)} genes could not be mapped to STRING IDs")
            print(f"Examples: {list(unmapped_genes)[:5]}")
            if len(unmapped_genes) > 5:
                print(f"... and {len(unmapped_genes)-5} more")
        
        # Step 2: Now query interactions
        print("Step 2: Querying network interactions...")
        
        network_url = "https://string-db.org/api/json/network"
        network_params = {
            "identifiers": "%0d".join(string_ids),
            "species": species,
            "caller_identity": "python_string_app",
            "required_score": int(score_threshold * 1000),
            "network_type": "physical"
        }
        
        # Add second shell parameter if requested
        if include_second_shell:
            network_params["add_nodes"] = additional_nodes
            
        network_response = requests.post(network_url, data=network_params)
        network_response.raise_for_status()
        network_data = network_response.json()
        
        if not network_data:
            print("No interactions found with alternative method")
            return pd.DataFrame(), set(gene_list) - unmapped_genes
        
        # Convert to DataFrame
        network_df = pd.DataFrame(network_data)
        
        # Step 3: Get link evidence details directly
        print("Step 3: Getting detailed evidence for each interaction...")
        
        # Create pairs of proteins for evidence query
        interaction_pairs = []
        for _, row in network_df.iterrows():
            protein_a = row['stringId_A']
            protein_b = row['stringId_B']
            interaction_pairs.append(f"{protein_a}.{protein_b}")
        
        # Query detailed evidence (in batches to avoid too large requests)
        all_evidence_data = []
        batch_size = 100
        
        for i in range(0, len(interaction_pairs), batch_size):
            batch = interaction_pairs[i:i+batch_size]
            
            evidence_url = "https://string-db.org/api/json/interaction_partners"
            evidence_params = {
                "identifiers": "%0d".join(batch),
                "species": species,
                "required_score": int(score_threshold * 1000),
                "limit": 5  # We only need the specific pair
            }
            
            try:
                evidence_response = requests.post(evidence_url, data=evidence_params)
                evidence_response.raise_for_status()
                evidence_data = evidence_response.json()
                all_evidence_data.extend(evidence_data)
            except:
                print(f"Warning: Could not get detailed evidence for batch {i//batch_size + 1}")
                continue
        
        # Process evidence data and merge with network data
        if all_evidence_data:
            evidence_df = pd.DataFrame(all_evidence_data)
            
            # Count interactions with experimental or database evidence
            if 'experimental' in evidence_df.columns and 'database' in evidence_df.columns:
                exp_db_count = sum((evidence_df['experimental'] > 0) | (evidence_df['database'] > 0))
                print(f"Found {exp_db_count} interactions with experimental or database evidence")
        
        # Add readable gene names if available
        if 'preferredName_A' not in network_df.columns and 'stringId_A' in network_df.columns:
            network_df['preferredName_A'] = network_df['stringId_A'].apply(
                lambda x: id_to_gene.get(x, x.split('.')[1])
            )
        
        if 'preferredName_B' not in network_df.columns and 'stringId_B' in network_df.columns:
            network_df['preferredName_B'] = network_df['stringId_B'].apply(
                lambda x: id_to_gene.get(x, x.split('.')[1])
            )
        
        print(f"Alternative method found {len(network_df)} total interactions")
        
        # Simple filtering for interactions with experimental or database evidence
        filtered_df = network_df.copy()
        
        # If we have evidence data, use it to filter
        if all_evidence_data:
            # Create a lookup for evidence scores
            evidence_lookup = {}
            for row in all_evidence_data:
                if 'stringId_A' in row and 'stringId_B' in row:
                    key = (row['stringId_A'], row['stringId_B'])
                    evidence_lookup[key] = {
                        'experimental': row.get('experimental', 0),
                        'database': row.get('database', 0)
                    }
            
            # Add evidence scores to network data
            for i, row in filtered_df.iterrows():
                key = (row['stringId_A'], row['stringId_B'])
                if key in evidence_lookup:
                    filtered_df.at[i, 'experimental_score'] = evidence_lookup[key]['experimental']
                    filtered_df.at[i, 'database_score'] = evidence_lookup[key]['database']
            
            # Filter for interactions with experimental or database evidence
            if 'experimental_score' in filtered_df.columns and 'database_score' in filtered_df.columns:
                mask = ((filtered_df['experimental_score'] > 0) | (filtered_df['database_score'] > 0))
                filtered_df = filtered_df[mask]
        
        print(f"After filtering: {len(filtered_df)} interactions with experimental or database evidence")
        
        # Identify which proteins are from the original set vs second shell
        mapped_genes = set(gene_list) - unmapped_genes
        filtered_df['A_original'] = filtered_df['preferredName_A'].isin(mapped_genes)
        filtered_df['B_original'] = filtered_df['preferredName_B'].isin(mapped_genes)
        filtered_df['interaction_type'] = 'second_shell'
        
        # Mark interactions between original proteins
        both_original = filtered_df['A_original'] & filtered_df['B_original']
        filtered_df.loc[both_original, 'interaction_type'] = 'original'
        
        # Get stats
        original_count = sum(filtered_df['interaction_type'] == 'original')
        second_shell_count = sum(filtered_df['interaction_type'] == 'second_shell')
        
        print(f"Original set interactions: {original_count}")
        print(f"Second shell interactions: {second_shell_count}")
        
        # Find genes with no interactions
        genes_with_interactions = set()
        for _, row in filtered_df.iterrows():
            if row['A_original']:
                genes_with_interactions.add(row['preferredName_A'])
            if row['B_original']:
                genes_with_interactions.add(row['preferredName_B'])
        
        genes_with_no_interactions = mapped_genes - genes_with_interactions
        genes_with_no_interactions.update(unmapped_genes)  # Add unmapped genes to no-interactions list
        
        print(f"\nGenes with no interactions: {len(genes_with_no_interactions)}")
        if genes_with_no_interactions:
            print(f"Examples: {list(genes_with_no_interactions)[:5]}")
            if len(genes_with_no_interactions) > 5:
                print(f"... and {len(genes_with_no_interactions)-5} more")
        
        return filtered_df, genes_with_no_interactions
        
    except requests.exceptions.RequestException as e:
        print(f"Error with alternative query: {e}")
        return pd.DataFrame(), set(gene_list)

def read_gene_list_from_file(file_path):
    """Read gene list from a comma-delimited file"""
    try:
        with open(file_path, 'r') as f:
            content = f.read().strip()
            # Split by comma and strip whitespace
            genes = [gene.strip() for gene in content.split(',')]
            # Remove empty strings
            genes = [gene for gene in genes if gene]
            return genes
    except Exception as e:
        print(f"Error reading gene list from file {file_path}: {str(e)}")
        sys.exit(1)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Query STRING database for protein-protein interactions')
    
    parser.add_argument('--species', type=int, default=9606,
                        help='NCBI taxonomy ID (default: 9606 for human)')
    
    parser.add_argument('--score_threshold', type=float, default=0.7,
                        help='Minimum interaction confidence score (0-1) (default: 0.7)')
    
    parser.add_argument('--evidence_contribution_threshold', type=float, default=0.2,
                        help='Minimum contribution of experimental/database evidence (0-1) (default: 0.2)')
    
    parser.add_argument('--required_evidence_types', type=str, default='experimental,database',
                        help='Comma-separated list of required evidence types (default: experimental,database)')
    
    parser.add_argument('--include_second_shell', type=str, default='True',
                        help='Whether to include second shell interactions (True/False) (default: True)')
    
    parser.add_argument('--additional_nodes', type=int, default=1000,
                        help='Number of additional proteins to include in second shell (default: 1000)')
    
    parser.add_argument('--input_list', type=str, required=True,
                        help='Path to comma-delimited file containing gene list')
    
    parser.add_argument('--output_prefix', type=str, default='string_output',
                        help='Prefix for output files (default: string_output)')
    
    parser.add_argument('--node_spacing', type=float, default=2.0,
                        help='Spacing factor for nodes in network visualization (default: 2.0)')
    
    args = parser.parse_args()
    
    # Process boolean argument
    args.include_second_shell = args.include_second_shell.lower() == 'true'
    
    # Process comma-separated evidence types
    args.required_evidence_types = [t.strip() for t in args.required_evidence_types.split(',')]
    
    return args

# Main function with command-line arguments
def query_string_physical_ppi(gene_list: List[str], 
                             species: int = 9606,
                             score_threshold: float = 0.7,
                             evidence_contribution_threshold: float = 0.2,
                             required_evidence_types: List[str] = ["experimental", "database"],
                             include_second_shell: bool = True,
                             additional_nodes: int = 1000,
                             output_prefix: str = "string_output",
                             node_spacing: float = 2.0) -> tuple:
    """
    Query STRING for physical protein interactions
    
    Args:
        gene_list: List of HUGO gene symbols
        species: NCBI taxonomy ID (default: 9606 for human)
        score_threshold: Minimum confidence score (0-1)
        evidence_contribution_threshold: Minimum contribution of experimental/database evidence
        required_evidence_types: List of evidence types to consider
        include_second_shell: Whether to include 2nd shell interactions
        additional_nodes: Number of additional proteins to include in 2nd shell
        output_prefix: Prefix for output files
        node_spacing: Spacing factor for nodes in network visualization
        
    Returns:
        Tuple of (DataFrame with interactions, Set of genes with no interactions)
    """
    fetcher = STRINGPhysicalInteractionFetcher(species=species)
    
    print(f"Querying STRING database for physical interactions among {len(gene_list)} genes")
    print(f"Species: {species}")
    print(f"Score threshold: {score_threshold}")
    print(f"Evidence contribution threshold: {evidence_contribution_threshold}")
    print(f"Required evidence types: {', '.join(required_evidence_types)}")
    print(f"Including second shell interactions: {include_second_shell}")
    if include_second_shell:
        print(f"Maximum additional nodes: {additional_nodes}")
    print(f"Node spacing factor: {node_spacing}")
        
    network_df, genes_with_no_interactions = fetcher.get_physical_interactions(
        gene_list, 
        score_threshold=score_threshold,
        evidence_contribution_threshold=evidence_contribution_threshold,
        required_evidence_types=required_evidence_types,
        include_second_shell=include_second_shell,
        additional_nodes=additional_nodes
    )
    
    # If no interactions found, try alternative method
    if network_df.empty:
        print("\nNo interactions found with primary method, trying alternative method...")
        network_df, genes_with_no_interactions = query_string_directly(
            gene_list, 
            species=species,
            score_threshold=score_threshold,
            include_second_shell=include_second_shell,
            additional_nodes=additional_nodes
        )
    
    if not network_df.empty:
        # Export interactions to CSV
        interactions_file = f"{output_prefix}_interactions.csv"
        fetcher.export_to_csv(network_df, interactions_file)
        
        # NEW: Export unique interactors list
        interactors_file = f"{output_prefix}_unique_interactors.csv"
        fetcher.export_unique_interactors(network_df, interactors_file)
        
        # Create visualization with adjusted node spacing
        visualization_file = f"{output_prefix}_network.png"
        fetcher.visualize_network(network_df, visualization_file, node_spacing=node_spacing)
    
    # Export genes with no interactions to a separate file
    if genes_with_no_interactions:
        no_interactions_file = f"{output_prefix}_no_interactions.csv"
        fetcher.export_no_interactions(genes_with_no_interactions, no_interactions_file)
    
    return network_df, genes_with_no_interactions

# Main script
if __name__ == "__main__":
    # Parse command line arguments
    args = parse_arguments()
    
    # Read gene list from file
    print(f"Reading gene list from file: {args.input_list}")
    gene_list = read_gene_list_from_file(args.input_list)
    print(f"Found {len(gene_list)} genes in input file")
    
    # Query interactions
    interactions, genes_with_no_interactions = query_string_physical_ppi(
        gene_list,
        species=args.species,
        score_threshold=args.score_threshold,
        evidence_contribution_threshold=args.evidence_contribution_threshold,
        required_evidence_types=args.required_evidence_types,
        include_second_shell=args.include_second_shell,
        additional_nodes=args.additional_nodes,
        output_prefix=args.output_prefix,
        node_spacing=args.node_spacing
    )
    
    # Display results summary
    if not interactions.empty:
        # Show statistics about the network
        unique_proteins = set(interactions['preferredName_A']).union(set(interactions['preferredName_B']))
        print(f"\nNetwork summary:")
        print(f"Total unique proteins: {len(unique_proteins)}")
        print(f"Total interactions: {len(interactions)}")
        
        # List output files
        output_files = [
            f"{args.output_prefix}_interactions.csv",
            f"{args.output_prefix}_unique_interactors.csv",  # NEW: unique interactors file
            f"{args.output_prefix}_network.png"
        ]
        
        if genes_with_no_interactions:
            output_files.append(f"{args.output_prefix}_no_interactions.csv")
        
        print("\nOutput files:")
        for file in output_files:
            if os.path.exists(file):
                file_size = os.path.getsize(file)
                print(f"- {file} ({file_size/1024:.1f} KB)")
    else:
        print("\nNo interactions found with any method. Try using different gene symbols " + 
              "or relaxing the filtering criteria.")
