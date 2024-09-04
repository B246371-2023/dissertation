import re
import numpy as np
import pandas as pd

def parse_data(file_path):
    results = []
    with open(file_path, 'r') as file:
        current_result = {}
        for line in file:
            if "Align" in line:
                match = re.search(r"Align (\S+)\.pdb \d+ with (\S+)\.pdb \d+", line)
                if match:
                    current_result['name1'] = match.group(1)
                    current_result['name2'] = match.group(2)
            elif "Twists" in line:
                match = re.search(r"Score (\d+\.\d+)", line)
                if match:
                    current_result['score'] = float(match.group(1))
            elif "P-value" in line:
                match = re.search(r"P-value (\S+)", line)
                if match:
                    current_result['p-value'] = float(match.group(1))
                    results.append(current_result)
                    current_result = {}
    return pd.DataFrame(results)



def initialize_similarity_matrix(proteins):
    size = len(proteins)
    sim_matrix = np.zeros((size, size))
    return sim_matrix

def fill_similarity_matrix_score(sim_matrix, proteins, df):
    protein_index = {protein: idx for idx, protein in enumerate(proteins)}
    for _, row in df.iterrows():
        if row['name1'] in protein_index and row['name2'] in protein_index:
            idx1 = protein_index[row['name1']]
            idx2 = protein_index[row['name2']]
            transformed_score = 1 - (row['score'] / 1000) #the formula that make tree more seperate
            # Ensure transformed_score is always positive; use a small positive number if it's not
            transformed_score = max(transformed_score, 1e-10)
            sim_matrix[idx1, idx2] = transformed_score
            sim_matrix[idx2, idx1] = transformed_score
            # sim_matrix[idx1, idx2] = row['score']
            # sim_matrix[idx2, idx1] = row['score']
    np.fill_diagonal(sim_matrix, 0)  # Diagonal elements are 0 (self-similarity)
    return sim_matrix

def fill_similarity_matrix_p_value(sim_matrix, proteins, df):
    protein_index = {protein: idx for idx, protein in enumerate(proteins)}
    for _, row in df.iterrows():
        if row['name1'] in protein_index and row['name2'] in protein_index:
            idx1 = protein_index[row['name1']]
            idx2 = protein_index[row['name2']]
            transformed_score = np.log10(row['p-value'] + 16 ) / 16
            # Ensure transformed_score is always positive; use a small positive number if it's not
            transformed_score = max(transformed_score, 0)
            sim_matrix[idx1, idx2] = transformed_score
            sim_matrix[idx2, idx1] = transformed_score
            # sim_matrix[idx1, idx2] = row['p-value']
            # sim_matrix[idx2, idx1] = row['p-value']
    np.fill_diagonal(sim_matrix, 0)  # Diagonal elements are 0 (self-similarity)
    return sim_matrix


def format_and_save_matrix(sim_matrix, proteins, file_path):
    """Formats the similarity matrix and saves it to a file."""
    num_proteins = len(proteins)
    with open(file_path, 'w') as file:
        # Write the number of proteins at the top of the file
        file.write(f"{num_proteins}\n")
        for i, protein in enumerate(proteins):
            # Format the protein name to be left-aligned with a fixed width for alignment
            formatted_name = f"{protein:<11}"
            # Format the similarity scores, ensuring they have a fixed number of decimal places
            # and are evenly spaced for readability
            values = '  '.join(f"{value:.12f}" for value in sim_matrix[i])
            # Write each line as 'protein_name scores...'
            file.write(f"{formatted_name} {values}\n")


def main():
    # Example of calling the function
    file_path = snakemake.input[0]
    df = parse_data(file_path)
    # Get a list of unique protein names
    proteins = pd.concat([df['name1'], df['name2']]).unique()
    sim_matrix = initialize_similarity_matrix(proteins)
    sim_matrix_score = fill_similarity_matrix_score(sim_matrix, proteins, df)
    sim_matrix = initialize_similarity_matrix(proteins)
    sim_matrix_p_value = fill_similarity_matrix_p_value(sim_matrix, proteins, df)
    output_score = snakemake.output[0]
    output_p_value = snakemake.output[1]
    # Example of saving the matrix
    format_and_save_matrix(sim_matrix_score, proteins, output_score)
    format_and_save_matrix(sim_matrix_p_value, proteins, output_p_value)

if __name__ == "__main__":
    main()

