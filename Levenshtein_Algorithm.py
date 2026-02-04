# #Levenshtein distance 
# #Levenshtein distance = minimum number of single‑character edits (insert, delete, substitute) 
# # needed to turn one string into another 
import numpy as np
def levenshtein_distance(source, target, show_matrix=True):
    m, n = len(source), len(target)
    matrix = np.zeros((m+1, n+1), dtype=int)
    scoring_mat = np.zeros((m+1, n+1), dtype=np.int_)
    
    for j in range(m+1): scoring_mat[j][0] = j
    for i in range(n+1): scoring_mat[0][i] = i
    
    for j in range(1, m+1):
        for i in range(1, n+1):
            cost = 0 if source[j-1] == target[i-1] else 1
            scoring_mat[j][i] = min(
                scoring_mat[j-1][i] + 1,
                scoring_mat[j][i-1] + 1,
                scoring_mat[j-1][i-1] + cost
            )
            matrix[i, j] = min(
                matrix[i-1, j] + 1,      # deletion
                matrix[i, j-1] + 1,      # insertion
                matrix[i-1, j-1] + cost  # substitution
            )
    
    # Display matrix
    if show_matrix:
        print("\nDP Matrix:")
        print(matrix)
    
    return matrix[m, n]


def parse_fasta(filename):
    seq = ""
    with open(filename) as f:           # Open & read
        for line in f:
            line = line.strip()
            if line.startswith(">"):    # Skip header
                continue
            seq += line.upper()         # Add sequence
    return seq

print("=== Simple FASTA Distance ===\n")

# 1. Get filenames interactively
f1 = input("File 1: ")
f2 = input("File 2: ")

# 2. Parse FASTA
print("\nParsing...")
seq1 = parse_fasta(f1)
seq2 = parse_fasta(f2)

# 3. Calculate & show
print(f"Seq1: {seq1}")
print(f"Seq2: {seq2}")
dist = levenshtein_distance(seq1, seq2)
print(f"Distance: {dist}")
