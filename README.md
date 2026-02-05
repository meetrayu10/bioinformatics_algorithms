
# Bioinformatics Basic Algorithms

A collection of fundamental algorithms for computational biology and bioinformatics.

## Algorithms Included
### Levenshtein Distance
Calculate edit distance between DNA/protein sequences from FASTA files.
**Features:**
- Dynamic programming implementation
- FASTA file parsing
- Displays DP matrix for educational purposes
- Time complexity: O(m×n)
- 
**##Binary search algorithm**
Traditional linear search operates at  O(n) complexity. By extracting and sorting k-mers (subsequences of length 
k), this implementation achieves 
O(logn) lookup time for exact pattern matching.
Workflow:
Indexing: Extract all k-mers from sequence (length = pattern length)
Sorting: Sort k-mers lexicographically (O(nlogn) preprocessing)
Search: Apply binary search to locate target patterns (O(logn) per query)

This approach is particularly efficient when:
Searching for multiple patterns in the same sequence (amortized preprocessing cost)
Analyzing large genomes (>1M bp) where repeated linear scans are prohibitive
Performing batch restriction enzyme site mapping 

#Genomic Binary Search 
A high-performance Python toolkit for searching and analyzing genomic sequences using binary search algorithms. Optimized for large-scale FASTA file processing with 
O(logn) search complexity.
