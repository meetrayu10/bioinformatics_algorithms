#GENOME DATA ANALYSIS USING BINARY SEARCH ALGORITHM 
#NEEDS sorted list 
#Time complexity: O(n) 
"""
Binary Search Algorithm for Genomic Data Analysis

WHY THIS APPROACH:
- Binary search is O(log n) vs linear O(n), making it much faster for large genomic datasets
- By sorting k-mers (sequence fragments) first, we can use binary search to find patterns
- This is ideal for genome-wide searches where the same operation repeats across millions of bases

APPLICATIONS:
- Finding start/stop codons in genes (ATG, TAA, TAG, TGA)
- Locating restriction enzyme sites for DNA cloning
- Identifying GC-rich regions (often regulatory elements)
- Pattern matching for promoters and enhancers
- Genomic sequence analysis and annotation

Efficiently searches for DNA/protein sequences in FASTA files
"""

from Bio import SeqIO
from Bio import Align
aligner = Align.PairwiseAligner()
from typing import List, Tuple, Dict


class GenomicBinarySearch:
    """Binary search implementation for genomic data analysis"""
    
    def __init__(self, fasta_file: str):
        """
        Initialize with FASTA file
        Args:
            fasta_file: Path to FASTA file
        """
        self.fasta_file = fasta_file
        self.sequences = {}
        self.sorted_kmers = []
        self.load_sequences()
    
    def load_sequences(self):
        """Load sequences from FASTA file"""
        try:
            for record in SeqIO.parse(self.fasta_file, "fasta"):
                self.sequences[record.id] = str(record.seq)
            print(f"✓ Loaded {len(self.sequences)} sequence(s) from FASTA file")
        except FileNotFoundError:
            print(f"Error: FASTA file '{self.fasta_file}' not found")
            return False
        except Exception as e:
            print(f"Error reading FASTA file: {e}")
            return False
        return True
    
    def binary_search(self, sorted_list: List[str], target: str) -> int:
        """
        Classic binary search algorithm
        Args:
            sorted_list: List of sorted items
            target: Item to search for
        Returns:
            Index of target if found, -1 otherwise
        """
        left, right = 0, len(sorted_list) - 1
        
        while left <= right:
            mid = (left + right) // 2
            mid_value = sorted_list[mid]
            
            if mid_value == target:
                return mid
            elif mid_value < target:
                left = mid + 1
            else:
                right = mid - 1
        
        return -1
    
    def search_substring_binary(self, seq_id: str, pattern: str) -> List[int]:
        """
        Search for pattern in sequence using sliding window + binary search
        Args:
            seq_id: Sequence identifier
            pattern: Pattern to search (e.g., 'ATG' for start codon)
        Returns:
            List of positions where pattern is found
        """
        if seq_id not in self.sequences:
            print(f"Sequence '{seq_id}' not found")
            return []
        
        sequence = self.sequences[seq_id]
        pattern_len = len(pattern)
        positions = []
        
        # Extract all k-mers of pattern length
        kmers = []
        for i in range(len(sequence) - pattern_len + 1):
            kmers.append((sequence[i:i + pattern_len], i))
        
        # Sort k-mers
        kmers_sorted = sorted(kmers, key=lambda x: x[0])
        sorted_sequences = [km[0] for km in kmers_sorted]
        
        # Binary search for pattern
        idx = self.binary_search(sorted_sequences, pattern)
        
        if idx != -1:
            # Find all occurrences
            for i in range(len(sorted_sequences)):
                if sorted_sequences[i] == pattern:
                    positions.append(kmers_sorted[i][1])
        
        return sorted(positions)
    
    def find_start_codons(self, seq_id: str) -> List[int]:
        """Find ATG start codons in sequence"""
        positions = self.search_substring_binary(seq_id, 'ATG')
        return positions
    
    def find_stop_codons(self, seq_id: str) -> Dict[str, List[int]]:
        """Find TAA, TAG, TGA stop codons"""
        stop_codons = {'TAA': [], 'TAG': [], 'TGA': []}
        for stop in stop_codons:
            positions = self.search_substring_binary(seq_id, stop)
            stop_codons[stop] = positions
        return stop_codons
    
    def find_gc_rich_regions(self, seq_id: str, window_size: int = 100, gc_threshold: float = 0.6) -> List[Tuple[int, int]]:
        """
        Find GC-rich regions using binary search on sorted GC content
        Args:
            seq_id: Sequence identifier
            window_size: Window size for calculating GC content
            gc_threshold: Threshold for GC content (0-1)
        Returns:
            List of (start, end) positions of GC-rich regions
        """
        if seq_id not in self.sequences:
            print(f"Sequence '{seq_id}' not found")
            return []
        
        sequence = self.sequences[seq_id]
        gc_regions = []
        gc_contents = []
        
        # Calculate GC content for each window
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            gc_count = window.count('G') + window.count('C')
            gc_percentage = gc_count / window_size
            gc_contents.append((gc_percentage, i))
        
        # Sort by GC content
        gc_sorted = sorted(gc_contents, key=lambda x: x[0])
        gc_values = [gc[0] for gc in gc_sorted]
        
        # Binary search to find threshold index
        threshold_idx = 0
        left, right = 0, len(gc_values) - 1
        while left <= right:
            mid = (left + right) // 2
            if gc_values[mid] < gc_threshold:
                threshold_idx = mid
                left = mid + 1
            else:
                right = mid - 1
        
        # Collect all regions above threshold
        for i in range(threshold_idx, len(gc_sorted)):
            gc_perc, position = gc_sorted[i]
            gc_regions.append((position, position + window_size, gc_perc))
        
        return gc_regions
    
    def search_restriction_sites(self, seq_id: str, enzyme_site: str) -> List[int]:
        """
        Search for restriction enzyme sites
        Args:
            seq_id: Sequence identifier
            enzyme_site: Restriction enzyme recognition site (e.g., 'GAATTC' for EcoRI)
        Returns:
            List of positions
        """
        return self.search_substring_binary(seq_id, enzyme_site)
    
    def analyze_sequence(self, seq_id: str) -> Dict:
        """
        Comprehensive genomic analysis of a sequence
        Args:
            seq_id: Sequence identifier
        Returns:
            Dictionary with analysis results
        """
        if seq_id not in self.sequences:
            return {"error": f"Sequence '{seq_id}' not found"}
        
        sequence = self.sequences[seq_id]
        
        results = {
            "sequence_id": seq_id,
            "length": len(sequence),
            "start_codons": self.find_start_codons(seq_id),
            "stop_codons": self.find_stop_codons(seq_id),
            "gc_rich_regions": self.find_gc_rich_regions(seq_id),
            "base_composition": {
                "A": sequence.count('A'),
                "T": sequence.count('T'),
                "G": sequence.count('G'),
                "C": sequence.count('C')
            }
        }
        
        return results
    
    def print_analysis(self, analysis: Dict):
        """Pretty print analysis results"""
        print("\n" + "="*60)
        print(f"GENOMIC ANALYSIS REPORT: {analysis['sequence_id']}")
        print("="*60)
        print(f"Sequence Length: {analysis['length']} bp")
        print(f"\nBase Composition:")
        for base, count in analysis['base_composition'].items():
            percentage = (count / analysis['length']) * 100
            print(f"  {base}: {count} ({percentage:.2f}%)")
        
        print(f"\nStart Codons (ATG): {len(analysis['start_codons'])} found at positions {analysis['start_codons'][:5]}" + 
              ("..." if len(analysis['start_codons']) > 5 else ""))
        
        print(f"\nStop Codons:")
        for stop, positions in analysis['stop_codons'].items():
            print(f"  {stop}: {len(positions)} found")
        
        print(f"\nGC-Rich Regions (>60% GC): {len(analysis['gc_rich_regions'])} found")
        for idx, (start, end, gc_perc) in enumerate(analysis['gc_rich_regions'][:5], 1):
            print(f"  Region {idx}: Position {start}-{end} (GC: {gc_perc*100:.2f}%)")
        
        print("="*60 + "\n")


def create_sample_fasta(filename: str):
    """Create a sample FASTA file for testing"""
    sample_fasta = """>sequence_1
ATGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTT
CCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTTAG
>sequence_2
ATGGGGTTTAAACCCGGGAAAATTTCCCGGGAAATTTCCCGGGAAATTTCCC
GGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTTAA
"""
    
    with open(filename, 'w') as f:
        f.write(sample_fasta)
    print(f"Sample FASTA file created: {filename}")


# Example usage
if __name__ == "__main__":
    # Create sample FASTA for testing
    sample_file = "sample_genomic.fasta"
    create_sample_fasta(sample_file)
    
    # Initialize genomic binary search
    print("\n" + "="*60)
    print("BINARY SEARCH FOR GENOMIC DATA ANALYSIS")
    print("="*60)
    
    genome = GenomicBinarySearch(sample_file)
    
    # Analyze each sequence
    for seq_id in genome.sequences.keys():
        analysis = genome.analyze_sequence(seq_id)
        genome.print_analysis(analysis)
    
    # Search for specific patterns
    print("\nCUSTOM PATTERN SEARCHES:")
    print("-" * 40)
    print("Searching for 'AAA' in sequence_1:")
    positions = genome.search_substring_binary('sequence_1', 'AAA')
    print(f"Found at positions: {positions}")
    
    print("\nSearching for 'GGG' in sequence_1:")
    positions = genome.search_substring_binary('sequence_1', 'GGG')
    print(f"Found at positions: {positions}")
    
    print("\nSearching for restriction site 'AAATTT' in sequence_1:")
    positions = genome.search_restriction_sites('sequence_1', 'AAATTT')
    print(f"Found at positions: {positions}")
