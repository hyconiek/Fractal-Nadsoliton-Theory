#!/usr/bin/env python3
import os
import sys

def decode_scd(file_path, n, d):
    """
    Decodes Meringer's binary shortcode (.scd) format.
    n: number of vertices
    d: regularity (degree)
    Returns: list of adjacency lists (1-based index)
    """
    m = n * d // 2  # Expected code length in bytes
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Shortcode file {file_path} not found.")
        
    with open(file_path, "rb") as f:
        data = f.read()
        
    graphs = []
    current_code = bytearray(m)
    pos = 0
    
    while pos < len(data):
        prefix_len = data[pos]
        pos += 1
        
        tail_len = m - prefix_len
        if pos + tail_len > len(data):
            break  # EOF / Padding
            
        tail = data[pos:pos+tail_len]
        pos += tail_len
        
        # Merge prefix and tail
        current_code[prefix_len:] = tail
        
        # Reconstruct adjacency list (1-based index)
        adj = {i: [] for i in range(1, n+1)}
        code_idx = 0
        
        for v in range(1, n+1):
            needed = d - len(adj[v])
            higher_neighbors = list(current_code[code_idx : code_idx + needed])
            code_idx += needed
            
            for w in higher_neighbors:
                adj[v].append(w)
                adj[w].append(v)
        
        # Sort neighbor lists for consistency
        decoded_adj = {v: sorted(adj[v]) for v in range(1, n+1)}
        graphs.append(decoded_adj)
        
    return graphs

if __name__ == "__main__":
    default_file = "16_4_4.scd"
    if len(sys.argv) > 1:
        default_file = sys.argv[1]
        
    if not os.path.exists(default_file):
        print(f"Usage: python3 decode_meringer.py [path_to_scd_file]")
        print(f"Default file '{default_file}' not found.")
        sys.exit(1)
        
    print(f"Decoding Meringer shortcode file: {default_file} (n=16, d=4)...")
    try:
        graphs = decode_scd(default_file, 16, 4)
        print(f"Successfully decoded {len(graphs)} graphs.")
        print("\nExample - Graph 1 adjacency list:")
        for v in sorted(graphs[0].keys()):
            print(f"  {v} : {' '.join(map(str, graphs[0][v]))}")
    except Exception as e:
        print(f"Error during decoding: {e}")
        sys.exit(1)
