"""
Prepare 2 TSV files with haplotype abundances. true.tsv predict.tsv
Usage: python distance.py xx.treefile true.tsv predict.tsv
"""
import sys
from ete3 import Tree

def read_abundance_file(tsv_file):
    """
    Reads TSV file with two columns: haplotype_name \t abundance
    Returns a dict {name: abundance}
    """
    d = {}
    with open(tsv_file) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                raise ValueError(f"Invalid line in {tsv_file}: {line}")
            name, freq = parts
            d[name] = float(freq)
    return d

def collect_leaves(t):
    # Return list of leaf names
    return [leaf.name for leaf in t.iter_leaves()]

def build_parent_dict(t):
    # Return parent dict {child_name: parent_name}
    parent = {}
    for node in t.traverse("preorder"):
        if node.is_root():
            continue
        for child in node.get_children():
            child_name = child.name if child.name else f"_internal_{id(child)}"
            parent[child_name] = node.name if node.name else "root"
    return parent

def compute_unifrac_emd(t, P, Q):
    """
    Simplified UniFrac-EMD:
    - t: ete3 Tree
    - P: ground truth {leaf: abundance}
    - Q: predicted {leaf: abundance}
    """
    # Assign abundance 0 to leaves missing in P or Q
    leaves = collect_leaves(t)
    for leaf in leaves:
        if leaf not in P:
            P[leaf] = 0.0
        if leaf not in Q:
            Q[leaf] = 0.0

    # Postorder traversal: compute branch contributions
    emd = 0.0
    for node in t.traverse("postorder"):
        if node.is_leaf():
            continue
        # sum of P and Q in subtree
        subtree_leaves = [l.name for l in node.iter_leaves()]
        sum_P = sum(P[l] for l in subtree_leaves)
        sum_Q = sum(Q[l] for l in subtree_leaves)
        imbalance = abs(sum_P - sum_Q)
        branch_length = node.dist if node.dist else 0.0
        emd += imbalance * branch_length
    return emd

def main():
    if len(sys.argv) != 4:
        print("Usage: python distance.py file.treefile truth.tsv predict.tsv")
        sys.exit(1)

    tree_file = sys.argv[1]
    truth_file = sys.argv[2]
    pred_file = sys.argv[3]

    # Load tree
    t = Tree(tree_file, format=1)
    leaves = collect_leaves(t)
    print("Leaves in tree:", leaves)

    # Load truth and prediction abundances
    P = read_abundance_file(truth_file)
    Q = read_abundance_file(pred_file)

    print("Labels in truth (P):", list(P.keys()))
    print("Labels in pred (Q):", list(Q.keys()))

    # Warn about missing leaves
    missing_in_tree = [x for x in set(list(P.keys()) + list(Q.keys())) if x not in leaves]
    if missing_in_tree:
        print("WARNING: following haplotypes not in tree leaves:", missing_in_tree)

    # Compute UniFrac EMD
    emd = compute_unifrac_emd(t, P, Q)
    print("\nWeighted UniFrac EMD:", emd)

if __name__ == "__main__":
    main()