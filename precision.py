from Bio import SeqIO
from Levenshtein import distance as edit_distance

def load_fasta(file):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(file, "fasta")}

def evaluate(truth_fasta, pred_fasta, threshold=85):
    truth = load_fasta(truth_fasta)
    preds = load_fasta(pred_fasta)

    TP, FP, FN = 0, 0, 0
    matched_truth = set()
    matched_pred = set()

    print(f"Using edit distance threshold = {threshold}\n")

    # Iterate over each truth haplotype
    for tid, tseq in truth.items():
        best_match = None
        best_dist = float("inf")
        for pid, pseq in preds.items():
            d = edit_distance(tseq, pseq)
            if d < best_dist:
                best_dist = d
                best_match = pid
        if best_dist <= threshold:
            TP += 1
            matched_truth.add(tid)
            matched_pred.add(best_match)
            print(f"[TP] Truth {tid} matched with Pred {best_match}, dist={best_dist}")
        else:
            FN += 1
            print(f"[FN] Truth {tid} no match (best={best_match}, dist={best_dist})")

    # Find FP
    for pid in preds.keys():
        if pid not in matched_pred:
            FP += 1
            print(f"[FP] Pred {pid} has no matching truth")

    print("\n=== Summary ===")
    print(f"TP = {TP}")
    print(f"FN = {FN}")
    print(f"FP = {FP}")

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    print(f"Precision = {precision:.3f}")
    print(f"Recall    = {recall:.3f}")

    return TP, FP, FN, precision, recall

if __name__ == "__main__":
    # change the file names
    # threshold could be 1% of your read length
    evaluate("truth.fasta", "predict.fasta", threshold=85)
