import argparse
import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve, auc, precision_score, recall_score, f1_score

def compute_metrics(real_df, pred_df):
    
    df = pd.merge(real_df, pred_df, on='Species_ID', how='outer').fillna(0)
    real = df['real']
    pred = df['predict']

    
    real_binary = (real > 0).astype(int)

    
    pred_binary = pred_df.set_index('Species_ID').reindex(df['Species_ID']).notnull().astype(int).values.flatten()


    
    precision = precision_score(real_binary, pred_binary, zero_division=0)
    recall = recall_score(real_binary, pred_binary, zero_division=0)
    f1 = f1_score(real_binary, pred_binary, zero_division=0)

    
    true_positive_idx = df.loc[(real_binary == 1) & (pred_binary == 1), 'Species_ID'].tolist()
    false_positive_idx = df.loc[(real_binary == 0) & (pred_binary == 1), 'Species_ID'].tolist()
    false_negative_idx = df.loc[(real_binary == 1) & (pred_binary == 0), 'Species_ID'].tolist()

    true_positive_count = len(true_positive_idx)
    false_positive_count = len(false_positive_idx)
    false_negative_count = len(false_negative_idx)

    
    try:
        prec_curve, rec_curve, _ = precision_recall_curve(real_binary, pred)
        aupr = auc(rec_curve, prec_curve)
    except ValueError:
        aupr = np.nan

    
    afe = np.abs(pred - real).mean()
    mask = real > 0
    rfe = (np.abs(pred[mask] - real[mask]) / real[mask]).mean() if mask.any() else np.nan

    
    l1 = np.sum(np.abs(pred - real))
    l2 = np.sqrt(np.sum((pred - real) ** 2))

    
    bc = np.sum(np.abs(pred - real)) / np.sum(pred + real + 1e-12)

    return {
        "Precision": precision,
        "Recall": recall,
        "F1_score": f1,
        "True_Positive_Count": true_positive_count,
        "False_Positive_Count": false_positive_count,
        "False_Negative_Count": false_negative_count,
        "True_Positive_Species_ID": true_positive_idx,
        "False_Positive_Species_ID": false_positive_idx,
        "False_Negative_Species_ID": false_negative_idx,
        "AUPR": aupr,
        "AFE": afe,
        "RFE": rfe,
        "L1_distance": l1,
        "L2_distance": l2,
        "BC_distance": bc
    }

def main(real_file, predict_file, output_file):
    real_df = pd.read_csv(real_file, sep='\t', header=None, names=['Species_ID', 'real'], dtype={0: str})
    pred_df = pd.read_csv(predict_file, sep='\t', header=None, names=['Species_ID', 'predict'], dtype={0: str})

    metrics = compute_metrics(real_df, pred_df)

    
    result_dict = {
        "Precision": metrics["Precision"],
        "Recall": metrics["Recall"],
        "F1_score": metrics["F1_score"],
        "True_Positive_Count": metrics["True_Positive_Count"],
        "False_Positive_Count": metrics["False_Positive_Count"],
        "False_Negative_Count": metrics["False_Negative_Count"],
        "AUPR": metrics["AUPR"],
        "AFE": metrics["AFE"],
        "RFE": metrics["RFE"],
        "L1_distance": metrics["L1_distance"],
        "L2_distance": metrics["L2_distance"],
        "BC_distance": metrics["BC_distance"]
    }

    
    result_dict["True_Positive_Species_ID"] = metrics["True_Positive_Species_ID"]
    result_dict["False_Positive_Species_ID"] = metrics["False_Positive_Species_ID"]
    result_dict["False_Negative_Species_ID"] = metrics["False_Negative_Species_ID"]

    
    result_df = pd.DataFrame([result_dict])
    result_df.to_csv(output_file, sep='\t', index=False)

    print(f"✅ Evaluation results saved to: {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate prediction against real species abundance.")
    parser.add_argument('--real', required=True, help="Path to real abundance file")
    parser.add_argument('--predict', required=True, help="Path to predicted abundance file")
    parser.add_argument('--output', required=True, help="Path to output result TSV file")
    args = parser.parse_args()

    main(args.real, args.predict, args.output)
