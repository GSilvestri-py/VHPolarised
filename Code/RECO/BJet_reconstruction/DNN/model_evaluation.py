import torch
import numpy as np
import pandas as pd
import configparser
import importlib
import pickle
import seaborn as sb
import matplotlib.pyplot as plt
from sklearn.metrics import (confusion_matrix, classification_report, roc_curve, auc)
import os
import shap



# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def recreate_model(config, test_df):

    training_output_path = config.get('saving', 'out_path')
    out_file_base        = config.get('saving', 'out_file_base')

    # --- import model ---
    library_name = config.get('model', 'library_name')
    model_name   = config.get('model', 'model_name')

    model_lib    = importlib.import_module(library_name)
    model_class  = getattr(model_lib, model_name)

    # --- get input dim - reflect what done in test_DNN
    X_test = torch.tensor(test_df.iloc[:, :-1].values, dtype=torch.float32)
    input_dim = X_test.shape[1]

    print(test_df)
    print(X_test)

    # --- hidden_units: = input_dim ---
    hidden_units = input_dim

    hidden_layers = config.getint('training', 'hidden_layers', fallback=3)

    out_dim      = config.getint('model', 'classes_num')
    dropout_prob = config.getfloat('training', 'dropout_prob')
    layer_norm   = config.getboolean('training', 'layer_norm')

    model = model_class(input_dim, hidden_units, hidden_layers, out_dim, dropout_prob, layer_norm)

    # --- load saved weights
    best_model_path = os.path.join(training_output_path, out_file_base + '_best_model_partial.pth')
    print(f" importing model saved in {best_model_path}")
    state_dict = torch.load(best_model_path, map_location="cpu")
    model.load_state_dict(state_dict)

    # --- print ---
    print("\n--- SHAPES NEL CHECKPOINT ---")
    for k, v in state_dict.items():
        if "weight" in k:
            print(k, v.shape)

    print("MODEL WEIGHT SAMPLE:", model.model[0].weight[0][:5])
    print("CHECKPOINT WEIGHT SAMPLE:", state_dict['model.0.weight'][0][:5])


    return model

# ------------------------------------------------------------------------------------------------------------------

def dnn_training_output(model, test_df, output_path):

    # --------- Prepare tensors ---------
    X = torch.tensor(test_df.iloc[:, :-1].values, dtype=torch.float32)
    y = torch.tensor(test_df.iloc[:, -1].values, dtype=torch.long)

    var_names = list(test_df.columns[:-1])

    print("SHAPE X =", X.shape)
    print("SHAPE y =", y.shape)

    device = next(model.parameters()).device
    X = X.to(device)
    y_np = y.numpy()   

    # --------- Forward pass ---------
    model.eval()
    with torch.no_grad():
        score = model(X)
        probs = torch.softmax(score, dim=1).cpu().numpy()

    y_prob_b    = probs[:, 1]
    #y_prob_nonb = probs[:, 0]
    y_pred      = np.argmax(probs, axis=1)

    # --------- Confusion matrix ---------

    cm = confusion_matrix(y_np, y_pred)
    sb.heatmap(cm, annot=True, fmt='.0f', cmap='Reds')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title('Confusion Matrix - DNN')
    plt.savefig(os.path.join(output_path, "DNN_confusion_matrix.png"), dpi=300)
    plt.close()

    cm_n = confusion_matrix(y_np, y_pred, normalize='true')
    sb.heatmap(cm_n, annot=True, cmap='Reds')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title('Confusion Matrix (normalized) - DNN')
    plt.savefig(os.path.join(output_path, "DNN_confusion_matrix_norm.png"), dpi=300)
    plt.close()

    # --------- Response score distribution ---------
    b_scores    = [p for (p, t) in zip(y_prob_b, y_np) if t == 1]
    nonb_scores = [p for (p, t) in zip(y_prob_b, y_np) if t == 0]

    plt.figure()
    sb.histplot(b_scores, bins=40, color='royalblue', stat='probability',
                binrange=(0,1), label='b jet', alpha=0.5)
    sb.histplot(nonb_scores, bins=40, color='firebrick', stat='probability',
                binrange=(0,1), label='non b jet', alpha=0.5)
    plt.title('DNN response score distribution')
    plt.xlabel('Score')
    plt.ylabel('Event fraction')
    plt.legend()
    plt.savefig(os.path.join(output_path, "DNN_response_score.png"), dpi=300)
    plt.close()

    # --------- ROC curve ---------
    fpr, tpr, _ = roc_curve(y_np, y_prob_b)
    auc_dnn = auc(fpr, tpr)

    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, label=f"AUC = {auc_dnn:.2f}", color="royalblue")
    plt.plot([0,1], [0,1], 'k--')
    plt.title('ROC curve - DNN')
    plt.xlabel('False Positive rate')
    plt.ylabel('True Positive rate')
    plt.legend()
    plt.savefig(os.path.join(output_path, "ROC_DNN.png"), dpi=300)
    plt.close()


    # --------- Classification report ---------

    print(classification_report(y_np, y_pred, target_names=["non b jet", "b jet"]))


    # --------- SHAP FEATURE IMPORTANCE (KernelExplainer) ---------
    print("------- SHAP FEATURE IMPORTANCE (KernelExplainer) -------")

    model_cpu = model.cpu()

    def f(X):
        X_t = torch.tensor(X, dtype=torch.float32)
        with torch.no_grad():
            out = model_cpu(X_t)
            prob = torch.softmax(out, dim=1)[:, 1]  #class 1
        return prob.numpy()

    #background sample
    X_cpu = X.cpu().numpy()
    background = X_cpu[np.random.choice(X_cpu.shape[0], size=min(100, len(X_cpu)), replace=False)]

    sample = X_cpu[np.random.choice(X_cpu.shape[0], size=min(300, len(X_cpu)), replace=False)]

    explainer = shap.KernelExplainer(f, background)
    shap_values = explainer.shap_values(sample)

    shap_class1 = shap_values

    plt.figure(figsize=(8,6))
    shap.summary_plot(shap_class1, sample, feature_names=var_names, show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "DNN_SHAP_summary.png"), dpi=300)
    plt.close()
    print(f"Saved: DNN_SHAP_summary.png in {output_path}")

    shap_abs = np.abs(shap_class1).mean(axis=0)
    imp_df = pd.DataFrame({
                            "feature": var_names,
                            "importance": shap_abs
                            }).sort_values("importance", ascending=False)

    plt.figure(figsize=(8,6))
    sb.barplot(data=imp_df, x="importance", y="feature", palette="viridis")
    plt.title("DNN Feature Importance (SHAP KernelExplainer)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "DNN_SHAP_barplot.png"), dpi=300)
    plt.close()

    print(f"Saved: DNN_SHAP_barplot.png in {output_path}")

    return

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

def plot_all_trends(out_path, out_file_base):

    trend_file = os.path.join(out_path, out_file_base + "_trend_partial.pickle")

    with open(trend_file, "rb") as f:
        train_log, valid_log = pickle.load(f)

    #from return of test_loop function -> change if return changes
    metric_names = ["Loss", "Accuracy", "Precision", "Recall", "F1"]

    num_vars = len(train_log[0])

    if num_vars != len(metric_names):
        print("Numero di metriche diverso da quello atteso")
        metric_names = [f"Var {i}" for i in range(num_vars)]

    #2x3 subplots
    rows, cols = 2, 3
    fig, axes = plt.subplots(rows, cols, figsize=(15, 8), sharex=True)

    axes = axes.flatten() 

    for i in range(num_vars):
        train_vals = [epoch[i] for epoch in train_log]
        valid_vals = [epoch[i] for epoch in valid_log]

        axes[i].plot(train_vals, label=f"Train {metric_names[i]}")
        axes[i].plot(valid_vals, label=f"Valid {metric_names[i]}")
        axes[i].set_title(metric_names[i])
        axes[i].grid(True)
        axes[i].legend()

    #remove empy subplots
    for j in range(num_vars, rows*cols):
        axes[j].axis("off")

    fig.suptitle("Training & Validation Trends", fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, "dnn_trends.png"), dpi=300)
    plt.close()

    return

    

# ------------------------------------------------------------------------------------------------------------------


def main():

    # read config file
    config = configparser.ConfigParser()
    config.read("configfile_path")

    # load test dataframe
    training_output_path = config.get('saving', 'out_path')
    out_file_base        = config.get('saving', 'out_file_base')
    file = os.path.join(training_output_path, "test_set_dataframe.pkl")
    test_df = pd.read_pickle(file)

    # re-create model
    model = recreate_model(config, test_df)

    #evaluate
    training_output_path = config.get('saving', 'out_path')
    dnn_training_output(model, test_df, training_output_path)

    plot_all_trends(training_output_path, out_file_base)

if __name__ == "__main__" : 
    main ()