import configparser
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
import xgboost as xgb
import os
import optuna

from sklearn.metrics import confusion_matrix, classification_report, roc_auc_score, roc_curve, auc
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score, log_loss
from sklearn.model_selection import train_test_split


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


def objective(trial, X_train, y_train, X_valid, y_valid):
    
    '''
    params = {
              "tree_method": "hist",
              "grow_policy": "lossguide",
              "max_leaves": trial.suggest_int("max_leaves", 64, 1024),
              "min_child_weight": trial.suggest_float("min_child_weight", 0.5, 2.0),
              "gamma": trial.suggest_float("gamma", 0.0, 1.0),
              "subsample": trial.suggest_float("subsample", 0.6, 1.0),
              "colsample_bytree": trial.suggest_float("colsample_bytree", 0.6, 1.0),
              "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.2),
              "n_estimators": 1200,
              "reg_alpha": trial.suggest_float("reg_alpha", 0.0, 1.0),
              "reg_lambda": trial.suggest_float("reg_lambda", 0.1, 5.0),
              "scale_pos_weight": trial.suggest_float("scale_pos_weight", 1.0, 5.0),
              "eval_metric": "logloss",
              "early_stopping_rounds": 50,
              }
    '''
    params = {
             "tree_method": "hist",
             "grow_policy": "lossguide",
             "max_leaves": trial.suggest_int("max_leaves", 700, 1200),
             "min_child_weight": trial.suggest_float("min_child_weight", 1.0, 2.5),
             "gamma": trial.suggest_float("gamma", 0.0, 0.4),
             "reg_alpha": trial.suggest_float("reg_alpha", 0.7, 1.2),
             "reg_lambda": trial.suggest_float("reg_lambda", 3.0, 6.0),
             "subsample": trial.suggest_float("subsample", 0.6, 0.9),
             "colsample_bytree": trial.suggest_float("colsample_bytree", 0.7, 1.0),
             "learning_rate": trial.suggest_float("learning_rate", 0.02, 0.06),
             "scale_pos_weight": trial.suggest_float("scale_pos_weight", 0.8, 1.5),
             "n_estimators": 3000,
             "eval_metric": "logloss",
             "early_stopping_rounds": 50,
             }


    model = xgb.XGBClassifier(**params)

    model.fit(
              X_train, y_train,
              eval_set=[(X_valid, y_valid)],
              verbose=False
              )

    y_prob = model.predict_proba(X_valid)[:, 1]
    loss = log_loss(y_valid, y_prob)

    return loss

def run_optuna(config, X_train, y_train, X_valid, y_valid, n_trials=50):

    storage = config["optuna"]["storage"]
    study_name = config["optuna"]["study_name"]

    study = optuna.create_study(study_name=study_name,
                                storage=storage,
                                load_if_exists=True,
                                direction="minimize"
                                )
    
    study.optimize(lambda trial: objective(trial, X_train, y_train, X_valid, y_valid), 
                   n_trials=n_trials,
                   show_progress_bar=True,
                   n_jobs=4)

    print(f"Optuna best trial = ")
    print(study.best_trial.params)

    return study.best_trial.params



def read_samples (config):

    #get info from config file
    tot_df = pd.read_pickle(config.get('samples', 'file'))

    tot_df["Tag"] = tot_df["Tag"].map({"b jet": 1, "non b jet": 0})     #assign 1 - b jet; 0 non b jet

    X_np = tot_df.iloc[:, 1:-1].values.astype(np.float32)
    y_np = tot_df.iloc[:, -1].values.astype(np.float32)

    
    #first split - x train % x temporary
    X_train, X_temp, y_train, y_temp = train_test_split(X_np, y_np,
                                                        test_size=0.30,
                                                        random_state=11,
                                                        stratify=y_np)
    
    #second split -> test (50% - 15% of OG) vs validation (50% - 15% of OG)
    X_valid, X_test, y_valid, y_test = train_test_split(X_temp, y_temp,
                                                        test_size=0.50,
                                                        random_state=11,
                                                        stratify = y_temp)
    
    #check shape
    print(f"features shape = {X_np.shape}")
    print(f"tags shape = {y_np.shape}")


    var_names = list(tot_df.columns[1:-1])
    print(f"List of input features = {var_names}")

    print('Training dataset with ' + str(len(X_train)) + ' elements generated')
    print('Validation dataset with ' + str(len(X_valid)) + ' elements generated')
    print('Test dataset with ' + str(len(X_test)) + ' elements generated')

  
    return var_names, X_train, X_valid, X_test, y_train, y_valid, y_test


def bdt_training(config, X_train, y_train, X_valid, y_valid, out_path):
    '''
    model = xgb.XGBClassifier(
                            tree_method="hist",
                            grow_policy="lossguide",
                            max_depth=0,              
                            max_leaves=512,
                            min_child_weight=1,
                            gamma=0,
                            subsample=1,
                            colsample_bytree=1,
                            learning_rate=0.05,
                            n_estimators=3000,
                            reg_alpha=0.01,
                            reg_lambda=2.0,
                            scale_pos_weight=2,
                            eval_metric="logloss",
                            early_stopping_rounds=50
                        )

                              
    print(f"MODEL TYPE = {type(model)}")

    #train model
    model.fit(X_train, y_train,
              eval_set = [(X_valid, y_valid)],
              verbose = True)
    '''
    print("Running optuna for best parameters research")
    best_params = run_optuna(config, X_train, y_train, X_valid, y_valid)

    best_params.update({
                        "tree_method": "hist",
                        "grow_policy": "lossguide",
                        "max_depth": 0,
                        "n_estimators": 3000,
                        "eval_metric": "logloss",
                        "early_stopping_rounds": 50
                        })

    model = xgb.XGBClassifier(**best_params)

    model.fit(
              X_train, y_train,
              eval_set=[(X_valid, y_valid)],
              verbose=True
             )

    
    #save model
    model.save_model(os.path.join(out_path, "best_bdt.json"))
    print(f"BDT best model saved as best_bdt.json")

    #print model parameters
    print(f"PARAMETERS = {model.get_params()}")

    return model


def bdt_training_output(model, X_test, y_test, out_path, var_names):

    y_pred = model.predict(X_test)
    probs = model.predict_proba(X_test)

    y_prob_b = probs[:, 1]
    y_prob_nonb = probs[:, 0]
    '''
    #creazione dizionario di probabilità
    BDT_prob_dict = {
               "b jet": y_prob_b,
               "non b jet"  : y_prob_nonb,
               "label"        : y_test
               }
   
    #creazione del dataframe
    df = pd.DataFrame(BDT_prob_dict)
        
    with open(os.path.join(out_path,"BDT_scores.pkl"), "wb") as f:
        pickle.dump(df, f)
        
    print("---- file created successfully ----")
    '''
    # --------- Confusion matrix ---------

    cm = confusion_matrix(y_test, y_pred)
    sb.heatmap(cm, annot = True, fmt='.0f', cmap='Reds')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title('Confusion Matrix - XGBOOST')
    plt.savefig(os.path.join(out_path, "BDT_confusion_matriz.png"), dpi=300)
    plt.close()

    cm_n = confusion_matrix(y_test, y_pred, normalize='true')
    sb.heatmap(cm_n, annot = True, cmap = 'Reds')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title('Confusion matrix - XGBOOST')
    plt.savefig(os.path.join(out_path, "BDT_confusion_matrix_norm.png"), dpi=300)
    plt.close()

    # --------- BDT response score distribution ---------
    
    b_scores    = [p for (p, t) in zip(y_prob_b, y_test) if t == 1]
    nonb_scores = [p for (p, t) in zip(y_prob_b, y_test) if t == 0]

    plt.figure()
    sb.set(style = 'whitegrid')
    sb.histplot(b_scores, bins=40, color='royalblue', stat='probability', binrange=(0,1), label = 'b jet', alpha = 0.5)
    sb.histplot(nonb_scores, bins=40, color='firebrick', stat='probability', binrange=(0,1), label='non b jet', alpha=0.5)
    plt.title('BDT response score distribution')
    plt.xlabel('Score')
    plt.ylabel('Event fraction')
    plt.legend()
    plt.savefig(os.path.join(out_path, "BDT_response_score.png"), dpi=300)
    plt.close()

    # --------------------------- BDT evaluation metrics ---------------------------

    print("-------------- BDT evaluation metrics --------------")
    print(classification_report(y_test, y_pred))
    print("AUC ROC:", roc_auc_score(y_test, y_prob_b))

    print("------- VARIABLES ROC CURVES -------")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for i, name in enumerate(var_names):

        ax          = axes[i]
        scores      = X_test[:, i]
        fpr, tpr, _ = roc_curve(y_test, scores)
        roc_auc     = auc(fpr, tpr)

        ax.plot(fpr, tpr, label = f"AUC = {roc_auc:.2f}", color="royalblue")

        ax.set_title(f"ROC curve for {name}")
        ax.set_xlabel("False Positive rate")
        ax.set_ylabel("True Positive rate")
        ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, f"ROC_variables.png"), dpi=300)
    plt.close()


    print("------- BDT roc curve -------")
    fpr, tpr, _ = roc_curve(y_test, y_prob_b)
    auc_bdt = auc(fpr, tpr)

    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, label = f"AUC = {auc_bdt:.2f}", color="royalblue")
    plt.title('ROC curve - XGBoost BDT')
    plt.xlabel('False Positive rate')
    plt.ylabel('True Positive rate')
    plt.legend()
    plt.savefig(os.path.join(out_path, "ROC_BDT.png"), dpi=300)
    plt.close()


    print("------- FEATURES GAIN IMPORTANCE -------")
    imp = model.get_booster().get_score(importance_type = 'gain')

    mapped_imp = {var_names[int(k[1:])]: v
                 for k, v in imp.items()}

    imp_df = pd.DataFrame({
        "feature": list(mapped_imp.keys()),
        "importance": list(mapped_imp.values())
        }).sort_values("importance", ascending = False)
    plt.figure(figsize=(8,6))
    sb.barplot(data=imp_df, x="importance", y="feature", palette = "viridis")
    plt.title("Feature importance (gain)")
    plt.tight_layout()
    plt.savefig(os.path.join(out_path,"BDT_feature_importance.png"), dpi=300)
    plt.close()


    return


def plot_bdt_metrics(model, X_valid, y_valid, X_test, y_test, out_path):

    result_dict = model.evals_result()
    logloss_list = result_dict['validation_0']['logloss']
    n_trees = model.best_iteration + 1                      #best_iteration starts from 0
    x_axis = list(range(1, n_trees+1))

    plot_log_list = logloss_list[:(n_trees)]
    minimum_log = np.min(plot_log_list)

    # ---------------------- logloss plot ----------------------

    print(len(plot_log_list))
    fig = plt.subplot()
    plt.plot(x_axis, plot_log_list, label=f"minimum = {minimum_log:.3f}", color='blue')
    plt.title("Logloss plot")
    plt.xlabel("N trees")
    plt.ylabel("Logloss")
    plt.legend()
    plt.savefig(os.path.join(out_path, "bdt_logloss_plotting.png"), dpi=300)
    plt.close()

    
    # ---------------------- valid vs test metrics ----------------------
    

    booster = model.get_booster()

    #Validation set
    precision_eval = []
    accuracy_eval  = []
    f1_eval        = []
    recall_eval    = []

    dvalid = xgb.DMatrix(X_valid)

    for i in x_axis:
        y_pred = booster.predict(dvalid, iteration_range=(0, i))
        y_pred = (y_pred > 0.5).astype(float)

        precision_eval.append(precision_score(y_valid, y_pred, zero_division=0))
        accuracy_eval.append(accuracy_score(y_valid, y_pred))
        f1_eval.append(f1_score(y_valid, y_pred))
        recall_eval.append(recall_score(y_valid, y_pred))

    #Test set
    precision_test = []
    accuracy_test  = []
    f1_test        = []
    recall_test    = []

    dtest = xgb.DMatrix(X_test)

    for i in x_axis:
        y_pred = booster.predict(dtest, iteration_range=(0, i))
        y_pred = (y_pred > 0.5).astype(float)

        precision_test.append(precision_score(y_test, y_pred, zero_division=0))
        accuracy_test.append(accuracy_score(y_test, y_pred))
        f1_test.append(f1_score(y_test, y_pred))
        recall_test.append(recall_score(y_test, y_pred))

    #plots
    metrics = [
        ("precision", precision_eval, precision_test),
        ("accuracy", accuracy_eval, accuracy_test),
        ("F1", f1_eval, f1_test),
        ("recall", recall_eval, recall_test),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for ax, (name, val_var, test_var) in zip(axes, metrics):
        ax.plot(x_axis, val_var, label="validation", color="blue")
        ax.plot(x_axis, test_var, label="testing", color="red")

        ax.set_title(f"{name} vs trees")
        ax.set_xlabel("N trees")
        ax.set_ylabel(name)
        ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(out_path, "bdt_metrics_plotting.png"), dpi=300)
    plt.close()


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


def main () :

    if (len (sys.argv) < 2):
        print ('usage:', sys.argv[0], 'input.cfg')
        sys.exit (1)

    # reading the config file
    # -----------------------

    config = configparser.ConfigParser ()
    config.read ('config.cfg')

    out_path  = config.get ('saving', 'out_path', fallback = None)
    os.makedirs(out_path, exist_ok=True)
    print("out_path:", out_path)
    
    print("Sezioni trovate:", config.sections())
    print("Contenuto saving:", dict(config['saving']) if 'saving' in config else "Non trovata")
    os.makedirs (out_path, exist_ok = True)     

    with open (out_path + 'config.cfg', 'w') as configfile :
        config.write (configfile)


    print('\n --- PREPARING SAMPLES ---\n')
    # ---------------------------------------

    # ------------------------  BDT ------------------------

    print("---- IMPLEMENTING BDT ----")

    var_names, X_train, X_valid, X_test, y_train, y_valid, y_test = read_samples(config)

    model = bdt_training(config, X_train, y_train, X_valid, y_valid, out_path)

    bdt_training_output(model, X_test, y_test, out_path, var_names)

    plot_bdt_metrics(model, X_valid, y_valid, X_test, y_test, out_path)
    
    
    # ------------------------ END of BDT impl ------------------------


   
if __name__ == "__main__" : 
    main ()
    




