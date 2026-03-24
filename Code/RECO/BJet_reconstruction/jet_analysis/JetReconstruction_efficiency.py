# READ SINGLE NANO FILE
import configparser
import importlib
import os
import ROOT
import torch
import pandas as pd
from tabulate import tabulate
import xgboost as xgb

import sys
sys.path.append("module_dir")


# ---------------------------------------------------------------------------------------------

class Jet:

    def __init__(self, pt, eta, phi, mass, flavb, flavcvb, flavcvl, flavqg):
        '''initialisation of element attributes'''
        self.pt      = pt
        self.eta     = eta
        self.phi     = phi
        self.mass    = mass
        self.flavb   = flavb 
        self.flavcvb  = flavcvb
        self.flavcvl = flavcvl
        self.flavqg  = flavqg

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- 

class Higgs:

    def __init__(self, pt, eta, phi, mass):
        ''' initialise higgs class'''
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.mass = mass
        self.children = []

    def add_child(self, pdgId):
        '''add a single id to self.children'''
        return self.children.append(pdgId)

    def has_bb(self):
        '''return True if at least two b index in self.children'''
        return sum(1 for id in self.children if abs(id) == 5) >= 2
    
# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

class Event:

    def __init__ (self, tree, entry):
        '''initialise Event class'''
        tree.GetEntry(entry)
        self.tree = tree
        self.jets = self._load_jets()

    
    def _find_higgs(self):
        '''find higgs in event and return Higgs class'''
        h_idx = None
        for j in range(int(self.tree.nGenPart)):
            if int(self.tree.GenPart_pdgId[j]) == 25:
                child_indices = [k for k in range(int(self.tree.nGenPart)) if int(self.tree.GenPart_genPartIdxMother[k]) == j]
                h_idx = j
                

        if h_idx is None:
            print(f"No higgs in event")
            return None
        
        higgs = Higgs(
                      pt = self.tree.GenPart_pt[h_idx],
                      eta = self.tree.GenPart_eta[h_idx],
                      phi=self.tree.GenPart_phi[h_idx],
                      mass=self.tree.GenPart_mass[h_idx]
                    )

        for k in child_indices:
            pdg = int(self.tree.GenPart_pdgId[k])
            higgs.add_child(pdg)

        return higgs


    def is_h_bb_event(self):
        '''return True if h > bb event'''
        higgs = self._find_higgs()
        if higgs is None:   return False
        return higgs.has_bb()


    def _load_jets(self):
        '''collect jets in the event'''
        jets = []
        for i in range(self.tree.nJet):
            jets.append(
                        Jet(
                            pt = self.tree.Jet_pt[i],
                            eta = self.tree.Jet_eta[i],
                            phi = self.tree.Jet_phi[i],
                            mass=self.tree.Jet_mass[i],
                            flavb=self.tree.Jet_btagDeepFlavB[i],
                            flavcvb=self.tree.Jet_btagDeepFlavCvB[i],
                            flavcvl=self.tree.Jet_btagDeepFlavCvL[i],
                            flavqg=self.tree.Jet_btagDeepFlavQG[i],
                        )
                        )
        return jets

    def true_bjet_indices(self):
        '''identifies true bjets indexes in event'''
        idx = []
        for i in range(self.tree.nJet):
            if abs(int(self.tree.Jet_partonFlavour[i])) == 5:
                idx.append(i)
        return idx

    def true_inv_mass(self):
        '''compute invariant mass between two jest class'''
        idx = self.true_bjet_indices()
        if len(idx) != 2:
            return None
        return self.jets[idx[0]].inv_mass(self.jets[idx[1]])


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

class Analyzer:

    def __init__(self, filename, model, dnn_model):
        '''initialize analyzer of event class'''
        self.file = ROOT.TFile.Open(filename)
        self.tree = self.file.Get("Events")
        self.model = model
        self.dnn_model = dnn_model

        self.counters = {
                            "DeepFlavB":   [0, 0, 0],
                            "DeepFlavCvL": [0, 0, 0],
                            "DeepFlavCvB": [0, 0, 0],
                            "JetPt":       [0, 0, 0],
                            "BDT":         [0, 0, 0],
                            "DNN":         [0, 0, 0]
                            }

        self.total_hbb_events = 0

    @staticmethod
    def top_2_index_order_list(values_list, reverse=True):
        ordered_list = sorted(values_list, key=lambda x: x[1], reverse=reverse)
        return ordered_list[:2]

    @staticmethod
    def counter_estimate(true_set, idx_top2, counters):
        pred_set = set(idx_top2)
        nmatches = len(true_set & pred_set)

        if nmatches == 2:
            counters[2] += 1
        elif nmatches == 1:
            counters[1] += 1
        else:
            counters[0] += 1

        return counters


    def run_analyzer(self):
        '''identify two most probable b jet in event and compare with true idxs -> return nested lists'''
        for i in range(self.tree.GetEntries()):

            event = Event(self.tree, i)

            is_hbb_event = event.is_h_bb_event()

            if not is_hbb_event:    continue

            true_set = set(event.true_bjet_indices())
            if len(true_set) != 2:
                continue

            FlavB    = []
            FlavCvB  = []
            FlavCvL  = []
            FlavQG   = []
            jet_pt   = []
            jet_eta  = []
            bdt_list = []
            dnn_list = []

            self.total_hbb_events += 1


            for (i, j) in enumerate(event.jets):
                FlavB.append((i, j.flavb))
                FlavCvB.append((i, j.flavcvb))
                FlavCvL.append((i, j.flavcvl))
                FlavQG.append((i, j.flavqg))
                jet_pt.append((i, j.pt))
                jet_eta.append((i, j.eta))

                my_X = [
                    FlavB[i][1],
                    FlavCvB[i][1],
                    FlavCvL[i][1],
                    FlavQG[i][1],
                    jet_pt[i][1],
                    jet_eta[i][1]
                    ]

                # bdt event input
                prob = float(self.model.predict_proba([my_X])[0, 1])
                bdt_list.append((i, prob))

                # dnn event input
                X_event = torch.tensor(my_X, dtype = torch.float32).unsqueeze(0)
                self.dnn_model.eval()
                with torch.no_grad():
                    score = self.dnn_model(X_event)
                    probs = torch.softmax(score, dim=1).numpy()
                
                dnn_list.append((i, float(probs[0][1])))


            


            idx_top2_B        = [i for (i, _) in self.top_2_index_order_list(FlavB)]
            idx_bottom2_CvB   = [i for (i, _) in self.top_2_index_order_list(FlavCvB, reverse=False)]
            idx_top2_CvL      = [i for (i, _) in self.top_2_index_order_list(FlavCvL)]
            idx_top2_pt       = [i for (i, _) in self.top_2_index_order_list(jet_pt)]
            idx_top2_bdt      = [i for (i, _) in self.top_2_index_order_list(bdt_list)]
            idx_top2_dnn      = [i for (i, _) in self.top_2_index_order_list(dnn_list)]

              
            # ---------------------------------------------------- update counters ----------------------------------------------------

            self.counters["DeepFlavB"]   = self.counter_estimate(true_set, idx_top2_B,   self.counters["DeepFlavB"])
            self.counters["DeepFlavCvL"] = self.counter_estimate(true_set, idx_top2_CvL, self.counters["DeepFlavCvL"])
            self.counters["DeepFlavCvB"] = self.counter_estimate(true_set, idx_bottom2_CvB, self.counters["DeepFlavCvB"])
            self.counters["JetPt"]       = self.counter_estimate(true_set, idx_top2_pt,  self.counters["JetPt"])
            self.counters["BDT"]         = self.counter_estimate(true_set, idx_top2_bdt, self.counters["BDT"])
            self.counters["DNN"]         = self.counter_estimate(true_set, idx_top2_dnn, self.counters["DNN"])

        results_data = [
                        ["DeepFlavB"]   + self.counters["DeepFlavB"],
                        ["DeepFlavCvL"] + self.counters["DeepFlavCvL"],
                        ["DeepFlavCvB"] + self.counters["DeepFlavCvB"],
                        ["JetPt"]       + self.counters["JetPt"],
                        ["BDT"]         + self.counters["BDT"],
                        ["DNN"]         + self.counters["DNN"],
                        ]

        eff_data = [
                    ["DeepFlavB",   100*self.counters["DeepFlavB"][2]/self.total_hbb_events,
                                    100*self.counters["DeepFlavB"][1]/self.total_hbb_events],
                    ["DeepFlavCvL", 100*self.counters["DeepFlavCvL"][2]/self.total_hbb_events,
                                    100*self.counters["DeepFlavCvL"][1]/self.total_hbb_events],
                    ["DeepFlavCvB", 100*self.counters["DeepFlavCvB"][2]/self.total_hbb_events,
                                    100*self.counters["DeepFlavCvB"][1]/self.total_hbb_events],
                    ["JetPt",       100*self.counters["JetPt"][2]/self.total_hbb_events,
                                    100*self.counters["JetPt"][1]/self.total_hbb_events],
                    ["BDT",         100*self.counters["BDT"][2]/self.total_hbb_events,
                                    100*self.counters["BDT"][1]/self.total_hbb_events],
                    ["DNN",         100*self.counters["DNN"][2]/self.total_hbb_events,
                                    100*self.counters["DNN"][1]/self.total_hbb_events],
                    ]

        print(f"total studied events = {self.total_hbb_events}")
        
        return results_data, eff_data

# -------------------------------------------------------------------------------------

def recreate_model(config):

    #training_output_path = config.get('saving', 'out_path')
    #out_file_base        = config.get('saving', 'out_file_base')

    # import model
    library_name = config.get('dnn_model', 'library_name')
    model_name   = config.get('dnn_model', 'model_name')

    model_lib    = importlib.import_module(library_name)
    model_class  = getattr(model_lib, model_name)

    '''
    # --- input_dim: deve essere identico al training ---
    X_test = torch.tensor(test_df.iloc[:, :-1].values, dtype=torch.float32)
    input_dim = X_test.shape[1]

    print(test_df)
    print(X_test)
    '''

    input_dim    = 6
    hidden_units = input_dim
    hidden_layers = config.getint('training', 'hidden_layers', fallback=3)
    out_dim      = config.getint('dnn_model', 'classes_num')
    dropout_prob = config.getfloat('training', 'dropout_prob')
    layer_norm   = config.getboolean('training', 'layer_norm')

    model = model_class(input_dim, hidden_units, hidden_layers, out_dim, dropout_prob, layer_norm)

    # load weights
    best_model_path = config.get('input_files', 'dnn_model_path')
    print(f" importing model saved in {best_model_path}")
    state_dict = torch.load(best_model_path, map_location="cpu")
    model.load_state_dict(state_dict)

    print("MODEL WEIGHT SAMPLE:", model.model[0].weight[0][:5])
    print("CHECKPOINT WEIGHT SAMPLE:", state_dict['model.0.weight'][0][:5])


    return model


# -----------------------------------------------------------------------------------------------


def main():

    filename = "nanofile_path"

    # -- config file ---
    config_file = configparser.ConfigParser()
    config_file.read("config_efficiency.cfg")

    # bdt
    bdt_model_path = config_file.get('input_files', 'bdt_model_path')
    bdt_model = xgb.XGBClassifier()
    bdt_model.load_model(bdt_model_path)

    # dnn - partial model -> set total model in config if existing
    dnn_model      = recreate_model(config_file)

    analyzer = Analyzer(filename, bdt_model, dnn_model)
    results_data, eff_data = analyzer.run_analyzer()

    # ------------------------------ print df ------------------------------

    df_results = pd.DataFrame(results_data,
                              columns=["Method", "0 matches", "1 match", "2 matches"]
                             )

    print("\n" + "="*58)
    print("{:^60}".format("b‑jet Identification Performance"))
    print("="*58 + "\n")

    print(tabulate(df_results.values.tolist(),
                headers=df_results.columns.tolist(),
                tablefmt="grid"))


    df_eff = pd.DataFrame(eff_data, 
                          columns=["Method", "Efficiency (2 matches) [%]", "Partial (1 match) [%]"]
                          )

    print("\n" + "="*58)
    print("{:^60}".format("b‑jet Tagging Efficiencies"))
    print("="*58 + "\n")

    print(tabulate(df_eff.values.tolist(), headers=df_eff.columns.tolist(), tablefmt="grid"))


if __name__ == "__main__":
    main()