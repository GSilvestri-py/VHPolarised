'''
plots of:
- 4 b tagging variables (true vs selected)
- ptj, m_jj
'''

from concurrent.futures import ProcessPoolExecutor          
import glob
import multiprocessing
import numpy as np
import os
import ROOT
from tqdm import tqdm
import vector
import xgboost as xgb

# ------------------------------------------------------------------------------------------------

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

    def vec4(self):
        '''Jet four vector from pt eta phi mass'''
        px = self.pt * np.cos(self.phi)
        py = self.pt * np.sin(self.phi)
        pz = self.pt * np.sinh(self.eta)
        p  = self.pt * np.cosh(self.eta)
        E  = np.sqrt(self.mass**2 + p**2)
        return vector.obj(px = px, py = py, pz =pz, E = E)
    
    def inv_mass(self, other):
        '''return invariant mass of self + other 4 vector'''
        return (self.vec4() + other.vec4()).M

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

class Higgs:

    def __init__(self, pt, eta, phi, mass):
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.mass = mass
        self.children = []

    def add_child(self, pdgId):
        return self.children.append(pdgId)

    def has_bb(self):
        return sum(1 for id in self.children if abs(id) == 5) >= 2

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------  

class Event:

    def __init__ (self, tree, entry):

        tree.GetEntry(entry)
        self.tree = tree
        self.jets = self._load_jets()

    
    def _find_higgs(self):

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

    def h_bb(self):

        higgs = self._find_higgs()
        if higgs is None:   return False
        return higgs.has_bb()


    def _load_jets(self):
        
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
        '''identifies true bjets in event'''
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

    def __init__(self, filename, model):

        self.file = ROOT.TFile.Open(filename)
        self.tree = self.file.Get("Events")
        self.model = model

        self.variables = {
                         "DeepFlavB": [], "DeepFlavCvB": [], "DeepFlavCvL": [],
                         "DeepFlavQG": [], "Jet_pt": [], "Jet_eta": [], "m_jj": []
                         }
        
        self.variables_true = {
                              "DeepFlavB": [], "DeepFlavCvB": [], "DeepFlavCvL": [],
                              "DeepFlavQG": [], "Jet_pt": [], "Jet_eta": [], "m_jj": []
                              }

    def run_analyzer(self):
        '''analyze the event file --> return variables & true variables dictionaries'''
        for i in range(self.tree.GetEntries()):
            event = Event(self.tree, i)

            is_hbb_event = event.h_bb()

            if not is_hbb_event:    continue

            # -------- true values --------

            true_mass = event.true_inv_mass()
            if true_mass is None:   continue    #skip to next event: no jj matches

            true_idxs = event.true_bjet_indices()

            for i in true_idxs:
                jet = event.jets[i]

                self.variables_true["DeepFlavB"].append(jet.flavb)
                self.variables_true["DeepFlavCvB"].append(jet.flavcvb)
                self.variables_true["DeepFlavCvL"].append(jet.flavcvl)
                self.variables_true["DeepFlavQG"].append(jet.flavqg)
                self.variables_true["Jet_pt"].append(jet.pt)
                self.variables_true["Jet_eta"].append(jet.eta)
            
            self.variables_true["m_jj"].append(true_mass)


            # -------- bdt selected --------

            scores = []
            for i, j in enumerate(event.jets):
                X = [[j.flavb, j.flavcvb, j.flavcvl, j.flavqg, j.pt, j.eta]]
                prob = self.model.predict_proba(X)[0, 1]    #first row, second element - class 1
                scores.append((i, prob))

            scores.sort(key=lambda x:x[1], reverse=True)
            j1, j2 = scores[0][0], scores[1][0]             #first element of first and second row

            jet1, jet2 = event.jets[j1], event.jets[j2]
            mass = jet1.inv_mass(jet2)
            self.variables["m_jj"].append(mass)

            for jet in [jet1, jet2]:
                self.variables["DeepFlavB"].append(jet.flavb)
                self.variables["DeepFlavCvB"].append(jet.flavcvb)
                self.variables["DeepFlavCvL"].append(jet.flavcvl)
                self.variables["DeepFlavQG"].append(jet.flavqg)
                self.variables["Jet_pt"].append(jet.pt)
                self.variables["Jet_eta"].append(jet.eta)

        return self.variables, self.variables_true

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

def plot_hist(variables, variables_true):

    for n, l in variables.items():
        print(f"\n {n} with length = {len(l)}")
    
    c = ROOT.TCanvas("c", "Variables after BDT", 1800, 1200)
    c.Divide(3, 2)

    hists = []  

    for i, name in enumerate(variables.keys(), start=1):

        val = variables[name]
        true_val = variables_true[name]

        if len(val) == 0 or len(true_val) == 0:
            print("empty list in dictionary")
            break  

        x_min = min(min(val), min(true_val))
        x_max = max(max(val), max(true_val))

        if x_max > 300:

            val      = [x for x in val if x < 300]
            true_val = [x for x in true_val if x < 300]

        h1 = ROOT.TH1F(f"{name}_bdt",  f"{name}; value; Normalized entries", 40, x_min, x_max)
        h2 = ROOT.TH1F(f"{name}_true", f"{name}; value; Normalized entries", 40, x_min, x_max)
        for v in val:
            h1.Fill(v)
        if h1.Integral() > 0:
            h1.Scale(1.0 / h1.Integral())
        for v in true_val:
            h2.Fill(v)
        if h2.Integral() > 0:
            h2.Scale(1.0 / h2.Integral())

        max_v = h1.GetMaximum()
        max_truev = h2.GetMaximum()

        ymax = max(max_v, max_truev) * 1.10
        h1.SetMaximum(ymax)
        h2.SetMaximum(ymax)

        h1.SetLineColor(ROOT.kBlue)
        h1.SetFillColorAlpha(ROOT.kBlue, 0.2)
        h2.SetLineColor(ROOT.kRed)
        h2.SetFillColorAlpha(ROOT.kRed, 0.10)
        c.cd(i)

        h1.Draw("HIST")
        h2.Draw("HIST SAME")

        ROOT.gStyle.SetOptStat(0)
        leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
        leg.SetBorderSize(1)
        leg.SetFillStyle(0)      
        leg.SetTextSize(0.035)
        leg.SetTextFont(42)
        leg.SetMargin(0.15)

        entry_bdt  = f"BDT jets"
        entry_true = f"True jets"
        leg.AddEntry(h1, entry_bdt,  "l")
        leg.AddEntry(h2, entry_true, "l")
        leg.Draw()
        
        hists.append(h1)  
        hists.append(h2)  

    c.SaveAs("allfiles_variables_plotting_after_bdt.root")
    c.SaveAs("allfiles_variables_plotting_after_bdt.png")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

def process_file(filename):
    #print(f"[Worker] Processing {filename}", flush=True)
    global MODEL
    analyzer = Analyzer(filename, MODEL)
    analyzer.run_analyzer()
    return analyzer.variables, analyzer.variables_true
    '''
    global DNN_MODEL
    
    '''

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

def merge_results(results):
    """Merge the list of (variables, variables_true) dictionaries 
       returned by process_file into two global dictionaries."""
    
    merged_vars = {k: [] for k in results[0][0].keys()}
    merged_true = {k: [] for k in results[0][1].keys()}

    for vars_i, true_i in results:
        for k in merged_vars:
            merged_vars[k].extend(vars_i[k])
            merged_true[k].extend(true_i[k])

    return merged_vars, merged_true

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

def init_worker(model_path):
    global MODEL
    MODEL = xgb.XGBClassifier(n_jobs=1)
    MODEL.load_model(model_path)

    '''
    global DNN_MODEL
    DNN_MODEL = recreate_model(config=configfile)
    '''
    
# ------------------------------------------------------------------------------------------------

def main():
    
    folder = "nanofiles_folder"
    files = sorted(glob.glob(os.path.join(folder, "*.root")))    
    
    model_path = "best_model_path"
    model = xgb.XGBClassifier()
    model.load_model(model_path)


    nblocks = 8

    with ProcessPoolExecutor(max_workers=nblocks,
                             mp_context=multiprocessing.get_context("spawn"),
                             initializer=init_worker,
                             initargs=(model_path,)
                             ) as pool:
        results = list(
                        tqdm(
                            pool.map(process_file, files),
                            total=len(files),
                            desc="Processing files"
                            )
                        )
        
    variables, variables_true = merge_results(results)
    plot_hist(variables, variables_true)

if __name__ == "__main__":
    main()