import awkward as ak
from concurrent.futures import ProcessPoolExecutor          #optimise file-loop
import glob
import numpy as np
import os
import pandas as pd
import ROOT
import sys
from tqdm import tqdm
import uproot                                               #avoid event-loop
import vector


def pdg_to_Name(pdgId):
    pdgNames = {
        11:"e-",  -11:"e+",    13:"mu-",  -13:"mu+",   15:"tau-",  -15:"tau+",
        12:"ve",  -12:"~ve",   14:"vmu",  -14:"~vmu",  16:"vtau",  -16:"~vtau",
        21:"g",    22:"γ",     23:"Z",     24:"W+",    -24:"W-",    25:"H",
        6:"t",    -6:"tbar",   5:"b",     -5:"bbar",   1:"d",      -1:"dbar",
        2:"u",    -2:"ubar",   3:"s",     -3:"sbar",   4:"c",      -4:"cbar"
    }
    return pdgNames.get(pdgId, f"PDG({pdgId})")


def build_4vec(pt, eta, phi, mass):
    return vector.obj(pt=pt, eta=eta, phi=phi, mass=mass)

def build_4vec_from_pt_eta_phi_m(pt, eta, phi, mass):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    p  = pt * np.cosh(eta)
    E  = np.sqrt(mass**2 + p**2)
    return vector.obj(px = px, py = py, pz =pz, E = E)

def expand_columns_vectors(df):
    vector_cols = []
    for col in df.columns:

        if df[col].apply(lambda v: hasattr(v, "pt")).all():
            vector_cols.append(col)
            df[f"{col}_pt"]   = df[col].apply(lambda v: v.pt)
            df[f"{col}_eta"]  = df[col].apply(lambda v: v.eta)
            df[f"{col}_phi"]  = df[col].apply(lambda v: v.phi)
            df[f"{col}_mass"] = df[col].apply(lambda v: v.mass)

    return df.drop(columns=vector_cols)


def is_b_hadron(pdgId):                     #return True if b hadron in event -> if b in event
    b_hadrons = [
                    511,   # B0
                    521,   # B+
                    531,   # Bs0
                    541,   # Bc+
                    5122,  # Lambda_b0
                    5132,  # Xi_b-
                    5232,  # Xi_b0
                    5332   # Omega_b-
                    ]
    return abs(pdgId) in b_hadrons


def mother_pdg(j, midx_evt, pdg_evt):
    m = int(midx_evt[j])
    return int(pdg_evt[m]) if m >= 0 else None


def get_VB_daughters(pdg_evt, midx_evt, pt_evt, eta_evt, phi_evt, mass_evt, pt_MET_evt, phi_MET_evt):

    daughters_VB = []
    lep_vec = antilep_vec = build_4vec(0,0,0,0)
    MET = build_4vec(0,0,0,0)

    event_type = "Z/W"

    for j in range(len(pdg_evt)):   #loop over Gen Particles

        pid = int(pdg_evt[j])
        mother = mother_pdg(j, midx_evt, pdg_evt)

        if pid == 23 and mother != 25:  #avoid VB from h
            event_type = "Z"
            child_indices = [k for k in range(len(pdg_evt)) if int(midx_evt[k]) == j]           #get daughters' indices
            lep_indices = [k for k in child_indices if abs(int(pdg_evt[k])) in [11, 13]]        #filter lepton indices

            if len(lep_indices) == 2:
                for k in lep_indices:
                    pid_k = int(pdg_evt[k])
                    daughters_VB.append(pdg_to_Name(pid_k))
                    vec = build_4vec(float(pt_evt[k]), float(eta_evt[k]), float(phi_evt[k]), float(mass_evt[k]))

                    if pid_k > 0: 
                        lep_vec = vec
                    else: 
                        antilep_vec = vec

        if abs(pid) == 24 and mother != 25:
            event_type = "W+/-"
            child_indices = [k for k in range(len(pdg_evt)) if int(midx_evt[k]) == j]               #get daughters' indices
            lep_indices = [k for k in child_indices if abs(int(pdg_evt[k])) in [11, 12, 13, 14]]    #filter lepton indices
            MET = build_4vec(pt_MET_evt, 0, phi_MET_evt, 0)

            if len(lep_indices) == 2:
                for k in lep_indices:
                    pid_k = int(pdg_evt[k])
                    daughters_VB.append(pdg_to_Name(pid_k))
                    vec = build_4vec(float(pt_evt[k]), float(eta_evt[k]), float(phi_evt[k]), float(mass_evt[k]))

                    if pid_k > 0: 
                        lep_vec = vec
                    else: 
                        antilep_vec = vec

    return lep_vec, antilep_vec, MET, event_type, daughters_VB, lep_indices



def clean_jets(lep_indices, pdg_evt, pt_evt, eta_evt, phi_evt, mass_evt,
                n_jets_evt, GJet_pt_evt, GJet_eta_evt, GJet_phi_evt, GJet_mass_evt):

    cleaned_list = []   #list with b jet candidate index

    for k in lep_indices:   
        pid = int(pdg_evt[k])   #get particle ID
        if abs(pid) in [11, 13]:
            lep_vec = build_4vec_from_pt_eta_phi_m(pt_evt[k], eta_evt[k], phi_evt[k], mass_evt[k])

            for jet_idx in range(n_jets_evt):
                jet_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[jet_idx], GJet_eta_evt[jet_idx], GJet_phi_evt[jet_idx], GJet_mass_evt[jet_idx])

                deltaR_lj = lep_vec.deltaR(jet_vec)
                
                if deltaR_lj > 0.1:
                    cleaned_list.append(jet_idx)  #save index if candidate b jet

    return cleaned_list

def get_b_jets(current_file, index_list, GJet_pt_evt, GJet_eta_evt, GJet_phi_evt, GJet_mass_evt):

    kin_info = []
    mass_jj_list = []
    pt_jj_dict = {}

    for jet_idx in index_list:
        j1_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[jet_idx], GJet_eta_evt[jet_idx], GJet_phi_evt[jet_idx], GJet_mass_evt[jet_idx])

        for jet2_idx in index_list:
            if jet2_idx <= jet_idx:    continue

            j2_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[jet2_idx], GJet_eta_evt[jet2_idx], GJet_phi_evt[jet2_idx], GJet_mass_evt[jet2_idx])

            vec_jj = j1_vec + j2_vec
            mass_jj = vec_jj.M
            pt_jj = GJet_pt_evt[jet_idx] + GJet_pt_evt[jet2_idx]
            pt_jj_dict[pt_jj] = (jet_idx, jet2_idx, mass_jj)

            tollerance = 20
            if abs(mass_jj - 125) < tollerance:
                kin_info.append(build_4vec(GJet_pt_evt[jet_idx], GJet_eta_evt[jet_idx], GJet_phi_evt[jet_idx], GJet_mass_evt[jet_idx]))
                kin_info.append(build_4vec(GJet_pt_evt[jet2_idx], GJet_eta_evt[jet2_idx], GJet_phi_evt[jet2_idx], GJet_mass_evt[jet2_idx]))
                mass_jj_list.append((jet_idx, jet2_idx, mass_jj))
    if not pt_jj_dict:
        print(f"Nessuna coppia trovata per evento in file {current_file}")
        return []

    if len(kin_info) == 0:  #no jets found but b in event -> get pair with higher pt
        max_value = max(list(pt_jj_dict.keys()))
        jet1_index, jet2_index, m = pt_jj_dict[max_value]
        kin_info.append(build_4vec(GJet_pt_evt[jet1_index], GJet_eta_evt[jet1_index], GJet_phi_evt[jet1_index], GJet_mass_evt[jet1_index]))
        kin_info.append(build_4vec(GJet_pt_evt[jet2_index], GJet_eta_evt[jet2_index], GJet_phi_evt[jet2_index], GJet_mass_evt[jet2_index]))
                
    if len(kin_info) > 2:   #multiple jet pairs found -> get pair with higher pt
        ptjj_selected = [GJet_pt_evt[j1] + GJet_pt_evt[j1] for (j1, j2, m) in mass_jj_list]
        max_idx = get_max_index(ptjj_selected)
        kin_info = kin_info[max_idx*2:(max_idx*2+2)]
    return kin_info



def get_max_index(list):
    max_v = max(list)
    return list.index(max_v)



def read_files(file_list, output_root):
    event_number = 0
    rows = []

    for full_path in tqdm(file_list, desc="Processing files", unit="file"):
        filename = os.path.basename(full_path)

        with uproot.open(full_path) as f:                                       
            if "Events" not in f:
                continue
            tree = f["Events"]

            Gen_pdg   = tree["GenPart_pdgId"].array(library="ak")
            Gen_midx  = tree["GenPart_genPartIdxMother"].array(library="ak")
            Gen_pt    = tree["GenPart_pt"].array(library="ak")
            Gen_eta   = tree["GenPart_eta"].array(library="ak")
            Gen_phi   = tree["GenPart_phi"].array(library="ak")
            Gen_mass  = tree["GenPart_mass"].array(library="ak")

            n_jets      = tree["nGenJet"].array(library="ak")
            GJet_pt     = tree["GenJet_pt"].array(library="ak")
            GJet_eta    = tree["GenJet_eta"].array(library="ak")
            GJet_phi    = tree["GenJet_phi"].array(library="ak")
            GJet_mass   = tree["GenJet_mass"].array(library="ak")
            GJet_h_flav = tree["GenJet_hadronFlavour"].array(library="ak")

            n_jets_vistau = tree["nGenVisTau"].array(library="ak")
            VisTau_pt     = tree["GenVisTau_pt"].array(library="ak")
            VisTau_eta    = tree["GenVisTau_eta"].array(library="ak")
            VisTau_phi    = tree["GenVisTau_phi"].array(library="ak")
            VisTau_mass   = tree["GenVisTau_mass"].array(library="ak")
            VisTau_mother = tree["GenVisTau_genPartIdxMother"].array(library="ak")
            VisTau_status = tree["GenVisTau_status"].array(library="ak") 

            pt_MET  = tree["MET_pt"].array(library="ak")
            phi_MET = tree["MET_phi"].array(library="ak")

            n_events = len(Gen_pdg)

            for i in range(n_events):
                event_number += 1
                lep_vec = antilep_vec = build_4vec(0,0,0,0)
                d1_vec  = d2_vec      = build_4vec(0,0,0,0)
                MET = build_4vec(0,0,0,0)
                decay      = "No Higgs"
                event_type = "Z/W"
                VB_decay   = ""

                daughters_VB = []
                pdg_evt   = Gen_pdg[i]
                midx_evt  = Gen_midx[i]
                pt_evt    = Gen_pt[i]
                eta_evt   = Gen_eta[i]
                phi_evt   = Gen_phi[i]
                mass_evt  = Gen_mass[i]

                #GenJet array i component
                n_jets_evt      = n_jets[i]
                GJet_pt_evt     = GJet_pt[i]
                GJet_eta_evt    = GJet_eta[i]
                GJet_phi_evt    = GJet_phi[i]
                GJet_mass_evt   = GJet_mass[i]
                GJet_h_flav_evt = GJet_h_flav[i]

                #GenVisTau array i component
                n_jets_vistau_evt = n_jets_vistau[i]
                VisTau_pt_evt     = VisTau_pt[i]
                VisTau_eta_evt    = VisTau_eta[i]
                VisTau_phi_evt    = VisTau_phi[i]
                VisTau_mass_evt   = VisTau_mass[i]
                VisTau_mother_evt = VisTau_mother[i]
                VisTau_status_evt = VisTau_status[i]

                #MET i component
                pt_MET_evt  = pt_MET[i]
                phi_MET_evt = phi_MET[i]
                

                is_b_inevent = False
                for j in range(len(pdg_evt)):
                    if is_b_hadron(int(pdg_evt[j])):
                        is_b_inevent = True

                #-------------------------------------------------- VB decay - GenPart --------------------------------------------------

                lep_vec, antilep_vec, MET, event_type, daughters_VB, lep_indices = get_VB_daughters(pdg_evt, midx_evt, pt_evt, eta_evt, phi_evt, mass_evt, pt_MET_evt, phi_MET_evt)
                
                if len(daughters_VB) == 0:

                    if event_type=="Z":
                        VB_decay = "v v~ / tau+ tau-"  
                    if event_type=="W+/-":
                        VB_decay = "tau v"

                else:
                    VB_decay = " + ".join(daughters_VB)

                #--------------------------------------------------- h decay - GenPart ---------------------------------------------------

                daughters, kin = [], []

                if n_jets_vistau_evt == 2:
                    daughters.append("tau+")
                    daughters.append("tau-")
                    kin.append(build_4vec(VisTau_pt_evt[0], VisTau_eta_evt[0], VisTau_phi_evt[0], VisTau_mass_evt[0]))
                    kin.append(build_4vec(VisTau_pt_evt[1], VisTau_pt_evt[1], VisTau_phi_evt[1], VisTau_mass_evt[1]))
                    
                elif is_b_inevent:
                    daughters.append("b")
                    daughters.append("bbar")

                    if n_jets_evt > 2:  #require at least two jets
                        clean_list = clean_jets(lep_indices, pdg_evt, pt_evt, eta_evt, phi_evt, mass_evt,
                                            n_jets_evt, GJet_pt_evt, GJet_eta_evt, GJet_phi_evt, GJet_mass_evt)

                        if len(clean_list) < 2:     #less than two jets remaining after clean up -> assume chance overlap
                            clean_list = np.arange(0, n_jets_evt)

                        kin = get_b_jets(filename, clean_list, GJet_pt_evt, GJet_eta_evt, GJet_phi_evt, GJet_mass_evt)


                    elif n_jets_evt == 2:   #avoid charged leptonic jets if only two jets in event
                        clean_list = np.arange(0, n_jets_evt)
                        kin = get_b_jets(filename, clean_list, GJet_pt_evt, GJet_eta_evt, GJet_phi_evt, GJet_mass_evt)

                    elif n_jets_evt < 2:     #only one jet or no jets in event
                        clean_list, kin = [], []    #no kin info found for this event
                        daughters = []              #empty list
                        daughters.append("?")       #add to avoid GenPart info later
                        daughters.append("?")
                    
                if len(kin) >= 2:
                    decay = " + ".join(daughters)
                    d1_vec = kin[0]
                    d2_vec = kin[1]

                for j in range(len(pdg_evt)):

                    if int(pdg_evt[j]) == 25:
                        child_indices = [k for k in range(len(pdg_evt)) if int(midx_evt[k]) == j]       #get daughters' indices

                        if len(child_indices) > 1:                                                      #pair decays
                            for k in sorted(child_indices):
                                if len(daughters) <= 1:
                                    daughters.append(pdg_to_Name(int(pdg_evt[k])))
                                    kin.append(build_4vec(pt_evt[k], eta_evt[k], phi_evt[k], mass_evt[k]))

                            if len(kin) > 1:
                                decay = " + ".join(daughters)
                                d1_vec = kin[0]
                                d2_vec = kin[1]

                #------------------------------------------------------- add to df ------------------------------------------------------

                row = {
                       "event": event_number-1,
                       "lep_vec": lep_vec,
                       "antilep_vec": antilep_vec,
                       "MET": MET,
                       "d1_vec": d1_vec,
                       "d2_vec": d2_vec,
                       "h_decay": decay,
                       "event_type": event_type,
                       "VB_decay": VB_decay
                      }
                rows.append(row)

    df = pd.DataFrame(rows)                             #create pandas df
    rdf = expand_columns_vectors(df)                    #expand vector obj columns
    root_df = ROOT.RDF.FromPandas(rdf)                  #create ROOT df
    root_df.Snapshot("Events", output_root)             #save ROOT df

    return df


def main():
    folder = "/gwpool/users/gisilvestri/VHPolar/nanofiles"
    files = sorted(glob.glob(os.path.join(folder, "*.root")))

    n_blocks = 20
    blocks = [files[i::n_blocks] for i in range(n_blocks)]

    output_dir = "/gwpool/users/gisilvestri/VHPolar/RECO/root_jet_df"
    os.makedirs(output_dir, exist_ok=True)   # crea la cartella se non esiste

    outputs = [os.path.join(output_dir, f"complete_jet_df_{i}.root") for i in range(n_blocks)]

    with ProcessPoolExecutor(max_workers=n_blocks) as pool:
        results = list(pool.map(read_files, blocks, outputs))

    print("Created ROOT df files:", outputs)

if __name__ == "__main__":
    main()