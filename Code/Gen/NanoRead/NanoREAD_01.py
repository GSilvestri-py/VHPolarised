import awkward as ak
from concurrent.futures import ProcessPoolExecutor          #optimise file-loop
import glob
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

            n_events = len(Gen_pdg)

            for i in range(n_events):
                event_number += 1
                lep_vec = antilep_vec = build_4vec(0,0,0,0)
                d1_vec  = d2_vec      = build_4vec(0,0,0,0)
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

                def mother_pdg(j):
                    m = int(midx_evt[j])
                    return int(pdg_evt[m]) if m >= 0 else None

                #-------------------------------------------------- VB decay - GenPart --------------------------------------------------

                for j in range(len(pdg_evt)):

                    pid = int(pdg_evt[j])
                    mother = mother_pdg(j)

                    if pid == 23 and mother != 25:
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

                        if len(lep_indices) == 2:
                            for k in lep_indices:
                                pid_k = int(pdg_evt[k])
                                daughters_VB.append(pdg_to_Name(pid_k))
                                vec = build_4vec(float(pt_evt[k]), float(eta_evt[k]), float(phi_evt[k]), float(mass_evt[k]))

                                if pid_k > 0: 
                                    lep_vec = vec
                                else: 
                                    antilep_vec = vec

                if len(daughters_VB) == 0:

                    if event_type=="Z":
                        VB_decay = "v v~ / tau+ tau-"  
                    if event_type=="W+/-":
                        VB_decay = "tau v"

                else:
                    VB_decay = " + ".join(daughters_VB)

                #--------------------------------------------------- h decay - GenPart ---------------------------------------------------

                daughters, kin = [], []

                for j in range(len(pdg_evt)):
                    if int(pdg_evt[j]) == 25:
                        child_indices = [k for k in range(len(pdg_evt)) if int(midx_evt[k]) == j]       #get daughters' indices

                        if len(child_indices) > 1:                                                      #pair decays
                            for k in sorted(child_indices)[:2]:
                                daughters.append(pdg_to_Name(int(pdg_evt[k])))
                                kin.append([float(pt_evt[k]), float(eta_evt[k]), float(phi_evt[k]), float(mass_evt[k])])
                            break

                if len(kin) >= 2:
                    decay = " + ".join(daughters)
                    d1_vec = build_4vec(*kin[0])
                    d2_vec = build_4vec(*kin[1])

                #------------------------------------------------------- add to df ------------------------------------------------------

                row = {
                       "event": event_number-1,
                       "lep_vec": lep_vec,
                       "antilep_vec": antilep_vec,
                       "d1_vec": d1_vec,
                       "d2_vec": d2_vec,
                       "h_decay": decay,
                       "event type": event_type,
                       "VB_decay": VB_decay
                      }
                rows.append(row)

    df = pd.DataFrame(rows)                             #create pandas df
    rdf = expand_columns_vectors(df)                    #expand vector obj columns
    root_df = ROOT.RDF.FromPandas(rdf)                  #create ROOT df
    root_df.Snapshot("Events", output_root)             #save ROOT df

    return df


def main():

    if len(sys.argv) < 2:
        print("use example: python3 script.py <folder_path>")
        sys.exit(1)

    folder = sys.argv[1] 
    files = sorted(glob.glob(os.path.join(folder, "*.root")))

    n_blocks = 20
    blocks = [files[i::n_blocks] for i in range(n_blocks)]
    outputs = [f"complete_df_{i}.root" for i in range(n_blocks)]

    with ProcessPoolExecutor(max_workers=n_blocks) as pool:
        results = list(pool.map(read_files, blocks, outputs))

    print("Created ROOT df files:", outputs)

if __name__ == "__main__":
    main()