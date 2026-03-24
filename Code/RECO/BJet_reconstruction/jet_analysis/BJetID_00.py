#Build pickle file with structure:
# - first column = index
# - last column = true tag
# - every other column = feature


from concurrent.futures import ProcessPoolExecutor          
import glob
import os
import pandas as pd
from tqdm import tqdm
import uproot                                               

# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

def read_files(file_list, output_path):
    event_number = 0
    
    rows = []
    njet = 0

    for full_path in tqdm(file_list, desc="Processing files", unit="file"):
        filename = os.path.basename(full_path)

        with uproot.open(full_path) as f:                                       
            if "Events" not in f:
                continue
            tree = f["Events"]

            Gen_pdg   = tree["GenPart_pdgId"].array(library="ak")
            Gen_midx  = tree["GenPart_genPartIdxMother"].array(library="ak")
            

            n_RECOjets = tree["nJet"].array(library="ak")
            Jet_pt     = tree["Jet_pt"].array(library="ak")
            Jet_eta    = tree["Jet_eta"].array(library="ak")
            Jet_p_flav = tree["Jet_partonFlavour"].array(library="ak")


            btagDeepFlavB   = tree["Jet_btagDeepFlavB"].array(library="ak")
            btagDeepFlavCvB = tree["Jet_btagDeepFlavCvB"].array(library="ak")
            btagDeepFlavCvL = tree["Jet_btagDeepFlavCvL"].array(library="ak")
            btagDeepFlavQG  = tree["Jet_btagDeepFlavQG"].array(library="ak")

            n_events = len(Gen_pdg)

            for i in range(n_events):
                event_number += 1
                
                pdg_evt   = Gen_pdg[i]
                midx_evt  = Gen_midx[i]
            
                #RECOjet array i component
                n_RECOjets_evt  = n_RECOjets[i]
                Jet_pt_evt      = Jet_pt[i]
                Jet_eta_evt     = Jet_eta[i]
                Jet_p_flav_evt  = Jet_p_flav[i]

                btagDeepFlavB_evt   = btagDeepFlavB[i]
                btagDeepFlavCvB_evt = btagDeepFlavCvB[i]
                btagDeepFlavCvL_evt = btagDeepFlavCvL[i]
                btagDeepFlavQG_evt  = btagDeepFlavQG[i]
                
                two_matches_bb = False
                b_child = 0

                #-------------------------------------------------- find bb events --------------------------------------------------

                for j in range(len(pdg_evt)):
                    if int(pdg_evt[j]) == 25:
                        #get daughters' indices
                        child_indices = [k for k in range(len(pdg_evt)) if int(midx_evt[k]) == j]

                        for n in child_indices:
                            if abs(pdg_evt[n]) == 5:
                                b_child += 1
                
                if b_child >= 2:

                    idx_list = []
                    for idx in range(n_RECOjets_evt):

                        if abs(int(Jet_p_flav_evt[idx])) == 5:
                            idx_list.append(idx)

                    if len(idx_list) == 2:
                        two_matches_bb = True           #require exactly two matches at parton flavour

                #--------------------------------------------------- h decay - GenPart ---------------------------------------------------
                
                #only bjets from 2 matches
                if two_matches_bb:
                    b_idx_list = []
                    
                    for idx in range(n_RECOjets_evt):

                        #require exactly two matches at parton flavour
                        if abs(int(Jet_p_flav_evt[idx])) == 5:
                            b_idx_list.append(idx)  
            
                    if len(b_idx_list) < 2 or len(b_idx_list) > 2:  continue
                    for idx in range(n_RECOjets_evt):

                        #get values
                        row = []
                        njet += 1

                        n_jets   = n_RECOjets_evt
                        DeepFlavB = btagDeepFlavB_evt[idx]
                        DeepFlavCvB = btagDeepFlavCvB_evt[idx]
                        DeepFlavCvL = btagDeepFlavCvL_evt[idx]
                        DeepFlavQG = btagDeepFlavQG_evt[idx]
                        pt_jet = Jet_pt_evt[idx]
                        eta_jet = Jet_eta_evt[idx]
                        abs_eta_jet = abs(eta_jet)

                        if abs(Jet_p_flav_evt[idx]) == 5:
                            tag = 'b jet'
                        else:
                            tag = 'non b jet'


                        row = {
                            "N jet": njet,
                            "N_jets_in_event": n_jets,
                            "FlavB": DeepFlavB,
                            "FlavCvB": DeepFlavCvB,
                            "FlavCvL": DeepFlavCvL,
                            "FlavQG": DeepFlavQG,
                            "Pt_jet": pt_jet,
                            "Eta_jet": abs_eta_jet,
                            "Tag": tag
                            }

                        rows.append(row)

    jet_df = pd.DataFrame(rows)
    jet_df.to_pickle(output_path)
    return jet_df


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------

        
def main():
    folder = "input_folder"
    files = sorted(glob.glob(os.path.join(folder, "*.root")))

    n_blocks = 20
    blocks = [files[i::n_blocks] for i in range(n_blocks)]

    output_dir = "output_folder"
    os.makedirs(output_dir, exist_ok=True)  

    outputs = [os.path.join(output_dir, f"b_vs_nonb_df_{i}.pkl") for i in range(n_blocks)]

    with ProcessPoolExecutor(max_workers=n_blocks) as pool:
        results = list(pool.map(read_files, blocks, outputs))

    # ----------------- Create single pkl file ----------------------

    pickles_path = output_dir
    dfs = []

    for file in os.listdir(pickles_path):
        if file.endswith(".pkl"):
            full_path = os.path.join(pickles_path, file)
            df = pd.read_pickle(full_path)
            dfs.append(df)

    tot_df = pd.concat(dfs, ignore_index=True)

    first_col = tot_df.columns[0]
    tot_df[first_col] = range(len(tot_df))

    print(f"Number of rows (len): {len(tot_df)}")
    print(f"DataFrame shape: {tot_df.shape}\n")

    print(tot_df.head(20))
    tot_df.to_pickle(os.path.join(output_dir, "total_b_vs_nonb_df.pkl"))

if __name__ == "__main__":
    main()