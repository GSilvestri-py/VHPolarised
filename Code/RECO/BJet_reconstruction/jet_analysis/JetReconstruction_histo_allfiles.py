import awkward as ak
from concurrent.futures import ProcessPoolExecutor          
import glob
import os
import pandas as pd
import ROOT
from tqdm import tqdm
import uproot                                             


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------


def get_jets_idx_lists(n_RECOjets_evt, Jet_p_flav_evt):
    '''
    gives b jet index list if exactly 2 matches -> otherwise empty list
    '''
    b_idx_list, non_b_idx_list = [], []
    for idx in range(n_RECOjets_evt):

        if abs(int(Jet_p_flav_evt[idx])) == 5:
            b_idx_list.append(idx)  #b idx list
        else:
            non_b_idx_list.append(idx)

    if len(b_idx_list) < 2 or len(b_idx_list) > 2:
        b_idx_list = []   #empty list if not 2 matches
        non_b_idx_list = []

    return b_idx_list, non_b_idx_list


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------


def build_variables_dicts_evt(btagDeepFlavB_evt, btagDeepFlavQG_evt, btagDeepFlavCvL_evt, btagDeepFlavCvB_evt, Jet_pt_evt, Jet_eta_evt, b_idx_list, non_b_idx_list):
    '''
    build dict with key = var name and value = array of values for b/non b jets
    '''
    var_dict = {
                "Jet_btagDeepFlavB": btagDeepFlavB_evt,
                "Jet_btagDeepFlavCvB": btagDeepFlavCvB_evt,
                "Jet_btagDeepFlavCvL": btagDeepFlavCvL_evt,
                "Jet_btagDeepFlavQG": btagDeepFlavQG_evt,
                "Jet_pt": Jet_pt_evt,
                "Jet_eta": Jet_eta_evt
            }

    result_b = {}
    result_nonb = {}

    for name, arr in var_dict.items():
        result_b[name] = [arr[i] for i in b_idx_list]
        result_nonb[name] = [arr[i] for i in non_b_idx_list]
    
    return result_b, result_nonb


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------


def read_files(file_list, output_root):
    event_number = 0
    rows_b, rows_nonb = [], []

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
            Jet_phi    = tree["Jet_phi"].array(library="ak")
            Jet_mass   = tree["Jet_mass"].array(library="ak")
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
                
                is_b_inevent = False
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
                    is_b_inevent = True

                #--------------------------------------------------- h decay - GenPart ---------------------------------------------------


                if is_b_inevent:
                    b_idx_list, non_b_idx_list = get_jets_idx_lists(n_RECOjets_evt, Jet_p_flav_evt)
                    result_b_dict, result_nonb_dict = build_variables_dicts_evt(btagDeepFlavB_evt, btagDeepFlavQG_evt, btagDeepFlavCvL_evt, btagDeepFlavCvB_evt, Jet_pt_evt, Jet_eta_evt, b_idx_list, non_b_idx_list)

                    row_b = {"event": i}
                    row_b.update(result_b_dict)

                    row_nonb = {"event": i}
                    row_nonb.update(result_nonb_dict)

                    rows_b.append(row_b)
                    rows_nonb.append(row_nonb)
                    
                #------------------------------------------------------- add to df ------------------------------------------------------

    df_b = pd.DataFrame(rows_b)                            
    df_nonb = pd.DataFrame(rows_nonb)                 

    return df_b, df_nonb


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------


def plot_global_distributions_root(df_b, df_nonb, nbins=50):

    variables = [col for col in df_b.columns if col != "event"]

    for var in variables:
        #values extraction
        all_b = df_b[var].explode().tolist()
        all_nonb = df_nonb[var].explode().tolist()

        if len(all_b) == 0 or len(all_nonb) == 0:
            print(f"Skipping {var}: empty data")
            continue

        xmin = min(all_b + all_nonb)
        xmax = max(all_b + all_nonb)

        if xmax > 300:
            all_b    = [x for x in all_b if x < 300]
            all_nonb = [x for x in all_nonb if x < 300]

            xmin = min(all_b + all_nonb)
            xmax = 300

        c = ROOT.TCanvas(f"c_{var}", f"{var} — global distribution", 800, 600)

        h_b = ROOT.TH1F(f"h_b_{var}", f"{var};{var};Counts", nbins, xmin, xmax)
        h_nonb = ROOT.TH1F(f"h_nonb_{var}", f"{var};{var};Counts", nbins, xmin, xmax)

        for v in all_b:
            h_b.Fill(v)
        for v in all_nonb:
            h_nonb.Fill(v)

        if h_b.Integral() > 0:
            h_b.Scale(1.0 / h_b.Integral())
        if h_nonb.Integral() > 0:
            h_nonb.Scale(1.0 / h_nonb.Integral())

        max_b = h_b.GetMaximum()
        max_nonb = h_nonb.GetMaximum()

        ymax = max(max_b, max_nonb) * 1.10
        h_b.SetMaximum(ymax)
        h_nonb.SetMaximum(ymax)

    
        h_b.SetLineColor(ROOT.kBlue)
        h_b.SetFillColorAlpha(ROOT.kBlue, 0.2)

        h_nonb.SetLineColor(ROOT.kRed)
        h_nonb.SetFillColorAlpha(ROOT.kRed, 0.2)

        h_b.Draw("HIST")
        h_nonb.Draw("HIST SAME")

        ROOT.gStyle.SetOptStat(0)
        leg = ROOT.TLegend(0.46, 0.70, 0.88, 0.88)
        leg.SetBorderSize(1)
        leg.SetFillStyle(0)   
        leg.SetTextSize(0.035)
        leg.SetTextFont(42)   


        entry_b = f"b jets     N={int(h_b.GetEntries())}, mean={h_b.GetMean():.3f}"
        entry_nonb = f"non-b jets N={int(h_nonb.GetEntries())}, mean={h_nonb.GetMean():.3f}"

        e1 = leg.AddEntry(h_b, entry_b, "l")
        e2 = leg.AddEntry(h_nonb, entry_nonb, "l")
        leg.SetMargin(0.15)

        leg.Draw()
        c.Update()
        c.SaveAs(f"{var}_histo.ROOT")


# ------- ------- ------- ------- ------- ------- ------- ------- ------- ------- -------


def main():
    
    folder = "nanofiles_folder"
    files = sorted(glob.glob(os.path.join(folder, "*.root")))

    n_blocks = 20
    blocks = [files[i::n_blocks] for i in range(n_blocks)]

    output_dir = "output_folder"
    os.makedirs(output_dir, exist_ok=True)   

    outputs = [os.path.join(output_dir, f"complete_jet_df_{i}.root") for i in range(n_blocks)]

    with ProcessPoolExecutor(max_workers=n_blocks) as pool:
        results = list(pool.map(read_files, blocks, outputs))

    df_b = pd.concat([df_b for (df_b, _) in results], ignore_index=True)
    df_nonb = pd.concat([df_nonb for (_, df_nonb) in results], ignore_index=True)

    plot_global_distributions_root(df_b, df_nonb)
    
if __name__ == "__main__":
    main()