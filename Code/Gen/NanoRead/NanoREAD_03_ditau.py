# Efficiency - h > tau+ + tau-
# Z : pt h > 40; eta h > 2.1
# W : pt h > 40; eta h > 2.1

import numpy as np
import pandas as pd
import ROOT
from tabulate import tabulate

def main():

    df_tot = ROOT.RDataFrame("df_ditau", "h_ditau.root")

    df_tot = (df_tot
              .Define("lep_vec", """
                  TLorentzVector v;
                  v.SetPtEtaPhiM(lep_vec_pt, lep_vec_eta, lep_vec_phi, lep_vec_mass);
                  return v;
              """)
              .Define("antilep_vec", """
                  TLorentzVector v;
                  v.SetPtEtaPhiM(antilep_vec_pt, antilep_vec_eta, antilep_vec_phi, antilep_vec_mass);
                  return v;
              """)
              .Define("vec_V", "lep_vec + antilep_vec")
              .Define("d1_vec", """
                  TLorentzVector v;
                  v.SetPtEtaPhiM(d1_vec_pt, d1_vec_eta, d1_vec_phi, d1_vec_mass);
                  return v;
              """)
              .Define("d2_vec", """
                  TLorentzVector v;
                  v.SetPtEtaPhiM(d2_vec_pt, d2_vec_eta, d2_vec_phi, d2_vec_mass);
                  return v;
              """)
              .Define("vec_h", "d1_vec + d2_vec")
             )

    #------------------------------------------ charged lepton selection ------------------------------------------
    
    df_tot = df_tot.Define("cut_e_mu", '(VB_decay != "v v~ / tau+ tau-") && (VB_decay != "tau v")')
    df_e_mu = df_tot.Filter("cut_e_mu")

    print(f"Fraction of e/mu events = {df_e_mu.Count().GetValue() / df_tot.Count().GetValue()}")

    #--------------------------------------------------------------------------------------------------------------

    #------------------------------------------------ ZH/WH Cutflow -----------------------------------------------
    
    df_ZH = df_e_mu.Filter("event_type == \"Z\"")
    df_WH = df_e_mu.Filter("event_type == \"W+/-\"")

    print(f"Fraction of ZH events = {df_ZH.Count().GetValue() / df_e_mu.Count().GetValue()}")
    print(f"Fraction of WH events = {df_WH.Count().GetValue() / df_e_mu.Count().GetValue()}")

    
    ##############################################
    # Define cuts
    ##############################################

    #ZH
    df_ZH = (df_ZH
        .Define("cut1_ZH", "vec_h.Pt() > 40")
        .Define("cut2_ZH", "vec_h.Eta() < 2.1"))

    # WH
    df_WH = (df_WH
        .Define("cut1_WH", "vec_h.Pt() > 40")
        .Define("cut2_WH", "vec_h.Eta() < 2.1"))

    
    ##############################################
    # Filter dataframe
    ##############################################

    df_ZH_f = (df_ZH
        .Filter("vec_h.Pt() > 40", "cut1_ZH")
        .Filter("vec_h.Eta() < 2.1", "cut2_ZH"))

    #WH
    df_WH_f = (df_WH
        .Filter("vec_h.Pt() > 40", "cut1_WH")
        .Filter("vec_h.Eta() < 2.1", "cut2_WH"))


    print("\n======================================== WH CutFlow Report ========================================")
    rep_wh = df_WH_f.Report()

    print("\n======================================== ZH CutFlow Report ========================================")
    rep_zh = df_ZH_f.Report()


    #------------------------------------------------ Print report -----------------------------------------------


    report_array_wh = []
    report_array_zh = []

    n_tot = df_WH.Count().GetValue()
    pass_prec = n_tot
    cut_list = ["cut1", "cut2"]

    i = 0
    for cut in rep_wh:
        cut_from_list = cut_list[i] + "_WH"
        name = cut.GetName()
        N_evt = cut.GetPass()
        part_eff = N_evt / pass_prec
        cumul_eff = N_evt / n_tot
        n_cut_only = df_WH.Filter(cut_from_list).Count().GetValue()
        tot_eff = n_cut_only / n_tot
        report_array_wh.append([
                                cut_list[i],
                                part_eff,      
                                cumul_eff,     
                                tot_eff        
                               ])
        pass_prec = N_evt
        i += 1

    n_tot = df_ZH.Count().GetValue()
    pass_prec = n_tot

    i = 0
    for cut in rep_zh:
        cut_from_list = cut_list[i] + "_ZH"
        name = cut.GetName()
        N_evt = cut.GetPass()
        part_eff = N_evt / pass_prec
        cumul_eff = N_evt / n_tot
        n_cut_only = df_ZH.Filter(cut_from_list).Count().GetValue()
        tot_eff = n_cut_only / n_tot
        report_array_zh.append([
                                cut_list[i],
                                part_eff,      
                                cumul_eff,     
                                tot_eff        
                               ])
        pass_prec = N_evt
        i += 1

    report_array_wh = np.array(report_array_wh, dtype=object)
    report_array_zh = np.array(report_array_zh, dtype=object)

    rows = []
    for i in range(len(report_array_wh)):
        rows.append(
                    [
                    report_array_zh[i][0],
                    report_array_zh[i][1],           #zh partial eff
                    report_array_zh[i][2],           #zh cumulative eff
                    report_array_zh[i][3],           #zh total eff
                    report_array_wh[i][1],           #wh partial eff
                    report_array_wh[i][2],           #wh cumulative eff
                    report_array_wh[i][3]            #wh total eff
                    ]            
                   )

    columns  = [
                "",
                "ZH partial eff", "ZH cumulative eff", "ZH total eff",
                "WH partial eff", "WH cumulative eff", "WH total eff"
                ]

    df_VH = pd.DataFrame(rows, columns = columns)

    print("\n" + "\u2550"*124)
    print("{:^124}".format("H > tau+ tau-   -  [kinematic cuts efficiency]"))
    print("\u2550"*124 + "\n")

    print(tabulate(df_VH.values.tolist(), headers=df_VH.columns.tolist(), colalign=["center"] * len(df_VH.columns), tablefmt="fancy_grid", floatfmt=".3f"))

if __name__ == "__main__":
    main()