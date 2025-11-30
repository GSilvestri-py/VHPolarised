# Efficiency - h > b + bar
# Z : pt Z; deltaphi Z jj
# W : pt W; pt j1; pt j2; pt jj; deltaphi W jj; deltaphi lepton MET

import ROOT

def main():
    df_tot = ROOT.RDataFrame("df_bbar", "h_bb.root")

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

    #ZH
    df_ZH = (df_ZH
        .Filter("vec_V.Pt() > 75", "cut1_ZH")
        .Filter("vec_V.DeltaPhi(vec_h) > 2.5", "cut2_ZH"))

    #WH
    df_WH = (df_WH
        .Define("charged_lep", """
            return (VB_decay == "e- + ~ve" || VB_decay == "ve~ + e-" || VB_decay == "mu- + ~vmu" || VB_decay == "vmu~ + mu-") ? lep_vec : antilep_vec;
        """)
        .Define("neutrino", """
            return (VB_decay == "e- + ~ve" || VB_decay == "ve~ + e-" || VB_decay == "mu- + ~vmu" || VB_decay == "vmu~ + mu-") ? antilep_vec : lep_vec;
        """)
        .Filter("vec_V.Pt() > 150", "cut1_WH")
        .Filter("(d1_vec.Pt() > 25) && (d2_vec.Pt() > 25)", "cut2_WH")
        .Filter("vec_h.Pt() > 100", "cut3_WH")
        .Filter("charged_lep.DeltaPhi(neutrino) < 2.0", "cut4_WH")
        .Filter("vec_V.DeltaPhi(vec_h) > 2.5", "cut5_WH"))

    '''
    cols = df_WH.AsNumpy(["VB_decay", "lep_vec", "antilep_vec", "charged_lep", "neutrino"])

    vb_decays    = cols["VB_decay"]
    lep_vecs     = cols["lep_vec"]
    antilep_vecs = cols["antilep_vec"]
    charged_leps = cols["charged_lep"]
    neutrinos    = cols["neutrino"]

    for i in range(10):
        vb   = vb_decays[i]
        lep  = lep_vecs[i]
        antl = antilep_vecs[i]
        chlp = charged_leps[i]
        nu   = neutrinos[i]

        print(f"Event {i} | VB_decay = {vb}")
        print(f"  lep_vec      -> pt={lep.Pt():.2f}, eta={lep.Eta():.2f}, phi={lep.Phi():.2f}, M={lep.M():.2f}")
        print(f"  antilep_vec  -> pt={antl.Pt():.2f}, eta={antl.Eta():.2f}, phi={antl.Phi():.2f}, M={antl.M():.2f}")
        print(f"  charged_lep  -> pt={chlp.Pt():.2f}, eta={chlp.Eta():.2f}, phi={chlp.Phi():.2f}, M={chlp.M():.2f}")
        print(f"  neutrino     -> pt={nu.Pt():.2f}, eta={nu.Eta():.2f}, phi={nu.Phi():.2f}, M={nu.M():.2f}")
        print("-"*70)
    '''

    print("======================================== WH CutFlow Report ========================================")
    df_WH.Report().Print()

    print("======================================== ZH CutFlow Report ========================================")
    df_ZH.Report().Print()

if __name__ == "__main__":
    main()