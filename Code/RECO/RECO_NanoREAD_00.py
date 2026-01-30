# READ SINGLE NANO FILE

import numpy as np
import ROOT
import pandas as pd
import vector

def pdg_to_Name(pdgId):
    #map PDGs to particle names; 
    #return PDG value if not in dictionary
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

def higgs_decay(n, tree):
    daughters_h = []
    child_indices = [k for k in range(len(tree.GenPart_pdgId)) if int(tree.GenPart_genPartIdxMother[k]) == n]

    if len(child_indices) > 1:
        daughters_h.extend(int(tree.GenPart_pdgId[k]) for k in child_indices)
    return daughters_h

def mother_pdg(j, tree):
    m = int(tree.GenPart_genPartIdxMother[j])
    return int(tree.GenPart_pdgId[m]) if m >= 0 else None

def process_outline(tree):
    #n_lhe_evt, lhe_pdg_evt, lhe_status_evt, pdg_evt, midx_evt):
    
    is_b_inevent = False
    father, daughters = [], []
    for j in range(tree.nLHEPart):
        pid = tree.LHEPart_pdgId[j]
        status = tree.LHEPart_status[j]

        if status == -1:
            father.append(pid)
        if status == 1:
            daughters.append(pid)
    
    daughters_h = []

    for n in range(len(tree.GenPart_pdgId)):                                               #loop over particles in GenPart

        pid = int(tree.GenPart_pdgId[n])
        mother = mother_pdg(n, tree)

        if pid == 25:
            daughters_h = higgs_decay(n, tree)

            if len(daughters_h) > 1:
                print(f"{father[0]} + {father[1]} --> " + " + ".join(map(str, daughters)) + ", 25 > " + " + ".join(map(str, daughters_h))  )
                print("\n")
    return


def get_z_daughters(j, t):

    daughters_VB = []
    lep_vec = antilep_vec = MET = build_4vec(0,0,0,0)

    child_indices = [k for k in range(int(t.nGenPart)) if int(t.GenPart_genPartIdxMother[k]) == j]      #get daughters' indices
    lep_indices = [k for k in child_indices if abs(int(t.GenPart_pdgId[k])) in [11, 13]]                #filter lepton indices

    if len(lep_indices) == 2:
        for k in lep_indices:
            pid = int(t.GenPart_pdgId[k])
            daughters_VB.append(pdg_to_Name(pid))
            vec = build_4vec(t.GenPart_pt[k], t.GenPart_eta[k], t.GenPart_phi[k], t.GenPart_mass[k])

            if pid > 0:
                lep_vec = vec
            elif pid < 0:
                antilep_vec = vec

    return lep_vec, antilep_vec, MET, daughters_VB, lep_indices


def get_w_daughters(j, t):

    daughters_VB = []
    lep_vec = antilep_vec = MET = build_4vec(0,0,0,0)

    child_indices = [k for k in range(int(t.nGenPart)) if int(t.GenPart_genPartIdxMother[k]) == j]      #get daughters' indices
    lep_indices = [k for k in child_indices if abs(int(t.GenPart_pdgId[k])) in [11, 13]]        #filter lepton indices
    MET = build_4vec(t.GenMET_pt, 0, t.GenMET_phi, 0)
    
    if len(lep_indices) == 2:

        for k in lep_indices:
            pid = int(t.GenPart_pdgId[k])
            daughters_VB.append(pdg_to_Name(pid))
            vec = build_4vec(t.GenPart_pt[k], t.GenPart_eta[k], t.GenPart_phi[k], t.GenPart_mass[k])

            if pid > 0:
                lep_vec = vec
            elif pid < 0:
                antilep_vec = vec
    
    return lep_vec, antilep_vec, MET, daughters_VB, lep_indices


def get_VB_daughters(t):
    event_type = "Z/W"

    for j in range(int(t.nGenPart)):
        pid = int(t.GenPart_pdgId[j])
        mother = int(t.GenPart_pdgId[t.GenPart_genPartIdxMother[j]]) if t.GenPart_genPartIdxMother[j] >= 0 else None

        if pid == 23 and mother != 25:
            event_type = "Z"
            lep_vec, antilep_vec, MET, daughters_VB, lep_indices = get_z_daughters(j, t)
            

        if abs(pid) == 24 and mother != 25:
            event_type = "W+/-"
            lep_vec, antilep_vec, MET, daughters_VB, lep_indices = get_w_daughters(j, t)

    return lep_vec, antilep_vec, MET, event_type, daughters_VB, lep_indices

'''
def clean_jets(lep_indices, t):
    cleaned_list = []   #list with b jet candidate index

    for k in lep_indices:   
        pid = int(t.GenPart_pdgId[k])   #get particle ID
        if abs(pid) in [11, 13]:
            lep_vec = build_4vec_from_pt_eta_phi_m(t.GenPart_pt[k], t.GenPart_eta[k], t.GenPart_phi[k], t.GenPart_mass[k])

            for jet_idx in range(t.nGenJet):
                jet_vec = build_4vec_from_pt_eta_phi_m(t.GenJet_pt[jet_idx], t.GenJet_eta[jet_idx], t.GenJet_phi[jet_idx], t.GenJet_mass[jet_idx])

                deltaR_lj = lep_vec.deltaR(jet_vec)
                if deltaR_lj < 0.1:
                    print(f"jet {jet_idx} is a charged lepton jet: ΔR = {deltaR_lj}")
                else:
                    cleaned_list.append(jet_idx)  #save index if candidate b jet

    return cleaned_list
'''

def get_b_jets(pair_inv_mass_list, t):
    '''
    - read partonflavour -> identify b jets
        - if less than two jets -> jet pt under threshold -> save only one
    - if t least two jets:
        - if two jets -> found b pair
        - if more:
            - check pt
    '''
    idx_list, kin_info = [], []
    for idx in range(t.nJet):

        if abs(int(t.Jet_hadronFlavour[idx])) == 5:
            idx_list.append(idx)  #b idx list

    if len(idx_list) < 2:
        print("Less than two jet for this event -> pt under threshold")
        kin_info.append(vector.obj(px = 0, py = 0, pz = 0, E = 0))
        kin_info.append(vector.obj(px = 0, py = 0, pz = 0, E = 0))

    if len(idx_list) == 2:
        j1 = build_4vec_from_pt_eta_phi_m(t.Jet_pt[idx_list[0]], t.Jet_eta[idx_list[0]], t.Jet_phi[idx_list[0]], t.Jet_mass[idx_list[0]])
        j2 = build_4vec_from_pt_eta_phi_m(t.Jet_pt[idx_list[1]], t.Jet_eta[idx_list[1]], t.Jet_phi[idx_list[1]], t.Jet_mass[idx_list[1]])
        kin_info.append(j1)
        kin_info.append(j2)
        pair_inv_mass_list.append((j1+j2).M)

    if len(idx_list) > 2:
        higher_pt = 0

        for id1 in idx_list:
            pt_id1 = t.Jet_pt[id1]

            for id2 in idx_list:
                if id2 != id1:
                    pt_id2 = t.Jet_pt[id2]
                    pt_sum = pt_id1 + pt_id2

                    if pt_sum > higher_pt:
                        higher_pt = pt_sum
                        pair = (id1, id2)
        
        idx1, idx2 = pair
        j1 = build_4vec_from_pt_eta_phi_m(t.Jet_pt[idx1], t.Jet_eta[idx1], t.Jet_phi[idx1], t.Jet_mass[idx1])
        j2 = build_4vec_from_pt_eta_phi_m(t.Jet_pt[idx2], t.Jet_eta[idx2], t.Jet_phi[idx2], t.Jet_mass[idx2])
        kin_info.append(j1)
        kin_info.append(j2)
        inv_m = (j1+j2).M
        pair_inv_mass_list.append(inv_m)

        print(f"pair with higher pt = {pair} and invariant mass = {inv_m}")

    return kin_info

def get_max_index(list):
    max_v = max(list)
    return list.index(max_v)


def read_file(filename, maxEvents=450):
    file = ROOT.TFile.Open(filename)
    t = file.Get("Events")

    #t.GetEntries()

    if not t:
        print("No Events tree in", filename)
        return pd.DataFrame()

    rows, inv_mass_list = [], []
    ditau_mass_list = []
    counter_empy_dict = 0
    counter_b_inevent = 0
    counter_0_jets, counter_1_jet, counter_2_jet, counter_2plus_jets = 0, 0, 0, 0

    for i in range(t.GetEntries()):

        print(f"\n================================================================================ Event {i} ================================================================================")

        t.GetEntry(i)

        d1_vec   = d2_vec = build_4vec(0,0,0,0)
        decay    = "No Higgs"
        VB_decay = ""

        is_b_inevent = False

        #------------------------------------------------------ LHE process -----------------------------------------------------

        print(f"PROCESS OUTLINE:")
        process_outline(t)
    
        for j in range(int(t.nGenPart)):
            if int(t.GenPart_pdgId[j]) == 25:
                #get daughters' indices
                child_indices = [k for k in range(int(t.nGenPart)) if int(t.GenPart_genPartIdxMother[k]) == j]

                for n in child_indices:
                    if abs(t.GenPart_pdgId[n]) == 5:
                        is_b_inevent = True

        #-------------------------------------------------- VB decay - GenPart --------------------------------------------------
 
        lep_vec, antilep_vec, MET, event_type, daughters_VB, lep_indices = get_VB_daughters(t)
        
        if len(daughters_VB) == 0:
            if event_type == "Z":
                VB_decay = "v v~ / tau+ tau-"
            if event_type == "W+/-":
                VB_decay = "tau v"
        else:
            VB_decay = " + ".join(daughters_VB)

        #--------------------------------------------------- h decay - GenPart ---------------------------------------------------

        daughters, kin = [], []
        bjet_idx_list = []
        
        if int(t.nGenVisTau) == 2:
            print(f"GenVisTau_genPartIdxMother = {t.GenVisTau_genPartIdxMother}")
            kin.append(build_4vec(t.GenVisTau_pt[0], t.GenVisTau_eta[0], t.GenVisTau_phi[0], t.GenVisTau_mass[0]))
            kin.append(build_4vec(t.GenVisTau_pt[1], t.GenVisTau_eta[1], t.GenVisTau_phi[1], t.GenVisTau_mass[1]))
            j1 = build_4vec_from_pt_eta_phi_m(t.GenVisTau_pt[0], t.GenVisTau_eta[0], t.GenVisTau_phi[0], t.GenVisTau_mass[0])
            j2 = build_4vec_from_pt_eta_phi_m(t.GenVisTau_pt[1], t.GenVisTau_eta[1], t.GenVisTau_phi[1], t.GenVisTau_mass[1])
            inv_mass = (j1+j2).M
            ditau_mass_list.append(inv_mass)
            daughters.append("tau+")
            daughters.append("tau-")

        if is_b_inevent:
            counter_b_inevent += 1
            daughters.append("b")
            daughters.append("bbar")

            print(f"n Jets = {int(t.nJet)}")
            print(f"t.Jet_partonFlavour = {list(t.Jet_hadronFlavour)}")
            print(f"t.Jet_pt = {list(t.Jet_pt)}")

            for j in range(t.nJet):
                print(f"idx = {j}")
                if abs(int(t.Jet_hadronFlavour[j])) == 5:
                    bjet_idx_list.append(j)

            if len(bjet_idx_list) == 0:
                counter_0_jets += 1
                print(f"no bjet in GenJet_hadronFlavour")

            if len(bjet_idx_list) == 1:
                counter_1_jet += 1
                print(f"one bjet found in event")
            
            if len(bjet_idx_list) == 2:
                counter_2_jet += 1

            if len(bjet_idx_list) > 2:
                counter_2plus_jets += 1
            
            kin = get_b_jets(inv_mass_list, t)
        
        if len(kin) >= 2:
            decay = " + ".join(daughters)
            d1_vec = kin[0]
            d2_vec = kin[1]

        for j in range(int(t.nGenPart)):
                        
            if int(t.GenPart_pdgId[j]) == 25:

                #get daughters' indices
                child_indices = [k for k in range(int(t.nGenPart)) if int(t.GenPart_genPartIdxMother[k]) == j]

                #work on particles pair decays
                if len(child_indices) > 1:
                    for k in sorted(child_indices):
                        if len(daughters) <= 1:         #skip if b+bbar or tau+ + tau-
                            daughters.append(pdg_to_Name(int(t.GenPart_pdgId[k])))
                            kin.append(build_4vec(t.GenPart_pt[k], t.GenPart_eta[k], t.GenPart_phi[k], t.GenPart_mass[k]))

                    
                    if len(kin) > 1:
                        decay = " + ".join(daughters)
                        d1_vec = kin[0]
                        d2_vec = kin[1]
                    print(f"daughters = {daughters} // GenPart kinematic info = {kin}")
                    
        #------------------------------------------------------- add to df ------------------------------------------------------

        row = {
               "event": i,
               "lep_vec": lep_vec,
               "antilep_vec": antilep_vec,
               "MET": MET,
               "d1_vec": d1_vec,
               "d2_vec": d2_vec,
               "h_decay": decay,
               "event_type": event_type,
               "VB_decay": VB_decay
              }

        print(f"\n---------------------------------------- event row ----------------------------------------")
        for key, value in row.items():
            print(f"{key}: {value}")
            
        rows.append(row)
    print(f"n 0 jets = {counter_0_jets}")
    print(f"n 1 jets = {counter_1_jet}")
    print(f"n 2 jets = {counter_2_jet}")
    print(f"2+ jets  = {counter_2plus_jets}")
    print(f"bb events = {counter_b_inevent}")

    return pd.DataFrame(rows), counter_empy_dict, inv_mass_list, ditau_mass_list


def expand_columns_vectors(df):
    vector_cols = []
    for col in df.columns:
        
        if df[col].apply(lambda v: hasattr(v, "pt")).all():         #check if vector.obj
            vector_cols.append(col)
            df[f"{col}_pt"]   = df[col].apply(lambda v: v.pt)
            df[f"{col}_eta"]  = df[col].apply(lambda v: v.eta)
            df[f"{col}_phi"]  = df[col].apply(lambda v: v.phi)
            df[f"{col}_mass"] = df[col].apply(lambda v: v.mass)

    df = df.drop(columns=vector_cols)                               #drop old columns

    return df


def main():

    #pt jets histogram    
    df, counter, mass_list, ditau_mass_list = read_file("nanoLatino_ZH_H_ZToLL_ZL__part108.root")
    
    print((ditau_mass_list))
    h = ROOT.TH1D("h", "b jet invariant mass", 30, min(mass_list), max(mass_list))
    for m in mass_list:
        h.Fill(m)
    h.GetXaxis().SetTitle("Invariant mass [GeV]")
    h.GetYaxis().SetTitle("Counts")
    canvas = ROOT.TCanvas("c", "c", 1400, 1000)
    h.Draw()
    canvas.SaveAs("RECO_bjet_inv_mass.png")

    df2 = ROOT.RDataFrame("Events", "nanoLatino_ZH_H_ZToLL_ZL__part108.root")
    df2 = df2.Define("reco_njets", "nJet")
    h2 = df2.Histo1D("reco_njets")
    c = ROOT.TCanvas("c2", "c2", 1400, 800)
    h2.Draw()
    c.Draw()
    c.SaveAs("recoNjets.png")
    '''
    h_tau = ROOT.TH1D("h", "tau jet invariant mass", 30, min(ditau_mass_list), max(ditau_mass_list))
    for m in ditau_mass_list:
        h_tau.Fill(m)
    h_tau.GetXaxis().SetTitle("Invariant mass [GeV]")
    h_tau.GetYaxis().SetTitle("Counts")
    c = ROOT.TCanvas("canvas", "canvas", 1400, 1000)
    h_tau.Draw()
    c.SaveAs("taujet_inv_mass.png")
    '''
    '''
    df = ROOT.RDataFrame("Events", "nanoLatino_ZH_H_ZToLL_ZL__part108.root")
    df = df.Define("pt_jets_reco", "Jet_pt")\
           .Define("eta_jets_reco", "Jet_eta")
    h = df.Histo1D( "pt_jets_reco")
    canvas = ROOT.TCanvas("canvas", "canvas", 1400, 800)
    h.GetXaxis().SetRangeUser(0, 200)  

    h.Draw()
    canvas.Draw()
    canvas.SaveAs("recojet_pt.png")

    h = df.Histo1D("eta_jets_reco")
    canvas = ROOT.TCanvas("c", "c", 1440, 1000)
    h.Draw()
    canvas.Draw()
    canvas.SaveAs("recojet_eta.png")
    '''

    
    df.to_pickle("events_withjets.pkl")

    print("\nShowing Pandas dataframe")
    print(df.head())

    rdf = expand_columns_vectors(df)
    print("Showing Expanded Pandas dataframe")
    print(rdf.head())

    #Pd dataframe --> ROOT df
    print(counter)
    root_df = ROOT.RDF.FromPandas(rdf)
    root_df.Snapshot("Events", "file_withjets.root")

if __name__ == "__main__":
    main()