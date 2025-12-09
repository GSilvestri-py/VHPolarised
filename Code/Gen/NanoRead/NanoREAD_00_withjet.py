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

def get_b_jets(index_list, t):

    kin_info = []

    for jet_idx in index_list:
        j1_vec = build_4vec_from_pt_eta_phi_m(t.GenJet_pt[jet_idx], t.GenJet_eta[jet_idx], t.GenJet_phi[jet_idx], t.GenJet_mass[jet_idx])

        for jet2_idx in index_list:
            if jet2_idx <= jet_idx:    continue

            print(f"jet pair indices = ({jet_idx}, {jet2_idx})")
            j2_vec = build_4vec_from_pt_eta_phi_m(t.GenJet_pt[jet2_idx], t.GenJet_eta[jet2_idx], t.GenJet_phi[jet2_idx], t.GenJet_mass[jet2_idx])

            vec_jj = j1_vec + j2_vec
            mass_jj = vec_jj.M
            print(f"jj invariant mass = {mass_jj}")

            tollerance = 15
            if abs(mass_jj - 125) < tollerance:
                print(f"> possible jet {jet_idx}-{jet2_idx} from h decay b pair\n")
                kin_info.append(build_4vec(t.GenJet_pt[jet_idx], t.GenJet_eta[jet_idx], t.GenJet_phi[jet_idx], t.GenJet_mass[jet_idx]))
                kin_info.append(build_4vec(t.GenJet_pt[jet2_idx], t.GenJet_eta[jet2_idx], t.GenJet_phi[jet2_idx], t.GenJet_mass[jet2_idx]))
    
    return kin_info


def read_file(filename, maxEvents=50):
    file = ROOT.TFile.Open(filename)
    t = file.Get("Events")

    if not t:
        print("No Events tree in", filename)
        return pd.DataFrame()

    rows = []
    for i in range(30):
        #t.GetEntries()):

        print(f"\n================================================================================ Event {i} ================================================================================")

        t.GetEntry(i)

        lep_vec = antilep_vec = build_4vec(0,0,0,0)
        d1_vec  = d2_vec = build_4vec(0,0,0,0)
        decay      = "No Higgs"
        event_type = "Z/W"
        VB_decay = ""

        is_b_inevent = False

        #------------------------------------------------------ LHE process -----------------------------------------------------

        print(f"PROCESS OUTLINE:")
        process_outline(t)

        for j in range(len(t.GenPart_pdgId)):
            if is_b_hadron(int(t.GenPart_pdgId[j])):
                is_b_inevent = True
                print(f">quark b in event")

        #-------------------------------------------------- VB decay - GenPart --------------------------------------------------

        daughters_VB = []

        for j in range(int(t.nGenPart)):
            pid = int(t.GenPart_pdgId[j])
            mother = int(t.GenPart_pdgId[t.GenPart_genPartIdxMother[j]]) if t.GenPart_genPartIdxMother[j] >= 0 else None

            if pid == 23 and mother != 25:
                event_type = "Z"
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

            if abs(pid) == 24 and mother != 25:
                event_type = "W+/-"
                child_indices = [k for k in range(int(t.nGenPart)) if int(t.GenPart_genPartIdxMother[k]) == j]      #get daughters' indices
                lep_indices = [k for k in child_indices if abs(int(t.GenPart_pdgId[k])) in [11, 12, 13, 14]]        #filter lepton indices

                if len(lep_indices) == 2:

                    for k in lep_indices:
                        pid = int(t.GenPart_pdgId[k])
                        daughters_VB.append(pdg_to_Name(pid))
                        vec = build_4vec(t.GenPart_pt[k], t.GenPart_eta[k], t.GenPart_phi[k], t.GenPart_mass[k])

                        if pid > 0:
                            lep_vec = vec
                        elif pid < 0:
                            antilep_vec = vec
       
        if len(daughters_VB) == 0:
            if event_type == "Z":
                VB_decay = "v v~ / tau+ tau-"
            if event_type == "W+/-":
                VB_decay = "tau v"
        else:
            VB_decay = " + ".join(daughters_VB)

        #--------------------------------------------------- h decay - GenPart ---------------------------------------------------

        daughters, kin = [], []
        
        if int(t.nGenVisTau) == 2:
            kin.append(build_4vec(t.GenVisTau_pt[0], t.GenVisTau_eta[0], t.GenVisTau_phi[0], t.GenVisTau_mass[0]))
            kin.append(build_4vec(t.GenVisTau_pt[1], t.GenVisTau_eta[1], t.GenVisTau_phi[1], t.GenVisTau_mass[1]))
            daughters.append("tau+")
            daughters.append("tau-")

        elif is_b_inevent:
            clean_list = clean_jets(lep_indices, t)
            daughters.append("b")
            daughters.append("bbar")
            kin = get_b_jets(clean_list, t)
          
        if len(kin) >= 2:
            decay = " + ".join(daughters)
            d1_vec = kin[0]
            d2_vec = kin[1]

        for j in range(int(t.nGenPart)):
                        
            if int(t.GenPart_pdgId[j]) == 25:

                #get daughters' indices
                child_indices = [k for k in range(int(t.nGenPart))
                                 if int(t.GenPart_genPartIdxMother[k]) == j]

                #work on particles pair decays
                if len(child_indices) > 1:
                    for k in sorted(child_indices):
                        if len(daughters) <= 1:         #skip if b+bbar or tau+ + tau-
                            daughters.append(pdg_to_Name(int(t.GenPart_pdgId[k])))
                            kin.append(build_4vec(t.GenPart_pt[k], t.GenPart_eta[k], t.GenPart_phi[k], t.GenPart_mass[k]))

                    decay = " + ".join(daughters)
                    if len(kin) > 1:
                        d1_vec = kin[0]
                        d2_vec = kin[1]
                    print(f"daughters = {daughters} // GenPart kinematic info = {kin}")
                    
        #------------------------------------------------------- add to df ------------------------------------------------------

        row = {
               "event": i,
               "lep_vec": lep_vec,
               "antilep_vec": antilep_vec,
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

    return pd.DataFrame(rows)


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
    df = read_file("nanoLatino_ZH_H_ZToLL_ZL__part0.root")
    df.to_pickle("events_withjets.pkl")

    print("\nShowing Pandas dataframe")
    print(df.head())

    rdf = expand_columns_vectors(df)
    print("Showing Expanded Pandas dataframe")
    print(rdf.head())

    #Pd dataframe --> ROOT df
    for col in df.columns:
        print(col, set(df[col].apply(type)))

    root_df = ROOT.RDF.FromPandas(rdf)
    root_df.Snapshot("Events", "file_withjets.root")

if __name__ == "__main__":
    main()