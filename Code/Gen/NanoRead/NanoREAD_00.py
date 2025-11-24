# READ SINGLE NANO FILE

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


def read_file(filename, maxEvents=50):
    file = ROOT.TFile.Open(filename)
    t = file.Get("Events")

    if not t:
        print("No Events tree in", filename)
        return pd.DataFrame()

    rows = []
    for i in range(t.GetEntries()):

        print(f"----------- Event {i} -----------")

        t.GetEntry(i)

        lep_vec = antilep_vec = build_4vec(0,0,0,0)
        d1_vec  = d2_vec = build_4vec(0,0,0,0)
        decay      = "No Higgs"
        event_type = "Z/W"
        VB_decay = ""

        #------------------------------------------------------ LHE process -----------------------------------------------------

        for j in range(int(t.nLHEPart)):
            pid = int(t.LHEPart_pdgId[j])       #get particle ID
            status = int(t.LHEPart_status[j])   #get particle status

            print(f" pid = {pid}, status = {status}")

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

        for j in range(int(t.nGenPart)):

            if int(t.GenPart_pdgId[j]) == 25:

                #get daughters' indices
                child_indices = [k for k in range(int(t.nGenPart))
                                 if int(t.GenPart_genPartIdxMother[k]) == j]

                #work on particles pair decays
                if len(child_indices) > 1:
                    for k in sorted(child_indices)[:2]:
                        daughters.append(pdg_to_Name(int(t.GenPart_pdgId[k])))
                        kin.append([t.GenPart_pt[k], t.GenPart_eta[k],
                                    t.GenPart_phi[k], t.GenPart_mass[k]])

                    print(f"daughters = {daughters} // kinematic info = {kin}")
                    break

        if len(kin) >= 2:
            decay = " + ".join(daughters)
            d1_vec = build_4vec(*kin[0])
            d2_vec = build_4vec(*kin[1])


        #------------------------------------------------------- add to df ------------------------------------------------------

        row = {
               "event": i,
               "lep_vec": lep_vec,
               "antilep_vec": antilep_vec,
               "d1_vec": d1_vec,
               "d2_vec": d2_vec,
               "h decay": decay,
               "event type": event_type,
               "VB decay": VB_decay
              }

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
    df.to_pickle("events.pkl")

    print("Showing Pandas dataframe")
    print(df.head())

    rdf = expand_columns_vectors(df)
    print("Showing Expanded Pandas dataframe")
    print(rdf.head())

    #Pd dataframe --> ROOT df
    root_df = ROOT.RDF.FromPandas(rdf)
    root_df.Snapshot("Events", "file.root")

if __name__ == "__main__":
    main()