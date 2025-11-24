import ROOT 

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
    vec = ROOT.TLorentzVector() 
    vec.SetPtEtaPhiM(pt, eta, phi, mass) 
    return vec


def read_file(filename, maxEvents=50):

    file = ROOT.TFile.Open(filename)    
    t = file.Get("Events")                      
    if not t: 
        print("No Events tree in", filename)
        return pd.DataFrame()

    n_entries = t.GetEntries()
    print(f"Number of events in file: {n_entries}")

    events_tree = ROOT.TTree("Events", "Final state particles")

    lep_4vec           = ROOT.TLorentzVector()
    antiL_4vec         = ROOT.TLorentzVector()
    h_daughter_1_4vec  = ROOT.TLorentzVector()
    h_daughter_2_4vec  = ROOT.TLorentzVector()
    event_type         = ROOT.std.string()
    VB_decay           = ROOT.std.string()
    decay_branch       = ROOT.std.string()

    events_tree.Branch("lepton",        lep_4vec)
    events_tree.Branch("anti_lepton",   antiL_4vec)
    events_tree.Branch("daughter_1",    h_daughter_1_4vec)
    events_tree.Branch("daughter_2",    h_daughter_2_4vec)
    events_tree.Branch("event_type",    event_type)
    events_tree.Branch("VB_decay",      VB_decay)
    events_tree.Branch("decay_branch",  decay_branch)

    for i in range(min(n_entries, maxEvents)):

        print(f"----------- event {i} -----------")

        t.GetEntry(i)

        lep_4vec.SetPtEtaPhiM(0,0,0,0)
        antiL_4vec.SetPtEtaPhiM(0,0,0,0)
        h_daughter_1_4vec.SetPtEtaPhiM(0,0,0,0)
        h_daughter_2_4vec.SetPtEtaPhiM(0,0,0,0)
        event_type.replace(0, len(event_type), "Z-W")
        VB_decay.replace(0, len(VB_decay), "")
        decay_branch.replace(0, len(decay_branch), "No Higgs")

        #LHE process 
        for j in range(int(t.nLHEPart)):
            pid = int(t.LHEPart_pdgId[j])       #get particle ID
            status = int(t.LHEPart_status[j])   #get particle status
            
            print(pid, status)

        #--------------------------------------------------- h decay - GenPart ---------------------------------------------------

        daughters, kin_info = [], []
        for j in range(int(t.nGenPart)):
            if int(t.GenPart_pdgId[j]) == 25:
                child_indices = [k for k in range(int(t.nGenPart)) if int(t.GenPart_genPartIdxMother[k]) == j]

                if len(child_indices) > 1:

                    for k in sorted(child_indices)[:2]:
                        daughters.append(pdg_to_Name(int(t.GenPart_pdgId[k])))
                        kin_info.append([t.GenPart_pt[k], t.GenPart_eta[k], t.GenPart_phi[k], t.GenPart_mass[k]])
                    print(f"daughters = {daughters} // kinematic info = {kin_info}")
                    
                    break

        if len(kin_info) >= 2:
            decay_branch.replace(0, len(decay_branch), " + ".join(daughters))
            h_daughter_1_4vec.SetPtEtaPhiM(*kin_info[0])
            h_daughter_2_4vec.SetPtEtaPhiM(*kin_info[1])

        events_tree.Fill()

    return events_tree
    

def main():
    
    tree = read_file('nanoLatino_ZH_H_ZToLL_ZL__part0.root', maxEvents=50) 
    rdf = ROOT.RDataFrame(tree)

if __name__ == "__main__":
    main()