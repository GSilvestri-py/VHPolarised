import awkward as ak
import uproot              
import numpy as np
import vector
import ROOT

def build_4vec_from_pt_eta_phi_m(pt, eta, phi, mass):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    p  = pt * np.cosh(eta)
    E  = np.sqrt(mass**2 + p**2)
    return vector.obj(px = px, py = py, pz =pz, E = E)

def compute_deltaR(vec1, vec2):

    deltaphi = vec1.deltaphi(vec2)
    deltaeta = abs(vec1.eta - vec2.eta)

    return np.sqrt( (deltaphi)**2 + (deltaeta)**2 )


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


def main():

    filename = "/gwpool/users/gisilvestri/VHPolar/RECO/nanoLatino_ZH_H_ZToLL_ZL__part0.root"

    with uproot.open(filename) as f:                                       

        tree = f["Events"]

        n_lhe      = tree["nLHEPart"].array(library="ak")
        lhe_pdg    = tree["LHEPart_pdgId"].array(library="ak")
        lhe_status = tree["LHEPart_status"].array(library="ak")

        Gen_pdg  = tree["GenPart_pdgId"].array(library="ak")
        Gen_midx = tree["GenPart_genPartIdxMother"].array(library="ak")
        Gen_pt   = tree["GenPart_pt"].array(library="ak")
        Gen_eta  = tree["GenPart_eta"].array(library="ak")
        Gen_phi  = tree["GenPart_phi"].array(library="ak")
        Gen_mass = tree["GenPart_mass"].array(library="ak")

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
        

        for i in range(130):                                                             #loop over events

            print(f"=============================================== EVENT {i} ===============================================")

            index_list = []                                                        
            e_counter = 0
            mu_counter = 0
            is_b_inevent = False

            #LHE array i component
            n_lhe_evt      = n_lhe[i]
            lhe_pdg_evt    = lhe_pdg[i]
            lhe_status_evt = lhe_status[i]

            #GenPart array i component
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

            print(f"\n------------------------------------------ process outline ------------------------------------------")
            father, daughters = [], []
            for j in range(n_lhe_evt):
                pid = lhe_pdg_evt[j]
                status = lhe_status_evt[j]

                if status == -1:
                    father.append(pid)
                if status == 1:
                    daughters.append(pid)
            
            
            def mother_pdg(j):
                m = int(midx_evt[j])
                return int(pdg_evt[m]) if m >= 0 else None

            daughters_h = []
            #clean e+/- mu+/- from GenJet
            for n in range(len(pdg_evt)):                                               #loop over particles in GenPart

                pid = int(pdg_evt[n])
                mother = mother_pdg(n)

                if abs(pid) == 15:
                    print(f"tau in event with mother {mother}")

                if is_b_hadron(abs(pid)):
                    is_b_inevent = True
                    print(f"quark b in event")

                if pid == 25:
                    child_indices = [k for k in range(len(pdg_evt)) if int(midx_evt[k]) == n]

                    if len(child_indices) > 1:
                        daughters_h.extend(int(pdg_evt[k]) for k in child_indices)
                        print(f"{father[0]} + {father[1]} --> " + " + ".join(map(str, daughters)) + ", 25 > " + " + ".join(map(str, daughters_h))  )
                        print("\n")


                if abs(pid) == 11 and mother in [23, 24, -24]:                          #electron/positron from Z/W
                    e_counter += 1
                    e_idx = [k for k in range(len(pdg_evt)) if pdg_evt[k] == pid]
                    e_vec = build_4vec_from_pt_eta_phi_m(pt_evt[e_idx[0]], eta_evt[e_idx[0]], phi_evt[e_idx[0]], mass_evt[e_idx[0]])

                    for j in range(n_jets_evt):                                         #loop over jets
                        jet_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[j], GJet_eta_evt[j], GJet_phi_evt[j], GJet_mass_evt[j])
                        deltaR_ej = compute_deltaR(e_vec, jet_vec)

                        if deltaR_ej < 0.1:
                            print(f"jet with index {j} is the e+/- jet: e j deltaR = {deltaR_ej}")
                            index_list.append(j)
                           


                if abs(pid) == 13 and mother in [23, 24, -24]:                          #muon from Z/W
                    mu_counter += 1
                    mu_idx = [k for k in range(len(pdg_evt)) if pdg_evt[k] == pid]
                    mu_vec = build_4vec_from_pt_eta_phi_m(pt_evt[mu_idx[0]], eta_evt[mu_idx[0]], 
                                                          phi_evt[mu_idx[0]], mass_evt[mu_idx[0]])

                    for j in range(n_jets_evt):                                         #loop over jets
                        jet_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[j], GJet_eta_evt[j], 
                                                               GJet_phi_evt[j], GJet_mass_evt[j])
                        deltaR_muj = compute_deltaR(mu_vec, jet_vec)

                        if deltaR_muj < 0.1:
                            print(f"jet with index {j} is the mu+/- jet: mu j deltaR = {deltaR_muj}")
                            index_list.append(j)
            
           
            #clean tau jets from GenJet
            for j_t in range(n_jets_vistau_evt):
                j_t_vec = build_4vec_from_pt_eta_phi_m(VisTau_pt_evt[j_t], VisTau_eta_evt[j_t],
                                                        VisTau_phi_evt[j_t], VisTau_mass_evt[j_t],)
                
                for j_g in range(n_jets_evt):
                    j_g_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[j_g], GJet_eta_evt[j_g], 
                                                            GJet_phi_evt[j_g], GJet_mass_evt[j_g])
                    deltaR_jj = j_t_vec.deltaR(j_g_vec)

                    if deltaR_jj < 0.1:
                        print(f"jet with index {j_g} is a tau jet: deltaR_jj = {deltaR_jj}")
                        index_list.append(j_g)


            print(f"\n------------------------------------------ before mask ------------------------------------------")
            print(f"electron, muon per event = {e_counter}, {mu_counter}")
            print(f"n_jets = {n_jets_evt}")
            print(f"GJet_pt = {GJet_pt_evt}")
            print(f"GJet_eta = {GJet_eta_evt}")
            print(f"GJet_phi = {GJet_phi_evt}")
            print(f"GJet_mass = {GJet_mass_evt}")
            print(f"GJet_h_flav = {GJet_h_flav_evt}")
            
            indices = ak.local_index(GJet_pt_evt)   
            mask = np.isin(indices, index_list)    
            
            GJet_pt_evt   = ak.where(mask, 0, GJet_pt_evt)
            GJet_eta_evt  = ak.where(mask, 0, GJet_eta_evt)
            GJet_phi_evt  = ak.where(mask, 0, GJet_phi_evt)
            GJet_mass_evt = ak.where(mask, 0, GJet_mass_evt)


            print(f"\n------------------------------------------ after mask ------------------------------------------")
            print(f"electron, muon per event = {e_counter}, {mu_counter}")
            print(f"n_jets = {n_jets_evt}")
            print(f"GJet_pt = {GJet_pt_evt}")
            print(f"GJet_eta = {GJet_eta_evt}")
            print(f"GJet_phi = {GJet_phi_evt}")
            print(f"GJet_mass = {GJet_mass_evt}")
            print(f"GJet_h_flav = {GJet_h_flav_evt}")

            print(f"n_jets_vistau_evt = {n_jets_vistau_evt}")
            print(f"VisTau_pt_evt = {VisTau_pt_evt}")
            print(f"VisTau_eta_evt = {VisTau_eta_evt}")
            print(f"VisTau_phi_evt = {VisTau_phi_evt}")
            print(f"VisTau_mass_evt = {VisTau_mass_evt}")
            print(f"VisTau_mother_evt = {VisTau_mother_evt}")
            print(f"VisTau_status_evt = {VisTau_status_evt}")
           
            print(f"\n------------------------------------------ b jets ------------------------------------------")
            
            #identify b quark jets
            jet_candidates = 0
            for e1 in range(n_jets_evt-1):                          #njets - 1 avoids duplicate combinations

                if GJet_pt_evt[e1] == 0:    continue                #avoid 'cleaned' jets
                
                j1_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[e1], GJet_eta_evt[e1], 
                                                               GJet_phi_evt[e1], GJet_mass_evt[e1])
                    
                for e2 in range(n_jets_evt):
                    if e2 <= e1:    continue                         #avoid duplicate pairs

                    print(f"jet pair indices = ({e1}, {e2})")
                    j2_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[e2], GJet_eta_evt[e2], 
                                                               GJet_phi_evt[e2], GJet_mass_evt[e2])

                    vec_jj = j1_vec + j2_vec
                    mass_jj = vec_jj.M
                    print(f"jj invariant mass = {mass_jj}")
                    tollerance = 15
                    if abs(mass_jj - 125) < tollerance:
                        jet_candidates += 1
                        print(f"> possible jet {e1}-{e2} from h decay b pair\n")

            b_from_h = False
            if jet_candidates > 0:
                b_from_h = True

            if b_from_h and not is_b_inevent:
                print(f"jets with expected invariant mass but no b in event")

            if is_b_inevent and not b_from_h:
                print(f"b in event but no candidate h > bb jets")

            print(f"\n------------------------------------------------------------------------------------------------------------------------------------------\n")

if __name__ == "__main__":
    main()

