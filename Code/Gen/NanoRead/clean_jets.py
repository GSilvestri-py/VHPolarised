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

def main():

    filename = "/gwpool/users/gisilvestri/VHPolar/RECO/nanoLatino_ZH_H_ZToLL_ZL__part0.root"

    with uproot.open(filename) as f:                                       

        tree = f["Events"]

        Gen_pdg   = tree["GenPart_pdgId"].array(library="ak")
        Gen_midx  = tree["GenPart_genPartIdxMother"].array(library="ak")
        Gen_pt    = tree["GenPart_pt"].array(library="ak")
        Gen_eta   = tree["GenPart_eta"].array(library="ak")
        Gen_phi   = tree["GenPart_phi"].array(library="ak")
        Gen_mass  = tree["GenPart_mass"].array(library="ak")

        n_jets      = tree["nGenJet"].array(library="ak")
        GJet_pt     = tree["GenJet_pt"].array(library="ak")
        GJet_eta    = tree["GenJet_eta"].array(library="ak")
        GJet_phi    = tree["GenJet_phi"].array(library="ak")
        GJet_mass   = tree["GenJet_mass"].array(library="ak")
        GJet_h_flav = tree["GenJet_hadronFlavour"].array(library="ak")


        for i in range(10):

            pdg_evt   = Gen_pdg[i]
            midx_evt  = Gen_midx[i]
            pt_evt    = Gen_pt[i]
            eta_evt   = Gen_eta[i]
            phi_evt   = Gen_phi[i]
            mass_evt  = Gen_mass[i]
            
            n_jets_evt      = n_jets[i]
            GJet_pt_evt     = GJet_pt[i]
            GJet_eta_evt    = GJet_eta[i]
            GJet_phi_evt    = GJet_phi[i]
            GJet_mass_evt   = GJet_mass[i]
            GJet_h_flav_evt = GJet_h_flav[i]

            def mother_pdg(j):
                m = int(midx_evt[j])
                return int(pdg_evt[m]) if m >= 0 else None

            for n in range(len(pdg_evt)):

                pid = int(pdg_evt[n])
                mother = mother_pdg(n)

                if abs(pid) == 11 and mother in [23, 24, -24]:      #electron/positron
                    e_idx = [k for k in range(len(pdg_evt)) if pdg_evt[k] == pid]
                    e_vec = build_4vec_from_pt_eta_phi_m(pt_evt[e_idx[0]], eta_evt[e_idx[0]], phi_evt[e_idx[0]], mass_evt[e_idx[0]])

                    for j in range(0, n_jets_evt-1):
                        jet_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[j], GJet_eta_evt[j], GJet_phi_evt[j], GJet_mass_evt[j])
                        deltaR_ej = compute_deltaR(e_vec, jet_vec)
                        print(f"e j deltaR = {deltaR_ej}")

                        if deltaR_ej < 0.1:
                            print(f"jet in position {j} is the e+/- jet")
                            


                if abs(pid) == 13 and mother in [23, 24, -24]:      #muon
                    mu_idx = [k for k in range(len(pdg_evt)) if pdg_evt[k] == pid]
                    mu_vec = build_4vec_from_pt_eta_phi_m(pt_evt[mu_idx[0]], eta_evt[mu_idx[0]], phi_evt[mu_idx[0]], mass_evt[mu_idx[0]])

                    for j in range(0, n_jets_evt-1):
                        jet_vec = build_4vec_from_pt_eta_phi_m(GJet_pt_evt[j], GJet_eta_evt[j], GJet_phi_evt[j], GJet_mass_evt[j])
                        deltaR_muj = compute_deltaR(mu_vec, jet_vec)
                        print(f"mu j deltaR = {deltaR_muj}")

                        if deltaR_muj < 0.1:
                            print(f"jet in position {j} is the mu+/- jet")
                            

            print(f"\n---------------------------------------------- event {i} ----------------------------------------------")
            print(f"n_jets = {n_jets_evt}")
            print(f"GJet_pt = {GJet_pt_evt}")
            print(f"GJet_eta = {GJet_eta_evt}")
            print(f"GJet_phi = {GJet_phi_evt}")
            print(f"GJet_mass = {GJet_mass_evt}")
            print(f"GJet_h_flav = {GJet_h_flav_evt}")

if __name__ == "__main__":
    main()

