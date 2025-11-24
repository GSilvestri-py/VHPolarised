#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>

// usage example: root -l -b -q 'scanGenDecay.C(folder/path)'
// first event of file.root


using namespace std;using namespace std;
namespace fs = std::filesystem;

void printGenDecay(const char *fname, int maxEvents = 1);

void scanGenDecay(const char *folderPath = "./") {
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (!entry.is_regular_file()) continue;
        string fname = entry.path().string();
        if (fname.find(".root") == string::npos) continue;

        cout << "\n==============================\n";
        cout << "Processing file: " << fname << endl;
        printGenDecay(fname.c_str(), 1); 
    }
}


//--------------------------------------------------
// Helper: map PDG IDs to human-readable particle names
//--------------------------------------------------
string pdgName(int pdgId) {
    static map<int, string> pdgNames = {
        {11, "e-"}, {-11, "e+"}, {13, "mu-"}, {-13, "mu+"}, {15, "tau-"}, {-15, "tau+"},
        {12, "ve"}, {-12, "~ve"}, {14, "vmu"}, {-14, "~vmu"}, {16, "vtau"}, {-16, "~vtau"},
        {21, "g"}, {22, "γ"}, {23, "Z"}, {24, "W+"}, {-24, "W-"}, {25, "H"},
        {6, "t"}, {-6, "tbar"}, {5, "b"}, {-5, "bbar"}, {1, "d"}, {-1, "dbar"},
        {2, "u"}, {-2, "ubar"}, {3, "s"}, {-3, "sbar"}, {4, "c"}, {-4, "cbar"},
    };
    auto it = pdgNames.find(pdgId);
    return (it != pdgNames.end()) ? it->second : "PDG(" + to_string(pdgId) + ")";
}

//--------------------------------------------------
// Main: Read GenPart info and print Higgs decays
//--------------------------------------------------
void printGenDecay(const char *fname = "NanoAOD.root", int maxEvents = 5) {
    cout << "Opening file: " << fname << endl;
    TFile *file = TFile::Open(fname);
    if (!file || file->IsZombie()) {
        cerr << "Error: could not open file " << fname << endl;
        return;
    }

    TTree *t = (TTree*)file->Get("Events");
    if (!t) {
        cerr << "Error: could not find TTree 'Events' in file " << fname << endl;
        return;
    }

    if (!t->GetBranch("nGenPart") || !t->GetBranch("GenPart_pdgId")) {
        cerr << "Error: required GenPart branches not found." << endl;
        return;
    }

    Int_t nGenPart;
    const int maxN = 2000;
    int   GenPart_pdgId[maxN];
    int   GenPart_status[maxN];
    int   GenPart_genPartIdxMother[maxN];
    float GenPart_pt[maxN];
    float GenPart_eta[maxN];
    float GenPart_phi[maxN];
    float GenPart_mass[maxN];

    t->SetBranchAddress("nGenPart", &nGenPart);
    t->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    t->SetBranchAddress("GenPart_status", GenPart_status);
    t->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    t->SetBranchAddress("GenPart_pt", GenPart_pt);
    t->SetBranchAddress("GenPart_eta", GenPart_eta);
    t->SetBranchAddress("GenPart_phi", GenPart_phi);
    t->SetBranchAddress("GenPart_mass", GenPart_mass);

    Long64_t nentries = t->GetEntries();
    int nToRead = min<int>(nentries, maxEvents);

    cout << "\n=== Gen Higgs Decay Summary ===\n";

    for (int i = 0; i < nToRead; ++i) {
        t->GetEntry(i);

        if (nGenPart < 0 || nGenPart > maxN) {
            cerr << "Warning: nGenPart = " << nGenPart << " is invalid → skipping event " << i << endl;
            continue;
        }

        int nSafe = std::min(nGenPart, maxN);
        cout << "\n--- Event " << i << " ---\n";

        // Find all Higgs bosons
        vector<int> higgsIndices;
        for (int j = 0; j < nSafe; ++j) {
            if (GenPart_pdgId[j] == 25) higgsIndices.push_back(j);
        }

        if (higgsIndices.empty()) {
            cout << "No Higgs boson in this event.\n";
            continue;
        }

        for (int idx : higgsIndices) {
            // Collect its direct daughters
            vector<int> daughters;
            for (int k = 0; k < nSafe; ++k) {
                if (GenPart_genPartIdxMother[k] == idx)
                    daughters.push_back(k);
            }

            if (daughters.empty()) {
                cout << "  No daughters found (possibly intermediate copy)\n";
                continue;
            }

            if (daughters.size() == 1) continue; // intermediate Higgs

            cout << "Higgs found: index " << idx
                 << ", pt=" << GenPart_pt[idx]
                 << ", eta=" << GenPart_eta[idx] << "\n";

            cout << "  H → ";
            for (size_t k = 0; k < daughters.size(); ++k) {
                int pdg = GenPart_pdgId[daughters[k]];
                cout << pdgName(pdg);
                if (k + 1 < daughters.size()) cout << " + ";
            }
            cout << endl;
        }
    }

    file->Close();
}