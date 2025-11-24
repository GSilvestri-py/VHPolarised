
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

// PDG to particle name
string pdgName(int pdgId) {
    static map<int, string> pdgNames = {
        {11,"e-"},{-11,"e+"},{13,"mu-"},{-13,"mu+"},{15,"tau-"},{-15,"tau+"},
        {12,"ve"},{-12,"~ve"},{14,"vmu"},{-14,"~vmu"},{16,"vtau"},{-16,"~vtau"},
        {21,"g"},{22,"γ"},{23,"Z"},{24,"W+"},{-24,"W-"},{25,"H"},
        {6,"t"},{-6,"tbar"},{5,"b"},{-5,"bbar"},{1,"d"},{-1,"dbar"},
        {2,"u"},{-2,"ubar"},{3,"s"},{-3,"sbar"},{4,"c"},{-4,"cbar"}
    };
    auto it = pdgNames.find(pdgId);
    return (it!=pdgNames.end()) ? it->second : "PDG("+to_string(pdgId)+")";
}

void printCombinedDecay(const char *fname="NanoAOD.root", int maxEvents=5) {
    TFile *file = TFile::Open(fname);
    if (!file || file->IsZombie()) { cerr<<"Error opening "<<fname<<endl; return; }
    TTree *t = (TTree*)file->Get("Events");
    if (!t) { cerr<<"No TTree 'Events' in "<<fname<<endl; return; }

    // --- LHE branches ---
    Int_t nLHEPart; const int maxLHE=20;
    float LHEPart_pt[maxLHE], LHEPart_eta[maxLHE], LHEPart_phi[maxLHE], LHEPart_mass[maxLHE];
    int   LHEPart_pdgId[maxLHE], LHEPart_status[maxLHE];
    t->SetBranchAddress("nLHEPart",&nLHEPart);
    t->SetBranchAddress("LHEPart_pt",LHEPart_pt);
    t->SetBranchAddress("LHEPart_eta",LHEPart_eta);
    t->SetBranchAddress("LHEPart_phi",LHEPart_phi);
    t->SetBranchAddress("LHEPart_mass",LHEPart_mass);
    t->SetBranchAddress("LHEPart_pdgId",LHEPart_pdgId);
    t->SetBranchAddress("LHEPart_status",LHEPart_status);

    // --- GEN branches ---
    Int_t nGenPart; const int maxGen=2000;
    int   GenPart_pdgId[maxGen], GenPart_status[maxGen], GenPart_genPartIdxMother[maxGen];
    float GenPart_pt[maxGen], GenPart_eta[maxGen], GenPart_phi[maxGen], GenPart_mass[maxGen];
    t->SetBranchAddress("nGenPart",&nGenPart);
    t->SetBranchAddress("GenPart_pdgId",GenPart_pdgId);
    t->SetBranchAddress("GenPart_status",GenPart_status);
    t->SetBranchAddress("GenPart_genPartIdxMother",GenPart_genPartIdxMother);
    t->SetBranchAddress("GenPart_pt",GenPart_pt);
    t->SetBranchAddress("GenPart_eta",GenPart_eta);
    t->SetBranchAddress("GenPart_phi",GenPart_phi);
    t->SetBranchAddress("GenPart_mass",GenPart_mass);

    Long64_t nentries = t->GetEntries();
    int nToRead = min<int>(nentries,maxEvents);

    for(int i=0;i<nToRead;++i){
        t->GetEntry(i);
        cout<<"\n================ Event "<<i<<" ================\n";

        // --- LHE process ---
        vector<int> incoming, outgoing;
        for(int j=0;j<nLHEPart;++j){
            int pdg=LHEPart_pdgId[j], status=LHEPart_status[j];
            if(status<0) incoming.push_back(j);
            else if(status>0) outgoing.push_back(j);
        }
        cout<<"LHE process: ";
        for(size_t k=0;k<incoming.size();++k){
            cout<<pdgName(LHEPart_pdgId[incoming[k]]);
            if(k+1<incoming.size()) cout<<" + ";
        }
        cout<<" -> ";
        for(size_t k=0;k<outgoing.size();++k){
            cout<<pdgName(LHEPart_pdgId[outgoing[k]]);
            if(k+1<outgoing.size()) cout<<" + ";
        }
        cout<<endl;

        // --- GEN Higgs decay ---
        vector<int> higgsIndices;
        int nSafe=min(nGenPart,maxGen);
        for(int j=0;j<nSafe;++j)
            if(GenPart_pdgId[j]==25) higgsIndices.push_back(j);

        if(higgsIndices.empty()){ cout<<"No Higgs boson at GEN level.\n"; continue; }

        for(int idx:higgsIndices){
            vector<int> daughters;
            for(int k=0;k<nSafe;++k)
                if(GenPart_genPartIdxMother[k]==idx) daughters.push_back(k);

            if(daughters.size()<=1) continue;
            cout<<"GEN Higgs decay: H → ";
            for(size_t k=0;k<daughters.size();++k){
                cout<<pdgName(GenPart_pdgId[daughters[k]]);
                if(k+1<daughters.size()) cout<<" + ";
            }
            cout<<endl;
        }
    }
    file->Close();
}

void scanDecay(const char *folderPath="./") {
    for(const auto& entry:fs::directory_iterator(folderPath)){
        if(!entry.is_regular_file()) continue;
        string fname=entry.path().string();
        if(fname.find(".root")==string::npos) continue;
        cout<<"\nProcessing file: "<<fname<<endl;
        printCombinedDecay(fname.c_str(),3);
    }
}

