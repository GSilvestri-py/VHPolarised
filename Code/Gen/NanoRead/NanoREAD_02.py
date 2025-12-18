# FILTER HIGGS DECAY

import glob
import os
import ROOT
import uproot

def main():

    folder = "/gwpool/users/gisilvestri/VHPolar/RECO/root_jet_df"
    files = sorted(glob.glob(os.path.join(folder, "*.root")))

    chain = ROOT.TChain("Events")
    for file in files:
        chain.Add(file)

    rdf = ROOT.RDataFrame(chain)
    print(f"Total events count: {rdf.Count().GetValue()}")

    print("\n------------ Dataframe columns ------------ ")
    cols = list(rdf.GetColumnNames())
    for i, col in enumerate(cols, start=1):
        print(f"{i:2d}. {col}")
    for col in rdf.GetColumnNames():
        print(col, rdf.GetColumnType(col))



    #define selections on h decay
    rdf = rdf.Define("H_bbar",  '(h_decay == "b + bbar") || (h_decay == "bbar + b")') \
            .Define("H_ditau",  '(h_decay == "tau+ + tau-") || (h_decay == "tau- + tau+")') \
            .Define("H_digamma",'h_decay == "γ + γ"')

    #filter on h decay
    df_bbar    = rdf.Filter("H_bbar")
    df_bbar.Snapshot("df_bbar", "h_bb.root")
    print("\n------------ Displaying b bbar df ------------ ")
    tree = uproot.open("h_bb.root")["df_bbar"]
    pdf = tree.arrays(library="pd")
    print(pdf.head())
    evt_fraction_bb = df_bbar.Count().GetValue() / rdf.Count().GetValue()

    df_ditau   = rdf.Filter("H_ditau") 
    print("\n------------ Displaying tau+ + tau- df ------------ ")
    df_ditau.Snapshot("df_ditau", "h_ditau.root")
    tree = uproot.open("h_ditau.root")["df_ditau"]
    pdf = tree.arrays(library="pd")
    print(pdf.head())

    evt_fraction_ditau = df_ditau.Count().GetValue() / rdf.Count().GetValue()

    df_digamma = rdf.Filter("H_digamma")
    print("\n------------ Displaying γ + γ df ------------ ")
    df_digamma.Snapshot("df_digamma", "df_digamma.root")
    tree = uproot.open("df_digamma.root")["df_digamma"]
    pdf = tree.arrays(library="pd")
    print(pdf.head())
    evt_fraction_digamma = df_digamma.Count().GetValue() / rdf.Count().GetValue()

    print(f"Total b + bbar events count = {df_bbar.Count().GetValue()}")
    print(f"Total tau+ + tau- events count = {df_ditau.Count().GetValue()}")
    print(f"Total γ + γ events count = {df_digamma.Count().GetValue()}")

    print(f"b + bbar event fraction = {evt_fraction_bb}")
    print(f"tau- + tau+ event fraction = {evt_fraction_ditau}")
    print(f"γ + γ event fraction = {evt_fraction_digamma}")


if __name__ == "__main__":
    main()