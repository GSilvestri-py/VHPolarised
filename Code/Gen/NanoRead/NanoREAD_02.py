# FILTER HIGGS DECAY

import glob
import os
import ROOT
import sys

def main():

    if len(sys.argv) < 2:
        print("use example: python3 script.py <origin_folder_path>")
        sys.exit(1)

    folder = sys.argv[1] 
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


    #define selections on h decay
    rdf = rdf.Define("H_bbar",  '(h_decay == "b + bbar") || (h_decay == "bbar + b")') \
            .Define("H_ditau",  '(h_decay == "tau+ + tau-") || (h_decay == "tau- + tau+")') \
            .Define("H_digamma",'h_decay == "γ + γ"')

    #filter on h decay
    df_bbar    = rdf.Filter("H_bbar")
    print("\n------------ Displaying b bbar df ------------ ")
    df_bbar.Display().Print()

    df_ditau   = rdf.Filter("H_ditau") 
    print("\n------------ Displaying tau+ + tau- df ------------ ")
    df_ditau.Display().Print()

    df_digamma = rdf.Filter("H_digamma")
    print("\n------------ Displaying γ + γ df ------------ ")
    df_digamma.Display().Print()

    print(f"Total b + bbar events count = {df_bbar.Count().GetValue()}")
    print(f"Total tau+ + tau- events count = {df_ditau.Count().GetValue()}")
    print(f"Total γ + γ events count = {df_digamma.Count().GetValue()}")


if __name__ == "__main__":
    main()