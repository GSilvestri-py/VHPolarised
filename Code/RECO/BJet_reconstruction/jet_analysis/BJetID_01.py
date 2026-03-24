import os
import pandas as pd


def main():
    
    pickles_path = "pickles_file_path"
    
    dfs = []

    for file in os.listdir(pickles_path):
        if file.endswith(".pkl"):
            full_path = os.path.join(pickles_path, file)
            df = pd.read_pickle(full_path)
            dfs.append(df)

    tot_df = pd.concat(dfs, ignore_index=True)

    first_col = tot_df.columns[0]
    tot_df[first_col] = range(len(tot_df))
    full_path = os.path.join(pickles_path, f"total_b_vs_nonb_df.pkl")
    tot_df.to_pickle(full_path)

    print(f"Number of rows (len): {len(tot_df)}")
    print(f"DataFrame shape: {tot_df.shape}\n")

    print(tot_df.head(20))
    
if __name__ == "__main__":
    main()