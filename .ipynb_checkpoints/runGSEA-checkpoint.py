import sys
import scanpy as sc
import pandas as pd
from cnmf import cNMF
import numpy as np
import tqdm as tqdm
import gseapy as gp
from tqdm.auto import tqdm

def summarize(all_gsea_results):
    top_pathways = []

    for program_name, results_df in all_gsea_results.items():
        if not results_df.empty:
            # Get the first row (the top-ranked pathway)
            top_hit = results_df.iloc[0]
            
            top_pathways.append({
                'Program': program_name,
                'Top_Pathway': top_hit['Term'],
                'NES': top_hit['NES'],         # Normalized Enrichment Score
                'FDR_q_val': top_hit['FDR q-val'] # Adjusted P-value
            })

    # Create a DataFrame from the list
    summary_df = pd.DataFrame(top_pathways)

    # Sort by the False Discovery Rate (most significant first)
    summary_df = summary_df.sort_values('FDR_q_val')

    return summary_df


def main():
    print("skibidi")
    cnmf_obj = cNMF(output_dir="cnmf_run", name="4050run")
    cnmf_obj.consensus(k=46, density_threshold=0.1)
    usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(K=46, density_threshold=0.1)
    print("1")

    my_program = 46
    all_gsea_results = {}
    gene_sets_db = "GO_Biolpogical_Process_2023"
    #gene_sets_db = "KEGG_2021_Human"
    #gene_sets_db = 'MSigDB_Hallmark_2020'


    program_ranked_list = spectra_scores[my_program].sort_values(ascending=False)
    prerank_df = pd.DataFrame(program_ranked_list)
    prerank_df.reset_index(inplace=True)
    prerank_df.columns = ['gene', 'score']

    print("2")

    for i in tqdm(range(1,my_program+1)):
        print(f"\n--- Running Program {i} of {my_program} ---")
        
        # 1. Get the ranked list for program 'i'
        program_ranked_list = spectra_scores[i].sort_values(ascending=False)
        
        # 2. Format it for gseapy
        prerank_df = pd.DataFrame(program_ranked_list)
        prerank_df.reset_index(inplace=True)
        prerank_df.columns = ['gene', 'score']
        
        # 3. Run Prerank
        try:
            prerank_results = gp.prerank(
                rnk=prerank_df,
                gene_sets=gene_sets_db,
                min_size=15,
                max_size=500,
                threads = -1,
                permutations = 1000,
                outdir=None,
                seed = 6,
                verbose = True,
            )
            
            # 4. Store the results in our dictionary
            if not prerank_results.res2d.empty:
                print(f"Program {i} - Top hit: {prerank_results.res2d.iloc[0]['Term']}")
                all_gsea_results[f'program_{i}'] = prerank_results.res2d
            else:
                print(f"Program {i} - No significant pathways found.")
                
        except Exception as e:
            print(f"Error processing Program {i}: {e}")
    
    print("Looping through results dictionary to build summary table...")

    print("\n\n✅ GSEA complete for all programs!")
    prerank_results.res2d.to_csv('GO_prerank_results_res2d.csv', index = False)
    print("successfully saved prerank_results!")
    summary_df = summarize(all_gsea_results)
    

    # --- And now you can save it ---
    print("Saving summary DataFrame to CSV...")
    summary_df.to_csv('GO_all_top_programs.csv', index=False)



main()