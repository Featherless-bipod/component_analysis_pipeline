import pandas as pd
import gseapy as gp
import scanpy as sc
from cnmf import cNMF
import os
import sys
from tqdm.auto import tqdm

print("--- Starting cNMF and ORA pipeline ---")


print("Loading cNMF results...")
cnmf_obj = cNMF(output_dir="cnmf_run", name="4050run") 
chosen_k = 46
density_threshold = 0.1

try:
    print(f"Running consensus for k={chosen_k}...")

    cnmf_obj.consensus(k=chosen_k, density_threshold=density_threshold)

    print("Loading final spectra scores...")
    usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
        K=chosen_k, 
        density_threshold=density_threshold
    )
except Exception as e:
    print(f"Error during cNMF consensus step: {e}")
    sys.exit(1) # Exit if cNMF fails

print("Successfully loaded cNMF programs.")

# --- 2. Run ORA (Enrichr) Loop ---
print(f"Starting ORA for all {chosen_k} programs...")
gene_sets_db = "KEGG_2021_Human"
#gene_sets_db = 'GO_Biological_Process_2025' 
#gene_sets_db = 'MSigDB_Hallmark_2020'
output_folder = "pathway_analysis_results"
os.makedirs(output_folder, exist_ok=True)
print(f"Results will be saved to: {output_folder}")

N_TOP_GENES = 300 
all_ora_results = {}

for i in tqdm(range(1, chosen_k+1 )):
    print(f"--- Running ORA on Program {i} of {chosen_k} ---")
    
    program_gene_list = top_genes[i].head(N_TOP_GENES).tolist()
    
    # 2. Run enrichr (this is the ORA function)
    try:
        enr_results = gp.enrichr(
            gene_list=program_gene_list,
            gene_sets=gene_sets_db,
            outdir=None, # Keep results in memory
            cutoff=0.05  # Only keep results with P-value < 0.05
        )
        
        if not enr_results.res2d.empty:
            print(f"Program {i} - Top hit: {enr_results.res2d.iloc[0]['Term']}")
            all_ora_results[f'program_{i}'] = enr_results.res2d
        else:
            print(f"Program {i} - No significant pathways found.")
            
    except Exception as e:
        print(f"Error processing Program {i}: {e}")


# --- 3. Save a Summary File ---
print("\n\n✅ ORA complete! Saving summary...")

top_pathways = []
for program_name, results_df in all_ora_results.items():
    if not results_df.empty:
        top_hit = results_df.iloc[0]
        top_pathways.append({
            'Program': program_name,
            'Top_Pathway': top_hit['Term'],
            'Adj_P-value': top_hit['Adjusted P-value'],
            'Combined_Score': top_hit['Combined Score'],
            'Genes': top_hit['Genes']
        })

summary_df = pd.DataFrame(top_pathways)
summary_df = summary_df.sort_values('Adj_P-value')

summary_file = os.path.join(output_folder, "ORA_summary_all_programs.csv")
summary_df.to_csv(summary_file, index=False)

print(f"Summary saved to {summary_file}")
print("--- Pipeline Complete ---")