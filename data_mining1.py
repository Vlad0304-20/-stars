import time
import pandas as pd
from Bio import Entrez

Entrez.email = "your.email@example.com"
search_term = "vaccinia"

try:
    with Entrez.esearch(db="pcassay", term=search_term, retmax=100) as handle:
        search_results = Entrez.read(handle)
        assay_ids = search_results.get("IdList", [])
except Exception:
    exit()

if not assay_ids:
    exit()

all_assays_data = []
chunk_size = 5 

for i in range(0, len(assay_ids), chunk_size):
    chunk = assay_ids[i:i + chunk_size]
    aid_list_str = ",".join(chunk)
    
    try:
        domain_parts = ["pubchem", "ncbi", "nlm", "nih", "gov"]
        base_host = "https://" + ".".join(domain_parts)
        pubchem_bulk_url = f"{base_host}/rest/pug/assay/aid/{aid_list_str}/concise/CSV"
        
        df_chunk = pd.read_csv(pubchem_bulk_url, low_memory=False)
        if not df_chunk.empty:
            all_assays_data.append(df_chunk)
    except Exception:
        pass
    
    time.sleep(1.5)

if not all_assays_data:
    exit()

df_raw_all = pd.concat(all_assays_data, ignore_index=True)
df_raw_all.columns = [str(col).upper() for col in df_raw_all.columns]

cid_col = 'PUBCHEM_CID'
aid_col = 'PUBCHEM_AID'
value_col = 'PUBCHEM_ACTIVITY_VALUE'

if cid_col not in df_raw_all.columns:
    cid_col = [col for col in df_raw_all.columns if 'CID' in col][0]
    aid_col = [col for col in df_raw_all.columns if 'AID' in col][0]
    value_col = [col for col in df_raw_all.columns if 'VALUE' in col][0]

extracted_rows = []
for idx, row in df_raw_all.iterrows():
    cid = row.get(cid_col)
    aid = row.get(aid_col)
    val = row.get(value_col)
    
    if pd.isna(cid) or cid == 0 or pd.isna(val): 
        continue
        
    try:
        float_val = float(val)
        row_str = str(row).upper()
        
        if 'CC50' in row_str or 'TC50' in row_str:
            standard_type = 'CC50'
        else:
            standard_type = 'IC50'
            
        extracted_rows.append({
            "molecule_id": int(cid),
            "assay_id": int(aid),
            "standard_type": standard_type,
            "standard_value": float_val
        })
    except ValueError:
        continue

df_activities = pd.DataFrame(extracted_rows)

if df_activities.empty:
    exit()

unique_mols = df_activities["molecule_id"].dropna().unique().tolist()
smiles_dfs = []
mol_chunk_size = 10 

for i in range(0, len(unique_mols), mol_chunk_size):
    chunk = [str(int(x)) for x in unique_mols[i:i + mol_chunk_size]]
    
    try:
        cid_list_str = ",".join(chunk)
        domain_parts = ["pubchem", "ncbi", "nlm", "nih", "gov"]
        base_host = "https://" + ".".join(domain_parts)
        pubchem_smiles_url = f"{base_host}/rest/pug/compound/cid/{cid_list_str}/property/CanonicalSMILES/CSV"
        
        df_chunk = pd.read_csv(pubchem_smiles_url)
        df_chunk = df_chunk.rename(columns={"CID": "molecule_id", "CanonicalSMILES": "SMILES"})
        smiles_dfs.append(df_chunk)
    except Exception:
        pass
        
    time.sleep(0.4)

df_smiles = pd.concat(smiles_dfs, ignore_index=True) if smiles_dfs else pd.DataFrame(columns=["molecule_id", "SMILES"])

unique_aids = df_activities["assay_id"].dropna().unique().tolist()
doi_mapping = {}

if unique_aids:
    aid_chunk_size = 20
    for i in range(0, len(unique_aids), aid_chunk_size):
        chunk_aids = unique_aids[i:i + aid_chunk_size]
        aid_str = ",".join([str(x) for x in chunk_aids])
        try:
            with Entrez.esummary(db="pcassay", id=aid_str) as handle:
                summaries = Entrez.read(handle)
                for summary in summaries:
                    aid_val = int(summary["Id"])
                    comment = summary.get("Comment", [])
                    pmid = None
                    for item in comment:
                        if "PMID" in item:
                            pmid = item.replace("PMID:", "").strip()
                            break
                    
                    if pmid and pmid.isdigit():
                        try:
                            with Entrez.esummary(db="pubmed", id=pmid) as p_handle:
                                p_summary = Entrez.read(p_handle)
                                article_ids = p_summary[0].get("ArticleIds", {})
                                doi_mapping[aid_val] = article_ids.get("doi", None)
                        except:
                            pass
            time.sleep(1.0)
        except Exception:
            pass

df_activities["doi_работа"] = df_activities["assay_id"].map(doi_mapping)

df_activities.columns = [str(col).upper() for col in df_activities.columns]
df_smiles.columns = [str(col).upper() for col in df_smiles.columns]

df_merged = pd.merge(df_activities, df_smiles, on="MOLECULE_ID", how="inner")
smiles_col = [col for col in df_merged.columns if 'SMILES' in col][0]

df_pivot = df_merged.pivot_table(
    index=['MOLECULE_ID', smiles_col],
    columns='STANDARD_TYPE',
    values='STANDARD_VALUE',
    aggfunc='mean'
).reset_index()

df_pivot = df_pivot.rename(columns={smiles_col: 'SMILES'})

if 'IC50' not in df_pivot.columns: df_pivot['IC50'] = None
if 'CC50' not in df_pivot.columns: df_pivot['CC50'] = None

df_pivot['IC50'] = pd.to_numeric(df_pivot['IC50'], errors='coerce')
df_pivot['CC50'] = pd.to_numeric(df_pivot['CC50'], errors='coerce')
df_pivot['SI'] = df_pivot['CC50'] / df_pivot['IC50']

if 'DOI_РАБОТА' in df_merged.columns:
    df_dois = df_merged.dropna(subset=['DOI_РАБОТА']).groupby('MOLECULE_ID')['DOI_РАБОТА'].first().reset_index()
    df_output = pd.merge(df_pivot, df_dois, on='MOLECULE_ID', how='left')
    df_output = df_output.rename(columns={'DOI_РАБОТА': 'doi_работа'})
else:
    df_output = df_pivot
    df_output['doi_работа'] = None

output_file = "collected_vaccinia_dataset2.csv"
df_output.to_csv(output_file, index=False)
