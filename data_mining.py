import time
import pandas as pd
from chembl_webresource_client.new_client import new_client
from Bio import Entrez
from io import StringIO

# === ЧАСТЬ 1: СБОР ИЗ BASE ChEMBL ===
assay = new_client.assay
assay_query = assay.search('vaccinia').only(['assay_chembl_id', 'description', 'assay_organism'])
assays_df = pd.DataFrame(list(assay_query))

df_chembl_out = pd.DataFrame()

if not assays_df.empty:
    assay_ids = assays_df['assay_chembl_id'].unique().tolist()
    activity = new_client.activity
    activities_list = []
    chunk_size = 500
    for i in range(0, len(assay_ids), chunk_size):
        chunk = assay_ids[i:i + chunk_size]
        query = activity.filter(
            assay_chembl_id__in=chunk,
            standard_type__in=['IC50', 'CC50']
        ).only([
            'molecule_chembl_id', 'standard_type', 'standard_value', 'standard_units', 'document_chembl_id', 'assay_chembl_id'
        ])
        activities_list.extend(list(query))
    df_raw = pd.DataFrame(activities_list)
else:
    df_raw = pd.DataFrame()

if not df_raw.empty:
    doc_ids = df_raw['document_chembl_id'].dropna().unique().tolist()
    documents = new_client.document.filter(document_chembl_id__in=doc_ids).only(['document_chembl_id', 'doi'])
    df_docs = pd.DataFrame(list(documents))

    df_merged = pd.merge(df_raw, df_docs, on='document_chembl_id', how='left')
    mol_ids = df_merged['molecule_chembl_id'].dropna().unique().tolist()

    mols_list = []
    chunk_size = 500
    for i in range(0, len(mol_ids), chunk_size):
        chunk = mol_ids[i:i + chunk_size]
        mols = new_client.molecule.filter(molecule_chembl_id__in=chunk).only(['molecule_chembl_id', 'molecule_structures'])
        mols_list.extend(list(mols))

    smiles_data = []
    for m in mols_list:
        chem_id = m.get('molecule_chembl_id')
        structs = m.get('molecule_structures')
        smiles = structs.get('canonical_smiles') if structs else None
        smiles_data.append({'molecule_chembl_id': chem_id, 'SMILES': smiles})

    df_smiles = pd.DataFrame(smiles_data)
    df_final = pd.merge(df_merged, df_smiles, on='molecule_chembl_id', how='left')

    assays_info = assays_df.rename(columns={'assay_organism': 'strain', 'description': 'culture'})
    df_final = pd.merge(df_final, assays_info, on='assay_chembl_id', how='left')
    df_final['standard_value'] = pd.to_numeric(df_final['standard_value'], errors='coerce')

    df_pivot = df_final.pivot_table(
        index=['molecule_chembl_id', 'SMILES'],
        columns='standard_type',
        values='standard_value',
        aggfunc='mean'
    ).reset_index()

    if 'IC50' not in df_pivot.columns: df_pivot['IC50'] = None
    if 'CC50' not in df_pivot.columns: df_pivot['CC50'] = None

    df_pivot['SI'] = df_pivot['CC50'] / df_pivot['IC50']

    df_strain = df_final.dropna(subset=['strain']).groupby('molecule_chembl_id')['strain'].first().reset_index()
    df_culture = df_final.dropna(subset=['culture']).groupby('molecule_chembl_id')['culture'].first().reset_index()
    df_dois = df_final.dropna(subset=['doi']).groupby('molecule_chembl_id')['doi'].first().reset_index()

    df_pivot = pd.merge(df_pivot, df_strain, on='molecule_chembl_id', how='left')
    df_pivot = pd.merge(df_pivot, df_culture, on='molecule_chembl_id', how='left')
    df_pivot = pd.merge(df_pivot, df_dois, on='molecule_chembl_id', how='left')

    df_pivot = df_pivot.rename(columns={'doi': 'doi_работа', 'molecule_chembl_id': 'molecule_id'})
    df_chembl_out = df_pivot.dropna(subset=['SI', 'SMILES'])
    df_chembl_out = df_chembl_out[['molecule_id', 'SMILES', 'strain', 'culture', 'CC50', 'IC50', 'SI', 'doi_работа']]

# === ЧАСТЬ 2: СБОР ИЗ PUBCHEM (BIOPYTHON) ===
Entrez.email = "suz29631@laoia.com"
search_term = "vaccinia"

try:
    with Entrez.esearch(db="pcassay", term=search_term, retmax=100) as handle:
        search_results = Entrez.read(handle)
        assay_ids = search_results.get("IdList", [])
except Exception:
    assay_ids = []

df_pubchem_out = pd.DataFrame()

if assay_ids:
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

    if all_assays_data:
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

        if not df_activities.empty:
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
            strain_mapping = {}
            culture_mapping = {}

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
                                name_text = str(summary.get("Name", ""))
                                description_text = " ".join(summary.get("Description", []))
                                strain_mapping[aid_val] = "Vaccinia virus" if "VACCINIA" in name_text.upper() else "Poxviridae"
                                
                                culture_found = "Vero cells"
                                if "HELA" in description_text.upper() or "HELA" in name_text.upper():
                                    culture_found = "HeLa cells"
                                elif "HEK293" in description_text.upper() or "HEK293" in name_text.upper():
                                    culture_found = "HEK293 cells"
                                elif "MDCK" in description_text.upper() or "MDCK" in name_text.upper():
                                    culture_found = "MDCK cells"
                                culture_mapping[aid_val] = culture_found

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
df_activities["strain"] = df_activities["assay_id"].map(strain_mapping)
df_activities["culture"] = df_activities["assay_id"].map(culture_mapping)
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
df_strain_final = df_merged.dropna(subset=['STRAIN']).groupby('MOLECULE_ID')['STRAIN'].first().reset_index()
df_culture_final = df_merged.dropna(subset=['CULTURE']).groupby('MOLECULE_ID')['CULTURE'].first().reset_index()

df_pivot = pd.merge(df_pivot, df_strain_final, on='MOLECULE_ID', how='left')
df_pivot = pd.merge(df_pivot, df_culture_final, on='MOLECULE_ID', how='left')

if 'DOI_РАБОТА' in df_merged.columns:
    df_dois = df_merged.dropna(subset=['DOI_РАБОТА']).groupby('MOLECULE_ID')['DOI_РАБОТА'].first().reset_index()
    df_output = pd.merge(df_pivot, df_dois, on='MOLECULE_ID', how='left')
    df_output = df_output.rename(columns={'DOI_РАБОТА': 'doi_работа', 'STRAIN': 'strain', 'CULTURE': 'culture', 'MOLECULE_ID': 'molecule_id'})
else:
    df_output = df_pivot
    df_output['doi_работа'] = None
    df_output = df_output.rename(columns={'STRAIN': 'strain', 'CULTURE': 'culture', 'MOLECULE_ID': 'molecule_id'})

df_output['molecule_id'] = "PUBCHEM_" + df_output['molecule_id'].astype(str)
df_pubchem_out = df_output.dropna(subset=['SI', 'SMILES'])
df_pubchem_out = df_pubchem_out[['molecule_id', 'SMILES', 'strain', 'culture', 'CC50', 'IC50', 'SI', 'doi_работа']]

# === ЧАСТЬ 3: СЛИЯНИЕ ДВУХ БАЗ ДАННЫХ ===
df_total_combined = pd.concat([df_chembl_out, df_pubchem_out], ignore_index=True)
df_final_clean = df_total_combined.drop_duplicates(subset=['SMILES'])

output_file = "collected_vaccinia_dataset.csv"
df_final_clean.to_csv(output_file, index=False)

