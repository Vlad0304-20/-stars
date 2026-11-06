import pandas as pd
from chembl_webresource_client.new_client import new_client

# Я думаю комментарии сделают мой код чуть понятнее
# Ищу все эксперементы и тесты где упоминаются vacciania оспа
assay = new_client.assay
assay_query = assay.search('vaccinia')
assays_df = pd.DataFrame(list(assay_query))

if not assays_df.empty:
    assay_ids = assays_df['assay_chembl_id'].unique().tolist()
    
    # Ищем запрашиваем idшки, небольшими объемами чтобы не перегружать сервера и тп,
    activity = new_client.activity
    activities_list = []
    chunk_size = 500
    for i in range(0, len(assay_ids), chunk_size):
        chunk = assay_ids[i:i + chunk_size]
        query = activity.filter(
            assay_chembl_id__in=chunk,
            standard_type__in=['IC50', 'CC50']
        ).only([
            'molecule_chembl_id', 'standard_type', 'standard_value', 'standard_units', 'document_chembl_id'
        ])
        activities_list.extend(list(query))
        
    df_raw = pd.DataFrame(activities_list)
else:
    df_raw = pd.DataFrame()

# id исследований и статей
if not df_raw.empty:
    doc_ids = df_raw['document_chembl_id'].dropna().unique().tolist()
    documents = new_client.document.filter(document_chembl_id__in=doc_ids).only(['document_chembl_id', 'doi'])
    df_docs = pd.DataFrame(list(documents))

    df_merged = pd.merge(df_raw, df_docs, on='document_chembl_id', how='left')

    # Склеиваем всё

    mol_ids = df_merged['molecule_chembl_id'].dropna().unique().tolist()

    # Качаем структуры молекул

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

    # Дальше лень описывать ну тут калькуляция происходит и финальный результат в табличке
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

    df_dois = df_final.dropna(subset=['doi']).groupby('molecule_chembl_id')['doi'].first().reset_index()
    df_pivot = pd.merge(df_pivot, df_dois, on='molecule_chembl_id', how='left')

    df_pivot = df_pivot.rename(columns={'doi': 'doi_работа'})
    df_output = df_pivot.dropna(subset=['SI', 'SMILES'])

    output_file = "collected_vaccinia_dataset.csv"
    df_output.to_csv(output_file, index=False)
