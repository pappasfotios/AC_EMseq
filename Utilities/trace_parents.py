import sys
import pandas as pd

def trace_parents(pedigree, individual_id, parent_type):
    generations = [individual_id]
    current_generation = 9
    while current_generation >= 0:
        if individual_id not in pedigree['Id'].values:
            break
        parent = pedigree[pedigree['Id'] == individual_id][f'{parent_type}'].values[0]
        generations.append(parent)
        individual_id = parent
        current_generation -= 1
    return generations

if len(sys.argv) != 2:
    print("Usage: python script.py <file_with_Tag_IDs>")
    sys.exit(1)

file_with_tag_ids = sys.argv[1]

pedigree_data = pd.read_csv('../ac_ped_updated_2023.txt', delimiter="\t")

pheno = pd.read_csv(file_with_tag_ids, delimiter="\t", header=None, names=['Tag_ID'])

ids_to_trace = [str(i) for i in pheno['Tag_ID']]

results_dam = []
results_sire = []

for individual_id in ids_to_trace:
    if individual_id in pedigree_data['Id'].values:
        dams_trace = trace_parents(pedigree_data, individual_id, 'Dam')
        sires_trace = trace_parents(pedigree_data, individual_id, 'Sire')
        results_dam.append({'ID': individual_id, 'Parents_Trace': dams_trace})
        results_sire.append({'ID': individual_id, 'Parents_Trace': sires_trace})
    else:
        print(f"ID {individual_id} not found in the pedigree data.")

results_dam_df = pd.DataFrame(results_dam)
results_sire_df = pd.DataFrame(results_sire)

num_generations = len(results_dam_df['Parents_Trace'].iloc[0])
columns = [f'Generation_{num_generations - i}' for i in range(num_generations)]
results_dam_df[columns] = pd.DataFrame(results_dam_df['Parents_Trace'].tolist(), columns=columns)
results_sire_df[columns] = pd.DataFrame(results_sire_df['Parents_Trace'].tolist(), columns=columns)

results_dam_df = results_dam_df.drop(columns=['Parents_Trace'])
results_sire_df = results_sire_df.drop(columns=['Parents_Trace'])

results_dam_df['ID'] = results_dam_df['ID'].astype(int)
results_sire_df['ID'] = results_sire_df['ID'].astype(int)

print(len(results_dam_df.index))
print(len(results_sire_df.index))
print(len(pheno.index))

merged_dam_df = results_dam_df.merge(pheno, left_on='ID', right_on='Tag_ID', how='inner')
merged_sire_df = results_sire_df.merge(pheno, left_on='ID', right_on='Tag_ID', how='inner')

merged_dam_df.iloc[:,1:11].to_csv(path_or_buf="./dam_traced.txt", sep="\t", index=False, na_rep="NA")
merged_sire_df.iloc[:,1:11].to_csv(path_or_buf="./sire_traced.txt", sep="\t", index=False, na_rep="NA")
