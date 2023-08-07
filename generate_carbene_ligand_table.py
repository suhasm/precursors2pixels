import pandas as pd
from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdqueries

halides_file_name = 'aromatic_halides_with_id.csv'
imidazoles_file_name = 'imidazoles_with_id.csv'

halide_df = pd.read_csv('input_data/' + halides_file_name)
imidazoles_df = pd.read_csv('input_data/' + imidazoles_file_name)

def smiles2mols(smiles_list):
    mols_list = []
    for smiles in smiles_list:
        mols_list.append(Chem.MolFromSmiles(smiles))
    return mols_list

def n_arylation(imidazole_smiles, halide_smiles):
    # Create RDKit molecule objects
    imidazole = Chem.MolFromSmiles(imidazole_smiles)
    halide = Chem.MolFromSmiles(halide_smiles)

    # Define n-arylation coupling reaction using reaction SMARTS
    # The carbene atom has to be non-aromatic. Otherwise you get a kekulize error.
    # Furthermore, it needs to be marked as having no Hs, or it will have an implicit H
    n_arylation_rxn_smarts = '[cH:1][c:2](-[I,Br,Cl,F:3])[c:4].[nX2:5][cH1:6][n:7]>>[c:1](-*):[cH0:2](-[N:5]-[CH0:6](->*)-[N:7]):[c:4].[I,Br,Cl,F:3]'
    n_arylation_rxn = AllChem.ReactionFromSmarts(n_arylation_rxn_smarts)

    # Perform the reaction
    product_set = n_arylation_rxn.RunReactants((halide, imidazole))

    # Collect the products and convert them to SMILES
    products = []
    for product in product_set:
        ligand = product[0]
        # Sanitize before adding. Errors indicate violations of the aromaticity model
        try:
            AllChem.SanitizeMol(ligand)
            Chem.MolToSmiles(ligand, isomericSmiles=False)
            products.append(ligand)
        except AllChem.KekulizeException:
            print(f"Could not sanitize combination of:\n{imidazole_smiles}\nwith\n{halide_smiles}")

    return set(products)

if __name__ == '__main__':

    result_data = []  # List to store the data

    for _, halide_row in tqdm(halide_df.iterrows(), total=len(halide_df)):
        for _, imidazole_row in imidazoles_df.iterrows():
            ligand_set = n_arylation(imidazole_row["imidazole_SMILES"], halide_row["halide_SMILES"])

            # If the set is empty, do nothing
            if not ligand_set:
                continue

            # Convert each ligand in the set to a row in the DataFrame and append to the list
            for i, ligand in enumerate(ligand_set):
                ligand_identifier = halide_row["halide_identifier"] + imidazole_row["imidazole_identifier"] + 'L' + hex(i)[2:]
                result_data.append({
                    "ligand_identifier": ligand_identifier,
                    "ligand_SMILES": ligand,
                    "halide_identifier": halide_row["halide_identifier"],
                    "halide_SMILES": halide_row["halide_SMILES"],
                    "imidazole_identifier": imidazole_row["imidazole_identifier"],
                    "imidazole_SMILES": imidazole_row["imidazole_SMILES"]
                })

    # Create the DataFrame using the list of data
    result_df = pd.DataFrame(result_data)

    result_df.to_csv("output_data/combinatorial_carbene_ligands.csv.gz", index=False)
