# precursors2pixels

## Combinatorial Organic LED Ligand Generator

Collection of cheminformatics scripts to perform SMARTS reactions between a set of organic precursors to generate ligands for Ir(III) based phosphorescent LEDs.
## Prerequisites

This script requires the following Python packages:

- pandas
- tqdm
- rdkit
- mols2grid

## Input Data

We demonstrate its use through carbene-based ligand generation. The script requires two CSV files containing the SMILES strings and identifiers of the halides and imidazoles to be used in the reaction. These files should be named `aromatic_halides_with_id.csv` and `imidazoles_with_id.csv` respectively, and should be placed in the `input_data/` directory.

The CSV files should have the following structure:

| halide_identifier | halide_SMILES |
|-------------------|---------------|
| id1               | SMILES1       |
| id2               | SMILES2       |
| ...               | ...           |

and

| imidazole_identifier | imidazole_SMILES |
|----------------------|------------------|
| id1                  | SMILES1          |
| id2                  | SMILES2          |
| ...                  | ...              |

## Output Data

The script writes the resulting carbene ligands to a CSV file named `combinatorial_carbene_ligands.csv.gz` in the `output_data/` directory. This file includes the identifiers and SMILES strings of the halides and imidazoles used to generate each ligand, as well as a unique identifier and the SMILES string for each ligand.

The output CSV file has the following structure:

| ligand_identifier | ligand_SMILES | halide_identifier | halide_SMILES | imidazole_identifier | imidazole_SMILES |
|-------------------|---------------|-------------------|---------------|----------------------|------------------|
| id1               | SMILES1       | id2               | SMILES2       | id3                  | SMILES3          |
| ...               | ...           | ...               | ...           | ...                  | ...              |

## Usage

You can run the script using Python 3 as follows:

```bash
python generate_carbene_ligand_table.py
