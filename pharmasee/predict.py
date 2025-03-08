# predict.py
import sys
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def mol_to_fp(smiles, radius=2, nBits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(nBits, dtype=int)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    arr = np.zeros((nBits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


names = ["lip", "sol"]
loaded_models = {name: joblib.load(f"saved_models/{name}.pkl") for name in names}


def main():
    if len(sys.argv) < 2:
        print("No SMILES string provided.")
        sys.exit(1)
    smiles = sys.argv[1]

# Load all models

    for name, model in loaded_models.items():
        preds = model.predict(smiles)



    print(f"Predicted ADMET properties for SMILES '{smiles}'")

if __name__ == '__main__':
    main()
