import sys
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

# Function to convert SMILES to Morgan Fingerprint
def mol_to_fp(smiles, radius=2, nBits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES string provided.")
        sys.exit(1)
    
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    arr = np.zeros((nBits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.reshape(1, -1)  # Reshape for sklearn model input

# Load all models
names = ["lip", "sol"]
loaded_models = {name: joblib.load(f"pharmasee/saved_models/{name}.pkl") for name in names}

def main():
    if len(sys.argv) < 2:
        print("No SMILES string provided.")
        sys.exit(1)

    smiles = sys.argv[1]
    
    # Convert SMILES to fingerprint
    fingerprint = mol_to_fp(smiles)

    predictions = {}
    
    for name, model in loaded_models.items():
        preds = model.predict(fingerprint)  # Pass the feature vector
        predictions[name] = preds[0]  # Assuming it's a single-value output

    # Print JSON-style output for easy parsing
    print(predictions)

if __name__ == '__main__':
    main()
