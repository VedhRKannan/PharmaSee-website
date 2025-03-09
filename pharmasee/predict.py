print("DEBUG: Script started")
import sys
import joblib
import numpy as np
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

print("DEBUG: Script started")  # ✅ See if Python runs at all

def mol_to_fp(smiles, radius=2, nBits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(json.dumps({"error": "Invalid SMILES string provided."}))
        sys.exit(1)

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    arr = np.zeros((nBits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.reshape(1, -1)  # Reshape for sklearn model input

# Try to load models, print errors if they fail
try:
    print("DEBUG: Loading models")  # ✅ See if model loading is slow
    names = ["lip", "sol"]
    loaded_models = {name: joblib.load(f"saved_models/{name}.pkl") for name in names}
    print("DEBUG: Models loaded successfully")  # ✅ Confirm success
except Exception as e:
    print(json.dumps({"error": f"Model loading failed: {str(e)}"}))
    sys.exit(1)

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"error": "No SMILES string provided."}))
        sys.exit(1)

    smiles = sys.argv[1]
    print(f"DEBUG: Received SMILES: {smiles}")  # ✅ Confirm input

    try:
        fingerprint = mol_to_fp(smiles)

        predictions = {}
        for name, model in loaded_models.items():
            preds = model.predict(fingerprint)
            predictions[name] = float(preds[0])

        print(json.dumps(predictions))  # ✅ Always return JSON
        sys.exit(0)  # ✅ Ensure script exits properly

    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)

if __name__ == '__main__':
    main()
