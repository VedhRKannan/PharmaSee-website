import sys
import joblib
import numpy as np
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def mol_to_fp(smiles, radius=2, nBits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(json.dumps({"error": "Invalid SMILES string provided."}))  # ✅ Ensure JSON format
        sys.exit(1)

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    arr = np.zeros((nBits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.reshape(1, -1)  # Reshape for sklearn model input

try:
    names = ["lip", "sol"]
    loaded_models = {name: joblib.load(f"saved_models/{name}.pkl") for name in names}
except Exception as e:
    print(json.dumps({"error": f"Model loading failed: {str(e)}"}))  # ✅ Ensure JSON error message
    sys.exit(1)

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"error": "No SMILES string provided."}))  # ✅ Ensure JSON format
        sys.exit(1)

    smiles = sys.argv[1]

    try:
        fingerprint = mol_to_fp(smiles)

        predictions = {}
        for name, model in loaded_models.items():
            preds = model.predict(fingerprint)
            predictions[name] = float(preds[0])  # ✅ Convert numpy float to standard float

        print(json.dumps(predictions))  # ✅ Ensure valid JSON output
        sys.exit(0)  # ✅ Exit successfully

    except Exception as e:
        print(json.dumps({"error": str(e)}))  # ✅ Catch all errors and print JSON
        sys.exit(1)

if __name__ == '__main__':
    main()
