from fastapi import FastAPI
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

app = FastAPI()

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

@app.post("/predict")
async def predict(smiles: str):
    preds = {name: model.predict([mol_to_fp(smiles)])[0] for name, model in loaded_models.items()}
    return preds
