import sys
import joblib
import numpy as np
import warnings
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

# Suppress all warnings
warnings.filterwarnings("ignore")

# Redirect stderr to null to suppress error messages
sys.stderr = open('/dev/null', 'w')

def mol_to_fp(smiles, radius=2, nBits=1024):
    """Convert SMILES to Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(nBits, dtype=int)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    arr = np.zeros((nBits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

# Load models without printing any warnings
try:
    names = ["lip", "sol"]
    loaded_models = {name: joblib.load(f"saved_models/{name}.pkl") for name in names}
except Exception:
    sys.exit(1)  # Exit silently if loading fails

def main():
    """Run predictions silently."""
    if len(sys.argv) < 2:
        sys.exit(1)  # Exit silently if no input

    smiles = sys.argv[1]
    
    try:
        preds = {name: model.predict([mol_to_fp(smiles)])[0] for name, model in loaded_models.items()}
        print(preds)  # Output should be JSON serializable
    except Exception:
        sys.exit(1)  # Exit silently on any prediction error

if __name__ == '__main__':
    main()
