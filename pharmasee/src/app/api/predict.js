// pages/api/predict.js
import { PythonShell } from "python-shell";

export default async function handler(req, res) {
  if (req.method !== "POST") {
    return res.status(405).json({ error: "Method not allowed" });
  }

  const { smiles } = req.body;
  if (!smiles) {
    return res.status(400).json({ error: "No SMILES string provided." });
  }

  let options = {
    mode: "text",
    // Change this path if needed. Use `which python3` to get the correct path.
    pythonPath: "/opt/homebrew/bin/python3", 
    args: [smiles],
  };

  PythonShell.run("predict.py", options, (err, results) => {
    if (err) {
      console.error("Python error:", err);
      return res.status(500).json({ error: err.message });
    }
    // results is an array of lines output by the python script
    return res.status(200).json({ prediction: results.join("\n") });
  });
}
