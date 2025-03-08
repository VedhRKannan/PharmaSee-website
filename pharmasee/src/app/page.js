"use client";

import { useState } from "react";

export default function Home() {
  const [smiles, setSmiles] = useState("");
  const [prediction, setPrediction] = useState(null);
  const [loading, setLoading] = useState(false);

  const examples = [
    "CCO",       // Ethanol
    "c1ccccc1",  // Benzene
    "CC(=O)O",   // Acetic acid
  ];

  const handleExampleClick = (example) => {
    setSmiles(example);
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    setPrediction(null);

    try {
      const response = await fetch("/api/predict", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles }),
      });
      const data = await response.json();
      setPrediction(data.prediction);
    } catch (error) {
      console.error("Error:", error);
      setPrediction("An error occurred.");
    }

    setLoading(false);
  };

  return (
    <main style={styles.main}>
      <h1 style={styles.title}>Pharmasee</h1>
      <p>Enter a SMILES string to predict ADMET properties:</p>
      <form onSubmit={handleSubmit} style={styles.form}>
        <input
          type="text"
          placeholder="Enter SMILES"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          style={styles.input}
        />
        <button type="submit" style={styles.button}>
          Predict
        </button>
      </form>
      {loading && <p>Loading...</p>}
      {prediction && (
        <div style={styles.result}>
          <h2>Prediction</h2>
          <pre>{prediction}</pre>
        </div>
      )}
      <h2>Examples</h2>
      <div style={styles.examples}>
        {examples.map((ex, idx) => (
          <button key={idx} onClick={() => handleExampleClick(ex)} style={styles.exampleButton}>
            {ex}
          </button>
        ))}
      </div>
    </main>
  );
}

const styles = {
  main: {
    maxWidth: "600px",
    margin: "0 auto",
    padding: "1rem",
    fontFamily: "Arial, sans-serif",
  },
  title: {
    textAlign: "center",
  },
  form: {
    display: "flex",
    flexWrap: "wrap",
    marginBottom: "1rem",
  },
  input: {
    flex: "1 1 auto",
    padding: "0.5rem",
    marginRight: "0.5rem",
    marginBottom: "0.5rem",
    minWidth: "0",
  },
  button: {
    padding: "0.5rem 1rem",
    marginBottom: "0.5rem",
  },
  result: {
    backgroundColor: "#f9f9f9",
    padding: "1rem",
    borderRadius: "5px",
  },
  examples: {
    display: "flex",
    flexWrap: "wrap",
    gap: "0.5rem",
  },
  exampleButton: {
    padding: "0.5rem 1rem",
    cursor: "pointer",
  },
};
