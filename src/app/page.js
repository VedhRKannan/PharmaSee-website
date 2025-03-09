'use client';

import { useState } from "react";

export default function Home() {
  const [smiles, setSmiles] = useState("");
  const [prediction, setPrediction] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const examples = ["CCO", "c1ccccc1", "CC(=O)O"];

  const handleExampleClick = (example) => {
    setSmiles(example);
  };

  const handleSubmit = async (event) => {
    event.preventDefault();
    setLoading(true);
    setError(null);
    setPrediction(null);

    try {
      const response = await fetch("/api/predict", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles }), // Ensure JSON format
      });

      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`HTTP error! Status: ${response.status}, Message: ${errorText}`);
      }

      const data = await response.json();
      console.log("API Response:", data); // ✅ Log full response for debugging
      setPrediction(data); // ✅ Store full JSON object instead of just `data.prediction`
    } catch (error) {
      console.error("Fetch error:", error);
      setError(`Failed to fetch prediction: ${error.message}`);
    }

    setLoading(false);
  };

  return (
    <div style={styles.container}>
      <h1>Pharmasee</h1>
      <p>Enter a SMILES string to predict ADMET properties:</p>
      <form onSubmit={handleSubmit} style={styles.form}>
        <input
          type="text"
          placeholder="Enter SMILES"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          style={styles.input}
        />
        <button type="submit" style={styles.button} disabled={loading}>
          {loading ? "Predicting..." : "Predict"}
        </button>
      </form>

      {loading && <p>Loading...</p>}
      {error && <p style={{ color: "red" }}>{error}</p>}

      {prediction && (
        <div style={styles.result}>
          <h2>Prediction</h2>
          <pre>{JSON.stringify(prediction, null, 2)}</pre> {/* ✅ Properly displays JSON */}
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
    </div>
  );
}

const styles = {
  container: {
    maxWidth: "600px",
    margin: "0 auto",
    padding: "1rem",
    fontFamily: "Arial, sans-serif",
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
  },
  button: {
    padding: "0.5rem 1rem",
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
