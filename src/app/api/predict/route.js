export const runtime = "nodejs"; // ✅ Ensure it's a function, not a static file

export async function POST(req) {
  try {
    console.log("API received a request...");

    const { smiles } = await req.json();
    if (!smiles) {
      console.log("Missing SMILES input.");
      return new Response(JSON.stringify({ error: "No SMILES provided" }), { status: 400 });
    }

    console.log("Running Python script with input:", smiles);

    const { spawn } = require("child_process");
    const pythonProcess = spawn("python3", ["predict.py", smiles]);

    return new Promise((resolve) => {
      let results = "";
      pythonProcess.stdout.on("data", (data) => {
        console.log("Python Output:", data.toString());
        results += data.toString();
      });
      pythonProcess.stderr.on("data", (err) => console.error("Python Error:", err.toString()));
      pythonProcess.on("close", () => {
        console.log("Sending Response:", results);
        resolve(new Response(results.trim(), { status: 200, headers: { "Content-Type": "application/json" } }));
      });
    });
  } catch (error) {
    console.error("API Server Error:", error);
    return new Response(JSON.stringify({ error: "Internal Server Error" }), { status: 500 });
  }
}
