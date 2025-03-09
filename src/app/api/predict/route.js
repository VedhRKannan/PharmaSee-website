export const runtime = "nodejs"; // ✅ Ensures it's a server function

export async function POST(req) {
  try {
    const { smiles } = await req.json();
    if (!smiles) {
      return new Response(JSON.stringify({ error: "No SMILES provided" }), { status: 400 });
    }

    console.log("DEBUG: Running Python script with input:", smiles);

    const { spawn } = require("child_process");
    const pythonProcess = spawn("python3", ["predict.py", smiles]);

    return new Promise((resolve) => {
      let results = "";
      pythonProcess.stdout.on("data", (data) => (results += data.toString()));
      pythonProcess.stderr.on("data", (err) => console.error("Python Error:", err.toString()));
      pythonProcess.on("close", () => {
        resolve(new Response(results.trim(), { status: 200, headers: { "Content-Type": "application/json" } }));
      });
    });
  } catch (error) {
    console.error("Server Error:", error);
    return new Response(JSON.stringify({ error: "Internal Server Error" }), { status: 500 });
  }
}
