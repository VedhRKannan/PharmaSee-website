export const runtime = "nodejs"; // ✅ Ensure Vercel treats this as a Node.js function

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
    return new Response(JSON.stringify({ error: "Server error" }), { status: 500 });
  }
}
