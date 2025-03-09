import { PythonShell } from "python-shell";

export const runtime = "nodejs"; // Required for child_process to work

export async function POST(req) {
  try {
    const { smiles } = await req.json();
    if (!smiles) {
      return new Response(JSON.stringify({ error: "No SMILES provided" }), { status: 400 });
    }

    let options = {
      mode: "text",
      pythonPath: "/Users/ramalingamkannan/coding/PharmaSee-website/pharmasee/venv/bin/python3", // Adjust if needed
      args: [smiles],
    };

    return new Promise((resolve) => {
      PythonShell.run("predict.py", options, (err, results) => {
        if (err) {
          console.error("Python error:", err);
          resolve(new Response(JSON.stringify({ error: "Internal server error" }), { status: 500 }));
        } else {
          try {
            const parsedOutput = JSON.parse(results.join("\n"));
            resolve(new Response(JSON.stringify(parsedOutput), { status: 200 }));
          } catch (jsonError) {
            console.error("JSON Parse Error:", jsonError);
            resolve(new Response(JSON.stringify({ error: "Invalid JSON response from Python script" }), { status: 500 }));
          }
        }
      });
    });
  } catch (error) {
    console.error("API Error:", error);
    return new Response(JSON.stringify({ error: "Server error" }), { status: 500 });
  }
}
