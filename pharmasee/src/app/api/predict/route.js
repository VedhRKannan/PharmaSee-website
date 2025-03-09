import { PythonShell } from "python-shell";

export const runtime = "nodejs"; // Required for child_process to work

export async function POST(req) {
  try {
    const { smiles } = await req.json();
    if (!smiles) {
      console.error("ERROR: No SMILES string provided.");
      return new Response(JSON.stringify({ error: "No SMILES provided" }), { status: 400 });
    }

    console.log("DEBUG: Running Python script with input:", smiles);

    let options = {
      mode: "text",
      pythonPath: "/Users/ramalingamkannan/coding/PharmaSee-website/pharmasee/venv/bin/python3", // ✅ Ensuring python3 is used
      args: [smiles],
    };

    return new Promise((resolve) => {
      const pyshell = new PythonShell("predict.py", options); // ✅ Correct PythonShell usage

      let results = [];
      let hasTimedOut = false;

      // Collect Python output
      pyshell.on("message", function (message) {
        if (!hasTimedOut) results.push(message);
      });

      // Handle errors
      pyshell.on("error", function (err) {
        console.error("Python error:", err);
        resolve(new Response(JSON.stringify({ error: "Internal server error" }), { status: 500 }));
      });

      // When Python script exits
      pyshell.on("close", function (code) {
        if (hasTimedOut) return; // Skip processing if timeout already triggered

        console.log("DEBUG: Python script exited with code", code);
        console.log("DEBUG: Raw Python output:", results);

        try {
          const parsedOutput = JSON.parse(results.join("\n")); // ✅ Ensure correct JSON parsing
          resolve(new Response(JSON.stringify(parsedOutput), { status: 200 }));
        } catch (jsonError) {
          console.error("JSON Parse Error:", jsonError);
          resolve(new Response(JSON.stringify({ error: "Invalid JSON response from Python script" }), { status: 500 }));
        }
      });

      // **Timeout if Python takes too long**
      setTimeout(() => {
        console.error("ERROR: Python script timed out.");
        hasTimedOut = true;
        pyshell.terminate();
        resolve(new Response(JSON.stringify({ error: "Python script timeout" }), { status: 500 }));
      }, 10000); // 10 seconds timeout
    });
  } catch (error) {
    console.error("API Error:", error);
    return new Response(JSON.stringify({ error: "Server error" }), { status: 500 });
  }
}
