"""Execute tutorial notebooks without running their Colab install cell.

Usage, from the repository root:
    uv run --group notebook python docs/tools/execute_notebooks.py docs/tutorials/*.ipynb
The install cell (any code cell whose source starts with "%pip install") is
blanked during execution, then restored with no outputs.
"""
import sys, time, nbformat
from nbclient import NotebookClient

for path in sys.argv[1:]:
    nb = nbformat.read(path, as_version=4)
    saved = {}
    for i, cell in enumerate(nb.cells):
        if cell.cell_type == "code" and cell.source.lstrip().startswith("%pip install"):
            saved[i] = cell.source
            cell.source = "pass"
    t = time.time()
    NotebookClient(nb, timeout=300, kernel_name="python3",
                   resources={"metadata": {"path": "docs/tutorials"}}).execute()
    n = 0
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        if i in saved:
            cell.source, cell.outputs, cell.execution_count = saved[i], [], None
            continue
        n += 1
        cell.execution_count = n
        for out in cell.outputs:
            if "execution_count" in out:
                out["execution_count"] = n
        cell.metadata = {}
        bad = [o for o in cell.outputs if o.output_type == "error" or o.get("name") == "stderr"]
        if bad:
            raise SystemExit(f"{path}: cell {i} has error/stderr output")
    nb.metadata = {k: v for k, v in nb.metadata.items() if k in ("kernelspec", "language_info")}
    nbformat.validate(nb)
    nbformat.write(nb, path)
    print(f"{path}: ok in {time.time() - t:.1f}s")
