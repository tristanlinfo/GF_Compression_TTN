#!/bin/bash
set -euo pipefail

DOT=results/ttn_graph.dot
PNG=results/ttn_graph.png
PDF=results/ttn_graph.pdf

mkdir -p results

if command -v dot >/dev/null 2>&1; then
  echo "Rendering PNG with dot (neato/dot)..."
  dot -Tpng "$DOT" -o "$PNG"
  echo "Rendered $PNG"
  dot -Tpdf "$DOT" -o "$PDF"
  echo "Rendered $PDF"
elif command -v sfdp >/dev/null 2>&1; then
  echo "Rendering PNG with sfdp..."
  sfdp -Tpng "$DOT" -o "$PNG"
  echo "Rendered $PNG"
  sfdp -Tpdf "$DOT" -o "$PDF"
  echo "Rendered $PDF"
else
  echo "Graphviz not found (dot/sfdp). Please install graphviz to render the DOT file." >&2
  exit 1
fi

#chmod +x src/render_ttn_graph.sh   
#./src/render_ttn_graph.sh 