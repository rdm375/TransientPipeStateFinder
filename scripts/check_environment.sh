#!/usr/bin/env bash
set -euo pipefail
SOLVER="${1:-../isothermal-pipe-sim}"
echo "Compiler: $(c++ --version | head -n 1)"
echo "CMake: $(cmake --version | head -n 1)"
if [[ ! -f "$SOLVER/CMakeLists.txt" ]]; then
  echo "Missing forward solver at $SOLVER" >&2
  exit 1
fi
git -C "$SOLVER" rev-parse HEAD
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPIPE_SIM_SOURCE_DIR="$SOLVER"
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
python3 -m unittest discover -s tests/python -v
