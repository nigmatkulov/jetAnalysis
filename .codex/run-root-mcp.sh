#!/bin/sh

set -eu

project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
root_mcp="$project_dir/.root-mcp-env/bin/root-mcp"

if [ ! -x "$root_mcp" ]; then
    echo "ROOT-MCP is not installed; follow the setup in README.md" >&2
    exit 1
fi

export MPLCONFIGDIR="$project_dir/.root-mcp-env/matplotlib"
mkdir -p "$MPLCONFIGDIR"

exec "$root_mcp" serve-stdio \
    --data-path "$project_dir" \
    --mode extended \
    --enable-root \
    --allowed-root "$project_dir"
