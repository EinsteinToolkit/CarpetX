#!/bin/bash

if [ -z "${CONTAINERLOCAL:-}" ]; then
  echo "✗ missing required env var: CONTAINERLOCAL" >&2
  exit 1
fi

if [ -z "${CONTAINERLOCALCACTUS:-}" ]; then
  echo "✗ missing required env var: CONTAINERLOCALCACTUS" >&2
  exit 1
fi

if [ -z "${CONTAINERLOCALMACHINE:-}" ]; then
  echo "✗ missing required env var: CONTAINERLOCALMACHINE" >&2
  exit 1
fi

log=$(mktemp)
trap 'rm -f "$log"' EXIT
if docker exec \
  -e CONTAINERLOCALCACTUS="$CONTAINERLOCALCACTUS" \
  -e CONTAINERLOCALMACHINE="$CONTAINERLOCALMACHINE" \
  "$CONTAINERLOCAL" zsh -c '
  cd "$CONTAINERLOCALCACTUS"
  ./simfactory/bin/sim build sim-carpetx -j8 \
    --machine "$CONTAINERLOCALMACHINE" \
    --thornlist=thornlists/carpetx.th
' > "$log" 2>&1; then
  echo "✓ build"
else
  cat "$log"
  echo "✗ build failed" >&2
  exit 1
fi
