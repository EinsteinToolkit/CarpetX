# How to access test logs and build artifacts

Both scripts run inside the running container named `$CONTAINERLOCAL`; build artifacts and testsuite logs live **inside the container**, not on the host. `$CONTAINERLOCAL` and `$CONTAINERLOCALCACTUS` are pre-set in the host shell:

```bash
docker exec -e CONTAINERLOCALCACTUS="$CONTAINERLOCALCACTUS" \
  "$CONTAINERLOCAL" zsh -c 'tail -80 "$CONTAINERLOCALCACTUS/<path>"'
```

Per-test logs: `$CONTAINERLOCALCACTUS/TEST/sim-carpetx/<Thorn>/<test>.log`.
Build artifacts: `$CONTAINERLOCALCACTUS/configs/sim-carpetx/`.

