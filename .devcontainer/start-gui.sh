#!/usr/bin/env bash
# Opens the desktop app on the codespace's virtual screen (port 6080 in the browser).
# Run it again to reopen the app after closing it.
export DISPLAY="${DISPLAY:-:1}"
cd "$(dirname "$0")/.."
for _ in $(seq 120); do
  [ -S /tmp/.X11-unix/X1 ] && break
  sleep 1
done
if pgrep -f "OTA_Hub_AntennaToolkit.py gui" >/dev/null; then
  echo "The app is already running - open port 6080."
else
  nohup python OTA_Hub_AntennaToolkit.py gui >/tmp/otahub-gui.log 2>&1 &
  echo "The app is starting on the desktop - open port 6080 (password: vscode)."
fi
