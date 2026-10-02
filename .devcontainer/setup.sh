#!/usr/bin/env bash
# Qt's X11 (xcb) platform libraries, then the toolkit with its GUI and 3-D extras.
set -euo pipefail
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  libegl1 libgl1 libfontconfig1 libdbus-1-3 libxkbcommon0 libxkbcommon-x11-0 \
  libxcb-cursor0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 \
  libxcb-render-util0 libxcb-shape0 libxcb-xinerama0 libxcb-xkb1 fonts-dejavu-core
python -m pip install --upgrade pip
python -m pip install -e ".[gui,dev]"
