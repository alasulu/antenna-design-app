#!/usr/bin/env python3
"""OTA Hub Antenna Toolkit — top-level launcher.

Run the command line interface:

    python OTA_Hub_AntennaToolkit.py list
    python OTA_Hub_AntennaToolkit.py show half_wave_dipole
    python OTA_Hub_AntennaToolkit.py synth rectangular_patch_inset --f0 2.4GHz --set eps_r=4.4
    python OTA_Hub_AntennaToolkit.py check

The graphical front end arrives in a later build session; when it does, this
launcher gains a ``gui`` subcommand and the CLI stays exactly as it is.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from otahub.cli.main import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
