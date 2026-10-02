# Trying the toolkit in a codespace

The desktop app runs on a virtual screen inside this codespace.

1. Open the **Ports** tab (next to Terminal) and click the globe icon on port
   **6080** (*Antenna app (desktop)*). It may already have opened in a new tab.
2. Click **Connect**; the password is `vscode`.
3. The app is on that desktop. If you closed it, reopen it from the terminal:

   ```bash
   bash .devcontainer/start-gui.sh
   ```

The command line works in the terminal here too:

```bash
python OTA_Hub_AntennaToolkit.py list
python OTA_Hub_AntennaToolkit.py synth rectangular_patch_inset --f0 2.4GHz --set eps_r=4.4 --set h=1.6mm
python OTA_Hub_AntennaToolkit.py export pyramidal_horn --f0 10GHz --set G_target=20 --format stl -o horn.stl
python OTA_Hub_AntennaToolkit.py check
```

Files you save (STL, CST/HFSS scripts) land in the codespace; right-click one in
the Explorer and choose **Download** to get it. Stop the codespace when you are
done (github.com/codespaces) so it does not use your monthly free hours.
