# Enhanced CSTR Saponification Simulator

## Overview
This repository contains a pilot-plant scale continuous stirred tank reactor (CSTR) saponification simulator built in MATLAB. The `wave.m` script implements an interactive 3D visualization of a 22,000 L reactor with live particle tracking, continuous feed plumbing, and theme-aware controls. The App Designer class `CSTR_Saponification_Modern` complements the visual model with data ingestion, kinetic studies, and export utilities for documenting NaOH and ethyl acetate neutralization studies.

The reaction tracked by both tools is `CH3COOC2H5 + NaOH -> CH3COONa + C2H5OH`.

## Key Features
- Real-time 3D reactor geometry with dual feed inlets, a discharge outlet, multiple Rushton turbines, baffles, and motor assembly rendered in `uiaxes`.
- Continuous feed management with enable toggles, feed-rate slider, residence-time controller, and particle coloration that signals conversion completion.
- Operations console with main power interlock, RPM knob and gauge (0-600 rpm), start/pause/stop/reset controls, emergency stop, and reactor purge.
- Fluid and reaction tuning via particle density, bulk viscosity, residence time, and target conversion sliders, letting you study hydrodynamics and conversion targets in a single run.
- Full interface theming: every panel, button, lamp, and axis responds to the dark/bright mode toggle for presentation-ready visuals.
- Data-centric App Designer workflow (`CSTR_Saponification_Modern`) featuring data upload, parameter entry, manual kinetic fitting, multi-level optimization, display toggles, and export buttons for plots and Excel summaries.
- Curated result packages (`Result 30c`, `Result 40c`, `result 50c`, `result 60c`) containing figures, spreadsheets, and manuscript-ready material grouped by temperature case studies.

## Repository Structure
- `wave.m`: Entry point for the enhanced 3D CSTR simulator (`Enhanced_CSTR_Fluid_Simulation` function).
- `CSTR_Saponification_Modern.m`: App Designer class that delivers a dashboard for simulation setup, plotting, and reporting.
- `CSTR_Saponification_Modern_needfix.m`: Work-in-progress variant that preserves debug scaffolding and legacy behaviour.
- `Result 30c` / `Result 40c` / `result 50c` / `result 60c`: Sample outputs (PNG, EMF, XLSX) produced by the App workflow for different jacket temperatures.
- `.vscode/`: Editor configuration used during development (optional).

## Requirements
- MATLAB R2022a or newer to access `uifigure`, `uiaxes`, modern UI components, and 3D rendering utilities.
- MATLAB App Designer toolbox installed to run `CSTR_Saponification_Modern`.
- Windows 10 or newer; GPU acceleration is optional but improves the 3D animation smoothness.

## Getting Started
1. Open MATLAB and add this folder to the path, for example `addpath('C:\Users\ASUS\Documents (3211)\Third portfolio')`.
2. Launch the 3D simulator by running `Enhanced_CSTR_Fluid_Simulation` in the Command Window. The reactor window and control panel will appear automatically.
3. Launch the App Designer workflow with:
   ```
   app = CSTR_Saponification_Modern;
   run(app);
   ```
   You can also open the file in App Designer and press Run.
4. Explore the `Result` folders for ready-made figures, datasets, and manuscript material that illustrate expected outputs.

## Operating the Enhanced CSTR Simulation
- Power: Toggle `MAIN POWER` to unlock the controls; the status lamp turns green when power is live.
- Mixing: Adjust the RPM knob (0-600 rpm). The gauge and numeric readout track current turbine speed.
- Feed: Enable continuous feed, fine-tune the feed-rate slider, or trigger `START FEED` / `CLEAR REACTOR` from the operations panel to manage particles in the tank.
- Residence and conversion: Use the residence time slider (5-120 s) and target conversion slider (0-100%) to control how long particles stay in the reactor and when they turn green.
- Visualization: Switch flow lines and particle visibility, and flip the dark/bright mode toggle for presentation or classroom settings.
- Safety controls: `PAUSE`, `STOP`, and `E-STOP` immediately halt animation; `RESET` restores the initial fluid state.

## Working with the App Designer Dashboard
- Import experimental or simulated CSV/XLSX datasets with the Upload button to populate concentration profiles.
- Set temperature, volume, flowrate, inlet concentration, and simulation horizon using the numeric inputs on the left pane.
- Use `Manual K`, `Optimize`, and `Optimize Ultra` to calibrate kinetics against uploaded data; progress indicators update in the status area.
- Toggle jacket and product displays to tailor the visualisation tabs, then export plots or full data tables with the provided buttons.
- Saved outputs are written into the corresponding `Result` folder, ready for documentation or thesis chapters.

## Troubleshooting and Tips
- Always run the scripts from MATLAB; GNU Octave and older MATLAB releases do not support the UI components used here.
- If you edit the App Designer class, call `clear classes` in MATLAB before re-running to avoid cached definitions.
- Reduce the impeller RPM or feed rate if frame rates drop; the particle system is intentionally dense for realism.
- When presenting, activate dark mode and hide flow lines for cleaner projector-friendly visuals.

---

Feel free to adapt the UI text and export templates to match your institutional branding or reporting standards.
