# MCP_SEY — secondary-electron avalanche in a single microchannel

Geant4 simulation of a **single channel** of a microchannel plate (MCP).
One electron is injected at the channel entrance, accelerated by the axial
field, and the secondary-electron avalanche is tracked trajectory by
trajectory until it leaves the plate.

## 1. Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Requires Geant4 ≥ 11.0 and CMake ≥ 3.16.

## 2. Run

```bash
./build/MCP_SEY macros/scan_wu.mac   # batch mode, multithreaded
./build/MCP_SEY                      # interactive mode, visualisation
```

With no argument the program runs `macros/vis.mac` and opens a 3D view,
in single-threaded mode.

> The output directory is created automatically if it does not exist. If it
> exists but is not writable, the run stops at once with a clear message

## 3. Geometry

Four volumes, built by `DetectorConstruction`:

| Volume    | Material   | Role                                    |
|-----------|------------|-----------------------------------------|
| `World`   | vacuum     | enclosing box                           |
| `MCP`     | lead glass | the plate, with one tilted pore cut out |
| `Channel` | vacuum     | the pore — **the only volume with field** |
| `Drift`   | vacuum     | 50 µm past the plate, field-free        |

Uniform electric field of magnitude `V/L` along `-z`, so the electron is
accelerated towards `+z`. The downstream face of `Drift` is the scoring plane.

## 4. Physics

`FTFP_BERT` extended by `FurmanPiviPhysics`, which attaches a discrete
secondary-emission process to the electron. Two interchangeable models:

| Model                     | Command                   | Status  |
|---------------------------|---------------------------|---------|
| Wu et al. 2008            | `/mcp/setWuModel true`    | default |
| Furman–Pivi / Peng et al. | `/mcp/setWuModel false`   | option  |

The active model is printed at the start of every run and stored in the
`EmissionModel` column of the `configuration` table.

> **The Furman–Pivi model falls 4 to 41 times below the curves published by
> Peng, even though it uses Peng's own parameters.**

## 5. Commands

### `/mcp/` — geometry and run control

| Command                       | Default  | When               |
|-------------------------------|----------|--------------------|
| `setPlateDiameter d mm`       | 0.5 mm   | before `initialize`|
| `setThickness l mm`           | 1.0 mm   | before / between runs |
| `setChannelDiameter d um`     | 15 µm    | before `initialize`|
| `setChannelAngle a deg`       | 8°       | before `initialize`|
| `setVoltage v volt`           | 2000 V   | before / between runs |
| `setMaxTracksPerEvent n`      | 100000   | before `initialize`|
| `setMaxGenerations n`         | 0        | before / between runs |
| `setWuModel b`                | true     | before / between runs |
| `setOutputFile path`          | mcp.root | before / between runs |
| `setTransitOutput b`          | false    | before `initialize`|

`0` means "no limit" for both caps.

### `/source/gauss/` — the injected electron

`muE`/`sigmaE` (keV), `muZ`/`sigmaZ` (mm), `muPz`/`sigmaPz` (keV), `zOrigin`.
A zero standard deviation gives a fixed value. A `muPz` larger than the total
momentum gives a purely axial `+z` shot.

> **Geometry** commands are only accepted **before** `/run/initialize`. Issued
> later, Geant4 answers *"command refused: illegal application state"*.

## 6. Example macro

```
/run/numberOfThreads 10
/random/setSeeds 12345 678

/mcp/setThickness 0.46 mm
/mcp/setChannelDiameter 10 um
/mcp/setChannelAngle 8 deg
/mcp/setVoltage 950 V
/mcp/setMaxTracksPerEvent 1000000

/run/initialize

/source/gauss/muE 300 eV
/source/gauss/sigmaE 0 eV
/source/gauss/muPz 1 MeV

/mcp/setOutputFile root_files/test.root
/run/beamOn 500
```

For a parameter scan, see the pair `macros/scan_wu.mac` and
`macros/scan_wu_point.mac`, driven by `/control/foreach`.

## 7. ROOT output

| Table           | Contents                                          |
|-----------------|---------------------------------------------------|
| `events`        | one row per injected electron (18 columns)        |
| `downstream`    | one row per electron at the scoring plane         |
| `configuration` | one row per run: geometry, caps, model, source    |

Plus two histograms, `impactEnergy` and `impactAngle`.

The column you want is `Gain`. `MeanTime_ns` is the time centroid of the
output bunch; its spread **between events** is the TTS.

> `TimeSpread_ns` is the internal width of the bunch
> within a single event.

> If `IsMultiplicationLimited` is non-zero, the track cap was hit and **the
> mean gain is underestimated**. Every analysis macro prints this fraction —
> read it first.

## 8. Analysis

```bash
root -l -b -q 'analysis/gain_voltage.C("root_files/scanWu_*.root","figures/gain")'
```

| Macro                | Purpose                                 |
|----------------------|-----------------------------------------|
| `gain_voltage.C`     | gain, transit time, TTS as FWHM         |
| `scan_thickness.C`   | plate-thickness scan                    |
| `scan_voltage.C`     | voltage scan, detailed version          |
| `scan_generations.C` | effect of the generation cap            |
| `scan_tranches.C`    | yield per depth slice                   |
| `transverse.C`       | spatial distribution at the exit        |
| `analyse.C`          | quick look at one file                  |

## 9. What the simulation does not do

- **No space charge.** The gain diverges where a real plate saturates, around
  10⁶ electrons per channel. Comparisons with measurements are only valid in
  the unsaturated regime.
- **One channel only.** No neighbours, no resistive coating, no recharge
  current.
- **No readout model.** The calculation stops at the scoring plane.

## 10. Layout

```
MCP_SEY.cc            entry point
include/ src/         one class per file, following the B1-B5 examples
macros/               run and scan macros
analysis/             ROOT analysis macros
```

The physics core is `FurmanPiviProcess` (detects the wall impact, creates the
secondaries) together with the two models `FurmanPiviModel` / `WuModel`.

## References

- M. A. Furman and M. T. F. Pivi, *Phys. Rev. ST Accel. Beams* **5** (2002) 124404
- L. Wu and C. A. Kruschwitz, *Rev. Sci. Instrum.* **79** (2008) 073104
- H. Peng et al., *Nucl. Instrum. Methods A* **1062** (2024) 169163
- F. Li et al., *Photonics* **9** (2022) 978
