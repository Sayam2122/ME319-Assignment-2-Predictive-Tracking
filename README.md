# ME319 Assignment 2 – Predictive Tracking and Elevation-Axis Control

**Team:** Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194)
**Department:** Mechanical Engineering, **Institute:** IIT Bombay
**Course:** ME319 Control Systems – Autumn 2026

> Student project for a course assignment. It has **not** been approved or evaluated by IIT Bombay.

## Overview
A simulation of a motorised **elevation tracking platform** that points at where a manoeuvring aerial object
*will be* 0.5 s ahead, using noisy 20 Hz range/elevation radar data that can drop out for 0.5 s.

Pipeline: radar (r, θ) → **EKF** (Singer acceleration model) → **open-loop prediction** to t+Tp →
elevation reference → **PI compensator** → geared DC-motor plant derived from first principles.

**Scope.** Predictive tracking **plus laser engagement**: `turret_controller` issues `fire_cmd` (0→1 edge) through a gate
(1 s recharge, estimated miss ≤ k_gate·R_eng, EKF uncertainty ≤ k_sig·R_eng). `simulation_runner` applies the assignment's
pulse rules and scores each pulse by beam-to-target distance (R_eng = 2 m) against the target position at t_fire + Tp
(assumption; the alternative reading, target at t_fire, is reported as well and is one switch, `P.eng.delay`).
Interception success, time and miss distance are reported for A–C and for random unseen cases.

## Repository Structure

```text
ME319-Assignment-2-Predictive-Tracking/
│
├── matlab/
│   ├── main.m
│   ├── parameters.m
│   ├── system_model.m
│   ├── target_motion_generator.m
│   ├── state_estimator.m
│   ├── target_predictor.m
│   ├── ca_singer_model.m
│   ├── controller_design.m
│   ├── turret_controller.m
│   ├── turret_plant.m
│   ├── simulation_runner.m
│   ├── performance_metrics.m
│   ├── step_metrics.m
│   ├── automated_checks.m
│   ├── design_study.m
│   ├── robustness_study.m
│   ├── random_scenario.m
│   ├── fire_tuning.m
│   ├── make_figures.m
│   └── write_results_csv.m
│
├── results/
│   ├── figures/
│   └── data/
│
├── paper/
│   ├── main.tex
│   ├── references.bib
│   ├── ai_appendix.tex
│   └── main.pdf
│
├── docs/
│   └── AI_usage.md
│
├── README.md
└── LICENSE
```text

### Folder Description

- **`matlab/`** – MATLAB source code for system modelling, target motion generation, state estimation, prediction, controller design, simulation and performance evaluation.
- **`results/`** – Generated simulation figures and result data in PNG, MAT and CSV formats.
- **`paper/`** – LaTeX source files, bibliography, AI usage appendix and compiled assignment report.
- **`docs/`** – Additional project documentation and AI usage details.
```
The three assignment interfaces are kept: `target_motion_generator`, `turret_controller`, `turret_plant`.
Extra optional outputs (ground truth, diagnostics) do not change the 3-output/2-output calls.

## Installation
* **Octave ≥ 7** + `control` package (`sudo apt install octave octave-control gnuplot`) – developed with Octave 8.4.0, or
* **MATLAB R2019b+** with Control System Toolbox.

## Run
```
cd matlab
octave --no-gui --eval "main"        % full run, ~8 min  (MATLAB: run main)
octave --no-gui --eval "main('quick')"  % skips the design study, 4 random cases
```
`main` runs the 10 self-checks, the controller analysis, the estimator design study, the firing-gate tuning,
scenarios A–C and the robustness study, and writes figures (`results/figures`) and CSV/MAT data (`results/data`).
It does **not** write LaTeX or build the paper. To compile the paper by hand:
`cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main` (needs TeX Live with IEEEtran).
Put your own repository URL in `\repourl` in `paper/main.tex`, and your own names in the `\author` block.
For the instructor's evaluation cases: build a scenario struct like `random_scenario` returns and call
`simulation_runner(cfg, parameters())`, then `performance_metrics(res, P).eng`.

## Mathematical model summary
* Plant: `Gp(s) = N·Kt / ( s[ Lm·Jeq s² + (Rm·Jeq + Lm·b) s + Rm·b + N²·Kt·Ke ] )`, `Jeq = Jturret + N²Jm = 3200 kg m²`
  → `1000/(320 s³ + 3350 s² + 1.0015e6 s)`; resonance ωn = 55.9 rad/s, ζ = 0.094.
* Estimator: per-axis Singer model `[p v a]`, EKF with `h = [√(x²+y²), atan2(y,x)]`, Joseph update; prediction `F(t−tk+Tp)·X̂k`.
* Controller: `Gc(s) = Kp (s+ωi)/s`, `Kp = ωc/Kv`, ωc = 5 rad/s, ωi = ωc/3.

## Results
The numbers quoted in the paper were copied by hand from `results/data/metrics_summary.csv` and `robustness_cases.csv`;
see the paper for the table and `results/figures/` for plots.

## Parameters
Every value in `matlab/parameters.m` is tagged `[ASSIGNMENT]`, `[DERIVED]`, `[DESIGN]` or `[ASSUMPTION]`.
Assumptions not in the assignment: balanced platform (no gravity), ideal encoder, rigid gears,
supply limit ±1500 V (never reached), platform initially aligned with the target.

## Reproducibility
Noise uses a Park–Miller + Box–Muller generator with fixed seed (`P.seed`), identical on MATLAB and Octave.
Tuning (`random_scenario` idx 101–116), robustness test (idx 201–210) and scenarios A–C are disjoint.

## Limitations
Simulation only; Octave-verified (not run in MATLAB here); 2-D idealised plant; smooth sinusoidal-acceleration manoeuvres; no engagement logic.

License: MIT.
