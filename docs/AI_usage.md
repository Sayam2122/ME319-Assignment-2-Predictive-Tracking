# AI and tool usage (ME319 Assignment 2)

**Team:** Sayam Midha (24B2240), Piyush Meena (24B2275), Monaal Kumar Arya (24B2156), Asmit Anand (24B2263), Shivam Gawande (24B2194). Full one-page appendix is in `paper/ai_appendix.tex`.

* **AI:** Anthropic Claude (Sonnet 5.5) via claude.ai with a sandboxed code-execution tool.
* **Software:** GNU Octave 8.4.0 + control 4.0.1, gnuplot, TeX Live (IEEEtran), Python (file patching only).
* **Team's work:** hand derivation of the model and design (plant, Singer/EKF, prediction, PI design), instructions to the AI, running and testing the simulations, review against hand calculations.
* **AI contributions:** Octave code, self-checks, simulations, figures, LaTeX set-up and report wording, README.
* **Accepted:** PI compensator (ωc = 5 rad/s), EKF + Singer model, open-loop prediction through loss.
* **Modified:** (q, α) re-selected by ensemble study; scenario-B periods lengthened; plotting/test bugs fixed.
* **Rejected:** notch filter, derivative action, `stepinfo`.
* **Added later (second Claude session):** firing gate + engagement scoring (`turret_controller.m`, `simulation_runner.m`, `fire_tuning.m`), PDF build step in `main.m`.
* **Validation:** 10 automated checks (`automated_checks.m`), NIS consistency, analytic Routh gain limit vs numerical margin,
  energy balance of the simulated plant, disjoint tuning/test sets.
