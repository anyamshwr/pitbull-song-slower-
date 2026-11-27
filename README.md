# Harmonic–Percussive Time-Scale Modification in MATLAB

Course project for **ECE 2003 Signal Representation and Processing** (Monsoon 2026, Shiv Nadar University).
We re-implement the method of

> J. Driedger, M. Müller and S. Ewert, “Improving Time-Scale Modification of Music Signals Using Harmonic-Percussive Separation,” *IEEE Signal Processing Letters*, vol. 21, no. 1, pp. 105–109, 2014.

## The idea

Time-scale modification (TSM) changes the tempo of audio without changing its pitch. Classic methods fail on one kind of sound each:

- the **phase vocoder** (long frames) keeps notes clean but smears drum hits;
- **overlap-add (OLA)** with tiny frames keeps drum hits sharp but makes notes warble.

The paper splits the music first and stretches each part with the method that suits it:

```
                ┌─► x_h (notes) ──► phase vocoder, N = 4096, phase locking ──► y_h ─┐
x ──► HPS ──────┤                                                                   ├─► y = y_h + y_p
   (N = 1024)   └─► x_p (drums) ──► OLA, N = 256 ───────────────────────────► y_p ─┘
```

## Repository layout

| Folder | Contents |
|---|---|
| `src/` | `hps.m` (separation), `tsm_pv.m` (phase vocoder), `tsm_ola.m` (OLA), `tsm_wsola.m` (WSOLA baseline), `tsm_hp.m` (full method), `load_audio.m`, `make_test_signals.m` |
| `tests/` | Acceptance tests for each block and for the integrated pipeline |
| `experiments/` | `run_pipeline.m` (quick demo), `exp_run_all.m` (all methods, runtimes, figures) |
| `audio/` | Synthetic test signals: a tone, clicks, and their mix |
| `report/` | Final report in IEEE format (LaTeX source, PDF, figures) |

## How to run

Plain MATLAB; no toolbox is needed for the included 22,050 Hz test signals. `load_audio` uses the Signal Processing Toolbox (`resample`) only when a file has a different sampling rate.

```matlab
setup                 % run once from the repository root
test_integration      % checks every block end to end
run_pipeline          % stretches synth_mix by 1.0, 1.2 and 1.8 into out/
```

Stretch your own audio:

```matlab
x = load_audio('audio/your_clip.wav');   % mono, 22,050 Hz
y = tsm_hp(x, 1.8);                      % 1.8x longer, same pitch
soundsc(y, 22050)
```

Every TSM function returns a column vector of exactly `round(alpha*length(x))` samples.

## Results on the test signals

| Check | Result |
|---|---|
| Separation: max \|x_h + x_p − x\| | 3.3 × 10⁻¹⁶ |
| Full pipeline at α = 1 (should equal input) | 4.1 × 10⁻¹² |
| Pitch after stretching by 1.8 | 440.0 Hz (unchanged) |
| Click positions after stretching by 1.8 | 7 of 7 within 1 ms |

Click sharpness at α = 1.8: the share of each click's high-passed energy that falls in its loudest 10 ms (higher is sharper).

| Original | **HP** | PV | WSOLA | OLA (whole signal) |
|---|---|---|---|---|
| 0.84 | **0.84** | 0.52 | 0.80 | 0.87, but notes destroyed |

HP keeps clicks as sharp as the original, the phase vocoder smears them, and WSOLA doubles some of them, as in the paper's Fig. 3.

## Team

| Member | Contribution |
|---|---|
| Anshika | Harmonic–percussive separation (`hps.m`) |
| Mridhula | Phase vocoder with identity phase locking (`tsm_pv.m`) |
| Anya Maheshwari | OLA, WSOLA baseline, integration, experiments |

Supervised by Prof. Vijay.

## Differences from the paper

Hop sizes and WSOLA parameters are our own choices, because the paper does not give them. We do not compare against the commercial élastique or the Nagel–Walther baselines, and our listening test is much smaller than the paper's.
