# forc-demo

A MATLAB demonstration of FORC (First Order Reversal Curve) analysis: a
small, working simulation-and-plotting pipeline you can point at any
single-particle magnetization solver, plus a growing set of worked
examples. Originally extracted from a project modeling interacting
Stoner-Wohlfarth particle packings (`meteor-dipoles`); this repo keeps
the FORC machinery and swaps the real physics for illustrative
placeholders you can follow end to end.

## What FORC analysis is

A FORC study repeats one experiment many times. Starting from
saturation, the field is swept down to some **reversal field** `Hr`
(this is the **major descending branch**), then swept back up to
saturation while the magnetization `M` is recorded — that upward sweep
is one **reversal curve**. Doing this from many different `Hr` builds a
family `M(Hr,H)`.

The **FORC distribution** is the mixed second derivative of that
family,

```
rho(Hr,H) = -1/2 * d2M / (dHr dH)
```

A single particle that switches irreversibly at one field contributes
a delta function to `rho`: the distribution reads as a density of
switching events. It's conventional to plot it in rotated coordinates,

```
Hc = (H - Hr) / 2      % coercivity coordinate
Hu = (H + Hr) / 2      % interaction coordinate
```

so that the spread along `Hc` reflects the spread of coercivities in
the sample, and the spread along `Hu` reflects the spread of
interaction (bias) fields between particles. `Hu = 0` means no bias.

## Repository structure

- **`FORCSimulation.m`** — the engine. A handle class that sweeps the
  major branch, picks reversal fields, measures every reversal curve,
  and computes the FORC distribution (Pike's method — see
  `forcDistribution.m`).
- **Solvers** — pluggable single-point magnetization models (see
  below). This is the piece each new demo replaces.
- **Demo scripts** — orchestrate one `FORCSimulation` run for one
  solver and write its results to a JSON file.
- **Plotting scripts** — read a demo's JSON and draw figures from it.
- **`forcDistribution.m`** — Pike's method fit, shared by the engine
  and every plotting script (one copy, not duplicated per file).
- **`createFigure.m`, `divergingColormap.m`, `savePlot.m`** — shared
  figure setup, a zero-centered blue-white-red colormap (rho changes
  sign, and the sign is the point), and a PDF exporter. Figures land
  in `plots/`, created on demand.

## The `SolveFcn` extension point

`FORCSimulation` never touches physics directly — every magnetization
value comes from one function handle, called with the signature

```matlab
[M,state] = SolveFcn(H,Hmax,Hmin,state)
```

`state` is opaque to the engine: it's just a solver's own memory
(a magnetization direction, an internal variable, whatever the model
needs), snapshotted after every step of the descending branch so a
reversal curve can restart mid-branch instead of resweeping from
saturation. Call your solver with `state=[]` to mean "start fresh at
saturation."

| Solver | Model | Used by |
|---|---|---|
| `simulateHysteresisCurve.m` | A single hysteron smoothed by a `tanh` saturation curve, lagged by a coercivity `Hc` in the direction the field last moved | `runFORCDemo.m` |
| `simulateIdealHysteron.m` | A single **ideal rectangular** hysteron: switches instantaneously at `+Hsw` (up) and `-Hsw` (down), no smoothing | `runFORCDemoIdealHysteron.m` |
| `simulateIdealHysteronPopulation.m` | Several independent ideal hysterons (wraps `simulateIdealHysteron.m` once per particle), averaged — no interaction between them yet | `runFORCDemoTwoHysterons.m` |
| `simulateLangevin.m` | Reversible equilibrium Langevin particle, evaluated from the applied field | `runFORCDemoLangevin.m` |

To write a new solver, follow either file as a template: keep the
signature, decide what `state` needs to remember, and let
`FORCSimulation` handle everything else.

## Running a demo

### Demo 1 — placeholder smoothed hysteron

```matlab
runFORCDemo            % writes forc_demo_output.json
plotFORCFamily
plotFORCDistribution
```

### Demo 2 — single ideal hysteron (one particle, one point)

```matlab
runFORCDemoIdealHysteron   % writes forc_demo_ideal_hysteron_output.json
plotFORCFamily("forc_demo_ideal_hysteron_output.json", ...
    "OutputPrefix","idealHysteron_","SolverLabel","ideal rectangular hysteron")
plotFORCDistribution("forc_demo_ideal_hysteron_output.json", ...
    "OutputPrefix","idealHysteron_","MarkerHcHu",[0.3 0])
```

This is the textbook validation case: an ideal hysteron's FORC
distribution is a true delta function at `(Hc,Hu) = (Hsw,0)` in the
continuum. The `MarkerHcHu` marker overlays that theoretical point on
the numerical result — but don't expect the numerical peak to land
*exactly* on it. See "A numerical surprise" below: even at perfect
grid alignment, Pike's fit spreads the delta over a small, predictable
cluster of cells rather than a single one, offset from the
theoretical point by half a grid step. That's a genuine, worthwhile
part of the demo, not an error to chase out.

### Demo 3 — two independent ideal hysterons (two particles, two points)

```matlab
runFORCDemoTwoHysterons   % writes forc_demo_two_hysterons_output.json
plotFORCFamily("forc_demo_two_hysterons_output.json", ...
    "OutputPrefix","twoHysterons_","SolverLabel","two non-interacting ideal hysterons")
plotFORCDistribution("forc_demo_two_hysterons_output.json", ...
    "OutputPrefix","twoHysterons_","MarkerHcHu",[0.3 0; 0.6 0])
```

Two particles with different switching fields (`Hsw = [0.3, 0.6]`),
measured together but not interacting: each one switches only on the
applied field, never on the other's state. The family plot shows a
three-level staircase (`-1`, `0`, `+1`) instead of the single
hysteron's two levels, since the two particles' switches no longer
coincide — but the distribution still shows two completely clean,
separate points, one per particle, with nothing connecting them.

That absence of coupling was checked rigorously, not just eyeballed:
because Pike's fit is linear in `M`, and this population's `M` is the
average of the two particles' individually-simulated `M`, the combined
distribution is *exactly* `0.5*(rho1 + rho2)` at every grid point
(verified to floating-point precision, `~1e-15`) — the two particles'
own smearing tails (see "A numerical surprise" above) can reach a
little way towards each other without ever really touching at this
separation, but there is no genuine cross-term. This is exactly what a
later *interacting* version of this demo would change: the two points
should start to shift, merge, or smear into a connecting ridge once
the particles can influence each other's switching field.

### Demo 4 — reversible Langevin particle

```matlab
runFORCDemoLangevin
plotFORCFamily("forc_demo_langevin_output.json", ...
    "OutputPrefix","langevin_","SolverLabel","equilibrium Langevin particle")
plotFORCDistribution("forc_demo_langevin_output.json", ...
    "OutputPrefix","langevin_")
animateFORCLangevin("forc_demo_langevin_output.json")
```

The normalized Langevin relation is `M = Ms*(coth(Alpha*H) - 1/(Alpha*H))`, evaluated with a series expansion near zero to avoid cancellation. Because this model is equilibrium and history-independent, all reversal curves overlap and the FORC distribution should be zero up to floating-point and fitting effects. A finite-relaxation version would need additional state and is a separate model.

## Visualizing results

Four plotting scripts, each reading a demo's JSON:

- **`plotFORCFamily`** — the raw family of reversal curves plus the
  major branch, i.e. the measurement the distribution is computed
  from.
- **`plotFORCDistribution`** — `rho` as a filled contour in
  `(Hc,Hu)`, diverging colormap, symmetric about zero.
- **`plotFORCSmoothingMaps`** — the same distribution at several
  smoothing factors `SF`, to see how the fit window trades resolution
  for noise.
- **`plotFORCSmoothingInfluence`** — a 4-panel diagnostic isolating
  true smoothing effects from the support-shrinkage artifact caused by
  a larger `SF` needing a bigger fit window (compares a fixed common
  support against each `SF`'s native one).

All four take `OutputFile` positionally and accept `OutputPrefix` to
namespace their output filenames — needed as soon as more than one
demo's figures are meant to coexist in `plots/` without overwriting
each other. `plotFORCFamily` also takes `SolverLabel` for its title,
and `plotFORCDistribution` takes `MarkerHcHu` to overlay theoretical
points.

## Animation

`animateFORCIdealHysteron.m` renders a two-panel video from either
ideal-hysteron demo's JSON (one particle or several): the left panel
traces the field sweeping down and each reversal curve being measured
back up, with thick dashed lines at each particle's `+-Hsw` throughout;
the right panel stays blank until every displayed curve is done (Pike's
fit needs the full measured neighborhood — there's no honest way to
show it "part way"), then reveals the FORC distribution with a matching
`Hc = Hsw` line and the theoretical point(s) marked per particle — so
it's visually unambiguous that the switching field and the peak's
coercivity are the same number, not just similarly placed.

```matlab
animateFORCIdealHysteron   % writes plots/idealHysteronAnimation.mp4

animateFORCIdealHysteron("forc_demo_two_hysterons_output.json", ...
    "OutputFile","twoHysteronsAnimation.mp4")
```

The frame rate (default 15) controls playback speed without changing
the animation's content — lower it further to slow things down more.
Pass `"AlsoWriteGif",true` for a GIF alongside the MP4 (better for
README/chat embedding; the MP4 is the one meant for slides). The final
frame is also saved as a static PDF, named after `OutputFile` (e.g.
`twoHysteronsAnimation_finalFrame.pdf`) so multiple demos' animations
don't overwrite each other's snapshot.

## Grid alignment — a gotcha that applies to any sharp-featured demo

The ideal hysteron's switching field `Hsw` must land exactly on both
the reversal-field lattice and the measurement lattice, or the delta
smears across many more grid cells than the (already unavoidable — see
below) minimum. Concretely: `Hsw` needs to be an exact multiple of
`Stride * (Hmax-Hmin)/(NumPoints-1)`. `simulateIdealHysteron.m` also
takes a floating-point `Tolerance` (default `1e-9`) around each
switching threshold, because `linspace`-generated grid points aren't
guaranteed to be bit-identical to a literal field value even when
they're mathematically the same point. Any future solver with a sharp
feature (a hard switching field, a discontinuity) will need the same
care.

### A numerical surprise: even perfect alignment doesn't give a single-cell delta

Running the default demo (`Hsw=0.3`, grid spacing `0.05`, any
`SmoothingFactor`) shows the peak of `rho` split evenly across a
**2x2 block of four grid cells**, each carrying the same magnitude,
centered at `(Hc,Hu) = (0.275, 0)` — offset from the theoretical
`(0.3, 0)` by exactly half a coercivity grid-step. This holds at
`SmoothingFactor = 1, 2` and `3` alike, so it isn't a smoothing
artifact in the usual sense (the kind `plotFORCSmoothingInfluence.m`
already documents).

The reason is geometric, not numerical error: the ideal hysteron's
switch is a genuine corner in `(Hr,H)` space (one boundary at
`Hr = -Hsw`, one at `H = +Hsw`), and when `Hsw` is chosen to land
exactly on the lattice, that corner sits at a **grid vertex shared by
four cells** rather than inside any one of them. Pike's fit is a
centered, symmetric estimator — the two grid points straddling each
boundary always see an equally strong response, because neither one's
fitting window is more "centered" on the corner than the other's. No
choice of on-grid `Hsw`, and no amount of extra smoothing, moves the
corner off that shared vertex. It's an intrinsic property of resolving
an on-lattice singularity with a centered fit, worth knowing before
reading too much into a FORC diagram's finest-scale structure.

## Requirements

MATLAB only — no toolboxes beyond base graphics, `VideoWriter`, and
`exportgraphics`.

## Roadmap

This is a planned series of demos, each swapping in a richer
`SolveFcn`:

1. **Single ideal hysteron** (done) — one particle, one point in FORC
   space.
2. **Two independent ideal hysterons** (done) — two particles, two
   points, exactly additive; sets up the contrast for interactions.
3. **Langevin-like (superparamagnetic) particle** — planned, not yet
   implemented.
4. **Interacting particles** — planned, not yet implemented; this is
   the physics this repo was originally extracted from. With step 2's
   exact-superposition result as a baseline, this is where the two
   points should start to shift, merge, or smear into a ridge.

## Provenance

The simulation engine and plotting suite were extracted from a
separate project (`meteor-dipoles`) that modeled interacting
Stoner-Wohlfarth particle packings. That project's actual physics
solver was stripped out; everything here works purely off the `H`/`M`
arrays a `SolveFcn` produces, which is what makes it safe to plug new
solvers in without touching the rest of the pipeline.
