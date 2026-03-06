# SU(2) Lattice Gauge Theory Simulation Framework

A Python/Numba-based simulation framework for **SU(2) lattice gauge theory** on reduced and nonstandard Euclidean lattice geometries.

This project began as an adaptation of a conventional square-lattice Monte Carlo code into a reduced **3D Euclidean \((t,x,y)\)** setting, and gradually evolved into a more flexible platform for studying how non-Abelian gauge observables depend on **lattice extent, loop size, loop orientation, boundary conditions, and spatial geometry**.

At the core of the simulation, each lattice link is assigned an **SU(2) matrix-valued gauge variable** \(U_\mu(x)\). Starting from either a **cold start** (identity links) or a **hot start** (random SU(2) matrices), the full gauge configuration is evolved through standard Monte Carlo update algorithms. The code currently implements three update schemes:

- **Heat Bath**
- **Metropolis**
- **Overrelaxation**

Each update is performed locally using the **staple contribution** surrounding a link, allowing the program to generate equilibrium gauge configurations under the Wilson action.

The purpose of this project is not limited to computing a single Wilson loop on a fixed rectangular lattice. Instead, it was developed as a customizable numerical framework for exploring:

- reduced-width SU(2) lattice systems,
- variable lattice extents,
- Wilson loops of different sizes,
- Wilson loops on different planes,
- zero-temperature scaling behavior,
- open boundary conditions,
- and finally, a **honeycomb lattice implementation** in the spatial plane using **brick-wall mapping**.

A major feature of the final version is that the code no longer only modifies measurement routines: it also adapts the **local staple calculus** to the geometry itself. In particular, the honeycomb implementation required rewriting the staple construction so that local update dynamics follow the connectivity of the new lattice rather than inheriting square-lattice plaquette rules unchanged.

---

## What This Code Does

This repository provides a custom SU(2) lattice simulation framework with the following capabilities:

- initialization of SU(2) gauge links with either hot start or cold start,
- Monte Carlo link updates using Heat Bath, Metropolis, and Overrelaxation,
- Wilson-loop / Wilson-action measurements on multiple planes,
- support for variable loop sizes such as \(1\times1\), \(1\times2\), \(2\times1\), \(2\times2\), and \(3\times3\),
- support for different lattice extents through local parameter control,
- support for reduced strip-like geometries and open boundary conditions,
- support for honeycomb / brick-wall-mapped spatial geometry,
- and geometry-aware staple construction in the final lattice implementation.

Overall, the code is designed as a computational platform for testing how gauge observables change under controlled modifications of geometry and simulation conditions, rather than as a single-purpose Wilson-loop calculator.

---

## Development History

### V1 — Linear Plaquette Chain Adaptation

The first stage of the project adapted an original square-lattice SU(2) code into a reduced **3D Euclidean \((t,x,y)\)** lattice geometry. Instead of working in a fully standard square-lattice setting, this version focused on a **linear plaquette chain** structure with:

- \(N_x \leq 20\),
- \(N_y = 2\),
- and large \(N_t\).

This stage established the computational backbone of the project: SU(2) link initialization, link-by-link Monte Carlo evolution, and the measurement of gauge observables on a reduced lattice.

---

### V2 — Zero-Temperature Limit

In the second stage, **\(N_t\)** was refactored from a global variable into a **local function argument**, making the temporal extent of the lattice easier to control and vary.

This allowed the code to investigate the **zero-temperature limit** \(N_t \gg N_s\). The main checks in this stage included:

- observing the \(1\times1\) Wilson action as a function of \(N_t\) at fixed \(\beta\),
- and observing the \(1\times W_t\) Wilson action as a function of temporal loop length at fixed large \(N_t\).

This version marked the transition from a fixed-geometry simulation into a framework capable of probing thermal scaling behavior.

---

### V3 — Orientation-Resolved Wilson Loops

In the third stage, **\(N_x\)** was also refactored from a global variable into a local function argument. At the same time, the Wilson-action measurement routine was generalized so that loop orientation could be chosen explicitly.

Instead of restricting measurements to one default loop construction, the code could now resolve Wilson loops on different planes, such as:

- \(x\)-\(y\),
- \(t\)-\(x\),
- \(t\)-\(y\).

This made the framework sensitive not only to loop size, but also to **directional and planar dependence** of the observables.

---

### V4 — Variable Wilson-Loop Size

The next stage generalized the Wilson-loop measurement itself by allowing the loop dimensions to vary.

Rather than measuring only a fixed plaquette size, the code was extended to support rectangular Wilson loops such as:

- \(1\times1\),
- \(1\times2\),
- \(2\times1\),
- \(2\times2\),
- \(3\times3\).

This step was physically important because it enabled comparisons across different loop areas and scales. It also supported consistency checks of the expected **area-law behavior**
\[
W(R,T) \sim \exp[-\sigma(\beta)RT],
\]
by comparing the string tension extracted from loops of different sizes.

---

### V5 — Open Boundary Conditions

In the fifth stage, the lattice was modified to include **open boundary conditions in the \(y\)-direction**.

This version fixed \(N_y = 2\) while removing periodic boundary identification along that axis, turning the geometry into a reduced strip-like system with an explicitly distinguished boundary direction.

This extension made the code suitable for studying how gauge observables behave when one spatial direction is no longer wrapped periodically, which is especially relevant for reduced-width lattice constructions.

---

### V6 — Honeycomb Lattice Implementation

The sixth stage introduced the most substantial geometric modification in the project: a **honeycomb lattice implementation** in the spatial \(x\)-\(y\) plane using **brick-wall mapping**.

This stage involved several distinct changes.

First, the spatial geometry itself was restructured. Instead of using the earlier square-derived spatial organization, the \(x\)-\(y\) plane was encoded through a brick-wall representation of the honeycomb lattice, allowing the simulation to preserve honeycomb connectivity within an array-based implementation.

Second, the measurement routines were specialized by plane. Separate Wilson-action functions were introduced for different directions, including:

- `WilsonAction12` for the \(x\)-\(y\) plane,
- `WilsonAction01` for the \(t\)-\(x\) plane,
- `WilsonAction02` for the \(t\)-\(y\) plane.

This was necessary because closed-loop constructions in the honeycomb geometry are no longer described by the same rectangular stepping rules used in earlier square-like versions.

Third, and most importantly, the **staple calculus itself was rewritten**. In earlier versions, the local staple surrounding a link followed square-lattice or strip-derived adjacency rules. In the honeycomb version, this is no longer sufficient, since the local neighborhood of a link depends on the new connectivity structure and, in practice, on the sublattice organization introduced by the brick-wall mapping.

As a result, the local update environment used by Heat Bath, Metropolis, and Overrelaxation was no longer inherited from the earlier geometry. Instead, the staple construction was adapted to the honeycomb lattice itself, so that local Monte Carlo updates were consistent with the new plaquette structure and local connectivity.

This final stage transformed the project from a reduced-lattice SU(2) simulation code into a more general **geometry-adaptive lattice gauge framework**.

---

## Physics Scope

This project is intended as a numerical platform for exploring how SU(2) gauge observables depend on controlled modifications of lattice structure and simulation setup.

Current physics directions explored in the code include:

- reduced-geometry SU(2) gauge simulations,
- zero-temperature behavior,
- loop-size dependence of Wilson observables,
- area-law consistency checks,
- open-boundary effects,
- and non-Abelian gauge dynamics on honeycomb-like spatial lattices.

---
