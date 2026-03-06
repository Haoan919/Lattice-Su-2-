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

### V6 — Honeycomb Lattice Implementation with Parallelogram Periodicity

The sixth stage introduced the largest geometric extension in the project: a full **honeycomb-lattice implementation** in the spatial \(x\)-\(y\) plane. Rather than treating the system as a simple square-derived lattice with modified measurements, this stage redefined both the spatial geometry and the local update structure of the simulation.

The honeycomb lattice was implemented through a **brick-wall embedding**, which allowed the non-rectangular connectivity of the honeycomb graph to be encoded inside an array-based lattice representation. In this construction, the three link directions were defined so that
- \(U(t,x,y,0)\) connects the temporal direction,
- \(U(t,x,y,1)\) connects \((t,x,y)\) to \((t,x+1,y)\),
- \(U(t,x,y,2)\) connects \((t,x,y)\) to \((t,x-1,y+1)\),

so that the effective plaquette structure in the embedded lattice reproduces a hexagonal plaquette through a skewed parallelogram representation rather than an ordinary rectangular one. The honeycomb plaquette itself was therefore represented by an ordered product of links adapted to this embedding, not by the standard square-lattice plaquette formula.

A key point of this stage is that the final implementation was **not based on a rectangular periodic construction**. Instead, the final method used a **parallelogram-based algorithm** to impose periodicity on the honeycomb lattice. This was essential because the natural periodic identification of a honeycomb tiling is not captured correctly by a rectangular fundamental domain. A rectangular identification distorts the true adjacency relations of the honeycomb lattice, breaks the intended matching of boundary links, and leads to incorrect closed-loop structure when Wilson loops or local plaquette environments cross the boundary. In contrast, the parallelogram construction preserves the translational structure of the embedded honeycomb tiling and gives the correct periodic continuation of links and plaquettes. In other words, the final results in this version were extracted using a **parallelogram fundamental-domain algorithm**, precisely because the ordinary rectangular construction was shown to be geometrically incorrect for periodic honeycomb simulations. 

This geometric change also required a full rewrite of the **Wilson-loop measurement routines** in the spatial plane. Separate functions were introduced for different planes, including:
- `WilsonAction12` for the \(x\)-\(y\) plane,
- `WilsonAction01` for the \(t\)-\(x\) plane,
- `WilsonAction02` for the \(t\)-\(y\) plane.

These were necessary because loop construction on the honeycomb lattice is no longer described by the rectangular stepping rules used in earlier square-like versions. In the \(x\)-\(y\) plane in particular, the Wilson loop had to be redefined as a path that follows the skewed link structure of the embedded honeycomb lattice. This means that the measurement layer was rebuilt so that gauge-invariant loops remain meaningful in the new geometry rather than simply reusing square-lattice formulas.
Just as importantly, the **staple calculus itself was rewritten**. In the earlier reduced strip and square-derived versions, the local staple surrounding a link could still be built from standard nearest-neighbor plaquette logic. That is no longer valid in the honeycomb setting, because the local neighborhood of a link depends on the honeycomb connectivity and on the sublattice structure introduced by the brick-wall embedding. As a result, the update kernel had to be made geometry-aware: different link directions require different staple constructions, and for some link types the surrounding local contribution is no longer a simple square-lattice three-link complement but a longer multi-link product adapted to the effective hexagonal plaquette structure. This ensured that Heat Bath, Metropolis, and Overrelaxation updates were driven by the correct local action in the new geometry rather than by a square-lattice approximation. 

Conceptually, this stage is important for two reasons. First, it extended the project from a reduced SU(2) lattice code into a genuinely **geometry-adaptive non-Abelian gauge simulation framework**. Second, it established and justified the use of a **parallelogram periodic algorithm** as the correct way to simulate periodic honeycomb lattices. In that sense, V6 was not just an implementation step; it also included a geometric result: the demonstration that the naive rectangular construction is wrong, while the parallelogram-based fundamental domain correctly preserves the periodic structure needed for honeycomb-lattice Monte Carlo evolution and Wilson-loop measurements. 

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
