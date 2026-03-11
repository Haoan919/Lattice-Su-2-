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
- and finally, a **honeycomb lattice implementation** in the spatial plane using a **brick-wall / parallelogram mapping**.

A major feature of the final version is that the code no longer only modifies measurement routines: it also adapts the **local staple calculus** to the geometry itself. In particular, the honeycomb implementation required rewriting the staple construction so that local update dynamics follow the connectivity of the new lattice rather than inheriting square-lattice plaquette rules unchanged.

---

## Overview

This repository is not just a single SU(2) Monte Carlo script, but a **progressively developed lattice gauge simulation framework**.

The project starts from a conventional reduced Euclidean lattice code and gradually extends it into a platform capable of handling:

- reduced-width strip-like systems,
- variable lattice extents,
- orientation-dependent Wilson loops,
- variable loop sizes,
- open boundary conditions,
- minimal transverse constructions,
- and finally a **geometry-aware honeycomb lattice simulation**.

The central goal is to study **how SU(2) gauge observables change when the lattice geometry itself is modified**, rather than only changing coupling parameters inside one fixed square-lattice setup.

In that sense, the repository documents a real **method-development path**: each version introduces a new geometric or algorithmic idea, while preserving the same underlying gauge-field evolution logic.

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

## Project Structure

This repository is organized as a sequence of related simulation stages rather than a single monolithic script.

- **Baseline reduced Euclidean implementation**  
  The project begins with a reduced **3D Euclidean \((t,x,y)\)** SU(2) lattice Monte Carlo framework, containing link initialization, Monte Carlo evolution, and basic Wilson-action measurements.

- **Thermal-scaling and extent-control versions**  
  Later versions promote lattice extents such as \(N_t\) and \(N_x\) to local function arguments, making the framework flexible enough to study zero-temperature limits and controlled scaling behavior.

- **Orientation-resolved measurement versions**  
  The code was extended to support Wilson loops on different planes, including \(x\)-\(y\), \(t\)-\(x\), and \(t\)-\(y\), allowing directional dependence of observables to be studied explicitly.

- **Variable Wilson-loop size versions**  
  Additional stages generalize the Wilson-loop size itself, allowing rectangular loops of multiple sizes rather than a single plaquette observable.

- **Finite-width and open-boundary geometries**  
  Several versions specialize to very small transverse width, especially \(N_y = 2\), creating strip-like systems in which one spatial direction is intentionally constrained. These versions are used to study open-boundary effects and reduced-width behavior.

- **Minimal transverse identification construction**  
  One stage introduces a special initialization-based field identification scheme that makes it possible to realize an effective \(y\)-length of 1 even though naive length-1 periodicity cannot be directly imposed in the original data structure.

- **Honeycomb geometry via brick-wall / parallelogram representation (`physics6`)**  
  The final stage represents a **honeycomb lattice** through a **parallelogram-based periodic construction** rather than a naive rectangular one. This is important because a honeycomb lattice is not naturally described by a standard 90-degree Cartesian \(x\)-\(y\) coordinate system. Its intrinsic spatial structure is skewed. The core difficulty is therefore computational as well as geometric: array-based simulation code naturally assumes orthogonal indexing, while the true honeycomb connectivity is oblique. The solution in this project is to embed the honeycomb structure into a vertically indexed array through a **brick-wall / parallelogram mapping**, which preserves the skewed adjacency relations inside an indexable lattice representation. In this construction, the middle connecting point is effectively removed, and the `%2` dependence in the code is used to encode the alternating sublattice pattern. This makes the final geometry equivalent to a honeycomb-type lattice rather than an ordinary rectangular grid.

Together, these components show that the repository is not just a collection of scripts, but a genuine **SU(2) lattice method-development pipeline**.

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

## V5 — Open Boundary Conditions and the `Ny = 2, Ly = 1` Construction

In the fifth stage, I introduced **open boundary conditions in the y-direction** and constructed a reduced transverse geometry with **Ny = 2** but **effective y-length = 1**.

The main difficulty was that although the y-direction could be restricted in size, a **true length-1 periodic identification could not be implemented directly** in the original lattice setup. A naive treatment would fail to realize the desired periodic structure at that minimal scale.

To resolve this, I introduced a mathematical construction at the level of the **initial field configuration**: the newly extended part of the lattice was assigned the **same field values as the original part**. This effectively enforced the required identification from the start, allowing the simulation to reproduce the intended length-1 behavior without relying on a naive periodic rule.

As a result, this version retained **two y-indexed layers** in the data structure, but the actual geometric extent in the y-direction was reduced to its minimal nontrivial form. This made the system a controlled intermediate case between the earlier reduced Euclidean lattices and the later honeycomb-lattice implementation, while also demonstrating a practical method for realizing otherwise inaccessible minimal periodic structure.

---

### V6 — Honeycomb Lattice Implementation with Parallelogram Periodicity

The sixth stage introduced the largest geometric extension in the project: a full **honeycomb-lattice implementation** in the spatial \(x\)-\(y\) plane. Rather than treating the system as a simple square-derived lattice with modified measurements, this stage redefined both the spatial geometry and the local update structure of the simulation.

The honeycomb lattice was implemented through a **brick-wall embedding** together with a **parallelogram fundamental-domain algorithm**. This was necessary because the honeycomb lattice is not naturally built from a standard 90-degree rectangular \(x\)-\(y\) coordinate system. Its true geometry is skewed: the natural translation directions are oblique rather than orthogonal. A naive rectangular treatment may be convenient for array indexing, but it does **not** preserve the real adjacency structure of the honeycomb tiling.

The key difficulty is that computers naturally interpret lattice coordinates through vertical and horizontal indexing. In other words, an array-based implementation wants to think in terms of rectangular coordinates, while the true honeycomb lattice lives on a skewed spatial structure. The central geometric idea of this stage was therefore to **transform that skewed geometry into a computationally manageable representation without destroying its connectivity**.

This was achieved by expressing the honeycomb lattice through a **parallelogram-based representation** inside a regular array structure. In effect, the skewed geometry was encoded into the indexing rules rather than approximated away. The construction is closely related to a brick-wall mapping, but the crucial point is that the final periodicity was imposed through a **parallelogram fundamental domain**, not through an ordinary rectangular one.

In this implementation, the three link directions were defined so that:

- \(U(t,x,y,0)\) connects the temporal direction,
- \(U(t,x,y,1)\) connects \((t,x,y)\) to \((t,x+1,y)\),
- \(U(t,x,y,2)\) connects \((t,x,y)\) to \((t,x-1,y+1)\).

This produces a skewed local connectivity pattern whose plaquette structure corresponds to the embedded honeycomb geometry. In addition, the middle connecting point of the naive rectangular picture is effectively removed, and the `%2` dependence in the code is used to encode the alternating sublattice structure. That is exactly what makes the geometry honeycomb-like rather than simply rectangular with modified labels.

A key point of this stage is that the final implementation was **not based on a rectangular periodic construction**. The ordinary rectangular approach is geometrically wrong for the honeycomb case because it forces the lattice into an orthogonal interpretation that does not match the true skewed translation structure. This leads to incorrect boundary matching, distorted adjacency relations, and wrong closed-loop behavior when Wilson loops or local plaquette environments cross the boundary.

By contrast, the **parallelogram construction** preserves the translational pattern of the embedded honeycomb tiling and provides the correct periodic continuation of links and plaquettes. In other words, the parallelogram approach is the mathematically correct way to let a computer simulate a lattice whose natural coordinate directions are slanted rather than perpendicular.

This geometric change also required a full rewrite of the **Wilson-loop measurement routines** in the spatial plane. Separate functions were introduced for different planes, including:

- `WilsonAction12` for the \(x\)-\(y\) plane,
- `WilsonAction01` for the \(t\)-\(x\) plane,
- `WilsonAction02` for the \(t\)-\(y\) plane.

These were necessary because loop construction on the honeycomb lattice is no longer described by the rectangular stepping rules used in earlier square-like versions. In the \(x\)-\(y\) plane in particular, the Wilson loop had to be redefined as a path that follows the skewed link structure of the embedded honeycomb lattice. This means that the measurement layer was rebuilt so that gauge-invariant loops remain meaningful in the new geometry rather than simply reusing square-lattice formulas.

Just as importantly, the **staple calculus itself was rewritten**. In the earlier reduced strip and square-derived versions, the local staple surrounding a link could still be built from standard nearest-neighbor plaquette logic. That is no longer valid in the honeycomb setting, because the local neighborhood of a link depends on the honeycomb connectivity and on the sublattice structure introduced by the brick-wall embedding. As a result, the update kernel had to be made geometry-aware: different link directions require different staple constructions, and for some link types the surrounding local contribution is no longer a simple square-lattice three-link complement but a longer multi-link product adapted to the effective hexagonal plaquette structure.

This ensured that Heat Bath, Metropolis, and Overrelaxation updates were driven by the correct local action in the new geometry rather than by a square-lattice approximation.

Conceptually, this stage is important for two reasons. First, it extended the project from a reduced SU(2) lattice code into a genuinely **geometry-adaptive non-Abelian gauge simulation framework**. Second, it established and justified the use of a **parallelogram periodic algorithm** as the correct way to simulate periodic honeycomb lattices.

In that sense, V6 was not just an implementation step; it also included a geometric result: the demonstration that the naive rectangular construction is wrong, while the parallelogram-based fundamental domain correctly preserves the periodic structure needed for honeycomb-lattice Monte Carlo evolution and Wilson-loop measurements.

---

## Validation

The framework was developed with validation through **algorithmic redundancy, observable comparison, and controlled geometric extension**.

1. **Multiple update algorithms were implemented for the same gauge system.**  
   The code includes Heat Bath, Metropolis, and Overrelaxation update schemes for the same SU(2) lattice variables. This makes it possible to compare equilibrium behavior across different local update dynamics rather than relying on a single sampler.

2. **Wilson loops were measured at multiple sizes and on multiple planes.**  
   The framework supports loop measurements such as \(1\times1\), \(1\times2\), \(2\times1\), \(2\times2\), and \(3\times3\), and also distinguishes between \(x\)-\(y\), \(t\)-\(x\), and \(t\)-\(y\) planes. This provides internal consistency checks across scale and orientation.

3. **The same Monte Carlo logic was ported across geometries.**  
   Instead of rewriting the project from scratch for each stage, the same basic gauge-update and observable structure was preserved while the geometry was progressively modified. This allows consistency to be checked through controlled extension.

4. **Parameter dependence was made explicit.**  
   Several stages promote \(N_t\), \(N_x\), loop dimensions, and other lattice parameters into local variables, making it possible to test whether observables change smoothly and sensibly under controlled variation.

5. **Geometric constraints were encoded directly into the simulation logic.**  
   In V5, the field-identification rule was imposed directly in the initial configuration to realize the minimal transverse construction. In V6, the staple calculus and Wilson-loop measurement rules were rewritten to follow the actual honeycomb connectivity. This makes the geometric assumptions explicit and testable rather than only conceptual.

Overall, the repository is validated not only by output values, but by the fact that **the same physical observables can be tracked consistently while the lattice geometry is systematically changed**.

---

## My Contributions

My contributions in this project include both **computational implementation** and **geometric / methodological extension**.

- I implemented SU(2) gauge-field Monte Carlo updates using **Heat Bath**, **Metropolis**, and **Overrelaxation** methods in reduced Euclidean lattice settings.
- I built Wilson-loop and Wilson-action measurement routines for several geometries, including different loop sizes and different planes.
- I refactored lattice extents such as \(N_t\) and \(N_x\) into locally controllable parameters, allowing systematic study of zero-temperature scaling and lattice-size dependence.
- I generalized the framework from fixed plaquette measurements to variable Wilson loops such as \(1\times1\), \(1\times2\), \(2\times1\), \(2\times2\), and \(3\times3\).
- I developed narrow-width transverse models with \(N_y = 2\) and used them to study open-boundary effects and minimal transverse geometry.
- I introduced an initialization-based field-identification rule that allows an effective \(y\)-length of 1 to be realized even though naive length-1 periodicity cannot be directly implemented in the original data structure.
- I implemented a constrained minimal-width geometry in which certain loop extensions are explicitly restricted so that the observable definitions remain consistent with the intended lattice structure.
- In the final stage, I constructed a **honeycomb lattice representation through a parallelogram / brick-wall periodic domain**, using `%2`-dependent structure and removal of the naive middle connection to reproduce true honeycomb connectivity inside an array-based simulation framework.
- I rewrote both the Wilson-loop measurement layer and the local staple calculus so that Monte Carlo updates follow the correct geometry of the honeycomb lattice instead of inheriting square-lattice approximations.

Taken together, these contributions reflect not only coding work, but a sustained effort to turn SU(2) lattice gauge theory into a flexible framework for exploring **nonstandard lattice geometry**.

---

## Key Results / Why It Matters

The main result of this repository is the construction of a **progressive SU(2) lattice gauge simulation framework** that remains computationally coherent while moving beyond the most standard square or rectangular periodic lattice setting.

The project shows that one can begin with a reduced Euclidean SU(2) Monte Carlo code and systematically extend it into reduced-width geometries, orientation-resolved Wilson-loop measurements, open-boundary systems, minimal transverse constructions, and finally a **honeycomb-type lattice represented through parallelogram periodicity**, without abandoning the core gauge-update logic.

A particularly important result is that the final honeycomb implementation does **not** simply rename a rectangular lattice. It addresses a genuine geometric problem: a honeycomb lattice is naturally skewed, while a computer array is naturally orthogonal. The solution developed here is to encode the slanted honeycomb geometry inside a regular array through a **parallelogram fundamental-domain representation**. This is significant because it preserves the true adjacency and periodic structure of the lattice instead of forcing it into a mathematically incorrect rectangular interpretation.

A second important result is methodological. The project does not treat geometric modification as a superficial change to measurement routines alone. In the final stage, the **staple calculus**, the **loop definitions**, and the **boundary identification rules** are all adapted to the geometry itself. This means the framework is not merely measuring different paths on the same old lattice; it is performing genuinely geometry-aware SU(2) Monte Carlo evolution.

More broadly, this repository matters because it documents a real development path from standard reduced-lattice simulation toward a **geometry-adaptive non-Abelian gauge framework**. It shows how gauge observables can be studied under controlled changes in dimensionality, extent, boundary conditions, and lattice connectivity. In particular, it includes a geometric conclusion of independent interest: for periodic honeycomb simulations, the naive rectangular construction is wrong, while the **parallelogram-based fundamental domain** correctly preserves the lattice structure needed for Wilson-loop measurements and local Monte Carlo updates.

---

## Physics Scope

This project is intended as a numerical platform for exploring how SU(2) gauge observables depend on controlled modifications of lattice structure and simulation setup.

Current physics directions explored in the code include:

- reduced-geometry SU(2) gauge simulations,
- zero-temperature behavior,
- loop-size dependence of Wilson observables,
- area-law consistency checks,
- open-boundary effects,
- minimal transverse-width constructions,
- and non-Abelian gauge dynamics on honeycomb-like spatial lattices.

---

## How to Run

This repository consists of several development-stage scripts. In general, each version can be run by executing the corresponding Python file after setting the desired lattice parameters inside the script.

Typical dependencies include:

- `numpy`
- `numba`
- `matplotlib`

A typical execution workflow is:

```bash
python your_script.py,
- area-law consistency checks,
- open-boundary effects,
- and non-Abelian gauge dynamics on honeycomb-like spatial lattices.

---
