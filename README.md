# Lattice Su(2)
## Version History

### V5: Honeycomb Lattice Implementation
*   **Refactoring:** Refactored `Ny` from a global variable to a local function argument.
*   **Geometry:** Transformed the lattice geometry in the x-y plane to a **Honeycomb Lattice** (implemented via Brick Wall mapping).
*   **Measurement:** Implemented specialized functions (`WilsonAction12`, `WilsonAction01`, and `WilsonAction02`) to calculate Wilson Actions on the x-y, t-x, and t-y planes, respectively.

### V4: Open Boundary Conditions
*   **Geometry:** Set `Ny=2` and removed periodic boundary conditions in the y-direction (Open Boundary Conditions).

### V3: Area Law Verification
*   **Refactoring:** Refactored `Nx` from a global variable to a local function argument.
*   **Functionality:** Updated the `WilsonAction` function to accept `nu` and `mu` parameters, allowing specification of the loop orientation.
*   **Physics Check:** Verified the Area Law $W(R,T) \approx \exp(-\sigma(\beta)RT)$ by comparing the string tension $\sigma(\beta)$ extracted from $1 \times 1$, $2 \times 2$, and $3 \times 3$ loops to ensure consistency.

### V2: Zero Temperature Limit
*   **Refactoring:** Refactored `Nt` from a global variable to a local function argument.
*   **Physics Check:** Investigated the zero-temperature limit ($N_t \gg N_s$) by:
    1.  Observing the $1 \times 1$ Wilson Action as a function of $N_t$ at fixed $\beta=2.5$.
    2.  Observing the $1 \times W_t$ Wilson Action as a function of temporal loop length $W_t$ with fixed $N_t=100$.

### V1: Linear Plaquette Chain Adaptation
*   **Initialization:** Adapted the original code from a 4-dimensional square lattice to a linear plaquette chain geometry (3-dimensional Euclidean: $x, y, t$).
*   **Configuration:** Setup initial dimensions with $N_x \le 20$, $N_y=2$, and large $N_t$ (starting at 100).
