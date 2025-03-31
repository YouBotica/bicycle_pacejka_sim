# Plan: Python Bicycle Simulation with Pacejka Tire Model

This document outlines the plan for creating a Python simulation of a dynamic bicycle model using Pacejka tire formulas based on provided parameters and a `.tir` file.

**User Provided Data:**

*   **Bicycle Parameters:**
    *   Mass (`m`): 800 kg
    *   Yaw Moment of Inertia (`Izz`): 1000 kg·m²
    *   Distance CoG to front axle (`a`): 1.6 m
    *   Distance CoG to rear axle (`b`): 1.2 m
    *   Wheelbase (`L = a + b`): 2.8 m
*   **Tire Data:** Coefficients to be extracted from `lf_tire.tir` (MF6.2 format).
*   **Libraries:** NumPy, SciPy, Matplotlib.
*   **Initial Scenario:** Constant radius turn.

**Phase 1: Foundation & Core Components**

1.  **Project Setup:**
    *   Create directory: `bicycle_pacejka_sim`
    *   Create files: `bicycle_params.py`, `pacejka_formulas.py`, `bicycle_model.py`, `simulation.py`, `tire_parser.py` (optional, for cleaner coefficient handling).
2.  **Define Bicycle Parameters (`bicycle_params.py`):**
    *   Store known parameters: `m=800`, `Izz=1000`, `a=1.6`, `b=1.2`, `L=2.8`.
    *   Include placeholders or typical values for CoG height (`h`), wheel radius (`R`), track width (`w`).
3.  **Process Pacejka Coefficients (`pacejka_formulas.py` or `tire_parser.py`):**
    *   Parse the `.tir` file content (already read).
    *   Extract coefficients for Fx (PCX1, PDX1...), Fy (PCY1, PDY1...), and Mz (QBZ1, QDZ1...).
    *   Store these coefficients in a structured way (e.g., Python dictionaries or a dedicated class).
4.  **Implement Pacejka Formula Functions (`pacejka_formulas.py`):**
    *   Create `calculate_Fx(slip_ratio, Fz, coefficients)`
    *   Create `calculate_Fy(slip_angle, Fz, coefficients)`
    *   Create `calculate_Mz(slip_angle, Fz, coefficients)`
    *   These functions will use the extracted coefficients and NumPy for calculations based on the MF6.2 standard indicated in the file.
5.  **Define State Variables:** Identify `u` (longitudinal velocity), `v` (lateral velocity), `r` (yaw rate), `X` (global X position), `Y` (global Y position), `psi` (yaw angle).

**Phase 2: Dynamic Model & Simulation**

6.  **Formulate Equations of Motion (`bicycle_model.py`):**
    *   Develop the function/class calculating state derivatives (`du/dt`, `dv/dt`, `dr/dt`, `dX/dt`, `dY/dt`, `dpsi/dt`).
    *   This involves calculating tire kinematics (slip angles, slip ratios), vertical loads (`Fz`), calling the Pacejka functions, and applying Newton-Euler equations using the parameters from `bicycle_params.py`.
7.  **Numerical Integration (`simulation.py`):**
    *   Implement the main simulation loop.
    *   Use `scipy.integrate.solve_ivp` to integrate the equations of motion over time.
8.  **Input Handling:**
    *   Define inputs for a constant radius turn scenario (e.g., constant steering angle derived from desired radius and speed, initial velocity).
9.  **Output & Visualization:**
    *   Store and plot key results (trajectory X vs Y, velocities `u`,`v`,`r` over time, tire forces Fx/Fy over time) using Matplotlib.

**Conceptual Flow Diagram:**

```mermaid
graph TD
    A[Start Simulation] --> B(Load Parameters & Coefficients);
    B --> C(Initialize State Vector [u, v, r, X, Y, psi]);
    C --> D{Simulation Loop (for each time step)};
    D --> E[Get Inputs (Steering δ, Torque Tx)];
    E --> F[Calculate Tire Kinematics (Slip Angles α, Slip Ratios κ)];
    F --> G[Calculate Vertical Loads (Fz_front, Fz_rear)];
    G --> H[Call Pacejka Functions (Fx, Fy, Mz for each tire)];
    H --> I[Calculate State Derivatives (du/dt, dv/dt, dr/dt, ...)];
    I --> J[Integrate Equations of Motion (Update State Vector)];
    J --> K{End Time Reached?};
    K -- No --> D;
    K -- Yes --> L[Store/Plot Results];
    L --> M[End Simulation];

    subgraph Bicycle Parameters
        P1[Mass, Inertia, Geometry]
    end

    subgraph Pacejka Coefficients
        P2[Extracted from lf_tire.tir]
    end

    subgraph Pacejka Formulas Module
        PF1[Inputs: α, κ, Fz, Coeffs] --> PF2{Pacejka Functions};
        PF2 --> PF3[Outputs: Fx, Fy, Mz];
    end

    B --> P1;
    B --> P2;
    F & G & P2 --> PF1;
    PF3 --> H;
```

**Next Steps:**

*   Proceed with implementation, likely starting with Phase 1, Step 1 (Project Setup) and Step 3 (Parsing Coefficients).