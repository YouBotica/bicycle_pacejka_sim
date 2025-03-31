"""
Main simulation script for the bicycle model using Pacejka tires.
Integrates the equations of motion and plots the results.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sys
import os

# Add the directory containing the modules to the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

try:
    from bicycle_model import bicycle_dynamics
    import bicycle_params as params
except ImportError as e:
    print(f"Error importing necessary modules: {e}")
    print("Ensure bicycle_params.py and bicycle_model.py are in the same directory.")
    sys.exit(1)

# --- Simulation Parameters ---
t_start = 0.0
t_end = 15.0  # Simulation duration (seconds)
# dt = 0.01   # Time step for saving output (solve_ivp adapts its internal step)
# t_eval = np.arange(t_start, t_end, dt) # Specific times to evaluate solution

# --- Initial Conditions ---
u0 = 15.0      # Initial longitudinal velocity (m/s) -> approx 54 km/h
v0 = 0.0       # Initial lateral velocity (m/s)
r0 = 0.0       # Initial yaw rate (rad/s)
X0 = 0.0       # Initial global X position (m)
Y0 = 0.0       # Initial global Y position (m)
psi0 = 0.0     # Initial yaw angle (rad)

initial_state = [u0, v0, r0, X0, Y0, psi0]

# --- Input Scenario: Constant Radius Turn ---
# Define target turn radius and calculate approximate steady-state steering angle
target_radius = 50.0 # meters
# Ackermann steering angle approximation (low speed)
# delta_approx = params.L / target_radius
# High-speed approximation requires considering understeer gradient (needs cornering stiffness)
# For now, let's use a fixed steering angle and see the resulting turn.
delta_steer_deg = 3.0 # degrees
delta_steer_rad = np.radians(delta_steer_deg)

# Drive force (set to zero for constant speed coasting turn, or add value)
Fx_drive = 0.0 # Newtons

# --- Run Simulation ---
print("Starting simulation...")
sol = solve_ivp(
    fun=lambda t, y: bicycle_dynamics(t, y, delta=delta_steer_rad, Fx_drive_rear=Fx_drive),
    t_span=[t_start, t_end],
    y0=initial_state,
    method='RK45',  # Common choice, others like 'LSODA' or 'BDF' available
    # t_eval=t_eval # Use t_eval if specific output times are needed
    dense_output=True # Allows evaluating solution at any point
)
print("Simulation finished.")

# Check if integration was successful
if not sol.success:
    print(f"Integration failed: {sol.message}")
    sys.exit(1)

# Extract results at desired time points
t_plot = np.linspace(t_start, t_end, 500) # 500 points for smooth plotting
states = sol.sol(t_plot) # Evaluate solution at t_plot points

u, v, r, X, Y, psi = states

# --- Plot Results ---
print("Plotting results...")

plt.figure(figsize=(12, 10))

# Plot 1: Trajectory (X vs Y)
plt.subplot(2, 2, 1)
plt.plot(Y, X, label='Trajectory') # Plot Y vs X to match typical vehicle coordinates
plt.xlabel("Y Position (m)")
plt.ylabel("X Position (m)")
plt.title("Vehicle Trajectory")
plt.grid(True)
plt.axis('equal')
plt.legend()

# Plot 2: Velocities vs Time
plt.subplot(2, 2, 2)
plt.plot(t_plot, u, label='Longitudinal Velocity (u)')
plt.plot(t_plot, v, label='Lateral Velocity (v)')
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.title("Velocities vs Time")
plt.grid(True)
plt.legend()

# Plot 3: Yaw Rate vs Time
plt.subplot(2, 2, 3)
plt.plot(t_plot, np.degrees(r), label='Yaw Rate (r)') # Convert to degrees/s for readability
plt.xlabel("Time (s)")
plt.ylabel("Yaw Rate (deg/s)")
plt.title("Yaw Rate vs Time")
plt.grid(True)
plt.legend()

# Plot 4: Yaw Angle vs Time
plt.subplot(2, 2, 4)
plt.plot(t_plot, np.degrees(psi), label='Yaw Angle (psi)') # Convert to degrees
plt.xlabel("Time (s)")
plt.ylabel("Yaw Angle (deg)")
plt.title("Yaw Angle vs Time")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

print("Script finished.")