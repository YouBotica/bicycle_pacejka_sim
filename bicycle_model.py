"""
Implements the dynamic bicycle model equations of motion.
Uses parameters from bicycle_params.py and tire forces from pacejka_formulas.py.
"""
import numpy as np
import sys
import os

# Add the directory containing the modules to the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

try:
    import bicycle_params as params
    import pacejka_formulas as pf
except ImportError as e:
    print(f"Error importing necessary modules: {e}")
    print("Ensure bicycle_params.py and pacejka_formulas.py are in the same directory.")
    sys.exit(1)

def bicycle_dynamics(t, state, delta, Fx_drive_rear=0.0):
    """
    Calculates the state derivatives for the bicycle model.

    Args:
        t (float): Current time (required by solve_ivp, but not used in this model).
        state (list or np.array): Current state vector [u, v, r, X, Y, psi]
            u: longitudinal velocity (m/s)
            v: lateral velocity (m/s)
            r: yaw rate (rad/s)
            X: global X position (m)
            Y: global Y position (m)
            psi: yaw angle (rad)
        delta (float): Steering angle of the front wheel (radians).
        Fx_drive_rear (float): Longitudinal drive/brake force applied at rear wheel center (N).
                                Positive for drive, negative for braking. Defaults to 0.

    Returns:
        np.array: Array of state derivatives [du_dt, dv_dt, dr_dt, dX_dt, dY_dt, dpsi_dt].
    """
    # Unpack state vector
    u, v, r, X, Y, psi = state

    # --- Calculate Vertical Loads (Fz) ---
    # Simple static load distribution (neglecting pitch/roll dynamics for now)
    Fz_f = params.m * params.g * params.b / params.L
    Fz_r = params.m * params.g * params.a / params.L
    # Add dynamic load transfer if needed (e.g., based on longitudinal/lateral acceleration)
    # Fz_f = Fz_f - params.m * ax * params.h / params.L # Example longitudinal transfer
    # Fz_r = Fz_r + params.m * ax * params.h / params.L # Example longitudinal transfer

    # --- Calculate Tire Kinematics ---
    # Slip Angles (alpha)
    # Avoid division by zero if u is very small
    if abs(u) < 0.1: # Threshold velocity to avoid instability at rest
        alpha_f = 0.0
        alpha_r = 0.0
    else:
        alpha_f = np.arctan((v + params.a * r) / u) - delta
        alpha_r = np.arctan((v - params.b * r) / u)

    # Slip Ratios (kappa) - Simplified: Assumes wheel speed matches vehicle speed adjusted by slip
    # For this version, we directly use the input drive force and assume kappa=0 for Fy/Mz calculations,
    # but calculate Fx based on the drive force for the rear and assume Fx_f=0 (no front drive/brake).
    # A more complex model would calculate kappa based on wheel torque and rotational dynamics.
    kappa_f = 0.0
    kappa_r = 0.0 # Placeholder for Fy/Mz calculations

    # --- Calculate Tire Forces using Pacejka Formulas ---
    # Assuming zero camber (gamma = 0) for simplicity
    gamma = 0.0

    # Need to calculate forces iteratively or make simplifying assumptions for combined slip
    # Approach: Calculate Fy first assuming kappa=0, then calculate Fx.

    # Front Tire
    Fy_f_calc = pf.calculate_Fy(alpha_f, Fz_f, gamma, kappa_f,
                                params.tire_coeffs['front']['LATERAL'],
                                params.tire_coeffs['front']['LONGITUDINAL'],
                                params.tire_coeffs['front']['ALIGNING'],
                                0.0) # Assume Fx_f = 0 for Fy calculation input
    Fx_f_calc = 0.0 # No drive/brake force on front wheel in this simple model

    # Rear Tire
    Fy_r_calc = pf.calculate_Fy(alpha_r, Fz_r, gamma, kappa_r,
                                params.tire_coeffs['rear']['LATERAL'],
                                params.tire_coeffs['rear']['LONGITUDINAL'],
                                params.tire_coeffs['rear']['ALIGNING'],
                                Fx_drive_rear) # Use actual drive force for Fy calculation input
    # Fx_r is primarily the drive force, but Pacejka Fx can model rolling resistance/induced drag
    # For simplicity, we'll use the input drive force directly in EoM,
    # but could calculate Pacejka Fx_r if needed (e.g., for tire limits)
    Fx_r_calc = Fx_drive_rear # Use the input drive force directly

    # Aligning Moments (can be calculated but not directly used in 3DOF model)
    # Mz_f = pf.calculate_Mz(alpha_f, Fz_f, gamma, kappa_f, params.tire_coeffs['front']['ALIGNING'], ...)
    # Mz_r = pf.calculate_Mz(alpha_r, Fz_r, gamma, kappa_r, params.tire_coeffs['rear']['ALIGNING'], ...)

    # --- Equations of Motion (Newton-Euler) ---
    # Sum of forces in longitudinal direction (x-body)
    du_dt = (Fx_f_calc * np.cos(delta) - Fy_f_calc * np.sin(delta) + Fx_r_calc) / params.m + v * r

    # Sum of forces in lateral direction (y-body)
    dv_dt = (Fx_f_calc * np.sin(delta) + Fy_f_calc * np.cos(delta) + Fy_r_calc) / params.m - u * r

    # Sum of moments about CoG (z-axis)
    dr_dt = ( (Fx_f_calc * np.sin(delta) + Fy_f_calc * np.cos(delta)) * params.a - Fy_r_calc * params.b ) / params.Izz

    # --- Kinematic Equations (Global Frame) ---
    dX_dt = u * np.cos(psi) - v * np.sin(psi)
    dY_dt = u * np.sin(psi) + v * np.cos(psi)
    dpsi_dt = r

    return np.array([du_dt, dv_dt, dr_dt, dX_dt, dY_dt, dpsi_dt])

# Example usage (for testing):
# initial_state = [10.0, 0.0, 0.0, 0.0, 0.0, 0.0] # u, v, r, X, Y, psi
# steering_angle = np.radians(5.0) # 5 degrees steering
# derivatives = bicycle_dynamics(0, initial_state, steering_angle)
# print(derivatives)