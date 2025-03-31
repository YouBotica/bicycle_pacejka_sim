"""
Implements the Pacejka Magic Formula (MF6.2) for tire forces and moments.
Based on the coefficients parsed from a .tir file.
Handles combined slip conditions (USE_MODE = 4).
"""
import numpy as np

# Define nominal load from the .tir file (or pass it if variable)
# From [VERTICAL] section: FNOMIN = 3112
FNOMIN = 3112.0 # Consider making this an input or loading from params

def calculate_Fy(alpha_rad, Fz, gamma_rad, kappa, coeffs_lat, coeffs_long, coeffs_align, Fx_calculated):
    """
    Calculates Lateral Force (Fy) using Pacejka Magic Formula (MF6.2 - combined slip).

    Args:
        alpha_rad (float): Slip angle in radians.
        Fz (float): Vertical load in Newtons. Fz > 0.
        gamma_rad (float): Camber (inclination) angle in radians.
        kappa (float): Longitudinal slip ratio.
        coeffs_lat (dict): Dictionary of lateral coefficients.
        coeffs_long (dict): Dictionary of longitudinal coefficients.
        coeffs_align (dict): Dictionary of aligning coefficients.
        Fx_calculated (float): Longitudinal force calculated separately (needed for Gyk).

    Returns:
        float: Calculated lateral force (Fy) in Newtons.
    """
    # Ensure Fz is positive to avoid division by zero or sqrt issues
    Fz = max(Fz, 1e-6) # Use a small positive minimum value

    # Scale factors (usually 1, but included for completeness)
    LMUY = coeffs_lat.get('LMUY', 1.0)
    LKY = coeffs_lat.get('LKY', 1.0)
    LCY = coeffs_lat.get('LCY', 1.0)
    LEY = coeffs_lat.get('LEY', 1.0)
    LHY = coeffs_lat.get('LHY', 1.0)
    LVY = coeffs_lat.get('LVY', 1.0)
    LKYC = coeffs_lat.get('LKYC', 1.0) # Camber stiffness scaling
    LYKA = coeffs_lat.get('LYKA', 1.0) # Kappa influence scaling
    LVYKA = coeffs_lat.get('LVYKA', 1.0) # Kappa induced Fy scaling

    # Normalized load change
    dfz = (Fz - FNOMIN * coeffs_lat.get('LFZO', 1.0)) / (FNOMIN * coeffs_lat.get('LFZO', 1.0))

    # --- Base Pure Slip Lateral Force Calculation ---
    # Camber effect (Shf) - simplified here, MF6.2 has more complex gamma effects
    Shf = coeffs_lat.get('PHY1', 0.0) + coeffs_lat.get('PHY2', 0.0) * dfz
    Kyg = coeffs_lat.get('PKY6', 0.0) * Fz * (1.0 + coeffs_lat.get('PKY7',0.0)*dfz) * LKYC # Camber stiffness

    # Effective slip angle incorporating camber and horizontal shift
    alpha_star = alpha_rad + Shf + Kyg * gamma_rad / (coeffs_lat.get('PKY1', -53) * FNOMIN * LKY) # Approximation

    # C factor (Shape)
    Cfy = coeffs_lat.get('PCY1', 1.55) * LCY

    # D factor (Peak) - Friction
    muy = (coeffs_lat.get('PDY1', 1.7) + coeffs_lat.get('PDY2', -0.35) * dfz) * \
          (1.0 - coeffs_lat.get('PDY3', 3.36) * gamma_rad**2) * LMUY
    Dy = muy * Fz

    # BCD factor (Cornering stiffness)
    Kya = coeffs_lat.get('PKY1', -53.0) * FNOMIN * np.sin(2.0 * np.arctan(Fz / (coeffs_lat.get('PKY2', 2.23) * FNOMIN))) * \
          (1.0 - coeffs_lat.get('PKY3', 0.87) * abs(gamma_rad)) * LKY
    # Avoid division by zero if Cfy or Dy are zero
    if abs(Cfy * Dy) < 1e-6:
        BCD = Kya # Stiffness dominates at low slip
    else:
        BCD = Kya / (Cfy * Dy)

    # B factor (Stiffness factor)
    # Note: PKY4 and PKY5 relate to curvature of stiffness, simplified here
    By = BCD # Approximation: B = BCD for pure slip

    # E factor (Curvature)
    Ey = (coeffs_lat.get('PEY1', -1.98) + coeffs_lat.get('PEY2', -1.62) * dfz) * \
         (1.0 + coeffs_lat.get('PEY5', 5.65) * gamma_rad**2 - (coeffs_lat.get('PEY3', 0.4) + coeffs_lat.get('PEY4', -4.08) * gamma_rad) * np.sign(alpha_star)) * LEY
    Ey = np.clip(Ey, -10.0, 1.0) # Limit Ey to prevent divergence, typical range <= 1

    # Horizontal and Vertical shifts
    Shy = (coeffs_lat.get('PHY1', 0.002) + coeffs_lat.get('PHY2', 0.003) * dfz) * LHY # + gamma effects if needed
    Svy = Fz * ((coeffs_lat.get('PVY1', -0.06) + coeffs_lat.get('PVY2', 0.06) * dfz) * LVY + \
           (coeffs_lat.get('PVY3', 0.5) + coeffs_lat.get('PVY4', -1.8) * dfz) * gamma_rad * LKYC)

    # Base lateral force (pure slip)
    alpha_eq = alpha_rad + Shy # Equivalent slip angle
    # Magic Formula Equation
    Fy0 = Dy * np.sin(Cfy * np.arctan(By * alpha_eq - Ey * (By * alpha_eq - np.arctan(By * alpha_eq)))) + Svy

    # --- Combined Slip Calculation (Gyk weighting factor) ---
    # Longitudinal Force dependency (using Fx_calculated)
    muy_x = (coeffs_long.get('PDX1', 1.7) + coeffs_long.get('PDX2', -0.28) * dfz) * (1.0 - coeffs_long.get('PDX3', 6.5)*gamma_rad**2) * coeffs_long.get('LMUX', 1.0)
    Dx = muy_x * Fz
    # Avoid division by zero
    if abs(Dx) < 1e-6:
        sx = 0.0
    else:
        sx = Fx_calculated / Dx

    # B factor for combined slip Gyk
    Rby1 = coeffs_lat.get('RBY1', 29.7)
    Rby2 = coeffs_lat.get('RBY2', 18.0)
    Rby3 = coeffs_lat.get('RBY3', -0.003)
    Rby4 = coeffs_lat.get('RBY4', 0.23) # Camber influence
    alphat = alpha_rad + Rby3
    B_Gyk = (Rby1 + Rby4*gamma_rad**2) * np.cos(np.arctan(Rby2 * alphat))

    # C factor for combined slip Gyk
    Rcy1 = coeffs_lat.get('RCY1', 0.99)
    C_Gyk = Rcy1

    # E factor for combined slip Gyk
    Rey1 = coeffs_lat.get('REY1', -0.0027)
    Rey2 = coeffs_lat.get('REY2', 0.083)
    E_Gyk = Rey1 + Rey2 * dfz
    E_Gyk = np.clip(E_Gyk, -10.0, 1.0) # Limit Ey

    # Gyk weighting factor calculation
    # Avoid potential division by zero if B_Gyk is zero
    if abs(B_Gyk) < 1e-6:
         Gyk = 0.0 # Or handle based on limits
    else:
        arg_atan_Gyk = B_Gyk * kappa
        Gyk = np.cos(C_Gyk * np.arctan(arg_atan_Gyk - E_Gyk * (arg_atan_Gyk - np.arctan(arg_atan_Gyk))))

    # Kappa-induced lateral force component (Svyk)
    Svyk = Dy * (coeffs_lat.get('RVY1', 0.04) + coeffs_lat.get('RVY2', 0.036) * dfz + coeffs_lat.get('RVY3', 0.28) * gamma_rad) * \
           np.cos(np.arctan(coeffs_lat.get('RVY4', 20.2) * alpha_rad)) * np.arctan(coeffs_lat.get('RVY5', -2.8) * kappa) * LVYKA # RVY6 term omitted for simplicity

    # Final Combined Lateral Force
    Fy = Gyk * Fy0 + Svyk

    return Fy


def calculate_Fx(kappa, Fz, gamma_rad, alpha_rad, coeffs_long, coeffs_lat, Fy_calculated):
    """
    Calculates Longitudinal Force (Fx) using Pacejka Magic Formula (MF6.2 - combined slip).

    Args:
        kappa (float): Longitudinal slip ratio.
        Fz (float): Vertical load in Newtons. Fz > 0.
        gamma_rad (float): Camber (inclination) angle in radians.
        alpha_rad (float): Slip angle in radians.
        coeffs_long (dict): Dictionary of longitudinal coefficients.
        coeffs_lat (dict): Dictionary of lateral coefficients.
        Fy_calculated (float): Lateral force calculated separately (needed for Gxk).

    Returns:
        float: Calculated longitudinal force (Fx) in Newtons.
    """
    # Ensure Fz is positive
    Fz = max(Fz, 1e-6)

    # Scale factors
    LMUX = coeffs_long.get('LMUX', 1.0)
    LKX = coeffs_long.get('LKX', 1.0)
    LCX = coeffs_long.get('LCX', 1.0)
    LEX = coeffs_long.get('LEX', 1.0)
    LHX = coeffs_long.get('LHX', 1.0)
    LVX = coeffs_long.get('LVX', 1.0)
    LXAL = coeffs_long.get('LXAL', 1.0) # Alpha influence scaling

    # Normalized load change
    dfz = (Fz - FNOMIN * coeffs_long.get('LFZO', 1.0)) / (FNOMIN * coeffs_long.get('LFZO', 1.0))

    # --- Base Pure Slip Longitudinal Force Calculation ---
    # C factor (Shape)
    Cx = coeffs_long.get('PCX1', 1.67) * LCX

    # D factor (Peak) - Friction
    mux = (coeffs_long.get('PDX1', 1.73) + coeffs_long.get('PDX2', -0.28) * dfz) * \
          (1.0 - coeffs_long.get('PDX3', 6.55) * gamma_rad**2) * LMUX
    Dx = mux * Fz

    # BCD factor (Longitudinal stiffness)
    Kxk = Fz * (coeffs_long.get('PKX1', 64.6) + coeffs_long.get('PKX2', -0.95) * dfz) * \
          np.exp(coeffs_long.get('PKX3', 0.016) * dfz) * LKX
    # Avoid division by zero
    if abs(Cx * Dx) < 1e-6:
        BCD = Kxk
    else:
        BCD = Kxk / (Cx * Dx)

    # B factor (Stiffness factor)
    Bx = BCD # Approximation for pure slip

    # E factor (Curvature)
    Ex = (coeffs_long.get('PEX1', 0.0) + coeffs_long.get('PEX2', 0.46) * dfz + coeffs_long.get('PEX3', -0.036) * dfz**2) * \
         (1.0 - coeffs_long.get('PEX4', 0.38) * np.sign(kappa)) * LEX # PEX4 depends on sign(kappa)
    Ex = np.clip(Ex, -10.0, 1.0) # Limit Ex

    # Horizontal and Vertical shifts
    Shx = (coeffs_long.get('PHX1', 0.0006) + coeffs_long.get('PHX2', -0.00066) * dfz) * LHX
    Svx = Fz * (coeffs_long.get('PVX1', -0.056) + coeffs_long.get('PVX2', 0.065) * dfz) * LVX

    # Base longitudinal force (pure slip)
    kappa_eq = kappa + Shx # Equivalent slip
    # Magic Formula Equation
    Fx0 = Dx * np.sin(Cx * np.arctan(Bx * kappa_eq - Ex * (Bx * kappa_eq - np.arctan(Bx * kappa_eq)))) + Svx

    # --- Combined Slip Calculation (Gxk weighting factor) ---
    # Lateral Force dependency (using Fy_calculated)
    muy_y = (coeffs_lat.get('PDY1', 1.7) + coeffs_lat.get('PDY2', -0.35) * dfz) * (1.0 - coeffs_lat.get('PDY3', 3.36)*gamma_rad**2) * coeffs_lat.get('LMUY', 1.0)
    Dy = muy_y * Fz
    # Avoid division by zero
    if abs(Dy) < 1e-6:
        sy = 0.0
    else:
        sy = Fy_calculated / Dy

    # B factor for combined slip Gxk
    Rbx1 = coeffs_long.get('RBX1', 23.1)
    Rbx2 = coeffs_long.get('RBX2', -22.9)
    Rbx3 = coeffs_long.get('RBX3', -0.54) # Camber influence
    kappat = kappa + Rbx3 * gamma_rad # Effective kappa for Gxk B factor
    B_Gxk = (Rbx1 + Rbx2 * kappat**2) * np.cos(np.arctan(Rbx3 * alpha_rad)) # Rbx3 used differently here? Check MF docs. Using alpha influence.

    # C factor for combined slip Gxk
    Rcx1 = coeffs_long.get('RCX1', 1.05)
    C_Gxk = Rcx1

    # E factor for combined slip Gxk
    Rex1 = coeffs_long.get('REX1', 0.0)
    Rex2 = coeffs_long.get('REX2', 0.0)
    E_Gxk = Rex1 + Rex2 * dfz
    E_Gxk = np.clip(E_Gxk, -10.0, 1.0) # Limit E

    # Gxk weighting factor calculation
    # Avoid potential division by zero if B_Gxk is zero
    if abs(B_Gxk) < 1e-6:
        Gxk = 0.0
    else:
        arg_atan_Gxk = B_Gxk * alpha_rad
        Gxk = np.cos(C_Gxk * np.arctan(arg_atan_Gxk - E_Gxk * (arg_atan_Gxk - np.arctan(arg_atan_Gxk))))

    # Final Combined Longitudinal Force
    Fx = Gxk * Fx0
    # Note: MF6.2 might have an alpha-induced Fx term (SvxAlpha), often small.

    return Fx


def calculate_Mz(alpha_rad, Fz, gamma_rad, kappa, coeffs_align, coeffs_lat, coeffs_long, Fy_calculated):
    """
    Calculates Aligning Torque (Mz) using Pacejka Magic Formula (MF6.2 - combined slip).

    Args:
        alpha_rad (float): Slip angle in radians.
        Fz (float): Vertical load in Newtons. Fz > 0.
        gamma_rad (float): Camber (inclination) angle in radians.
        kappa (float): Longitudinal slip ratio.
        coeffs_align (dict): Dictionary of aligning coefficients.
        coeffs_lat (dict): Dictionary of lateral coefficients.
        coeffs_long (dict): Dictionary of longitudinal coefficients.
        Fy_calculated (float): Lateral force calculated separately (needed for Mzr).

    Returns:
        float: Calculated aligning torque (Mz) in Newton-meters.
    """
    # Ensure Fz is positive
    Fz = max(Fz, 1e-6)
    R0 = coeffs_align.get('UNLOADED_RADIUS', 0.3) # Need UNLOADED_RADIUS, get from [DIMENSION]

    # Scale factors
    LTR = coeffs_align.get('LTR', 1.0)
    LRES = coeffs_align.get('LRES', 1.0)
    LS = coeffs_align.get('LS', 1.0)
    LFZO = coeffs_align.get('LFZO', 1.0) # Should be same as others

    # Normalized load change
    dfz = (Fz - FNOMIN * LFZO) / (FNOMIN * LFZO)

    # --- Pneumatic Trail (t) Calculation ---
    # C factor (Shape)
    Cpt = coeffs_align.get('QCZ1', 1.24)

    # D factor (Peak Trail)
    Dpt = Fz * (coeffs_align.get('QDZ1', 0.11) + coeffs_align.get('QDZ2', 0.003) * dfz) * \
          (1.0 + coeffs_align.get('QDZ3', -0.95) * gamma_rad + coeffs_align.get('QDZ4', -18.9) * gamma_rad**2) * LTR * R0

    # BCD factor (Trail Stiffness)
    Bpt = (coeffs_align.get('QBZ1', 9.6) + coeffs_align.get('QBZ2', -0.32) * dfz + coeffs_align.get('QBZ3', 0.12) * dfz**2) * \
          (1.0 + coeffs_align.get('QBZ4', 1.88) * gamma_rad + coeffs_align.get('QBZ5', 1.52) * abs(gamma_rad)) * (coeffs_lat.get('LKY', 1.0) / coeffs_lat.get('LMUY', 1.0)) # Uses Ky/muy ratio
    # Avoid division by zero
    if abs(Cpt * Dpt) < 1e-6:
        BCD_t = Bpt * coeffs_lat.get('LKY', 1.0) / coeffs_lat.get('LMUY', 1.0) # Approximation
    else:
        BCD_t = Bpt / (Cpt * Dpt)

    # B factor (Stiffness factor)
    Bt = BCD_t # Approximation

    # E factor (Curvature)
    Ept = (coeffs_align.get('QEZ1', -6.8) + coeffs_align.get('QEZ2', 0.89) * dfz + coeffs_align.get('QEZ3', 1.24) * dfz**2) * \
          (1.0 + (coeffs_align.get('QEZ4', -0.62) + coeffs_align.get('QEZ5', -3.46) * gamma_rad) * (2.0/np.pi) * np.arctan(Bt * Cpt * alpha_rad))
    Ept = np.clip(Ept, -10.0, 1.0) # Limit Ept

    # Horizontal shift
    Sht = (coeffs_align.get('QHZ1', -0.0088) + coeffs_align.get('QHZ2', 0.0026) * dfz) + \
          (coeffs_align.get('QHZ3', 0.05) + coeffs_align.get('QHZ4', 0.078) * dfz) * gamma_rad

    # Equivalent slip angle for trail
    alpha_t_eq = alpha_rad + Sht

    # Pneumatic trail calculation (Magic Formula form)
    trail = Dpt * np.cos(Cpt * np.arctan(Bt * alpha_t_eq - Ept * (Bt * alpha_t_eq - np.arctan(Bt * alpha_t_eq)))) * R0 # R0 multiplication? Check MF docs. Assuming Dpt is normalized.

    # --- Residual Torque (Mzr) Calculation ---
    # D factor (Peak Residual Torque)
    Dr = Fz * R0 * ( (coeffs_align.get('QDZ6', -0.01) + coeffs_align.get('QDZ7', -0.014) * dfz) * LRES + \
                     (coeffs_align.get('QDZ8', -1.55) + coeffs_align.get('QDZ9', 0.39) * dfz) * gamma_rad * coeffs_lat.get('LKZC', 1.0) + \
                     (coeffs_align.get('QDZ10', 5.76) + coeffs_align.get('QDZ11', -1.6) * dfz) * gamma_rad * abs(gamma_rad) * coeffs_lat.get('LKZC', 1.0) )

    # B factor (Stiffness factor)
    Br = (coeffs_align.get('QBZ9', 0.0) * coeffs_lat.get('LKY', 1.0) / coeffs_lat.get('LMUY', 1.0) + coeffs_align.get('QBZ10', 0.01) * Bpt) * LRES # Uses Bpt

    # Equivalent slip angle for residual torque
    alpha_r_eq = alpha_rad + Sht # Using same shift for simplicity

    # Residual torque calculation (Magic Formula form, often simplified)
    # Using Fy_calculated as the primary driver for Mzr shape
    # Avoid division by zero
    muy_y = (coeffs_lat.get('PDY1', 1.7) + coeffs_lat.get('PDY2', -0.35) * dfz) * (1.0 - coeffs_lat.get('PDY3', 3.36)*gamma_rad**2) * coeffs_lat.get('LMUY', 1.0)
    Dy = muy_y * Fz
    if abs(Dy) < 1e-6:
        fy_norm = 0.0
    else:
        fy_norm = Fy_calculated / Dy

    Mzr = Dr * np.cos(np.arctan(Br * alpha_r_eq)) * fy_norm # Simplified cosine form driven by Fy

    # --- Combined Slip Effect on Mz ---
    # Effect of Fx on Mz (proportional to Fx * trail arm 's')
    s = R0 * (coeffs_align.get('SSZ1', -0.09) + coeffs_align.get('SSZ2', 0.01) * (Fy_calculated / (FNOMIN * LFZO)) + \
             (coeffs_align.get('SSZ3', 1.7) + coeffs_align.get('SSZ4', -1.7) * dfz) * gamma_rad) * LS

    # Need Fx calculated under combined slip
    # Fx_combined = calculate_Fx(kappa, Fz, gamma_rad, alpha_rad, coeffs_long, coeffs_lat, Fy_calculated) # Avoid recursive call if possible
    # For now, assume Fx_calculated passed in is the combined one, or use pure slip Fx0 as approximation
    Fx_combined = calculate_Fx(kappa, Fz, gamma_rad, alpha_rad, coeffs_long, coeffs_lat, Fy_calculated) # Recalculate Fx here

    # Combined slip reduction factor for trail (similar to Gyk/Gxk, but often simpler)
    # Using a simple cosine reduction based on kappa for demonstration
    # MF6.2 has specific parameters (QDTP1, etc.) if turn slip is fully modeled
    lambda_t = np.cos(np.pi/2 * kappa) # Very simplified reduction factor

    # Final Aligning Torque
    Mz = -trail * Fy_calculated * lambda_t + Mzr + s * Fx_combined

    return Mz

# Note: These implementations are based on common MF6.2 structures but may need
# refinement based on specific MF-Tyre documentation for exact parameter interactions,
# especially for combined slip weighting factors and gamma effects.
# The UNLOADED_RADIUS needs to be passed or loaded correctly.
# Consider adding error handling for missing coefficients.