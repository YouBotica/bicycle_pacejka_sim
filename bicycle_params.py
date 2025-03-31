"""
Defines the parameters for the bicycle model, including loading tire coefficients.
"""
import sys
import os

# Add the directory containing tire_parser to the Python path
# This assumes bicycle_params.py and tire_parser.py are in the same directory
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

try:
    from tire_parser import parse_tir_coeffs
except ImportError as e:
    print(f"Error importing tire_parser: {e}")
    print("Please ensure tire_parser.py is in the same directory or accessible in PYTHONPATH.")
    sys.exit(1)

# --- Bicycle Parameters ---
# Known Parameters (from user input)
m = 800.0      # Total mass (kg)
Izz = 1000.0   # Yaw moment of inertia (kg*m^2)
a = 1.6        # Distance from CoG to front axle (m)
b = 1.2        # Distance from CoG to rear axle (m)
L = a + b      # Wheelbase (m)

# Placeholder/Typical Parameters (to be refined if necessary)
h = 0.5        # Height of CoG above ground (m) - Placeholder
w = 1.5        # Track width (m) - Placeholder (assuming car-like, adjust if needed)

# Gravitational acceleration
g = 9.81       # m/s^2

# --- Tire Data Loading ---

# Content from lf_tire.tir (obtained via read_file)
# It's generally better to read this from the file, but embedding for now.
tir_file_content = """
[MDI_HEADER]
FILE_TYPE        = 'tir'
FILE_VERSION     =  3.0
FILE_FORMAT      = 'ASCII'
! : TIRE_VERSION :  MF62
! : USE_MODE :      4
! : COMPATIBILITY : TNO MF-Tyre MF-Swift 6.2
! : COMMENT :       Stackpole Engineering Services, Inc.
! : COMMENT :       Created By : Michael Arbogast
! : COMMENT :       Date Created : July 30, 2021
! : COMMENT :       Customer : Bridgestone Americas
! : COMMENT :       Manufacturer : Firestone
! : COMMENT :       Construction : IAC LF, B4112T
! : COMMENT :       DOT : XXXX
! : COMMENT :       Tyre Size : 10.0/23.5R15
! : COMMENT :       Rim Width (in) : 10
! : COMMENT :       Infl. Pressure (psi) : 29.09
! : COMMENT :       Test Speed (mph) : 149.8
! : COMMENT :       Position : Left
!
! Copyright Stackpole Engineering Services, Inc. 2021
!
! USE_MODE specifies the type of calculation performed:
!       0: Fz only, no Magic Formula evaluation
!       1: Fx,My only
!       2: Fy,Mx,Mz only
!       3: Fx,Fy,Mx,My,Mz uncombined force/moment calculation
!       4: Fx,Fy,Mx,My,Mz combined force/moment calculation
!       5: Fx,Fy,Mx,My,Mz combined force/moment calculation + turnslip
!      +0: steady state behaviour
!     +10: including relaxation behaviour
!     +20: including relaxation behaviour (nonlinear)
!     +30: including rigid ring dynamics
!    +100: smooth road contact
!    +200: smooth road contact (circular cross section, motorcycles)
!    +400: road contact for 2D roads (using travelled distance)
!    +500: road contact for 3D roads
!
$-----------------------------------------------------------------UNITS
[UNITS]
LENGTH                   = 'meter'
FORCE                    = 'newton'
ANGLE                    = 'radians'
MASS                     = 'kg'
TIME                     = 'second'
$-----------------------------------------------------------------MODEL
[MODEL]
FITTYP                   =  62                      $Magic Formula version number
PROPERTY_FILE_FORMAT     = 'MF62'                   $Tyre model selection (ADAMS only)
USE_MODE                 =  4                       $Tyre use mode switch (ADAMS only)
TYRESIDE                 = 'Left'                   $Position of tyre during measurements
LONGVL                   =  66.9666                 $Reference speed
VXLOW                    =  1                       $Lower boundary velocity in slip calculation
$-----------------------------------------------------------------DIMENSIONS
[DIMENSION]
UNLOADED_RADIUS          =  0.3004                  $Free tyre radius
WIDTH                    =  0.283464                $Nominal section width of the tyre
RIM_RADIUS               =  0.1905                  $Nominal rim radius
RIM_WIDTH                =  0.254                   $Rim width
ASPECT_RATIO             =  0.388                   $Nominal aspect ratio
$-----------------------------------------------------------------OPERATING_CONDITIONS
[OPERATING_CONDITIONS]
INFLPRES                 =  200573                  $Tyre inflation pressure
NOMPRES                  =  200573                  $Nominal inflation pressure used in MF equations
$-----------------------------------------------------------------SHAPE
[SHAPE]
{radial width}
 1.0    0.0
 1.0    0.4
 1.0    0.9
 0.9    1.0
$-----------------------------------------------------------------INERTIA
[INERTIA]
MASS                     =  6.439                   $Tyre Mass
IXX                      =  0                       $Tyre diametral moment of inertia
IYY                      =  0                       $Tyre polar moment of inertia
$-----------------------------------------------------------------VERTICAL
[VERTICAL]
FNOMIN                   =  3112                    $Nominal wheel load
VERTICAL_STIFFNESS       =  293478                  $Tyre vertical stiffness
VERTICAL_DAMPING         =  50                      $Tyre vertical damping
MC_CONTOUR_A             =  0                       $Motorcycle contour ellipse A
MC_CONTOUR_B             =  0                       $Motorcycle contour ellipse B
BREFF                    =  1.67769                 $Low load stiffness of effective rolling radius
DREFF                    =  0.105132                $Peak value of effective rolling radius
FREFF                    =  0.0763774               $High load stiffness of effective rolling radius
Q_RE0                    =  0.994334                $Ratio of free tyre radius with nominal tyre radius
Q_V1                     =  0.00542853              $Tyre radius increase with speed
Q_V2                     =  0.0867414               $Vertical stiffness increase with speed
Q_FZ1                    =  28.1283                 $Linear term in load vs. deflection
Q_FZ2                    =  2.78052                 $Quadratic term in load vs. deflection
Q_FCX                    =  0.147144                $Longitudinal force influence on vertical stiffness
Q_FCY                    =  0.171698                $Lateral force influence on vertical stiffness
Q_FCY2                   =  0                       $Explicit load dependency for including the lateral force influence on vertical stiffness
Q_CAM                    =  0                       $Stiffness reduction due to camber
Q_CAM1                   =  0                       $Linear load dependent camber angle influence on vertical stiffness
Q_CAM2                   =  0                       $Quadratic load dependent camber angle influence on vertical stiffness
Q_CAM3                   =  0                       $Linear load and camber angle dependent reduction on vertical stiffness
Q_FYS1                   =  0                       $Combined camber angle and side slip angle effect on vertical stiffness (constant)
Q_FYS2                   =  0                       $Combined camber angle and side slip angle linear effect on vertical stiffness
Q_FYS3                   =  0                       $Combined camber angle and side slip angle quadratic effect on vertical stiffness
PFZ1                     =  0.896152                $Pressure effect on vertical stiffness
BOTTOM_OFFST             =  0.01                    $Distance to rim when bottoming starts to occur
BOTTOM_STIFF             =  2.93478e+06             $Vertical stiffness of bottomed tyre
$-----------------------------------------------------------------STRUCTURAL
[STRUCTURAL]
LONGITUDINAL_STIFFNESS   =  0                       $Tyre overall longitudinal stiffness
LATERAL_STIFFNESS        =  0                       $Tyre overall lateral stiffness
YAW_STIFFNESS            =  0                       $Tyre overall yaw stiffness
DAMP_RESIDUAL            =  0.002                   $Residual damping (proportional to stiffness)
DAMP_VLOW                =  0.001                   $Additional low speed damping (proportional to stiffness)
PCFX1                    =  0                       $Tyre overall longitudinal stiffness vertical deflection dependency linear term
PCFX2                    =  0                       $Tyre overall longitudinal stiffness vertical deflection dependency quadratic
PCFX3                    =  0                       $Tyre overall longitudinal stiffness pressure dependency
PCFY1                    =  0                       $Tyre overall lateral stiffness vertical deflection dependency linear term
PCFY2                    =  0                       $Tyre overall lateral stiffness vertical deflection dependency quadratic
PCFY3                    =  0                       $Tyre overall lateral stiffness pressure dependency
PCMZ1                    =  0                       $Tyre overall yaw stiffness pressure dependency
$-----------------------------------------------------------------INFLATION_PRESSURE_RANGE
[INFLATION_PRESSURE_RANGE]
PRESMIN                  =  172715                  $Minimum allowed inflation pressure
PRESMAX                  =  214409                  $Maximum allowed inflation pressure
$-----------------------------------------------------------------VERTICAL_FORCE_RANGE
[VERTICAL_FORCE_RANGE]
FZMIN                    =  100                     $Minimum allowed wheel load
FZMAX                    =  8000                    $Maximum allowed wheel load
$-----------------------------------------------------------------LONGITUDINAL_SLIP_RANGE
[LONG_SLIP_RANGE]
KPUMIN                   = -1.00                    $Minimum valid wheel slip
KPUMAX                   =  1.00                    $Maximum valid wheel slip
$-----------------------------------------------------------------SLIP_ANGLE_RANGE
[SLIP_ANGLE_RANGE]
ALPMIN                   = -0.523599                $Minimum valid slip angle
ALPMAX                   =  0.523599                $Maximum valid slip angle
$-----------------------------------------------------------------INCLINATION_ANGLE_RANGE
[INCLINATION_ANGLE_RANGE]
CAMMIN                   = -0.087266                $Minimum valid camber angle
CAMMAX                   =  0.087266                $Maximum valid camber angle
$-----------------------------------------------------------------SCALING
[SCALING_COEFFICIENTS]
LFZO                     =  1                       $Scale factor of nominal (rated) load
LCX                      =  1                       $Scale factor of Fx shape factor
LMUX                     =  1                       $Scale factor of Fx peak friction coefficient
LEX                      =  1                       $Scale factor of Fx curvature factor
LKX                      =  1                       $Scale factor of Fx slip stiffness
LHX                      =  1                       $Scale factor of Fx horizontal shift
LVX                      =  1                       $Scale factor of Fx vertical shift
LCY                      =  1                       $Scale factor of Fy shape factor
LMUY                     =  1                       $Scale factor of Fy peak friction coefficient
LEY                      =  1                       $Scale factor of Fy curvature factor
LKY                      =  1                       $Scale factor of cornering stiffness
LKYC                     =  1                       $Scale factor of camber stiffness
LKZC                     =  1                       $Scale factor of camber moment stiffness
LHY                      =  1                       $Scale factor of Fy horizontal shift
LVY                      =  1                       $Scale factor of Fy vertical shift
LTR                      =  1                       $Scale factor of Peak of pneumatic trail
LRES                     =  1                       $Scale factor for offset of residual torque
LXAL                     =  1                       $Scale factor of alpha influence on Fx
LYKA                     =  1                       $Scale factor of alpha influence on Fy
LVYKA                    =  1                       $Scale factor of kappa induced Fy
LS                       =  1                       $Scale factor of Moment arm of Fx
LMX                      =  1                       $Scale factor of overturning couple
LVMX                     =  1                       $Scale factor of Mx vertical shift
LMY                      =  1                       $Scale factor of rolling resistance torque
LMP                      =  1                       $Scale factor of parking moment
$-----------------------------------------------------------------LONGITUDINAL_FORCE
[LONGITUDINAL_COEFFICIENTS]
PCX1                     =  1.66628                 $Shape factor Cfx for longitudinal force
PDX1                     =  1.72645                 $Longitudinal friction Mux at Fznom
PDX2                     = -0.278479                $Variation of friction Mux with load
PDX3                     =  6.54553                 $Variation of friction Mux with camber
PEX1                     = -5.04761e-06             $Longitudinal curvature Efx at Fznom
PEX2                     =  0.463966                $Variation of curvature Efx with load
PEX3                     = -0.0359398               $Variation of curvature Efx with load squared
PEX4                     =  0.38592                 $Factor in curvature Efx while driving
PKX1                     =  64.6315                 $Longitudinal slip stiffness Kfx/Fz at Fznom
PKX2                     = -0.95099                 $Variation of slip stiffness Kfx/Fz with load
PKX3                     =  0.015721                $Exponent in slip stiffness Kfx/Fz with load
PHX1                     =  0.00060979              $Horizontal shift Shx at Fznom
PHX2                     = -0.000659181             $Variation of shift Shx with load
PVX1                     = -0.0560991               $Vertical shift Svx/Fz at Fznom
PVX2                     =  0.0646569               $Variation of shift Svx/Fz with load
RBX1                     =  23.0974                 $Slope factor for combined slip Fx reduction
RBX2                     = -22.9129                 $Variation of slope Fx reduction with kappa
RBX3                     = -0.5445                  $Influence of camber on stiffness for Fx combined
RCX1                     =  1.05471                 $Shape factor for combined slip Fx reduction
REX1                     =  3.97084e-08             $Curvature factor of combined Fx
REX2                     =  5.73425e-09             $Curvature factor of combined Fx with load
RHX1                     =  0.00655663              $Shift factor for combined slip Fx reduction
PPX1                     =  0                       $Linear pressure effect on slip stiffness
PPX2                     =  0                       $Quadratic pressure effect on slip stiffness
PPX3                     =  0                       $Linear pressure effect on longitudinal friction
PPX4                     =  0                       $Quadratic pressure effect on longitudinal friction
$-----------------------------------------------------------------OVERTURNING_MOMENT
[OVERTURNING_COEFFICIENTS]
QSX1                     = -0.00469662              $Overturning moment offset
QSX2                     = -5.13093                 $Camber induced overturning couple
QSX3                     =  0.117227                $Fy induced overturning couple
QSX4                     =  0.764013                $Mixed load, lateral force, and camber on Mx
QSX5                     = -0.569854                $Load effect on Mx with lateral force and camber
QSX6                     = -3.0208                  $B-factor of load with Mx
QSX7                     = -4.74982                 $Camber with load on Mx
QSX8                     = -5.9049                  $Lateral force with load on Mx
QSX9                     =  0.0425242               $B-factor of lateral force with load on Mx
QSX10                    = -2.94964                 $Vertical force with camber on Mx
QSX11                    =  13.9701                 $B-factor of vertical force with camber on Mx
QSX12                    =  0                       $Camber squared induced overturning moment
QSX13                    =  0                       $Lateral force induced overturning moment
QSX14                    =  0                       $Lateral force induced overturning moment with camber
PPMX1                    = -0.192994                $Influence of inflation pressure on overturning moment
$-----------------------------------------------------------------LATERAL_FORCE
[LATERAL_COEFFICIENTS]
PCY1                     =  1.5496                  $Shape factor Cfy for lateral forces
PDY1                     =  1.70839                 $Lateral friction Muy
PDY2                     = -0.354062                $Variation of friction Muy with load
PDY3                     =  3.35762                 $Variation of friction Muy with squared camber
PEY1                     = -1.98642                 $Lateral curvature Efy at Fznom
PEY2                     = -1.62092                 $Variation of curvature Efy with load
PEY3                     =  0.408939                $Zero order camber dependency of curvature Efy
PEY4                     = -4.08114                 $Variation of curvature Efy with camber
PEY5                     =  5.64917                 $Camber curvature Efc
PKY1                     = -53.0946                 $Maximum value of stiffness Kfy/Fznom
PKY2                     =  2.23333                 $Load at which Kfy reaches maximum value
PKY3                     =  0.872111                $Variation of Kfy/Fznom with camber
PKY4                     =  2.10492                 $Curvature of stiffness Kfy
PKY5                     =  0.95099                 $Peak stiffness variation with camber squared
PKY6                     = -2.40799                 $Camber stiffness factor
PKY7                     = -2.70746                 $Load dependency of camber stiffness factor
PHY1                     =  0.00234094              $Horizontal shift Shy at Fznom
PHY2                     =  0.00326645              $Variation of shift Shy with load
PVY1                     = -0.0625501               $Vertical shift in Svy/Fz at Fznom
PVY2                     =  0.0628882               $Variation of shift Svy/Fz with load
PVY3                     =  0.50157                 $Variation of shift Svy/Fz with camber
PVY4                     = -1.80795                 $Variation of shift Svy/Fz with camber and load
RBY1                     =  29.6888                 $Slope factor for combined Fy reduction
RBY2                     =  17.9759                 $Variation of slope Fy reduction with alpha
RBY3                     = -0.00317577              $Shift term for alpha in slope Fy reduction
RBY4                     =  0.229929                $Influence of camber on stiffness of Fy combined
RCY1                     =  0.993821                $Shape factor for combined Fy reduction
REY1                     = -0.00271301              $Curvature factor of combined Fy
REY2                     =  0.0831273               $Curvature factor of combined Fy with load
RHY1                     =  0.00988171              $Shift factor for combined Fy reduction
RHY2                     =  0.00695088              $Shift factor for combined Fy reduction with load
RVY1                     =  0.042693                $Kappa induced side force Svyk/Muy*Fz at Fznom
RVY2                     =  0.0364021               $Variation of Svyk/Muy*Fz with load
RVY3                     =  0.284909                $Variation of Svyk/Muy*Fz with camber
RVY4                     =  20.1816                 $Variation of Svyk/Muy*Fz with alpha
RVY5                     = -2.79814                 $Variation of Svyk/Muy*Fz with kappa
RVY6                     = -20.8171                 $Variation of Svyk/Muy*Fz with atan(kappa)
PPY1                     =  0.716732                $Pressure effect on cornering stiffness magnitude
PPY2                     =  1.64013                 $Pressure effect on location of cornering stiffness peak
PPY3                     = -0.0572745               $Linear pressure effect on lateral friction
PPY4                     = -0.543181                $Quadratic pressure effect on lateral friction
PPY5                     = -0.996956                $Influence of inflation pressure on camber stiffness
$-----------------------------------------------------------------ROLLING_MOMENT
[ROLLING_COEFFICIENTS]
QSY1                     =  0.0291029               $Rolling resistance torque coefficient
QSY2                     =  0                       $Rolling resistance torque depending on Fx
QSY3                     =  0.0106797               $Rolling resistance torque depending on speed
QSY4                     =  0.00150037              $Rolling resistance torque depending on speed ^4
QSY5                     =  0                       $Rolling resistance torque depending on camber squared
QSY6                     =  0                       $Rolling resistance torque depending on load and camber squared
QSY7                     =  0.293107                $Rolling resistance torque coefficient load dependency
QSY8                     = -1.27727                 $Rolling resistance torque coefficient pressure dependency
$-----------------------------------------------------------------ALIGNING_TORQUE
[ALIGNING_COEFFICIENTS]
QBZ1                     =  9.58435                 $Trail slope factor for trail Bpt at Fznom
QBZ2                     = -0.320338                $Variation of slope Bpt with load
QBZ3                     =  0.117051                $Variation of slope Bpt with load squared
QBZ4                     =  1.88053                 $Variation of slope Bpt with camber
QBZ5                     =  1.52129                 $Variation of slope Bpt with absolute camber
QBZ9                     =  0                       $Slope factor Br of residual torque Mzr
QBZ10                    =  0.0106191               $Slope factor Br of residual torque Mzr
QCZ1                     =  1.24211                 $Shape factor Cpt for pneumatic trail
QDZ1                     =  0.109243                $Peak trail Dpt" = Dpt*(Fz/Fznom*R0)
QDZ2                     =  0.00324279              $Variation of peak Dpt with load
QDZ3                     = -0.956087                $Variation of peak Dpt with camber
QDZ4                     = -18.9248                 $Variation of peak Dpt with camber squared
QDZ6                     = -0.00999806              $Peak residual torque Dmr = Dmr/(Fz*R0)
QDZ7                     = -0.0136136               $Variation of peak factor Dmr with load
QDZ8                     = -1.5519                  $Variation of peak factor Dmr with camber
QDZ9                     =  0.392793                $Variation of peak factor Dmr with camber and load
QDZ10                    =  5.75951                 $Variation of peak factor Dmr with camber squared
QDZ11                    = -1.59703                 $Variation of Dmr with camber squared and load
QEZ1                     = -6.82198                 $Trail curvature Ept at Fznom
QEZ2                     =  0.891157                $Variation of curvature Ept with load
QEZ3                     =  1.23965                 $Variation of curvature Ept with load squared
QEZ4                     = -0.621374                $Variation of curvature Ept with sign of Alpha-t
QEZ5                     = -3.45744                 $Variation of Ept with camber and sign Alpha-t
QHZ1                     = -0.0088485               $Trail horizontal shift Sht at Fznom
QHZ2                     =  0.00259685              $Variation of shift Sht with load
QHZ3                     =  0.0518341               $Variation of shift Sht with camber
QHZ4                     =  0.0782043               $Variation of shift Sht with camber and load
SSZ1                     = -0.0916292               $Nominal value of s/R0: effect of Fx on Mz
SSZ2                     =  0.0115969               $Variation of distance s/R0 with Fy/Fznom
SSZ3                     =  1.69749                 $Variation of distance s/R0 with camber
SSZ4                     = -1.7076                  $Variation of distance s/R0 with load and camber
PPZ1                     =  1.09916                 $Linear pressure effect on pneumatic trail
PPZ2                     =  0.41332                 $Influence of inflation pressure on residual aligning torque
$-----------------------------------------------------------------TURN-SLIP_PARAMETERS
[TURNSLIP_COEFFICIENTS]
PDXP1                    =  0                       $Peak Fx reduction due to spin parameter
PDXP2                    =  0                       $Peak Fx reduction due to spin with load parameter
PDXP3                    =  0                       $Peak Fx reduction due to spin with kappa parameter
PKYP1                    =  0                       $Cornering stiffness reduction due to spin
PDYP1                    =  0                       $Peak Fy reduction due to spin parameter
PDYP2                    =  0                       $Peak Fy reduction due to spin with varying load parameter
PDYP3                    =  0                       $Peak Fy reduction due to spin with alpha parameter
PDYP4                    =  0                       $Peak Fy reduction with square root of spin parameter
PHYP1                    =  0                       $Fy -alpha curve lateral shift limitation
PHYP2                    =  0                       $Fy -alpha curve maximum lateral shift parameter
PHYP3                    =  0                       $Fy -alpha curve maximum lateral shift varying with load parameter
PHYP4                    =  0                       $Fy -alpha curve maximum lateral shift parameter
PECP1                    =  0                       $Camber w.r.t. spin reduction factor parameter in camber stiffness
PECP2                    =  0                       $Camber w.r.t. spin reduction factor varying with load parameter in camber stiffness
QDTP1                    =  0                       $Pneumatic trail reduction factor due to turn slip parameter
QCRP1                    =  0                       $Turning moment at constant turning with zero forward speed parameter
QCRP2                    =  0                       $Turn slip moment (at alpha=90deg) parameter for increase with spin
QBRP1                    =  0                       $Residual (spin) torque reduction factor parameter due to side slip
QDRP1                    =  0                       $Turn slip moment peak magnitude parameter
"""

# Parse the coefficients
tire_coeffs_raw = parse_tir_coeffs(tir_file_content)

if tire_coeffs_raw is None:
    print("Error: Failed to parse tire coefficients. Exiting.")
    sys.exit(1)

# Extract key parameters from parsed data
# Add DIMENSION and VERTICAL sections to parser if needed, or parse manually here
def find_param(section, key, default=None):
    """Helper to find parameters in specific sections"""
    in_section = False
    for line in tir_file_content.splitlines():
        line = line.strip()
        if not line or line.startswith('$') or line.startswith('!'):
            continue
        if line.upper() == f"[{section.upper()}]":
            in_section = True
            continue
        if line.startswith('['): # Start of another section
            in_section = False
            continue
        if in_section:
            parts = line.split('=')
            if len(parts) == 2 and parts[0].strip().upper() == key.upper():
                try:
                    # Extract numeric value before potential comment
                    value_part = parts[1].split('$')[0].strip()
                    return float(value_part)
                except ValueError:
                    return default
    return default

FNOMIN = find_param('VERTICAL', 'FNOMIN', 3112.0) # Default from previous observation
UNLOADED_RADIUS = find_param('DIMENSION', 'UNLOADED_RADIUS', 0.30) # Default placeholder

# Update bicycle radius parameter
R = UNLOADED_RADIUS

# Structure coefficients for front/rear (assuming same tire for now)
# The parser returns a flat structure, let's nest it slightly
tire_coeffs = {
    'FNOMIN': FNOMIN,
    'UNLOADED_RADIUS': UNLOADED_RADIUS,
    'front': tire_coeffs_raw, # Contains LONGITUDINAL, LATERAL, ALIGNING keys
    'rear': tire_coeffs_raw   # Using the same coefficients for the rear tire
}

# Clean up temporary variables if desired
del tir_file_content
del tire_coeffs_raw
del find_param
del current_dir
del sys
del os

# Optional: Print loaded parameters for verification
# print(f"Loaded Bicycle Parameters:")
# print(f"  m    = {m} kg")
# print(f"  Izz  = {Izz} kg*m^2")
# print(f"  a    = {a} m")
# print(f"  b    = {b} m")
# print(f"  L    = {L} m")
# print(f"  h    = {h} m (placeholder)")
# print(f"  R    = {R} m (from .tir)")
# print(f"  w    = {w} m (placeholder)")
# print(f"  g    = {g} m/s^2")
# print(f"Loaded Tire Parameters:")
# print(f"  FNOMIN = {tire_coeffs['FNOMIN']} N")
# print(f"  UNLOADED_RADIUS = {tire_coeffs['UNLOADED_RADIUS']} m")
# print(f"  Front Lateral PCY1 = {tire_coeffs['front']['LATERAL'].get('PCY1', 'Not Found')}")
# print(f"  Rear Longitudinal PCX1 = {tire_coeffs['rear']['LONGITUDINAL'].get('PCX1', 'Not Found')}")