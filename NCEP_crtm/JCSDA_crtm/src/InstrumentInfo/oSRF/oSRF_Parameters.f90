MODULE oSRF_Parameters

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds, ONLY: fp
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &


  ! Planck coefficient parameters
  ! ...Dimension
  INTEGER, PUBLIC, PARAMETER :: N_PLANCK_COEFFS = 2
  ! ...First Planck function constant (C1) scale factors. Units of C1 are W.m^2.
  !    Length scaling: To convert to W/(m^2.cm^-4) requires a scaling of m->cm,
  !                    which is 100, to the fourth power, which is 1.0e+08.
  !    Power scaling:  To convert to mW.m^2 requires a scaling of 1000.
  REAL(fp), PUBLIC, PARAMETER :: C1_LENGTH_SCALE_FACTOR = 1.0e+08_fp
  REAL(fp), PUBLIC, PARAMETER :: C1_POWER_SCALE_FACTOR  = 1.0e+03_fp
  REAL(fp), PUBLIC, PARAMETER :: C1_SCALE_FACTOR = C1_LENGTH_SCALE_FACTOR * C1_POWER_SCALE_FACTOR
  ! ...Second Planck function constant (C2) scale factor. Units of C2 are K.m,
  !    So to convert to K.cm, a scaling of 100 is applied.
  REAL(fp), PUBLIC, PARAMETER :: C2_SCALE_FACTOR = 100.0_fp


  ! Polychromatic coefficient parameters
  ! ...Dimension. Currently, default is a linear fit, y = a + b*x
  INTEGER, PUBLIC, PARAMETER :: N_POLY_COEFFS = 2
  ! ...Temperature settings
  REAL(fp), PUBLIC, PARAMETER :: D_TEMP   = 5.0_fp    ! ** Changing any of these values
  REAL(fp), PUBLIC, PARAMETER :: MIN_TEMP = 150.0_fp  ! ** will generate different numbers
  REAL(fp), PUBLIC, PARAMETER :: MAX_TEMP = 340.0_fp  ! ** for the polychromatic coeffs.
  INTEGER , PUBLIC, PARAMETER :: N_TEMPS  = INT((MAX_TEMP-MIN_TEMP)/D_TEMP + 1.5_fp)

  
  ! Tolerance values for comparisons
  REAL(fp), PUBLIC, PARAMETER :: INTEGRATE_TOLERANCE    = 1.0e-11_fp
  REAL(fp), PUBLIC, PARAMETER :: CONVOLVE_TOLERANCE     = 1.0e-11_fp
  REAL(fp), PUBLIC, PARAMETER :: EFFECTIVE_T_TOLERANCE  = 1.0e-11_fp
  REAL(fp), PUBLIC, PARAMETER :: POLY_COEFF_TOLERANCE   = 1.0e-11_fp
  REAL(fp), PUBLIC, PARAMETER :: F0_TOLERANCE           = 1.0e-11_fp
  REAL(fp), PUBLIC, PARAMETER :: PLANCK_COEFF_TOLERANCE = 1.0e-18_fp
  
  ! Flag specific parameters. Include file auto-generated by gen_flag_procedures.rb
  INCLUDE 'oSRF_Flag_Parameters.inc'

END MODULE oSRF_Parameters
