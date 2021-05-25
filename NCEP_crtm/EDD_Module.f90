MODULE EDD_Module

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE Type_Kinds,      ONLY: fp
  USE Message_Handler, ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters, ONLY: SET, ZERO, ONE, TWO, THREE, FOUR, ONEpointFIVE, &
                             PI, POINT_5, &
                             MAX_N_LAYERS, MAX_N_ANGLES, MAX_N_STOKES, &
                             DEGREES_TO_RADIANS, &
                             SCATTERING_ALBEDO_THRESHOLD, &
                             OPTICAL_DEPTH_THRESHOLD
  USE CRTM_Utility
  USE RTV_Define     , ONLY: RTV_type, &
                             DELTA_OPTICAL_DEPTH
  USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type

  ! Disable all implicit typing
  IMPLICIT NONE

  ! --------------------
  ! Default visibilities
  ! --------------------
  ! Everything private by default
  PRIVATE
  ! Procedures
  PUBLIC :: CRTM_EDD
  PUBLIC :: CRTM_EDD_Version

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*),  PARAMETER :: MODULE_VERSION_ID = &
  '$Id$'

  ! Version Id for the module
!  CHARACTER(*),  PARAMETER :: MODULE_VERSION_ID = &
!     '$Id$'

CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

    SUBROUTINE CRTM_EDD ( n_Layers, & ! Input, Number of atmospheric layers
                                 w, & ! Input, Layer single scattering albedo
                                 g, & ! Input, Layer asymmetry parameter
                              T_OD, & ! Input, Layer optical depth
                 cosmic_background, & ! Input, Cosmic background radiance
                        emissivity, & ! Input, Surface emissivity
                      reflectivity, & ! Input, Surface reflectivity
                   Index_Sat_Angle, & ! Input, Satellite angle index 
                               Atm, & ! Input, Atmospheric profiles
                               RTV)   ! IN/Output, Upward radiance and others

! ------------------------------------------------------------------------- !
!                                                                           !
! FUNCTION:                                                                 !
!                                                                           !
!   This subroutine calculates unpolarized IR/MW radiance at the top of     !
!   a scattering atmosphere using the Delta-Eddington approximation         !
!                                                                           ! 
!   Modified and adapted by Tom Greenwald    tomg@ssec.wisc.edu             !
!                                                                           !
!   EDDRT was originally written by K. F. Evans  Frank.Evans@colorado.edu   !
!                                                                           !
!   EDDRT computes the upwelling radiance for a plane-parallel              !
!   thermally emitting atmosphere using Eddington's second approximation.   !
!   The exitting radiance is calculated for a specified observation angle.  !
!   The only sources of radiation are internal thermal emission and         !
!   incident thermal radiation, so the radiation is azimuthally symmetric.  !
!   The atmosphere is divided up into a number of homogeneous layers each   !
!   having a Planck function linear with height.  The subroutine was        !
!   coded with the Rayleigh-Jeans approximation in mind so the Planck       !
!   function values are called temperatures. The height, Planck function,   !
!   extinction, single scattering albedo, and asymmetry factor are          !
!   specified for each layer.  The thermal radiation from the ground and    !
!   from above are input.  There are two types of ground emissivity: one for!
!   the diffuse calculation and one for the particular observation angle.   !
!                                                                           !
!   The model works by calculating the reflection, transmission, and        !      
!   thermal source terms for each layer from the input properties. A        !
!   tri-diagonal matrix solver is then used to compute the internal         !
!   diffuse radiation at each layer from the applied boundary conditions.   !
!   The diffuse Eddington radiation incident on each layer is then used     !
!   to find the scattering integral source term in the radiative transfer   !
!   equation.  With the known scattering and thermal source terms the       !
!   radiative transfer equation can be analytically integrated.  The        !
!   radiation exiting at the observation angle is summed over the layers,   !
!   first from the top down to the surface, then reflected, and summed      !
!   from the bottom to the top of the atmosphere.                           !
! ------------------------------------------------------------------------- !
    IMPLICIT NONE

    !------------
    !  Arguments
    !------------

    INTEGER, INTENT( IN ) :: n_Layers
    REAL (fp), INTENT( IN ), DIMENSION( : ) ::  w, g, T_OD
    REAL (fp), INTENT( IN ) ::  cosmic_background
    REAL (fp), INTENT( IN ), DIMENSION( : ) ::  emissivity
    REAL (fp), INTENT( IN ), DIMENSION( :, : ) ::  reflectivity
    INTEGER, INTENT( IN ) :: Index_Sat_Angle
    TYPE(CRTM_Atmosphere_type), INTENT( IN ) :: Atm
    TYPE(RTV_type), INTENT( INOUT ) :: RTV


   ! -------------- internal variables --------------------------------- !

    INTEGER, PARAMETER :: MAXN = 2 * MAX_N_LAYERS + 1
    REAL (fp), PARAMETER :: RADIANCE_RATIO_THRESH = 1.000_fp
    LOGICAL, PARAMETER :: SPECULAR = .TRUE.    
    INTEGER :: i, l, n
    REAL (fp) :: mu, asy, omega, ssa_tmp
    REAL (fp) :: surface_emissivity, surface_reflectivity 
    REAL (fp) :: UPRAD, DOWNRAD, RAD_RATIO, UPRAD_OLD
    REAL (fp) :: R, T, CP, CM, term1, term2, term3, term4
    REAL (fp) :: RADP1P, RADP1M, RADP2P, RADP2M, RADH1P, RADH2M, DELRAD
    REAL (fp), DIMENSION(MAX_N_LAYERS) :: opd, LAMBDA, X1, X2
    REAL (fp), DIMENSION(MAX_N_LAYERS) :: EXLP, EXLM, EX
    REAL (fp), DIMENSION(MAX_N_LAYERS) :: D, CPLUS, CMINUS
    REAL (fp), DIMENSION(MAX_N_LAYERS) :: D0, D1, DPLUS, DMINUS
    REAL (fp), DIMENSION(MAX_N_LAYERS) :: V, B0MU, B1, T1, T2
    REAL (fp), DIMENSION(MAX_N_LAYERS) :: REFLECT, TRANS
    REAL (fp), DIMENSION(MAXN) :: LOWER, UPPER, DIAG, RHS
    REAL (fp), DIMENSION(2,MAX_N_LAYERS) :: FLUX, SOURCE
!    REAL (fp), DIMENSION(0:n_Layers) :: total_opt, temp

    ! Cosine of satellite observation zenith angle
    mu = RTV%COS_Angle( Index_Sat_Angle )

    surface_emissivity = emissivity( Index_Sat_Angle )
    surface_reflectivity = reflectivity( Index_Sat_Angle, Index_Sat_Angle ) 

    ! Dump input to file
!    write(*,*) 'n_Layers: ', n_layers
!    write(56,*) n_layers
!    write(*,*) 'Index_Sat_angle: ', Index_sat_angle
!    write(*,*) 'cosmic_background: ', cosmic_background
!    write(56,*) cosmic_background
!    write(*,*) 'emissivity: ', surface_emissivity
!    write(56,*) surface_emissivity
!    write(*,*) 'reflectivity: ', surface_reflectivity
!    write(56,*) surface_reflectivity
!    write(*,*) 'Planck surface: ', RTV%Planck_surface
!    write(56,*) RTV%Planck_surface
!    write(*,*) 'ilev    opd    omega    asymf'
!    do i=1,n_layers
!      write(56,'(i3,5e13.5)') i, T_OD(i), w(i), g(i), Atm%Temperature(i), &
!                                Atm%Absorber(i,1)
!      write(*,'(i3,5e13.5)') i, T_OD(i), w(i), g(i), Atm%Temperature(i), &
!                                Atm%Absorber(i,1)
!    enddo
!    write(*,*) 'ilev   Planck'
!    do i=1,n_layers+1
!      write(56,*) i, Atm%Level_Pressure(i-1), RTV%Planck_atmosphere(i-1) 
!      write(*,*) i, Atm%Level_Pressure(i-1), RTV%Planck_atmosphere(i-1) 
!    enddo

!    temp(0) = Atm%Temperature(1)-0.3
!    do i=1,n_layers
!      temp(i) = Atm%Temperature(i)
!    enddo

    i = 1
    DO l = 1, N_Layers

      ! QC - Do I really need this?
!      IF ( w( l ) .GE. ONE ) THEN
!        ssa_tmp = 0.999
!      ELSE
!        ssa_tmp = w( l )
!      END IF

      ! Apply delta-scaling
      opd( l ) = ( ONE - w( l ) * g( l ) ** 2 ) * T_OD( l )
      omega = ( ONE - g( l ) ** 2 ) * w( l ) / ( ONE - w( l ) * g( l ) ** 2 )
      asy = g( l ) / ( ONE + g( l ) )
!      opd( l ) = ( ONE - ssa_tmp * g( l ) ** 2 ) * T_OD( l )
!      omega = ( ONE - g( l ) ** 2 ) * ssa_tmp / ( ONE - ssa_tmp * g( l ) ** 2 )
!      asy = g( l ) / ( ONE + g( l ) )
!      opd( l ) = T_OD( l )
!      omega = ssa_tmp
!      asy = g( l )

      ! Compute the reflection and transmission
      LAMBDA( l ) = SQRT( THREE * ( ONE - omega ) * ( ONE - omega * asy ) )
      R = ( ONE - omega * ( FOUR - THREE * asy ) ) / FOUR
      T = ( 7.0_fp - omega * ( FOUR + THREE * asy ) ) / FOUR
      X1( l ) = -R
      X2( l ) = LAMBDA( l ) + T
      EXLP( l ) = EXP( MIN( LAMBDA( l ) * opd( l ), 75.0_fp ) )
      EXLM( l ) = ONE / EXLP( l )
      EX( l ) = EXP( -opd( l ) / mu )
      D( l ) = X2( l )**2 * EXLP( l ) - X1( l )**2 * EXLM( l )
      D0( l ) = omega * ( T + LAMBDA( l ) - R )
      D1( l ) = ONEpointFIVE * omega * asy * ( T + LAMBDA( l ) + R ) * mu
      TRANS( l ) = TWO * LAMBDA(l) / ( X2( l ) * EXLP( l ) + ( LAMBDA( l ) - T ) * EXLM( l ) )
      REFLECT( l ) = X1( l ) * ( EXLP( l ) - EXLM( l ) ) * TRANS( l ) / ( TWO * LAMBDA( l ) )

      ! Calculate thermal source terms
      V( l ) = ( RTV%Planck_Atmosphere( l ) - RTV%Planck_Atmosphere( l - 1 ) ) / &
               ( THREE * ( ONE - omega * asy ) * opd( l ) )
!      V( l ) = ( temp( l ) - temp( l - 1 ) ) / &
!               ( THREE * ( ONE - omega * asy ) * opd( l ) )
      B1( l ) = ( RTV%Planck_Atmosphere( l ) - RTV%Planck_Atmosphere( l - 1 ) ) / opd( l )
!      B1( l ) = ( temp( l ) - temp( l - 1 ) ) / opd( l )
      B0MU( l ) = THREE * omega * asy * V( l ) * mu
      T1( l ) = POINT_5 * RTV%Planck_Atmosphere( l - 1 )
      T2( l ) = POINT_5 * RTV%Planck_Atmosphere( l )
!      T1( l ) = POINT_5 * temp( l - 1 )
!      T2( l ) = POINT_5 * temp( l )
      RADP1P = -V( l ) + T1( l )
      RADP2M =  V( l ) + T2( l )
      CP  =  ( X1( l ) * EXLM( l ) * RADP1P - X2( l ) * RADP2M ) / D( l )
      CM = ( -X2( l ) * EXLP( l ) * RADP1P + X1( l ) * RADP2M ) / D( l )
      RADP2P = -V( l ) + T2( l )
      RADP1M =  V( l ) + T1( l )
      SOURCE( 1, l ) = X1( l ) * CP * EXLP( l ) + X2( l ) * CM * EXLM( l ) + RADP2P
      SOURCE( 2, l ) = X2( l ) * CP + X1( l ) * CM + RADP1M
      DIAG( i ) = -REFLECT( l )
      DIAG( i + 1) = -REFLECT( l )
      LOWER( i ) = ONE
      LOWER( i + 1 ) = -TRANS( l )
      UPPER( i ) = -TRANS( l )
      UPPER( i + 1 ) = ONE
      RHS( i ) = SOURCE( 2, l )
      RHS( i + 1 ) = SOURCE( 1, l )
      i = i + 2
!write(*,*) l, lambda(l), d0(l),d1(l)
!write(*,*) l, EXLP( l ), EXLM( l ), EX( l )
    END DO
 
    ! Setup for and call the tri-diagonal matrix solver
    N = 2 * N_Layers + 1
    RHS( 1 ) = RHS( 1 ) + REFLECT( 1 ) * cosmic_background / TWO
    RHS( 2 ) = RHS( 2 ) + TRANS( 1 ) * cosmic_background / TWO
!    RHS( 1 ) = RHS( 1 ) + REFLECT( 1 ) * 2.7 / TWO
!    RHS( 2 ) = RHS( 2 ) + TRANS( 1 ) * 2.7 / TWO
    RHS( N ) = RTV%Planck_surface * surface_emissivity / TWO  ! Diffuse emissivity?
!    RHS( N ) = 304.368 * surface_emissivity / TWO  ! Diffuse emissivity?
    DIAG( 1 ) = ONE
    LOWER( 2 ) = ZERO
!    DIAG( N ) = -( ONE - surface_emissivity )
    DIAG( N ) = -surface_reflectivity
    LOWER( N ) = ONE
!    do i=1,n
!      print *, lower(i),diag(i),upper(i),rhs(i)
!    enddo
    CALL TRIDAG(LOWER, DIAG, UPPER, RHS, FLUX( 2, 1 ), N)
 
   ! FLUX is the Eddington fluxes at layer interfaces.
    ! FLUX(1,L) is upwelling, FLUX(2,L) is downwelling,
    !    L=1 is top, L=NUML+1 is bottom 
    FLUX( 1, 1 ) = FLUX( 2, 1 )
    FLUX( 2, 1 ) = cosmic_background / TWO
!    FLUX( 2, 1 ) = 2.7 / TWO

!print *, 'FLUX: ', flux(1,1),flux(2,1)
 
    ! Compute the homogeneous constants C+ and C-
    DO l = 1, N_Layers
      RADH1P = FLUX( 2, l ) + V( l ) - T1( l )
      RADH2M = FLUX( 1, l + 1 ) - V( l ) - T2( l )
      CPLUS( l )  = ( X2( l ) * RADH2M - X1( l ) * EXLM( l ) * RADH1P ) / D( l )
      CMINUS( l ) = ( X2( l ) * EXLP( l ) * RADH1P - X1( l ) * RADH2M ) / D( l )
!      print *, l, cplus(l), cminus(l)
!      print *, l, x1(l), x2(l),radh1p,radh2m,d(l)
!      print *, l, flux(2,l), flux(1,l+1)
    END DO
 
    !  Using the Eddington radiances, calculate the Eddington
    !  radiative transfer integral to get the downwelling radiance
    DOWNRAD = cosmic_background
!    DOWNRAD = 2.7
    DO l = 1, N_Layers
      DPLUS( l )  = CPLUS( l ) * ( D0( l ) - D1( l ) )
      DMINUS( l ) = CMINUS( l ) * ( D0( l ) + D1( l ) )
!      DELRAD = ( temp( l ) - B0MU( l ) ) * ( ONE - EX( l ) ) &
      DELRAD = ( RTV%Planck_Atmosphere( l ) - B0MU( l ) ) * ( ONE - EX( l ) ) &
               + B1( l ) * ( opd( l ) - mu * ( ONE - EX( l ) ) ) &
               + DPLUS( l ) / ( ONE + LAMBDA( l ) * mu ) * ( EXLP( l ) - EX( l ) ) &
               + DMINUS( l ) / ( ONE - LAMBDA( l ) * mu ) * ( EXLM( l ) - EX( l ) )
!      term1 = ( temp( l ) - B0MU( l ) ) * ( ONE - EX( l ) )
!      term2 = B1( l ) * ( opd( l ) - mu * ( ONE - EX( l ) ) )
!      term3 = DPLUS( l ) / ( ONE + LAMBDA( l ) * mu ) * ( EXLP( l ) - EX( l ) )
!      term4 = DMINUS( l ) / ( ONE - LAMBDA( l ) * mu ) * ( EXLM( l ) - EX( l ) )
!      DELRAD = term1 + term2 + term3 + term4
      DOWNRAD = DOWNRAD * EX( l ) + DELRAD
! write(*,*) 'DOWNRAD: ', l, term1, term2, term3, term4, DOWNRAD
! write(*,*) 'DOWNRAD: ', l, term1, term2, term3, term4, DOWNRAD
! write(*,*) 'DOWNRAD: ', l, B1(l), opd(l), ex(l), DOWNRAD
! print *, l, DOWNRAD, dplus(l),dminus(l),delrad
! print *, l, cplus(l),cminus(l),d0(l),d1(l)
    END DO
 
    ! The upwelling radiance at the surface is the surface emission
    ! plus the reflected downwelling radiance if specular reflection
    ! or the reflected downwelling flux if Lambertian refelection.
    IF ( SPECULAR ) THEN
!      UPRAD = surface_emissivity * RTV%Planck_Surface + ( ONE - surface_emissivity ) * DOWNRAD
      UPRAD = surface_emissivity * RTV%Planck_Surface + surface_reflectivity * DOWNRAD
!      UPRAD = surface_emissivity * 304.368 + ( ONE - surface_emissivity ) * DOWNRAD
    ELSE
      ! Lambertian
!      UPRAD = surface_emissivity * RTV%Planck_Surface + ( ONE - surface_emissivity ) * TWO &
!                * FLUX( 2, N_Layers + 1 )
      UPRAD = surface_emissivity * RTV%Planck_Surface + surface_reflectivity * TWO &
                * FLUX( 2, N_Layers + 1 )
!      UPRAD = surface_emissivity * 304.368 + ( ONE - surface_emissivity ) * TWO &
!                * FLUX( 2, N_Layers + 1 )
    END IF

!print *, 'UPRAD BC: ', UPRAD

    !  Using the Eddington radiances, calculate the Eddington
    !  radiative transfer integral to get the upwelling radiance
    DO l = N_Layers, 1, -1
      DPLUS( l )  =  CPLUS( l ) * ( D0( l ) + D1( l ) )
      DMINUS( l ) = CMINUS( l ) * ( D0( l ) - D1( l ) )
!      DELRAD = (temp( l ) + B0MU( l ) ) * ( ONE - EX( l ) ) &
      DELRAD = (RTV%Planck_Atmosphere( l ) + B0MU( l ) ) * ( ONE - EX( l ) ) &
              - B1( l ) * ( ( opd( l ) + mu ) * EX( l ) - mu ) &
              - DPLUS( l ) / ( ONE - LAMBDA( l ) * mu ) * ( EXLP( l ) * EX( l ) - ONE ) &
              - DMINUS( l ) / ( ONE + LAMBDA( l ) * mu ) * ( EXLM( l ) * EX( l ) - ONE )
      UPRAD = UPRAD * EX( l ) + DELRAD
! print *, l, UPRAD, dplus(l),dminus(l),delrad
    END DO

!     l = N_Layers
!     UPRAD_OLD = UPRAD
!     RAD_RATIO = 100.0_fp
!     DO WHILE ( ( l > 0 ) .AND. ( RAD_RATIO > RADIANCE_RATIO_THRESH ) ) 
!      DPLUS( l )  =  CPLUS( l ) * ( D0( l ) + D1( l ) )
!      DMINUS( l ) = CMINUS( l ) * ( D0( l ) - D1( l ) )
!      DELRAD = (RTV%Planck_Atmosphere( l ) + B0MU( l ) ) * ( ONE - EX( l ) ) &
!              - B1( l ) * ( ( opd( l ) + mu ) * EX( l ) - mu ) &
!              - DPLUS( l ) / ( ONE - LAMBDA( l ) * mu ) * ( EXLP( l ) * EX( l ) - ONE ) &
!              - DMINUS( l ) / ( ONE + LAMBDA( l ) * mu ) * ( EXLM( l ) * EX( l ) - ONE )
!      UPRAD = UPRAD_OLD * EX( l ) + DELRAD
!      RAD_RATIO = UPRAD / UPRAD_OLD
!print *, l, uprad, rad_ratio
!      UPRAD_OLD = UPRAD
!      l = l - 1
!    END DO
                                     
    RTV%s_Rad_UP  = UPRAD
!print *, 'TB: ', uprad, downrad

    END SUBROUTINE CRTM_EDD

    SUBROUTINE TRIDAG( A, B, C, R, U, N )

      IMPLICIT NONE

      !------------
      !  Arguments
      !------------

      INTEGER :: N
      REAL (fp) ::  A( N ), B( N ), C( N ), R( N ), U( N )

      ! Internal variables

      INTEGER, PARAMETER :: NMAX = 2 * MAX_N_LAYERS + 1
      INTEGER :: J
      REAL (fp) :: GAM( NMAX ), BET
 
      IF ( B( 1 ) .EQ. ZERO ) STOP 'TRIDAG FAILURE!'
      BET = B( 1 )
      U( 1 ) = R( 1 ) / BET
!print *, 'BET: ', BET
      DO J = 2, N
        GAM( J ) = C( J- 1 ) / BET
        BET = B( J ) - A( J ) * GAM( J )
        IF ( BET .EQ. ZERO ) STOP 'TRIDAG FAILURE!'
        U( J ) = ( R( J ) - A( J ) * U( J - 1 ) ) / BET
!print *, j, U(J),a(j),r(j),u(j-1),bet
      END DO
      DO J = N - 1, 1, -1
        U( J ) = U( J ) - GAM( J + 1 ) * U( J + 1 )
      END DO

      RETURN

    END SUBROUTINE TRIDAG
 
!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_EDD_Version
!
! PURPOSE:
!       Subroutine to return the module version information.
!
! CALLING SEQUENCE:
!       CALL CRTM_EDD_Version( Id )
!
! OUTPUT ARGUMENTS:
!       Id:            Character string containing the version Id information
!                      for the module.
!                      UNITS:      N/A
!                      TYPE:       CHARACTER(*)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_EDD_Version( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE CRTM_EDD_Version

END MODULE EDD_Module

