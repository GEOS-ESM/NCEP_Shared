!
! P2S_Module
!
! Module containing the Polarized 2-Stream (P2S) radiative
! transfer solution procedures used in the CRTM.
!
!
! CREATION HISTORY:
!       Written by:     Tom Greenwald, CIMSS/SSEC; tom.greenwald@ssec.wisc.edu
!


MODULE P2S_Module

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE Type_Kinds,      ONLY: fp
  USE Message_Handler, ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters, ONLY: SET, ZERO, ONE, TWO, FOUR, PI, &
                             MAX_N_LAYERS, MAX_N_ANGLES, MAX_N_LEGENDRE_TERMS, &
                             DEGREES_TO_RADIANS, &
                             SCATTERING_ALBEDO_THRESHOLD, &
                             OPTICAL_DEPTH_THRESHOLD
  USE CRTM_Utility
  USE RTV_Define     , ONLY: RTV_type, &
                             DELTA_OPTICAL_DEPTH, &
                             MAX_N_QUAD, MAX_N_LEG, MAX_N_STOKES

  ! Disable all implicit typing
  IMPLICIT NONE

  ! --------------------
  ! Default visibilities
  ! --------------------
  ! Everything private by default
  PRIVATE
  ! Procedures
  PUBLIC :: CRTM_P2S
  PUBLIC :: CRTM_P2S_Version

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*),  PARAMETER :: MODULE_VERSION_ID = &
!     '$Id$'
     '$Id$'

CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

    SUBROUTINE CRTM_P2S ( n_Layers, & ! Input, Number of atmospheric layers
                                 w, & ! Input, Layer single scattering albedo
                                 g, & ! Input, Layer asymmetry parameter
                              T_OD, & ! Input, Layer optical depth
                 cosmic_background, & ! Input, Cosmic background radiance
                        emissivity, & ! Input, Surface emissivity
                      reflectivity, & ! Input, Surface reflectivity
                   Index_Sat_Angle, & ! Input, Satellite angle index 
                               RTV)   ! IN/Output, Upward radiance and others

! ------------------------------------------------------------------------- !
!                                                                           !
! FUNCTION:                                                                 !
!                                                                           !
!   This subroutine calculates polarized (I,Q) IR/MW radiance at the top of !
!   a scattering atmosphere using the two-stream approximation              !
!                                                                           ! 
!   Modified and adapted by Tom Greenwald    tomg@ssec.wisc.edu             !
!                                                                           !
!   References:                                                             !
!                                                                           !
!      Liu, Q., and F. Weng, 2002: A microwave polarimetric two-stream      !
!          radiative transfer model. J. Atmos. Sci., 59, 2396-2402.         !
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
    TYPE(RTV_type), INTENT( INOUT ) :: RTV


   ! -------------- internal variables --------------------------------- !

    INTEGER, PARAMETER :: select = 1     ! 0: Standard, 1: Improved, 2: Polarimetric
    INTEGER, PARAMETER :: n_Stokes = 1   ! Number of Stokes parameters actually computed
    INTEGER, PARAMETER :: n_Leg = 12     ! Number of expansion terms actually used
    INTEGER, PARAMETER :: n_Quad = 6
    INTEGER i
    !  Gaussian points and weights for 6 quadrature nodes
    REAL (fp), DIMENSION( n_Quad ) :: an6 = [ 1.00, 0.944899272, 0.819279322, 0.632876153, 0.399530941, 0.136552933 ]
    REAL (fp), DIMENSION( n_Quad ) :: anw6 = [ 0.015151515, 0.091684517, 0.157974706, 0.212508418, 0.251275603, 0.271405241 ]
    REAL (fp) :: mu
    REAL (fp) :: surface_emissivity( MAX_N_STOKES ), surface_reflectivity( MAX_N_STOKES ) 

    ! Precomputed Legendre polynomials for Gaussian nodes
    REAL (fp) :: ple( MAX_N_QUAD, MAX_N_LEG ) = [ &
      1.000000000E+00,  1.000000000E+00,  1.000000000E+00,  1.000000000E+00,  1.000000000E+00,  1.000000000E+00,&
      1.000000000E+00,  9.448992610E-01,  8.192793131E-01,  6.328761578E-01,  3.995309472E-01,  1.365529299E-01,&
      1.000000000E+00,  8.392518759E-01,  5.068279505E-01,  1.007983685E-01, -2.605625391E-01, -4.720299542E-01,&
      1.000000000E+00,  6.917479634E-01,  1.458698511E-01, -3.155959547E-01, -4.398585856E-01, -1.984637380E-01,&
      1.000000000E+00,  5.144174099E-01, -1.709816903E-01, -4.251317680E-01, -1.121180654E-01,  3.065960705E-01,&
      1.000000000E+00,  3.215323985E-01, -3.688430488E-01, -2.318236083E-01,  2.712565064E-01,  2.341308594E-01,&
      1.000000000E+00,  1.283143312E-01, -4.115220010E-01,  8.529780060E-02,  2.921198905E-01, -1.968827695E-01,&
      1.000000000E+00, -5.043155700E-02, -3.099872470E-01,  2.989600003E-01, -1.575669274E-02, -2.506127357E-01,&
      1.000000000E+00, -2.016239315E-01, -1.161047220E-01,  2.801231444E-01, -2.674085498E-01,  1.081063598E-01,&
      1.000000000E+00, -3.150323331E-01,  9.586896002E-02,  6.912615895E-02, -1.877991408E-01,  2.506510913E-01,&
      1.000000000E+00, -3.841187358E-01,  2.537268102E-01, -1.689890772E-01,  9.810771793E-02, -3.226415440E-02,&
      1.000000000E+00, -4.065181911E-01,  3.096951246E-01, -2.670176327E-01,  2.455572635E-01, -2.362756282E-01,&
      1.000000000E+00, -3.841187656E-01,  2.537268102E-01, -1.689891070E-01,  9.810774773E-02, -3.226410970E-02,&
      1.000000000E+00, -3.227400780E-01,  1.138836071E-01,  4.080633074E-02, -1.512892395E-01,  2.096279711E-01,&
      1.000000000E+00, -2.314493656E-01, -5.566297099E-02,  2.067244947E-01, -2.076720297E-01,  8.516549319E-02,&
      1.000000000E+00, -1.215888336E-01, -1.944581717E-01,  2.148540318E-01, -1.920808107E-02, -1.731688827E-01,&
      1.000000000E+00, -5.614042282E-03, -2.564898431E-01,  6.964926422E-02,  1.798237115E-01, -1.256581694E-01,&
      1.000000000E+00,  1.041391790E-01, -2.248931974E-01, -1.166497469E-01,  1.575422883E-01,  1.296738386E-01,&
      1.000000000E+00,  1.966374964E-01, -1.160241514E-01, -2.093281597E-01, -4.744429141E-02,  1.531081051E-01,&
      1.000000000E+00,  2.631679773E-01,  2.794728801E-02, -1.474747509E-01, -1.861638576E-01, -8.213456720E-02 ]


    ! Cosine of satellite observation zenith angle
    mu = RTV%COS_Angle( Index_Sat_Angle )

    IF ( select .GT. 0 ) THEN
      !-----------------------------------------------------
      !  Compute scattering phase function expansion 
      !  terms assuming a Henyey-Greenstein phase function
      !-----------------------------------------------------

      CALL scat_phase_HG( &
                        n_Layers,   & ! Input, Number of atmospheric layers
                        n_Leg,      & ! Input, Number of Legendre expansion terms to use
                        g,          & ! Input, Layer asymmetry parameter
                        n_Stokes,   & ! Input, Number of Stokes parameters to compute
                        RTV         ) ! IN/Output, Scattering phase matrix

      !--------------------------------------------------
      !  Compute polarization difference qq and
      !  forward and backward scattering phase functions
      !--------------------------------------------------

      CALL phaseq( &
                 n_Layers,       & ! Input, Number of atmospheric layers
                 n_Leg,          & ! Input, Number of Legendre polynomial terms
                 n_Quad,         & ! Input, Number of quadrature nodes
                 n_Stokes,       & ! Input, Number of Stokes parameters to compute
                 mu,             & ! Input, Cosine of observation angle
                 ple,            & ! Input, Legendre polynomials at quadrature nodes
                 RTV             ) ! Input/output, Polarization difference, scattering phase matrix,
                                   !               forward/backward phase functions
    END IF

    !--------------------------------------------------
    ! Compute reflection and transmission matrices
    !--------------------------------------------------

    CALL t_r( &
           n_Layers,           & ! Input, Number of atmospheric layers
           n_Quad,             &  ! Input, Number of quadrature nodes
           mu,                 & ! Input, Cosine of observation angle
           w,                  & ! Input, Layer single scattering albedo
           g,                  & ! Input, Layer asymmetry parameter
           T_OD,               & ! Input, Layer optical depth
           select,            & ! Input, Selector for 2-stream type
           an6,                & ! Input, Quadrature nodes 
           anw6,               & ! Input, Quadrature weights
           RTV                 ) ! Input/Output, qq, pff, pbb, transmission, reflection, qq1

    !--------------------------
    ! Do radiative transfer
    !--------------------------

    IF ( n_Stokes .EQ. 1 ) THEN
      surface_emissivity( 1 ) = emissivity( Index_Sat_Angle )
      surface_reflectivity( 1 ) = reflectivity( Index_Sat_Angle, Index_Sat_Angle ) 
    ELSE
      ! Need message handler here
    END IF

    CALL two_RT_model( &
                   n_Layers,             & ! Input, 
                   surface_emissivity,   & ! Input, 
                   surface_reflectivity, & ! Input, 
                   cosmic_background,    & ! Input, 
                   n_Stokes,             & ! Input,
                   RTV                   ) ! Input/Output, 


    END SUBROUTINE CRTM_P2S


    SUBROUTINE two_RT_model( n_Layers,       & ! Input, Number of atmospheric layers
                             emissivity,     & ! Input, Surface emissivity
                             reflectivity,   & ! Input, Surface reflectivity
                             tsky,           & ! Input
                             n_Stokes,       & ! Input
                             RTV             ) ! Input/Output


    IMPLICIT NONE

    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    !------------
    !  Arguments
    !------------

    INTEGER, INTENT( IN ) :: n_Layers
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: emissivity
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: reflectivity
    REAL (fp), INTENT( IN ) :: tsky
    INTEGER, INTENT( IN ) :: n_Stokes
    TYPE(RTV_type), INTENT( INOUT ) :: RTV


    !-----------------
    ! Local variables
    !-----------------

    INTEGER i, is
    REAL (fp) :: Rdo( MAX_N_LAYERS )
    REAL (fp ) :: Tdo( MAX_N_LAYERS )
    REAL (fp) :: Sdo( MAX_N_STOKES, MAX_N_LAYERS )
    REAL (fp) :: Source( MAX_N_STOKES, MAX_N_LAYERS )
    REAL (fp) :: Tbup( MAX_N_LAYERS )
    REAL (fp) :: Tbdo( MAX_N_LAYERS )
    REAL (fp) :: r
    REAL (fp) :: a
    REAL (fp) :: b
    REAL (fp) :: c( MAX_N_LAYERS )
    REAL (fp) :: aa
    REAL (fp) :: bb
    REAL (fp) :: sign
    REAL (fp) :: den


!    write(11,*) emissivity(1), reflectivity(1), tsky

    DO is = 1, n_Stokes
      IF ( n_Stokes .LT. 3 ) THEN
        sdo( is, 1 ) = 0.0

        rdo( 1 ) = RTV%reflection( 1 )

        IF ( is .EQ. 1 ) THEN
          sign = ONE
        ELSE
          sign = -ONE
        END IF
        source( is, 1 ) = ( ONE - RTV%transmission( 1 ) - RTV%reflection( 1 ) ) * RTV%Planck_Atmosphere( 1 ) &
                           + sign * RTV%qq( 1 ) * RTV%PLanck_Atmosphere( 1 )

        tdo( 1 ) = 1.0
        i_nlayer_loop: DO i = 2, n_Layers
          r = rdo( i - 1 ) * RTV%reflection( i )

          source( is, i ) = ( ONE - RTV%transmission( i ) - RTV%reflection( i ) ) * RTV%Planck_Atmosphere( i ) &
                             + sign * RTV%qq( i ) * RTV%Planck_Atmosphere( i )

          rdo( i ) = RTV%reflection( i ) + RTV%transmission(i) * rdo( i - 1 ) * RTV%transmission( i ) / ( ONE - r )
! write(11,'(i3,4e12.4)') i, RTV%Reflection( i ), Rdo( i - 1 ), RTV%transmission(i), r
          tdo( i ) = RTV%transmission( i ) * tdo( i - 1 ) / ( ONE - r )

          sdo( is, i ) = sdo( is, i - 1 ) * RTV%transmission( i ) / ( ONE - r ) + source( is, i ) &
                        + RTV%Transmission( i ) * rdo( i - 1 ) / ( ONE - r ) * source( is, i )
        END DO i_nlayer_loop

        !--------------------
        !   at the surface
        !--------------------

        a = emissivity( is ) * RTV%Planck_Surface
        b = Rdo( n_Layers ) * a + Sdo( is, n_Layers ) + Tdo( n_Layers ) * Tsky
!        r = ( 1.0 - emissivity( is ) ) * Rdo( n_Layers )
!  Replaced 1-e with reflectivity - TG 6/22/2017
        r = reflectivity( is ) * Rdo( n_Layers )

!        Tbup( n_Layers + 1 ) = a + ( 1.0 - emissivity( is ) ) * b / ( 1.0 - r )
!  Replaced 1-e with reflectivity - TG 6/22/2017
        Tbup( n_Layers + 1 ) = a + reflectivity( is ) * b / ( ONE - r )
        Tbdo( n_Layers + 1 ) = Sdo( is, n_Layers ) + Rdo( n_Layers ) * Tbup( n_Layers + 1 ) &
                              + Tdo( n_Layers ) * Tsky

        !---------------------------------------------------
        ! 	 from surface to the top of the atmosphere
        !---------------------------------------------------
        i_nlayer_loop_reverse: DO i = n_Layers, 1, -1

          IF ( i .GT. 1 ) THEN

            !-----------------------------------------------
            !  calculating upward atmospheric contribution
            !-----------------------------------------------
            den = ONE - RTV%Reflection( i ) * Rdo( i - 1 )
            aa = RTV%Transmission( i ) * r / den
            bb = (Source( is, i ) + RTV%Reflection( i ) * Sdo( is, i - 1 ) ) / den
            r = aa + bb
            a = RTV%Transmission( i ) * Tbup( i + 1 ) / den
            b = (Source( is, i ) + RTV%Reflection( i ) * Sdo( is, i - 1 ) ) / den
            IF ( is .EQ. 1 ) c( i ) = Tdo( i - 1 ) * RTV%Reflection( i ) * Tsky / den
            Tbup( i ) = a + b + c( i )

          ELSE

            Tbup( i ) = RTV%Transmission( i ) * Tbup( i + 1 ) + Source( is, i ) + RTV%Reflection( i ) * Tsky

          ENDIF

! write(11,'(i3,6e12.4)') i, a, b, c(i), den, RTV%Reflection( i ), Rdo( i - 1 )

        END DO i_nlayer_loop_reverse

        ! Note: Currently set up just for n_Stokes=1
        RTV%s_Rad_UP = Tbup( 1 )

      END IF

    END DO


    END SUBROUTINE two_RT_model 



    SUBROUTINE p_two( w0,       &   ! Input, Single scattering albedo for this layer
                      gg,       &   ! Input, Asymmetry parameter for this layer
                      cu,       &   ! Input, Cosine of observation angle
                      q,        &   ! Input, Polarization difference for this layer
                      tau,      &   ! Input, Optical depth for this layer
                      t1,       &   ! Output, Polarimetric transmission for this layer
                      r1        )   ! Output, Polarimetric reflection for this layer
                      
    IMPLICIT NONE

    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    !------------
    !  Arguments
    !------------

    REAL (fp), INTENT( IN ) :: w0
    REAL (fp), INTENT( IN ) :: gg
    REAL (fp), INTENT( IN ) :: cu
    REAL (fp), INTENT( IN ) :: q
    REAL (fp), INTENT( IN ) :: tau
    REAL (fp), INTENT( INOUT ), DIMENSION( : ) :: t1
    REAL (fp), INTENT( INOUT ), DIMENSION( : ) :: r1 

    !-----------------
    ! Local variables
    !-----------------

    REAL (fp) :: amda1
    REAL (fp) :: amda2
    REAL (fp) :: a
    REAL (fp) :: b
    REAL (fp) :: c
    REAL (fp) :: d
    REAL (fp) :: a11
    REAL (fp) :: a12
    REAL (fp) :: a21
    REAL (fp) :: a22
    REAL (fp) :: det
    REAL (fp) :: a8
    REAL (fp) :: b8
    REAL (fp) :: t11
    REAL (fp) :: d11
    REAL (fp) :: t12
    REAL (fp) :: d12
    REAL (fp) :: r11
    REAL (fp) :: r12
    REAL (fp) :: u11
    REAL (fp) :: u12

    amda1 = ( ONE - w0 ) * ( ONE - w0 * gg ) - w0 * q * ( ONE - w0 * gg ) * ( ONE - w0 )
    amda1 = dsqrt( amda1 )
    amda2 = ( ONE - w0 ) * ( ONE - w0 * gg ) + w0 * q * ( ONE - w0 * gg ) * ( ONE - w0 )
    amda2 = dsqrt( amda2 )

    a = ( ( ONE + ( ONE - w0 * gg ) / amda1 ) ** 2 ) * dexp( amda1 * tau / cu )
    b = ( ( ONE + ( ONE - w0 * gg ) / amda2 ) ** 2 ) * dexp( amda2 * tau / cu )
    c = ( ( ONE - ( ONE - w0 * gg ) / amda1 ) ** 2 ) * dexp( -amda1 * tau / cu )
    d = ( ( ONE - ( ONE - w0 * gg ) / amda2 ) ** 2 ) * dexp( -amda2 * tau / cu )
    a11 = ( a + b - ( c + d ) ) / TWO
    a12 = ( ( a - b ) - ( c - d ) ) / TWO
    a21 = a12
    a22 = a11
    det = a11 * a22 - a12 * a21
    a8 = ( ( ONE + ( ONE - w0 * gg ) / amda1 ) ** 2 )
    b8 = ( ( ONE + ( ONE - w0 * gg ) / amda2 ) ** 2 )
    t11 = ( ( a8 + b8 ) * a22 - a21 * ( a8 - b8 ) ) / det / TWO
    t12 = ( -a12 * ( a8 + b8 ) + a11 * ( a8 - b8 ) ) / det / TWO  

    a8 = ( ( ONE - ( ONE - w0 * gg ) / amda1 ) ** 2 )
    b8 = ( ( ONE - ( ONE - w0 * gg ) / amda2 ) ** 2 )

    d11 = -( ( a8 + b8 ) * a22 - a21 * ( a8 - b8 ) ) / det / 2.0 
    d12 = -( -a12 * ( a8 + b8 ) + a11 * ( a8 - b8 ) ) / det / 2.0      
    a8 = ( ONE + ( ONE - w0 * gg ) / amda1 ) * ( ONE - ( ONE - w0 * gg ) / amda1 ) * &
         dexp( -amda1 * tau / cu )
    b8 = ( ONE + ( ONE - w0 * gg ) / amda2 ) * ( ONE - ( ONE - w0 * gg ) / amda2 ) * &
         dexp( -amda2 * tau / cu )
    r11 = ( ( a8 + b8 ) * a22 - a21 * ( a8 - b8 ) ) / det / TWO 
    r12 = ( -a12 * ( a8 + b8 ) + a11 * ( a8 - b8 ) ) / det / TWO   

    a8 = ( ONE + ( ONE - w0 * gg ) / amda1 ) * ( ONE - ( ONE - w0 * gg ) / amda1 ) * &
         dexp( amda1 * tau / cu )
    b8 = ( ONE + ( ONE - w0 * gg ) / amda2 ) * ( ONE - ( ONE - w0 * gg ) / amda2 ) * &
         dexp( amda2 * tau / cu )
    u11 = -( ( a8 + b8 ) * a22 - a21 * ( a8 - b8 ) ) / det / TWO
    u12 = -( -a12 * ( a8 + b8 ) + a11 * ( a8 - b8 ) ) / det / TWO 
  
    t1( 1 ) = t11 + d11
    t1( 2 ) = t12 + d12
    t1( 4 ) = t1( 1 )
    t1( 3 ) = t1( 2 )
    r1( 1 ) = r11 + u11
    r1( 2 ) = r12 + u12
    r1( 3 ) = r1( 2 )
    r1( 4 ) = r1( 1 )

    END SUBROUTINE p_two


    SUBROUTINE t_r( n_Layers,       &  ! Input, Number of atmospheric layers
                    n_Quad,         &  ! Input, Number of quadrature nodes
                    mu,             &  ! Input, Cosine of observation zenith angle
                    omega,          &  ! Input, Layer single scattering albedo
                    asy_g,          &  ! Input, Layer asymmetry parameter
                    tau,            &  ! Input, Layer optical depth
                    select,         &  ! Input, Selector for 2-stream type
                    an6,            &  ! Input, Quadrature nodes
                    anw6,           &  ! Input, Quadrature weights
                    RTV             )  ! Output, qq, qq1, pff, pbb, transmission, reflection

    IMPLICIT NONE

    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    !------------
    !  Arguments
    !------------
    
    INTEGER, INTENT( IN ) :: n_Layers, n_Quad
    REAL (fp), INTENT( IN ) :: mu
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: omega
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: asy_g
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: tau
    INTEGER, INTENT( IN ) :: select
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: an6
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: anw6
    TYPE(RTV_type), INTENT( INOUT ) :: RTV

    !-----------------
    ! Local variables
    !-----------------

    INTEGER i
    INTEGER j
    REAL (fp) :: a
    REAL (fp) :: b
    REAL (fp) :: c
    REAL (fp) :: d
    REAL (fp) :: k
    REAL (fp) :: t1
    REAL (fp) :: r1
    REAL (fp) :: ssa, asm, opd, ssa_tmp
    REAL (fp) :: expdif
    REAL (fp) :: t11( MAX_N_STOKES )
    REAL (fp) :: r11( MAX_N_STOKES )        
    REAL (fp) :: pf( MAX_N_QUAD )
    REAL (fp) :: pb( MAX_N_QUAD )
      
! select=0 using two-stream method; select=1 using improved two-stream method
! select=2 using polarized two-stream method

    i_nlayer_loop: DO i = 1, n_Layers

      RTV%qq1( i ) = ZERO

      ! QC 
      IF ( omega(i) .GE. ONE ) THEN
        ssa_tmp = 0.999
      ELSE
        ssa_tmp = omega(i)
      END IF

      ! Apply delta-scaling
      opd = ( ONE - ssa_tmp * asy_g( i ) ** 2 ) * tau( i )
      ssa = ( ONE - asy_g( i ) ** 2 ) * ssa_tmp / &
            ( ONE - asy_g( i ) ** 2 * ssa_tmp )
      asm = asy_g( i ) / ( ONE + asy_g( i ) )

      IF ( select .EQ. 0 ) THEN

        !--------------------------------
        !  Standard two-stream method
        !--------------------------------
        a = ( ONE - ssa )
        b = ( ONE - ssa * asm )
        k = SQRT( a * b )
        expdif = EXP( k * opd / mu ) - EXP( -k * opd / mu )
        c = EXP( k * opd / mu ) + EXP( -k * opd / mu ) + ( SQRT( a / b ) &
            + SQRT( b / a ) ) / TWO * expdif
        d = ( SQRT( b / a ) - SQRT( a / b ) ) * expdif / FOUR
        RTV%transmission( i ) = TWO / c
        RTV%reflection( i ) = RTV%transmission( i ) * d
!write(11,'(i3,6e12.4)') i, ssa, asy_g(i), a, b, c, d
      ELSE IF (select .EQ. 1) THEN

        !--------------------------------
        !  Improved two-stream method
        !--------------------------------
        CALL s_c_two(an6,        & ! Input
                     anw6,       & ! Input
                     n_Quad,     & ! Input
                     RTV%pff2( :, i ), & ! Input 
                     RTV%pbb2( :, i ), & ! Input
!                     omega( i ), & ! Input
                     ssa, & ! Input
                     asy_g( i ), & ! Input
                     tau( i ),   & ! Input
                     mu,         & ! Input
                     r1,         & ! Output
                     t1          ) ! Output
        RTV%transmission( i ) = t1
        RTV%reflection( i ) = r1
      ELSE

        !----------------------------------
        !  Polarimetric two-stream method
        !----------------------------------      
!        CALL p_two( omega( i ), & ! Input
        CALL p_two( ssa, & ! Input
                    asy_g( i ), & ! Input
                    mu,         & ! Input 
                    RTV%qq( i ),& ! Input
                    tau( i ),   & ! Input
                    t11,        & ! Output
                    r11         ) ! Output
        RTV%qq1( i ) = ( t11( 2 ) + r11( 2 ) )
        RTV%transmission( i ) = t11( 1 )
        RTV%reflection( i ) = r11( 1 )
        CALL s_c_two(an6,        & ! Input
                     anw6,       & ! Input
                     n_Quad,     & ! Input
                     RTV%pff2( :, i ), & ! Input 
                     RTV%pbb2( :, i ), & ! Input
!                     omega( i ), & ! Input
                     ssa, & ! Input
                     asy_g( i ), & ! Input
                     tau( i ),   & ! Input
                     mu,         & ! Input
                     r1,         & ! Output
                     t1          ) ! Output
        !
        ! Why are these arrays overwritten?
        !
        RTV%transmission( i ) = t1
        RTV%reflection( i ) = r1
     ENDIF
      
    END DO i_nlayer_loop

    END SUBROUTINE t_r


    SUBROUTINE s_c_two ( angx,          &   ! Input, Quadrature points
                         angw,          &   ! Input, Quadrature weights
                         n_Quad,        &   ! Input, Number of quadrature nodes
                         pf,            &   ! Input, Forward phase function for this layer
                         pb,            &   ! Input, Backward phase function for this layer
                         w0,            &   ! Input, Single scattering albedo for this layer
                         gg,            &   ! Input, Asymmetry parameter for this layer
                         tau,           &   ! Input, Optical depth for this layer
                         mu,            &   ! Input, Cosine of observation zenith angle
                         r1,            &   ! Output, Reflection for this layer
                         t1             )   ! Output, Transmission for this layer

    IMPLICIT NONE

    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    !------------
    !  Arguments
    !------------
    
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: angx
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: angw
    INTEGER, INTENT( IN ) :: n_Quad
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: pf
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: pb
    REAL (fp), INTENT( IN ) :: w0
    REAL (fp), INTENT( IN ) :: gg
    REAL (fp), INTENT( IN ) :: tau
    REAL (fp), INTENT( IN ) :: mu 
    REAL (fp), INTENT( OUT ) :: r1
    REAL (fp), INTENT( OUT ) :: t1
    
    !-----------------
    ! Local variables
    !-----------------

    INTEGER i, j
    REAL (fp) :: ak
    REAL (fp) :: a1
    REAL (fp) :: a11
    REAL (fp) :: a2
    REAL (fp) :: cu
    REAL (fp) :: det
    REAL (fp) :: Ra1( MAX_N_QUAD )
    REAL (fp) :: Rb1( MAX_N_QUAD )
    REAL (fp) :: Ta1( MAX_N_QUAD )
    REAL (fp) :: Tb1( MAX_N_QUAD )


! //   corrected two-stream model
    
    cu = mu
    IF ( tau .GT. 0.0001 ) THEN
      ak = dsqrt( ( ONE - w0 ) * ( ONE - w0 * gg ) )
      i_quad_loop: DO i = 1, n_Quad
        a1 = ( ( ONE - w0 + ak ) + ( ONE - w0 - ak ) * dexp( -tau * ak / angx( i ) ) )
        det = ( ( ONE - w0 - ak ) ** 2 ) * dexp( -ak * tau / angx( i ) ) - ( ( ONE - w0 + ak ) ** 2 ) * &
              dexp( ak * tau / angx( i ) )
        a2 = ( -( ONE - w0 - ak ) - ( ONE - w0 + ak ) * dexp( tau * ak / angx( i ) ) )
        a1 = a1 / det
        a2 = a2 / det
        a11 = ONE - cu * ak / angx( i )
        IF ( ABS( a11 ) .GT. 0.000001 ) THEN
          Ta1( i ) = ( ONE - w0 + ak ) * a1 * ( ONE - dexp( -tau * ( ONE / cu - ak / angx( i ) ) ) ) / &
                     ( ONE - cu * ak / angx( i ) )
          Tb1( i ) = -( ONE - w0 - ak ) * a1 * ( ONE - dexp( -tau * ( -ak / angx( i ) + ONE / cu ) ) ) / &
                      ( ONE - cu * ak / angx( i ) )
        ELSE
          Ta1( i ) = ( ONE - w0 + ak ) * a1 * tau / cu
          Tb1( i ) = -( ONE - w0 - ak ) * a1 * tau / cu
        ENDIF
        Ra1( i ) = ( ONE - w0 - ak ) * a2 * ( ONE - dexp( -tau * ( ak / angx( i ) + ONE / cu ) ) ) / &
                   ( ONE + cu * ak / angx( i ) )
        Rb1( i ) = -( ONE - w0 + ak ) * a2 * ( ONE - dexp( -tau * ( ONE / cu + ak / angx( i ) ) ) ) / &
                    ( ONE + cu * ak / angx( i ) )
      END DO i_quad_loop
      a1 = ZERO
      a2 = ZERO
      j_quad_loop: DO j = 1, n_Quad
        a1 = a1 + w0 / TWO * ( pf( j ) * Ta1( j ) + pb( j ) * Tb1( j ) ) * angw( j )
        a2 = a2 + w0 / TWO * (pf( j ) * Ra1( j ) + pb( j ) * Rb1( j ) ) * angw( j )
      END DO j_quad_loop
      t1 = dexp( -tau / cu ) - a1
      r1 = -a2
    ELSE
      t1 = dexp( -tau / cu )
      r1 = 0.0
    ENDIF

    END SUBROUTINE s_c_two


!----------------------------------------------------------
!     Calculation of polarization difference qq and
!     forward and backward scattering phase function.
!----------------------------------------------------------

    SUBROUTINE phaseq( n_Layers,    &   ! Input
                       n_Leg,       &   ! Input
                       n_Quad,      &   ! Input
                       n_Stokes,    &   ! Input
                       mu,          &   ! Input
                       ple,         &   ! Input
                       RTV          ) ! Input/output, Polarization difference, forward/backward phase functions

    IMPLICIT NONE

    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    !------------
    !  Arguments
    !------------

    INTEGER, INTENT( IN ) :: n_Layers, n_Leg, n_Quad, n_Stokes
    REAL (fp), INTENT( IN ) :: mu
    REAL (fp), INTENT( IN ), DIMENSION( :, : ) :: ple
    TYPE(RTV_type), INTENT( INOUT ) :: RTV
    
    !-----------------
    ! Local variables
    !-----------------

    INTEGER i, j, k
    REAL (fp) ::  x
    REAL (fp) ::  pp1( MAX_N_LEG )
    REAL (fp) ::  pp2( MAX_N_LEG )

    !--------------------------------------------
    !  Compute polynomials at observation angle
    !--------------------------------------------
    CALL pleg( mu, n_Leg, pp1)
    CALL pleg( -mu, n_Leg, pp2)
 
    i_layer_loop: DO i = 1, n_Layers

      RTV%qq( i ) = ZERO

      k_quad_loop: DO k = 1, n_Quad

        RTV%pff2( k, i ) = ZERO
        RTV%pbb2( k, i ) = ZERO

      END DO k_quad_loop

      IF ( n_Stokes .EQ. 1 ) THEN

        k_quad2_loop: DO k = 1, n_Quad

          j_legp2_loop: DO j = 1, n_Leg 

!****** Is this correct???
            RTV%pff2( k, i ) = RTV%pff2( k, i ) + RTV%scat_phase( i, j, 1 ) * pp1( j ) * ple( k, j ) / TWO
            RTV%pbb2( k, i ) = RTV%pbb2( k, i ) + RTV%scat_phase( i, j, 1 ) * pp2( j ) * ple( k, j ) / TWO

          END DO j_legp2_loop

        END DO k_quad2_loop

      ELSE

        j_legp_loop: DO j = 1, n_Leg

    !--------------------------------------------
    !  Polarization difference of phase elements
    !--------------------------------------------
          RTV%qq( i ) = RTV%qq( i ) + ( RTV%scat_phase( i, j, 2 ) - RTV%scat_phase( i, j, 1 ) ) * pp1( j ) / 16.0_fp

        END DO j_legp_loop

        k_quad3_loop: DO k = 1, n_Quad

          j_legp3_loop: DO j = 1, n_Leg 

            RTV%pff2( k, i ) = RTV%pff2( k, i ) + ( RTV%scat_phase( i, j, 2 ) + RTV%scat_phase( i, j, 1 ) ) * &
                               pp1( j ) * ple( k, j ) / TWO
            RTV%pbb2( k, i ) = RTV%pbb2( k, i ) + ( RTV%scat_Phase( i, j, 2 ) + RTV%scat_phase( i, j, 1 ) ) * &
                               pp2( j ) * ple( k, j ) / TWO

          END DO j_legp3_loop

        END DO k_quad3_loop

      END IF      

    END DO i_layer_loop


    END SUBROUTINE phaseq


!---------------------------
!   Legendre polynomials
!---------------------------

    SUBROUTINE pleg( x,      &   ! Input
                     nleg,   &   ! Input
                     p       )   ! Output

    IMPLICIT NONE

    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    !------------
    !  Arguments
    !------------
    REAL (fp), INTENT( IN ) :: x
    INTEGER, INTENT( IN )   :: nleg
    REAL (fp), INTENT( IN OUT ), DIMENSION( : ) :: p
    
    INTEGER i

    P( 1 ) = ONE
    P( 2 ) = x
    i_nleg_loop: DO i = 3, nleg
      p( i ) = ( ( 2 * i - 3 ) * x * p( i - 1 ) - ( i - 2 ) * p( i - 2 ) ) / ( i - 1 )
    END DO i_nleg_loop
    
    END SUBROUTINE pleg

    SUBROUTINE scat_phase_HG( &
                        n_Layers,   & ! Input, Number of atmospheric layers
                        n_Leg,      & ! Input, Number of Legendre expansion terms to use
                        g,          & ! Input, Asymmetry parameter
                        n_Stokes,   & ! Input, Number of Stokes parameters to use
                        RTV         ) ! Input/Output, Scattering phase matrix

    IMPLICIT NONE

    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    !------------
    !  Arguments
    !------------
    INTEGER, INTENT( IN ) :: n_Layers, n_Leg
    REAL (fp), INTENT( IN ), DIMENSION( : ) :: g
    INTEGER, INTENT( IN ) :: n_Stokes
    TYPE(RTV_type), INTENT( INOUT ) :: RTV

    !-----------------
    ! Local variables
    !-----------------

    INTEGER :: i, j, is

!-------------------------------------------------------------------------------------
! Set up Legendre polynomial expansion coefficients for scattering phase matrix using
! the Henyey-Greenstein phase function
!-------------------------------------------------------------------------------------

    DO is = 1, n_Stokes
      DO i = 1, n_Layers
        IF ( g( i ) .EQ. 0 ) THEN
          RTV%scat_phase( i, 1 : n_Leg, is ) = ZERO
        ELSE
          RTV%scat_phase( i, 1, is ) = ONE
          DO j = 2, n_Leg
! Is this correct?
            RTV%scat_phase( i, j, is ) = ( 2 * ( j - 1 ) + 1 ) * g( i ) ** ( j - 1 )
          END DO
        END IF
      END DO
    END DO

    END SUBROUTINE scat_phase_HG

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_P2S_Version
!
! PURPOSE:
!       Subroutine to return the module version information.
!
! CALLING SEQUENCE:
!       CALL CRTM_P2S_Version( Id )
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

  SUBROUTINE CRTM_P2S_Version( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE CRTM_P2s_Version


END MODULE P2S_Module
