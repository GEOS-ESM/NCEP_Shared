C-----------------------------------------------------------------------
      SUBROUTINE SPFFTE(IMAX,INCW,INCG,KMAX,W,G,IDIR,AFFT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:  SPFFTE     PERFORM MULTIPLE FAST FOURIER TRANSFORMS
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-02-20
C
C ABSTRACT: THIS SUBPROGRAM PERFORMS MULTIPLE FAST FOURIER TRANSFORMS
C           BETWEEN COMPLEX AMPLITUDES IN FOURIER SPACE AND REAL VALUES
C           IN CYCLIC PHYSICAL SPACE.
C           SUBPROGRAM SPFFTE MUST BE INVOKED FIRST WITH IDIR=0
C           TO INITIALIZE TRIGONEMETRIC DATA.  USE SUBPROGRAM SPFFT1
C           TO PERFORM AN FFT WITHOUT PREVIOUS INITIALIZATION.
C           THIS VERSION INVOKES A GENERIC FFT.
C
C PROGRAM HISTORY LOG:
C 1998-12-18  IREDELL
C 2007-04-26  R.YANG  use a general FFT (SGI library)
C
C USAGE:    CALL SPFFTE(IMAX,INCW,INCG,KMAX,W,G,IDIR,AFFT)
C
C   INPUT ARGUMENT LIST:
C     IMAX     - INTEGER NUMBER OF VALUES IN THE CYCLIC PHYSICAL SPACE
C                (SEE LIMITATIONS ON IMAX IN REMARKS BELOW.)
C     INCW     - INTEGER FIRST DIMENSION OF THE COMPLEX AMPLITUDE ARRAY
C                (INCW >= IMAX/2+1)
C     INCG     - INTEGER FIRST DIMENSION OF THE REAL VALUE ARRAY
C                (INCG >= IMAX)
C     KMAX     - INTEGER NUMBER OF TRANSFORMS TO PERFORM
C     W        - COMPLEX(INCW,KMAX) COMPLEX AMPLITUDES IF IDIR>0
C     G        - REAL(INCG,KMAX) REAL VALUES IF IDIR<0
C     IDIR     - INTEGER DIRECTION FLAG
C                IDIR=0 TO INITIALIZE TRIGONOMETRIC DATA
C                IDIR>0 TO TRANSFORM FROM FOURIER TO PHYSICAL SPACE
C                IDIR<0 TO TRANSFORM FROM PHYSICAL TO FOURIER SPACE
C     AFFT       REAL(8) (25+2*IMAX) AUXILIARY ARRAY IF IDIR<>0
C
C   OUTPUT ARGUMENT LIST:
C     W        - COMPLEX(INCW,KMAX) COMPLEX AMPLITUDES IF IDIR<0
C     G        - REAL(INCG,KMAX) REAL VALUES IF IDIR>0
C     AFFT       REAL(8) (25+2*IMAX) AUXILIARY ARRAY IF IDIR=0
C
C SUBPROGRAMS CALLED:
C   SCRFT        IBM ESSL COMPLEX TO REAL FOURIER TRANSFORM
C   DCRFT        IBM ESSL COMPLEX TO REAL FOURIER TRANSFORM
C   SRCFT        IBM ESSL REAL TO COMPLEX FOURIER TRANSFORM
C   DRCFT        IBM ESSL REAL TO COMPLEX FOURIER TRANSFORM
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C
C REMARKS:
C   THE RESTRICTIONS ON IMAX ARE THAT IT MUST BE A MULTIPLE
C   OF 1 TO 25 FACTORS OF TWO, UP TO 2 FACTORS OF THREE,
C   AND UP TO 1 FACTOR OF FIVE, SEVEN AND ELEVEN.
C
C   IF IDIR=0, THEN W AND G NEED NOT CONTAIN ANY VALID DATA.
C   THE OTHER PARAMETERS MUST BE SUPPLIED AND CANNOT CHANGE
C   IN SUCCEEDING CALLS UNTIL THE NEXT TIME IT IS CALLED WITH IDIR=0.
C
C   THIS SUBPROGRAM IS THREAD-SAFE.
C
C$$$
        IMPLICIT NONE
        INTEGER,INTENT(IN):: IMAX,INCW,INCG,KMAX,IDIR
        REAL,INTENT(INOUT):: W(2*INCW,KMAX)
        REAL,INTENT(INOUT):: G(INCG,KMAX)
        INTEGER:: I,K

        REAL(8) :: AFFT(*) ! expects IMAX+256 by sgi_fft()

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c local variable
        integer j,ico,ico2,imaxd2
        REAL(8) :: dummy(imax)
        REAL(8) :: field(imax)
        COMPLEX(8) :: spcoef(imax/2 +1)
c-----------------------------------

c       write (6,*) 'SUBROUTINE SPFFTE IS CALLED, ldafft=',ldafft

        imaxd2 = imax/2 + 1
        SELECT CASE(IDIR)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C 1.  INITIALIZATION.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CASE(0)
          CALL sgi_fft (0,imax,1,dummy,1,dummy,1,afft)
c check afft

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C 2.  FOURIER TO PHYSICAL TRANSFORM.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CASE(1:)
c        input:   w is input
c        output: g is output field
c   arrays mapping: w is zonal wave coefficient, mapping w to spcoef
c they both have the same number of coefficient
          
          do k=1,kmax
c---------------------------------------
c  do transformation on each latitude
c---------------------------------------
             do j =1,imaxd2
              ico = (j-1)*2 +1
              ico2 = ico + 1
              spcoef(j) = cmplx(w(ico,k),w(ico2,k))
             enddo
          call sgi_fft (1,imax,1,field,imax,spcoef,imaxd2,afft)

c   mapping the field back to the NCEP arrays
            DO i=1,IMAX
              G(i,K)=field(i)
            ENDDO
          ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C 3.  PHYSICAL TO FOURIER TRANSFORM.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CASE(:-1)
c       write (6,*) 'in spffte: transfer from fld-sp'
c  input: g (imax,kmax)
c  output: w (imax+2,kmax)
          do k=1,kmax
             do i=1,imax
             field(i)=g(i,k)
             enddo
c  
c       write (6,*) 'in spffte do physical to foure transform,input fld'
c       write (6,*)  (field (i),i=1,200)
   
          call sgi_fft (-1,imax,1,field,imax,spcoef,imaxd2,afft)
c  mapping spcoef to NCEP arrays
            do i=1,imaxd2
              ico= (i-1)*2 + 1
              ico2= ico + 1 
              w(ico,k) = real(spcoef(i))
              w(ico2,k) = imag(spcoef(i))
            enddo
c check whether the w(2,k) and w(imax+2,k) =0.0 
c             write (6,*) 'in after sgi_fft, coeff. '

c  do the same assignment as in the NCEP spffte code for these two elements:
            W(2,K)=0.
            W(IMAX+2,K)=0.
          enddo
        END SELECT
      END SUBROUTINE  SPFFTE
