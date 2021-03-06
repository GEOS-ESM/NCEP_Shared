C-----------------------------------------------------------------------
      SUBROUTINE SPLAT(IDRT,JMAX,SLAT,WLAT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:  SPLAT      COMPUTE LATITUDE FUNCTIONS
C   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-02-20
C
C ABSTRACT: COMPUTES COSINES OF COLATITUDE AND GAUSSIAN WEIGHTS
C           FOR ONE OF THE FOLLOWING SPECIFIC GLOBAL SETS OF LATITUDES.
C             GAUSSIAN LATITUDES (IDRT=4)
C             EQUALLY-SPACED LATITUDES INCLUDING POLES (IDRT=0)
C             EQUALLY-SPACED LATITUDES EXCLUDING POLES (IDRT=256)
C           THE GAUSSIAN LATITUDES ARE LOCATED AT THE ZEROES OF THE
C           LEGENDRE POLYNOMIAL OF THE GIVEN ORDER.  THESE LATITUDES
C           ARE EFFICIENT FOR REVERSIBLE TRANSFORMS FROM SPECTRAL SPACE.
C           (ABOUT TWICE AS MANY EQUALLY-SPACED LATITUDES ARE NEEDED.)
C           THE WEIGHTS FOR THE EQUALLY-SPACED LATITUDES ARE BASED ON
C           ELLSAESSER (JAM,1966).  (NO WEIGHT IS GIVEN THE POLE POINT.)
C           NOTE THAT WHEN ANALYZING GRID TO SPECTRAL IN LATITUDE PAIRS,
C           IF AN EQUATOR POINT EXISTS, ITS WEIGHT SHOULD BE HALVED.
C           THIS VERSION INVOKES THE IBM ESSL MATRIX SOLVER.
C
C PROGRAM HISTORY LOG:
C   96-02-20  IREDELL
C   97-10-20  IREDELL  ADJUST PRECISION
C   98-06-11  IREDELL  GENERALIZE PRECISION USING FORTRAN 90 INTRINSIC
C 1998-12-03  IREDELL  GENERALIZE PRECISION FURTHER
C 1998-12-03  IREDELL  USE BLAS CALLS
C 2007-04-26  YANG     FOR EQUALLY-SPACED GRID EXCLUDING POLES, CHANGE 
c                      DIMENSION OF A FEW ARRAYS
C 2008-04-08  Todling  Replaced errexit() with exit()
C
C USAGE:    CALL SPLAT(IDRT,JMAX,SLAT,WLAT)
C
C   INPUT ARGUMENT LIST:
C     IDRT     - INTEGER GRID IDENTIFIER
C                (IDRT=4 FOR GAUSSIAN GRID,
C                 IDRT=0 FOR EQUALLY-SPACED GRID INCLUDING POLES,
C                 IDRT=256 FOR EQUALLY-SPACED GRID EXCLUDING POLES)
C     JMAX     - INTEGER NUMBER OF LATITUDES.
C
C   OUTPUT ARGUMENT LIST:
C     SLAT     - REAL (JMAX) SINES OF LATITUDE.
C     WLAT     - REAL (JMAX) GAUSSIAN WEIGHTS.
C
C SUBPROGRAMS CALLED:
C   DGEF         MATRIX FACTORIZATION
C   DGES         MATRIX SOLVER
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C
C$$$
      REAL SLAT(JMAX),WLAT(JMAX)
      INTEGER,PARAMETER:: KD=SELECTED_REAL_KIND(15,45)
      REAL(KIND=KD):: PK(JMAX/2),PKM1(JMAX/2),PKM2(JMAX/2)
      REAL(KIND=KD):: SLATD(JMAX/2),SP,SPMAX,EPS=10.*EPSILON(SP)
      PARAMETER(JZ=50)
      REAL BZ(JZ)
      DATA BZ        / 2.4048255577,  5.5200781103,
     $  8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679,
     $ 21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684,
     $ 33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132,
     $ 46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550,
     $ 58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299,
     $ 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711,
     $ 84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819,
     $ 96.6052679510, 99.7468198587, 102.888374254, 106.029930916,
     $ 109.171489649, 112.313050280, 115.454612653, 118.596176630,
     $ 121.737742088, 124.879308913, 128.020877005, 131.162446275,
     $ 134.304016638, 137.445588020, 140.587160352, 143.728733573,
     $ 146.870307625, 150.011882457, 153.153458019, 156.295034268 /
      REAL(8) DLT
      REAL(8) AWORK((JMAX+1)/2,((JMAX+1)/2)),BWORK(((JMAX+1)/2))
      INTEGER(4):: NRHS,LDA,LDB,INFO
      INTEGER(4):: JHE,JHO,J0
      INTEGER(4) IPVT((JMAX+1)/2)
      PARAMETER(PI=3.14159265358979,C=(1.-(2./PI)**2)*0.25)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GAUSSIAN LATITUDES
      J0=0
      IF(IDRT.EQ.4) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        R=1./SQRT((JMAX+0.5)**2+C)
        DO J=1,MIN(JH,JZ)
          SLATD(J)=COS(BZ(J)*R)
        ENDDO
        DO J=JZ+1,JH
          SLATD(J)=COS((BZ(JZ)+(J-JZ)*PI)*R)
        ENDDO
        SPMAX=1.
        DO WHILE(SPMAX.GT.EPS)
          SPMAX=0.
          DO J=1,JH
            PKM1(J)=1.
            PK(J)=SLATD(J)
          ENDDO
          DO N=2,JMAX
            DO J=1,JH
              PKM2(J)=PKM1(J)
              PKM1(J)=PK(J)
              PK(J)=((2*N-1)*SLATD(J)*PKM1(J)-(N-1)*PKM2(J))/N
            ENDDO
          ENDDO
          DO J=1,JH
            SP=PK(J)*(1.-SLATD(J)**2)/(JMAX*(PKM1(J)-SLATD(J)*PK(J)))
            SLATD(J)=SLATD(J)-SP
            SPMAX=MAX(SPMAX,ABS(SP))
          ENDDO
        ENDDO
CDIR$ IVDEP
        DO J=1,JH
          SLAT(J)=SLATD(J)
          WLAT(J)=(2.*(1.-SLATD(J)**2))/(JMAX*PKM1(J))**2
          SLAT(JMAX+1-J)=-SLAT(J)
          WLAT(JMAX+1-J)=WLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          SLAT(JHE)=0.
          WLAT(JHE)=2./JMAX**2
          DO N=2,JMAX,2
            WLAT(JHE)=WLAT(JHE)*N**2/(N-1)**2
          ENDDO
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  EQUALLY-SPACED LATITUDES INCLUDING POLES
      ELSEIF(IDRT.EQ.0) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE-1
        DLT=PI/(JMAX-1)
        SLAT(1)=1.
        DO J=2,JH
          SLAT(J)=COS((J-1)*DLT)
        ENDDO
        DO JS=1,JHO
          DO J=1,JHO
            AWORK(JS,J)=COS(2*(JS-1)*J*DLT)
          ENDDO
        ENDDO
        DO JS=1,JHO
          BWORK(JS)=-1./(4*(JS-1)**2-1)
        ENDDO
!cjfe   CALL DGEF(AWORK,JHE,JHO,IPVT)
!cjfe   CALL DGES(AWORK,JHE,JHO,IPVT,BWORK,J0)
	NRHS=1
	LDA=size(AWORK,1)
	LDB=size(BWORK,1)
        CALL DGESV(JHO,NRHS,AWORK,LDA,IPVT,BWORK,LDB,INFO)
	  if(INFO/=0) call die_at_dgesv_messages_(INFO)
        WLAT(1)=0.
        DO J=1,JHO
          WLAT(J+1)=BWORK(J)
        ENDDO
CDIR$ IVDEP
        DO J=1,JH
          SLAT(JMAX+1-J)=-SLAT(J)
          WLAT(JMAX+1-J)=WLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          SLAT(JHE)=0.
          WLAT(JHE)=2.*WLAT(JHE)
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  EQUALLY-SPACED LATITUDES EXCLUDING POLES
      ELSEIF(IDRT.EQ.256) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE
        DLT=PI/(JMAX+1)
        DO J=1,JH
          SLAT(J)=COS(J*DLT)
        ENDDO
        DO JS=1,JHO
          DO J=1,JHO
            AWORK(JS,J)=COS(2*(JS-1)*J*DLT)
          ENDDO
        ENDDO
        DO JS=1,JHO
          BWORK(JS)=-1./(4*(JS-1)**2-1)
        ENDDO
!cjfe   CALL DGEF(AWORK1,JHE,JHO,IPVT)
!cjfe   CALL DGES(AWORK1,JHE,JHO,IPVT,BWORK,J0)
	NRHS=1
	LDA=size(AWORK,1)
	LDB=size(BWORK,1)
        CALL DGESV(JHO,NRHS,AWORK,LDA,IPVT,BWORK,LDB,INFO)
	  if(INFO/=0) call die_at_dgesv_messages_(INFO)
        DO J=1,JHO
          WLAT(J)=BWORK(J)
        ENDDO
CDIR$ IVDEP
        DO J=1,JH
          SLAT(JMAX+1-J)=-SLAT(J)
          WLAT(JMAX+1-J)=WLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          SLAT(JHE)=0.
          WLAT(JHE)=2.*WLAT(JHE)
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
	contains
	subroutine die_at_dgesv_messages_(INFO)
	integer(4),intent(in) :: INFO

	integer(4),parameter :: IEXIT=2

	if(INFO<0) then
	  write(*,*) "NCEP_sp/splat(): >>>> ERROR <<<< the ",
     &		-INFO,"-th argument had an illegal value"
	endif

	if(INFO>0) then
	  write(*,*) "NCEP_sp/splat(): >>>> ERROR <<<< U(",
     &		INFO,",",INFO,") is singular"
	endif

	call exit(IEXIT)
	end subroutine die_at_dgesv_messages_
      END
