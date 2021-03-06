      SUBROUTINE W3CNVTOV (IBUFTN,IFLDUN,stnid,INSTR,KINDX)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    W3CNVTOV    CONVERT RTOVS RPT -NMCEDS TO IW3UNPBF FMT
C   PRGMMR: KEYSER           ORG: NP22       DATE: 1998-09-21
C
C ABSTRACT: CONVERTS AN RTOVS REPORT ORIGINALLY IN UNPACKED NMCEDS
C   FORMAT TO UNPACKED IW3UNPBF FORMAT (AS DESCRIBED IN W3LIB
C   IW3UNPBF).  THE UNPACKED NMCEDS FORMAT IS FILLED ONLY WITH THOSE
C   VALUES NEEDED FOR RTOVS PROCESSING BY THE PREPDATA PROGRAM.
C
C PROGRAM HISTORY LOG:
C 1998-02-17  D. A. KEYSER -- ORIGINAL AUTHOR (ADAPTED FROM W3LIB
C        ROUTINE W3FI43).
C 1998-06-15  D. A. KEYSER -- ADAPTED FOR USE ONLY WITH RTOVS DATA
C        (AFTER TOVS DEMISE) (I.E., NO PARTLY-CLOUDY PATH AVAILABLE,
C        ONLY ESSENTIAL DATA IN UNPACKED NMCEDS FORMAT); WRITES OUT
C        ONLY CATEGORY 1 IW3UNPBF DATA SINCE THIS IS ALL THAT IS
C        PROCESSED IN PREPDATA; OTHERWISE STREAMLINED
C 1998-09-21  D .A. KEYSER -- SUBROUTINE NOW Y2K AND FORTRAN 90
C        COMPLIANT
C
C USAGE:    CALL W3CNVTOV(IBUFTN,IFLDUN,STNID,INSTR,KINDX)
C   INPUT ARGUMENT LIST:
C     IBUFTN   - ADDRESS HOLDING A SINGLE RTOVS REPORT (140 INTEGER
C              - WORDS) IN UNPACKED NMCEDS FORMAT (THE UNPACKED NMCEDS
C              - FORMAT IS FILLED ONLY WITH THOSE VALUES NEEDED FOR
C              - RTOVS PROCESSING BY THE PREPDATA PROGRAM)
C     INSTR    - INDICATOR FOR RETRIEVAL PATH (EITHER 1 FOR CLEAR OR
C              - 3 FOR CLOUDY)
C     KINDX    - INTEGER  1-5 DIGIT NUMBER USED TO GENERATE FIRST
C              - 5 CHARACTERS OF STATION ID (USUALLY JUST A REPORT
C              - COUNTER INDEX EXCEPT FIRST NUMBER MAY BE NADIR
C              - PROXIMITY INDICATOR -- SEE PREPDATA PROGRAM)
C
C   OUTPUT ARGUMENT LIST:
C     IFLDUN   - INTEGER 273-WORD ARRAY HOLDING A SINGLE RTOVS
C              - REPORT IN UNPACKED IW3UNPBF FORMAT (SEE W3LIB
C              - IW3UNPBF) (NOTE: DOES NOT INCLUDE STATION ID)
C     STNID    - CHARACTER*8 SINGLE REPORT STATION IDENTIFICATION (UP
C              - TO 8 CHARACTERS, LEFT-JUSTIFIED)
C
C
C REMARKS: MUST BE CALLED AFTER CALL TO W3FA07 WHICH FILLS IN VALUES
C          IN COMMON BLOCK /FA07AA/.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM-SP, CRAY, SGI
C
C$$$
 
      REAL  RDATA(273),PMAND(20)
 
      INTEGER  IDATA(273),IFLDUN(273),IBUFTN(140)
 
      CHARACTER*1  CSUFX(3,4)
      CHARACTER*8  STNID
 
      COMMON/FA07AA/TM(20),Z(20),MLVLS
 
      EQUIVALENCE (IDATA,RDATA)
 
      SAVE
 
      DATA  CSUFX /'W','?','Y',    'E','?','G',
     $             'S','?','U',    'A','?','C'/
      DATA  XMSG/99999./,IMSG/99999/
      DATA  PMAND/10000.,8500.,7000.,5000.,4000.,3000.,2500.,2000.,
     $ 1500.,1000.,700.,500.,300.,200.,100.,70.,50.,30.,20.,10./
 

C  INITIALIZE ALL CATEGORY TYPES AND NUMBER OF LEVELS TO ZERO
 
      IDATA(13:42) = 0
 
C  ALLOWS 21 LEVELS FOR CAT. 1 - THIS IS THE ONLY CATEGORY THAT IS
C   PROCESSED  (SET ALL WORDS IN CAT. 1 TO MISSING)
 
      RDATA(43:273) = XMSG
 
C  SET Q.M.'S TO 2 FOR CATEGORY 1 (DEFAULT)
 
      RDATA(49:269:11) = 2.0
      RDATA(50:270:11) = 2.0
      RDATA(51:271:11) = 2.0
      RDATA(52:272:11) = 2.0
      RDATA(53:273:11) = 2.0
 
C  FILL IN IW3UNPBF REPORT HEADER
 
      RDATA(1)  = IBUFTN(5)/100.
      RDATA(2)  = IBUFTN(6)/100.
      IF(IBUFTN(6).LT.0)  RDATA(2) = 360. + IBUFTN(6)/100.
      IF(RDATA(2).EQ.360.0)  RDATA(2) = 0.0
      RDATA(3)  = 0.
      IHR  = MOD(IBUFTN(3),256)
      IB4  = IBUFTN(4)/256
      XMIN = IB4/60.
      RDATA(4)  = IHR + XMIN
      IDATA(5)  = IMSG
      IDATA(6)  = IMSG
      RDATA(7)  = IBUFTN(8)
      IDATA(8)  = IMSG
      IDATA(9)  = 61
      RDATA(10) = 0.
      RDATA(11) = XMSG
      IDATA(12) = IMSG
 
C  STN. ID: POS. 1-5 FROM 'KINDX', POS. 6 FROM CHAR. BASED ON SATELLITE
C   NUMBER AND RETRIEVAL PATH
 
C       POSITION 6 CHARACTER POSSIBILITIES ARE:
C            ODD  SATELLITE NUMBERS 3, 7, 11, 15, ETC.:  A, C
C            ODD  SATELLITE NUMBERS 1, 5,  9, 13, ETC.:  E, G
C            EVEN SATELLITE NUMBERS 2, 6, 10, 14, ETC.:  S, U
C            EVEN SATELLITE NUMBERS 4, 8, 12, 16, ETC.:  W, Y
C    WHERE: CHARACTERS  A, E, S, W  ARE FOR CLEAR PATH (DEFAULT)
C           CHARACTERS  C, G, U, Y  ARE FOR CLOUDY (MICROWAVE) PATH
 
      MODSAT = MOD(IBUFTN(1),4) + 1
 
C  STATION IDENTIFICATION IN "STNID" (8 CHARACTERS)
 
      WRITE(STNID,1)  KINDX,CSUFX(INSTR,MODSAT)
    1 FORMAT(I5.5,A1,'  ')
 
C  FILL IN IW3UNPBF CATEGORY 1 (MANDATORY LEVEL DATA)
 
      IDATA(13) = MLVLS + 1
      IDATA(14) = 43
      K = 43
 
      DO I = 1,MLVLS
         IF(I.EQ.2)  THEN
            RDATA(K) = 9250.
            K = K + 11
         END IF
         RDATA(K) = PMAND(I)
         RDATA(K+1) = Z(I) + 0.5
         IF(TM(I).LT.10273.)  RDATA(K+2) = (TM(I) - 273.16) * 10.
         K = K + 11
      ENDDO
 
C  COPY IW3UNPBF FIELD (IDATA) TO IFLDUN FOR TRANSFER OUT OF SUBROUTINE
 
      IFLDUN = IDATA
 
      RETURN
      END
