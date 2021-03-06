C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:  W3ERSUNB      DECODES SINGLE ERS RPT FROM BUFR MESSAGES
C   PRGMMR: KEYSER           ORG: NP22        DATE: 1999-01-20
C
C ABSTRACT: THIS SUBROUTINE DECODES A SINGLE ERS SCATTEROMETER REPORT
C   FROM BUFR MESSAGES IN A SEQUENTIAL DATA FILE.  REPORT IS RETURNED
C   IN THE FORMAT DESCRIBED IN THE REMARKS 4.  PROCESSES ONLY EDITION
C   2 OR GREATER MESSAGES.
C
C PROGRAM HISTORY LOG:
C 1998-02-17  KEYSER -- ORIGINAL AUTHOR
C 1998-06-15  D. A. KEYSER -- REDEFINED UNITS FOR UNPACKED WORDS 1
C             (LATITUDE), 2 (LONGITUDE), 4 (OBS. TIME) AND 11
C             (RECEIPT TIME) - ALL TO CONFORM WITH UNPACKED IW3UNPBF
C             FORMAT AND TO STREAMLINE PROCESSING IN PREPDATA PROGRAM
C 1998-09-21  D. A. KEYSER -- SUBROUTINE NOW Y2K AND FORTRAN 90
C             COMPLIANT
C 1998-10-09  D. A. KEYSER -- ADDED SECONDS TO DECODED REPORT DATE
C             (ONLY USED IN PRINTS OF REPORT DATE)
C 1999-01-20 D. A. KEYSER -- INCORPORATED BOB KISTLER'S CHANGES NEEDED
C             TO PORT THE CODE TO THE IBM SP
C
C USAGE:    CALL W3ERSUNB(IDATE,IHE,IHL,INDTA,INTBB,INTBD,RDATA,STNID,
C          $              IRET)
C   INPUT ARGUMENT LIST:
C     IDATE    - 4-WORD ARRAY HOLDING "CENTRAL" DATE TO PROCESS
C              - (YYYY, MM, DD, HH)
C     IHE      - NUMBER OF WHOLE HOURS RELATIVE TO "IDATE" FOR DATE OF
C              - EARLIEST BUFR MESSAGE THAT IS TO BE DECODED; EARLIEST
C              - DATE IS "IDATE" + "IHE" HOURS (IF "IHE" IS POSITIVE,
C              - LATEST MESSAGE DATE IS AFTER "IDATE"; IF "IHE" IS
C              - NEGATIVE LATEST MESSAGE DATE IS PRIOR TO "IDATE")
C              - EXAMPLE: IF IHE=1, THEN EARLIEST DATE IS 1-HR AFTER
C              - IDATE; IF IHE=-3, THEN EARLIEST DATE IS 3-HR PRIOR
C              - TO IDATE
C     IHL      - NUMBER OF WHOLE HOURS RELATIVE TO "IDATE" FOR DATE OF
C              - LATEST BUFR MESSAGE THAT IS TO BE DECODED; LATEST
C              - DATE IS "IDATE" + ("IHL" HOURS PLUS 59 MIN) IF "IHL"
C              - IS POSITIVE (LATEST MESSAGE DATE IS AFTER "IDATE"),
C              - AND "IDATE" + ("IHL"+1 HOURS MINUS 1 MIN) IF "IHL"
C              - IS NEGATIVE (LATEST MESSAGE DATE IS PRIOR TO "IDATE")
C              - EXAMPLE: IF IHL=3, THEN LATEST DATE IS 3-HR 59-MIN
C              - AFTER IDATE; IF IHL=-2, THEN LATEST DATE IS 1-HR 1-MIN
C              - PRIOR TO IDATE
C     INDTA    - FORTRAN UNIT NUMBER FOR INPUT DATA FILE
C     INTBB    - FORTRAN UNIT NUMBER FOR BUFR TABLE B FILE
C     INTBD    - FORTRAN UNIT NUMBER FOR BUFR TABLE D FILE
C     IRET     - CONTROLS DEGREE OF UNIT 6 PRINTOUT (.GE. 0 -LIMITED
C              - PRINTOUT; = -1 SOME ADDITIONAL DIAGNOSTIC PRINTOUT;
C              = .LT. -1 -EXTENSIVE PRINTOUT) (SEE REMARKS 3.)
C
C   OUTPUT ARGUMENT LIST:
C     RDATA    - SINGLE REPORT RETURNED THE FORMAT DESCRIBED IN THE
C              - REMARKS 4 SECTION OF THIS DOCBLOCK (MUST BE
C              - DIMENSIONED TO AT LEAST 14 WORDS BY CALLING PROGRAM)
C              - (NOTE: DOES NOT INCLUDE STATION ID)
C     STNID    - CHARACTER*8 SINGLE REPORT STATION IDENTIFICATION (UP
C              - TO 8 CHARACTERS, LEFT-JUSTIFIED)
C     IRET     - RETURN CODE AS FOLLOWS:
C       IRET = 0 ---> REPORT SUCCESSFULLY RETURNED
C       IRET > 0 ---> NO REPORT RETURNED DUE TO:
C            = 1 ---> ALL REPORTS READ IN, END
C            = 2 ---> EITHER LAT AND/OR LON DESCRIPTOR NOT FOUND, OR
C                      LAT AND/OR LON DATA MISSING
C            = 3 ---> EITHER SOME/ALL DESCRIPTORS FOR DATE INFORMATION
C                      NOT FOUND, OR SOME/ALL DATE INFORMATION MISSING
C
C   INPUT FILES:
C     UNIT AA  - (WHERE AA IS INDTA ABOVE) FILE HOLDING THE DATA IN
C              - THE FORM OF SEQUENTIAL RECORDS CONTAINING BUFR MSGS
C     UNIT BB  - (WHERE BB IS INTBB ABOVE) BUFR TABLE B
C     UNIT CC  - (WHERE CC IS INTBD ABOVE) BUFR TABLE D
C
C   OUTPUT FILES:
C     UNIT 06  - PRINTOUT
C
C   SUBPROGRAMS CALLED:
C     UNIQUE        - ERSUNB01 ERSUNB02 ERSUNB03 ERSUNB04
C     LIBRARY:
C       COMMON      - EXIT
C       W3LIB       - W3FI04   W3FI88   W3MOVDAT W3DIFDAT
C                   - GBYTE    GBYTES
C
C REMARKS: 1) A CONDITION CODE (STOP) OF 15 WILL OCCUR IF THE INPUT
C     DATES FOR START AND/OR STOP TIME ARE SPECIFIED INCORRECTLY.
C          2) A CONDITION CODE (STOP) OF 22 WILL OCCUR IF THE
C     CHARACTERS ON THIS MACHINE ARE NEITHER ASCII NOR EBCDIC.
C          3) THE INPUT ARGUMENT "IRET" SHOULD BE SET PRIOR TO EACH
C     CALL TO THIS SUBROUTINE.
C
C   ***************************************************************
C          4)
C    BELOW IS THE FORMAT OF AN UNPACKED REPORT IN OUTPUT ARRAY RDATA
C     (EACH WORD REPRESENTS A FULL-WORD ACCORDING TO THE MACHINE)
C     (NOTE: DOES NOT INCLUDE STATION IDENTIFICATION - SEE OUTPUT
C                      ARGUMENT "STNID" ABOVE)
C   ***************************************************************
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                FORMAT FOR ERS SCATTEROMETER REPORTS
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                           HEADER
C   WORD   CONTENT                   UNIT                 FORMAT
C   ----   ----------------------    -------------------  ---------
C     1    LATITUDE                  DEGREES (N+,S-)      REAL
C     2    LONGITUDE                 DEGREES EAST         REAL
C     3    NOT USED                  ZEROED
C     4    OBSERVATION TIME          HOURS (UTC)          REAL
C     5    YEAR/MONTH                YYYYMM               INTEGER
C     6    DAY/HOUR                  DDHH                 INTEGER
C     7    STATION ELEVATION         10 METERS (CONSTANT) REAL
C     8    INSTRUMENT TYPE           99 (MISSING)         INTEGER
C     9    REPORT TYPE               581 (CONSTANT)       INTEGER
C    10    NOT USED                  ZEROED
C    11    RECEIPT TIME              HOURS (UTC)          REAL
C    12    NOT USED                  ZEROED
C                                    LEFT-JUSTIFIED
C    13    HORIZ. WIND DIRECTION     DEGREES              REAL
C    14    HORIZ. WIND SPEED         0.1 METERS/REC       REAL
C
C  MISSING VALUES ARE:
C         99999. FOR REAL
C         99999  FOR INTEGER (EXCEPT FOR WORDS 5 AND 6, WHERE
C                 MISSING IS 999999)
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM-SP, CRAY, SGI
C
C$$$
      SUBROUTINE W3ERSUNB(IDATE,IHE,IHL,INDTA,INTBB,INTBD,RDATA,STNID,
     $ IRET)
C PARAMETER "MAXD" IS THE NUMBER OF EXPANDED DESCRIPTORS
C  NEEDED TO CONTAIN A SINGLE REPORT (MAXIMUM = 25 FOR ERS
C  SCATTEROMETER REPORTS)
C PARAMETER "MAXR" IS THE NUMBER OF SUBSETS IN A SINGLE BUFR MESSAGE
C  (MAXIMUM = 500 FOR ERS SCATTEROMETER REPORTS)
C PARAMETER "MAXL" IS THE SIZE NEEDED FOR THE ILOC AND SC ARRAYS FOR A
C  SINGLE REPORT (MAXIMUM = 25 FOR ERS SCATTEROMETER REPORTS)
      PARAMETER (MAXD = 25, MAXR = 500, MAXL =  25)
      INTEGER  IDATE(4),LSDATE(4),IDATA(14),IPTR(45)
      INTEGER JDATE(8),K1DATE(8),L1DATE(8)
      CHARACTER*8  STNID
      DIMENSION RINC(5)
      REAL  RDATA(*),RDATX(14)
      COMMON /ERSUBB/IST,KDATE(8),LDATE(8),MSTACK(2,MAXD),
     $ KDATA(MAXR,MAXD),IPRINT
      COMMON /ERSUCC/INDEX,ILOC(MAXL),SC(MAXL)
cTEMP#############
c note:  currently the date in sec. 1 appears to be garbage, so we must
c        date check each individual report here
      COMMON/TMPTIM/IYEAR,IMNTH,IDAY,IHOUR,IMIN,ISEC
cTEMP#############
C
      SAVE
C
      EQUIVALENCE(RDATX,IDATA)
      DATA ITM/0/,INDTAL/-99/,KOUNT/0/
cTEMP#############
c note:  currently the date in sec. 1 appears to be garbage, so we must
c        date check each individual report here
      DATA ITCNT/0/
cTEMP#############
      IPRINT = 0

      write(*,*) 'Sorry missing date routines, cannot run'
      iret = 1
      RETURN
      END
