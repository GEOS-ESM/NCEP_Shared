C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:  W3UNPKB7      DECODES SINGLE REPORT FROM BUFR MESSAGES
C   PRGMMR: KEYSER           ORG: NP22        DATE: 1999-01-20
C
C ABSTRACT: THIS SUBROUTINE DECODES A SINGLE REPORT FROM BUFR MESSAGES
C   IN A JBUFR-TYPE DATA FILE.  CURRENTLY WIND PROFILER, NEXRAD (VAD)
C   WIND AND GOES SOUNDING/RADIANCE DATA TYPES ARE VALID.  REPORT IS
C   RETURNED IN QUASI-IW3UNPBF UNPACKED FORMAT (SEE REMARKS 4.). ALSO,
C   INFORMATION ABOUT THE INPUT DATA SET ITSELF (NAME, CENTER DATE,
C   DUMP TIME) IS RETURNED
C
C PROGRAM HISTORY LOG:
C 1998-02-17  D. A. KEYSER -- ORIGINAL AUTHOR (BASED ON W3LIB ROUTINE
C             W3UNPK77)
C 1998-06-14  D. A. KEYSER -- MODIFIED WIND PROFILER CAT. 11 (HEIGHT,
C             HORIZ. SIGNIFICANCE, VERT. SIGNIFICANCE) AND VAD WIND
C             HEADER (STATION ID) PROCESSING TO ACCOUNT FOR UPDATES
C             TO BUFRTABLE MNEMONICS IN /dcom; CHANGED CHAR. 6 OF GOES
C             STNID TO BE UNIQUE FOR TWO DIFFERENT EVEN OR ODD
C             SATELLITE ID'S (EVERY OTHER EVEN OR ODD SAT. ID NOW GETS
C             SAME CHAR. 6 TAG)
C 1998-06-15  D. A. KEYSER -- REDEFINED UNITS FOR UNPACKED WORDS 1
C             (LATITUDE), 2 (LONGITUDE), 4 (OBS. TIME) AND 11
C             (RECEIPT TIME) - ALL TO CONFORM WITH UNPACKED IW3UNPBF
C             FORMAT AND TO STREAMLINE PROCESSING IN PREPDATA PROGRAM
C 1998-09-21  D. A. KEYSER -- SUBROUTINE NOW Y2K AND FORTRAN 90
C             COMPLIANT
C 1998-10-09  D. A. KEYSER -- CORRECTED ERROR IN RETURNING GOES
C             CAT. 8 DATA WHEN GREATER THAN 9 "LEVELS" ARE PRESENT
C 1999-01-20 D. A. KEYSER -- INCORPORATED BOB KISTLER'S CHANGES NEEDED
C             TO PORT THE CODE TO THE IBM SP
C
C USAGE:    CALL W3UNPKB7(IDATE,IHE,IHL,LUNIT,RDATA,STNID,DSNAME,
C                         IDSDAT,IDSDMP,IRET)
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
C     LUNIT    - FORTRAN UNIT NUMBER FOR INPUT DATA FILE
C     IRET     - CONTROLS DEGREE OF UNIT 6 PRINTOUT (.GE. 0 -LIMITED
C              - PRINTOUT; = -1 SOME ADDITIONAL DIAGNOSTIC PRINTOUT;
C              = .LT. -1 -EXTENSIVE PRINTOUT) (SEE REMARKS 3.)
C
C   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
C     RDATA    - SINGLE REPORT RETURNED AN A QUASI-IW3UNPBF
C              - UNPACKED FORMAT (SEE REMARKS 4.) (MINIMUM SIZE IS
C              - 720 WORDS) (NOTE: DOES NOT INCLUDE STATION ID)
C     STNID    - CHARACTER*8 SINGLE REPORT STATION IDENTIFICATION (UP
C              - TO 8 CHARACTERS, LEFT-JUSTIFIED) (SEE % IN DOCBLOCK
C              - FOR INFO ON GOES SOUNDING STNID)
C     DSNAME   - CHARACTER*8 DATA SET NAME (SAME FOR ALL REPORTS IN
C              - A COMMON INPUT DATA SET - SEE REMARKS FOR IRET=1)
C     IDSDAT   - INTEGER DATA SET CENTER DATE IN FORM YYYYMMDDHH (SAME
C              - FOR ALL REPORTS IN A COMMON INPUT DATA SET - SEE
C              - REMARKS FOR IRET=1)
C     IDSDMP   - INTEGER DATA SET DUMP TIME IN FORM YYYYMMDDHHMM (SAME
C              - FOR ALL REPORTS IN A COMMON INPUT DATA SET - SEE
C              - REMARKS FOR IRET=1)
C     IRET     - RETURN CODE AS FOLLOWS:
C            =  0 OBSERVATION READ AND UNPACKED INTO OUTPUT ARGUMENT
C                 LOCATIONS (SEE ABOVE).  SEE REMARKS FOR CONTENTS.
C                 NEXT CALL TO W3UNPKB7 WILL RETURN NEXT OBSERVATION
C                 IN DATA SET.
C            =  1 INFORMATION ABOUT THE BUFR DATASET IS RETURNED IN
C                 THE OUTPUT ARGUMENTS DSNAME, IDSDAT, IDSDMP (SEE
C                 OUTPUT ARGUMENT LIST ABOVE)
C
C                 THIS SHOULD ALWAYS OCCUR AFTER THE FIRST CALL TO
C                 THIS SUBROUTINE FOR A NEW UNIT NUMBER (I.E., A NEW
C                 DATA SET) OR IF THE INPUT DATE/TIME OR RANGE IN TIME
C                 HAS BEEN CHANGED FROM LAST CALL (I.E., A DATA SET
C                 IS ABOUT TO BE READ FROM THE TOP).  NO REPORT IS
C                 UNPACKED AT THIS POINT, AND ONLY DSNAME, IDSDAT, AND
C                 IDSDMP CONTAIN INFORMATION.  ALL SUBSEQUENT CALLS TO
C                 W3UNPKB7 SHOULD RETURN THE OBSERVATIONS IN THIS DATA
C                 SET, SEQUENTIALLY, (IRET=0) UNTIL THE END OF FILE IS
C                 ENCOUNTERED (IRET=2).  THE VALUES STORED IN DSNAME,
C                 IDSDAT, AND IDSDMP WILL CONTINUE TO BE RETURNED ALONG
C                 WITH EACH REPORT WHEN IRET = 0.
C            =  2 FOR NORMAL END-OF-FILE ENCOUNTERED.
C            =  3 LAT AND/OR LON DATA MISSING -- NO REPORT RETURNED.
C            =  4 SOME/ALL DATE INFORMATION MISSING -- NO REPORT
C                 RETURNED.
C            =  5 NO DATA LEVELS PROCESSED (ALL LEVELS ARE MISSING) --
C                 NO REPORT RETURNED.
C            =  6 NUMBER OF LEVELS IN REPORT HEADER IS NOT 1 -- NO
C                 REPORT RETURNED.
C            =  7 NUMBER OF LEVELS IN ANOTHER SINGLE LEVEL SEQUENCE IS
C                 NOT 1 -- NO REPORT RETURNED.
C
C   INPUT FILES:
C     UNIT AA  - (WHERE AA IS LUNIT ABOVE) FILE HOLDING THE DATA
C              - IN THE FORM OF BUFR MESSAGES
C
C   OUTPUT FILES:
C     UNIT 06  - STANDARD OUTPUT PRINT
C
C   SUBPROGRAMS CALLED:
C     UNIQUE        - UNPKB701 UNPKB702 UNPKB703 UNPKB704 UNPKB705
C                   - UNPKB706 UNPKB707 UNPKB708 UNPKB709
C     LIBRARY:
C       COMMON      - EXIT
C       W3LIB       - W3FI04   W3MOVDAT W3DIFDAT
C       BUFRLIB     - DATELEN  DUMPBF   OPENBF   READMG   UFBCNT
C                   - READSB   UFBINT   CLOSBF
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
C                FORMAT FOR WIND PROFILER REPORTS
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                           HEADER
C   WORD   CONTENT                   UNIT                 FORMAT
C   ----   ----------------------    -------------------  ---------
C     1    LATITUDE                  DEGREES (N+,S-)      REAL
C     2    LONGITUDE                 DEGREES EAST         REAL
C     3    TIME SIGNIFICANCE (BUFR CODE TABLE "0 08 021") INTEGER
C     4    OBSERVATION TIME          HOURS (UTC)          REAL
C     5    YEAR/MONTH                YYYYMM               INTEGER
C     6    DAY/HOUR                  DDHH                 INTEGER
C     7    STATION ELEVATION         METERS               REAL
C     8    SUBMODE/EDITION NO.       (SM X 10) + ED. NO.  INTEGER
C                                    (ED. NO.=2, CONSTANT)
C     9    REPORT TYPE               71 (CONSTANT)        INTEGER
C    10    AVERAGING TIME            MINUTES              INTEGER
C                                 (NEGATIVE MEANS PRIOR TO OBS. TIME)
C    11    RECEIPT TIME              HOURS (UTC)          REAL
C    12    NOT USED                  ZEROED
C
C 13-34    ZEROED OUT - NOT USED                          INTEGER
C    35    CATEGORY 10, NO. LEVELS   COUNT                INTEGER
C    36    CATEGORY 10, DATA INDEX   COUNT                INTEGER
C    37    CATEGORY 11, NO. LEVELS   COUNT                INTEGER
C    38    CATEGORY 11, DATA INDEX   COUNT                INTEGER
C 39-42    ZEROED OUT - NOT USED                          INTEGER
C
C 43-END   UNPACKED DATA GROUPS      (FOLLOWS)            REAL
C
C   CATEGORY 10 - WIND PROFILER SFC DATA (EACH LEVEL, SEE WORD 35 ABOVE)
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    SEA-LEVEL PRESSURE   0.1 MILLIBARS       REAL
C(SEE *)2    STATION PRESSURE     0.1 MILLIBARS       REAL
C       3    HORIZ. WIND DIR.     DEGREES             REAL
C       4    HORIZ. WIND SPEED    0.1 METERS/SEC      REAL
C       5    AIR TEMPERATURE      0.1 DEGREES K       REAL
C       6    RELATIVE HUMIDITY    PERCENT             REAL
C       7    RAINFALL RATE        0.0000001 M/SEC     REAL
C
C   CATEGORY 11 - WIND PROFILER UPPER-AIR DATA (FIRST LEVEL IS SURFACE)
C                 (EACH LEVEL, SEE WORD 37 ABOVE)
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    HEIGHT ABOVE SEA-LVL METERS              REAL
C       2    HORIZ. WIND DIR.     DEGREES             REAL
C       3    HORIZ. WIND SPEED    0.1 METERS/SEC      REAL
C       4    QUALITY CODE         (SEE %)             INTEGER
C       5    VERT. WIND COMP. (W) 0.01 METERS/SEC     REAL
C       6    HORIZ. CONSENSUS NO. (SEE $)             INTEGER
C       7    VERT.  CONSENSUS NO. (SEE $)             INTEGER
C       8    SPECTRAL PEAK POWER  DB                  REAL
C       9    HORIZ. WIND SPEED    0.1 METERS/SEC      REAL
C            STANDARD DEVIATION
C      10    VERT. WIND COMPONENT 0.1 METERS/SEC      REAL
C            STANDARD DEVIATION
C      11    MODE                 (SEE #)             INTEGER
C
C  *-  ALWAYS MISSING
C  %-  0 - MEDIAN AND SHEAR CHECKS BOTH PASSED
C      2 - MEDIAN AND SHEAR CHECK RESULTS INCONCLUSIVE
C      4 - MEDIAN CHECK PASSED; SHEAR CHECK FAILED
C      8 - MEDIAN CHECK FAILED; SHEAR CHECK PASSED
C     12 - MEDIAN AND SHEAR CHECKS BOTH FAILED
C  $-  NO. OF INDIVIDUAL 6-MINUTE AVERAGE MEASUREMENTS THAT WERE
C      INCLUDED IN FINAL ESTIMATE OF AVERAGED WIND (RANGE: 0, 2-10)
C      (BASED ON A ONE-HOUR AVERAGE)
C  #-  1 - DATA FROM LOW MODE
C      2 - DATA FROM HIGH MODE
C      3 - MISSING
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C               FORMAT FOR GOES SOUNDING/RADIANCE REPORTS
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                           HEADER
C   WORD   CONTENT                   UNIT                 FORMAT
C   ----   ----------------------    -------------------  ---------
C     1    LATITUDE                  DEGREES (N+,S-)      REAL
C     2    LONGITUDE                 DEGREES EAST         REAL
C     3    FIELD OF VIEW NUMBER      NUMERIC              INTEGER
C     4    OBSERVATION TIME          HOURS (UTC)          REAL
C     5    YEAR/MONTH                YYYYMM               INTEGER
C     6    DAY/HOUR                  DDHH                 INTEGER
C     7    STATION ELEVATION         METERS               REAL
C     8    PROCESS. TECHNIQUE        (=21-CLEAR;          INTEGER
C                                   =23-CLOUD-CORRECTED)
C     9    REPORT TYPE               61 (CONSTANT)        INTEGER
C    10    QUALITY FLAG      (BUFR CODE TABLE "0 33 002") INTEGER
C    11    RECEIPT TIME              HOURS (UTC)          REAL
C    12    NOT USED                  ZEROED
C
C 13-26    ZEROED OUT - NOT USED
C    27    CATEGORY 08, NO. LEVELS   COUNT                INTEGER
C    28    CATEGORY 08, DATA INDEX   COUNT                INTEGER
C 29-38    ZEROED OUT - NOT USED
C    39    CATEGORY 12, NO. LEVELS   COUNT                INTEGER
C    40    CATEGORY 12, DATA INDEX   COUNT                INTEGER
C    41    CATEGORY 13, NO. LEVELS   COUNT                INTEGER
C    42    CATEGORY 13, DATA INDEX   COUNT                INTEGER
C
C 43-END   UNPACKED DATA GROUPS      (FOLLOWS)            REAL
C
C   CATEGORY 12 - SATELLITE SOUNDING LEVEL DATA (FIRST LEVEL IS SURFACE;
C                 EACH LEVEL, SEE 39 ABOVE)
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    PRESSURE             0.1 MILLIBARS       REAL
C       2    GEOPOTENTIAL         METERS              REAL
C       3    TEMPERATURE          0.1 DEGREES C       REAL
C       4    DEWPOINT TEMPERATURE 0.1 DEGREES C       REAL
C       5    NOT USED             SET TO MISSING      REAL
C       6    NOT USED             SET TO MISSING      REAL
C       7    Q.M. FOR GEOPOT      (SEE &)             REAL
C       8    Q.M. FOR TEMPERATURE (SEE &)             REAL
C       9    Q.M. FOR DEWPT TEMP  (SEE &)             REAL
C
C   CATEGORY 13 - SATELLITE RADIANCE "LEVEL" DATA (EACH "LEVEL", SEE
C                 41 ABOVE)
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    CHANNEL NUMBER       NUMERIC             INTEGER
C       2    BRIGHTNESS TEMP.     0.01 DEG. KELVIN    REAL
C       3    Q.M. FOR BTEMP       (SEE &)             REAL
C
C   CATEGORY 08 - ADDITIONAL (MISCELLANEOUS) DATA (EACH LEVEL, SEE @
C                 BELOW)
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    VARIABLE             (SEE @)             REAL
C       2    CODE FIGURE          (SEE @)             REAL
C       3    Q.M. FOR THE DATUM   (SEE &)             REAL
C       4    NOT USED                                 REAL
C
C  %-  SIXTH CHARACTER OF STATION ID (stnid) IS A TAGGED AS FOLLOWS:
C          "I" - GOES-EVEN-1 (252, 256, ...) SAT. , CLEAR COLUMN  RETR.
C          "J" - GOES-EVEN-1 (252, 256, ...) SAT. , CLD-CORRECTED RETR.

C          "L" - GOES-ODD-1  (253, 257, ...) SAT. , CLEAR COLUMN  RETR.
C          "M" - GOES-ODD-1  (253, 257, ...) SAT. , CLD-CORRECTED RETR.

C          "O" - GOES-EVEN-2 (254, 258, ...) SAT. , CLEAR COLUMN  RETR.
C          "P" - GOES-EVEN-2 (254, 258, ...) SAT. , CLD-CORRECTED RETR.

C          "Q" - GOES-ODD-2  (251, 255, ...) SAT. , CLEAR COLUMN  RETR.
C          "R" - GOES-ODD-2  (251, 255, ...) SAT. , CLD-CORRECTED RETR.

C          "?" - EITHER SATELLITE AND/OR RETRIEVAL TYPE UNKNOWN

C  &-  2.0 - INDICATES DATA NOT SUSPECT
C      3.0 - INDICATES DATA ARE SUSPECT
C     13.0 - INDICATES DATA ARE BAD
C  @-  NUMBER OF "LEVELS" FROM WORD 27.  MAXIMUM IS 12, AND ARE ORDERED
C      AS FOLLOWS (IF A DATUM ARE MISSING THAT LEVEL NOT STORED)
C           1 - LIFTED INDEX ---------- .01 DEG. KELVIN -- C. FIG. 250.
C           2 - TOTAL PRECIP. WATER  -- .01 MILLIMETERS -- C. FIG. 251.
C           3 - 1. TO .9 SIGMA P.WATER- .01 MILLIMETERS -- C. FIG. 252.
C           4 - .9 TO .7 SIGMA P.WATER- .01 MILLIMETERS -- C. FIG. 253.
C           5 - .7 TO .3 SIGMA P.WATER- .01 MILLIMETERS -- C. FIG. 254.
C           6 -  SKIN TEMPERATURE ----- .01 DEG. KELVIN -- C. FIG. 255.
C           7 -  CLOUD TOP TEMPERATURE- .01 DEG. KELVIN -- C. FIG. 256.
C           8 -  CLOUD TOP PRESSURE --- .1 MILLIBARS ----- C. FIG. 257.
C           9 -  CLOUD AMOUNT (BUFR TBL. C.T. 0-20-011) -- C. FIG. 258.
C          10 -  INSTR. DATA USED IN PROC.
C                             (BUFR TBL. C.T. 0-02-021) -- C. FIG. 259.
C          11 -  SOLAR ZENITH ANGLE --- .01 DEGREE ------- C. FIG. 260.
C          12 -  SAT. ZENITH ANGLE ---- .01 DEGREE ------- C. FIG. 261.
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                FORMAT FOR NEXRAD (VAD) WIND REPORTS
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C                           HEADER
C   WORD   CONTENT                   UNIT                 FORMAT
C   ----   ----------------------    -------------------  ---------
C     1    LATITUDE                  DEGREES              REAL
C     2    LONGITUDE                 DEGREES EAST         REAL
C     3    ** RESERVED **            SET TO 99999         INTEGER
C     4    OBSERVATION TIME          HOURS (UTC)          REAL
C     5    YEAR/MONTH                YYYYMM               INTEGER
C     6    DAY/HOUR                  DDHH                 INTEGER
C     7    STATION ELEVATION         METERS               REAL
C     8    ** RESERVED **            SET TO 99999         INTEGER
C
C     9    REPORT TYPE               72 (CONSTANT)        INTEGER
C    10    ** RESERVED **            SET TO 99999         INTEGER
C    11    RECEIPT TIME              HOURS (UTC)          REAL
C    12    NOT USED                  ZEROED
C
C 13-18    ZEROED OUT - NOT USED                          INTEGER
C    19    CATEGORY 04, NO. LEVELS   COUNT                INTEGER
C    20    CATEGORY 04, DATA INDEX   COUNT                INTEGER
C 21-42    ZEROED OUT - NOT USED                          INTEGER
C
C 43-END   UNPACKED DATA GROUPS      (FOLLOWS)            REAL
C
C   CATEGORY 04 - UPPER-AIR WINDS-BY-HEIGHT DATA(FIRST LEVEL IS SURFACE)
C                 (EACH LEVEL, SEE WORD 19 ABOVE)
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    HEIGHT ABOVE SEA-LVL METERS              REAL
C       2    HORIZ. WIND DIR.     DEGREES             REAL
C       3    HORIZ. WIND SPEED    0.1 METERS/SEC      REAL
C       4    Q.M. FOR HEIGHT      (HARDWIRED = 2.0)   REAL
C       5    Q.M. FOR WIND        (SEE %)             REAL
C
C  %-  "CONFIDENCE LEVEL" WHICH IS RELATED TO THE ROOT-MEAN-SQUARE
C       VECTOR ERROR FOR THE HORIZONTAL WIND.  IT IS DEFINED AS
C       FOLLOWS:
C                1.0  = RMS OF  1.9 KNOTS
C                2.0  = RMS OF  3.9 KNOTS
C                3.0  = RMS OF  5.8 KNOTS
C                4.0  = RMS OF  7.8 KNOTS (DEFAULT)
C                5.0  = RMS OF  9.7 KNOTS
C                6.0  = RMS OF 11.7 KNOTS
C                7.0  = RMS  > 13.6 KNOTS
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C   FOR ALL REPORT TYPES, MISSING VALUES ARE:
C                       99999. FOR REAL
C                       99999  FOR INTEGER (except for words 5 and 6,
C                               where missing is 999999)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM-SP, CRAY, SGI
C
C$$$
      SUBROUTINE W3UNPKB7(IDATE,IHE,IHL,LUNIT,RDATA,STNID,DSNAME,IDSDAT,
     $ IDSDMP,IRET)
      CHARACTER*4  CBUFR
      INTEGER  IDATE(4),LSDATE(4),IDATA(720)
      INTEGER JDATE(8)
      INTEGER(8) IDSDMP
      CHARACTER*8  STNID,SUBSET,DSNAME
      DIMENSION  RINC(5)
      REAL  RDATA(*),RDATX(720)
      COMMON /PKB7BB/KDATE(8),LDATE(8),IPRINT
      COMMON /PKB7CC/INDEX
      COMMON /PKB7DD/LSHE,LSHL,ICDATE(5),IDDATE(5)
      COMMON /PKB7FF/IFOV(3),KNTSAT(250:260)
 
      SAVE
 
      EQUIVALENCE (RDATX,IDATA)
      DATA ITM/0/,LUNITL/-99/,KOUNT/0/
      IPRINT = 0
      IF(IRET.LT.0)  IPRINT = IABS(IRET)
      write(*,*) 'Sorry cannot run this - missing a few routines.'
      IRET =  2
      return
      END