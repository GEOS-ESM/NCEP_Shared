C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    IW3UNPBF    UNPACKS REPORT FROM JBUFR INTO SPEC. FMT
C   PRGMMR: KEYSER           ORG: NP22       DATE: 1999-01-20
C
C ABSTRACT: READS AND UNPACKS ONE REPORT FROM INPUT JBUFR FILE INTO
C   SPECIFIED FORMAT.  FUNCTION RETURNS THE UNPACKED REPORT TO THE
C   CALLING PROGRAM IN THE ARRAY 'OBS' AS WELL AS IN VARIABLES 'STNID',
C   'CRES1', AND 'CRES2'.  IT ALSO RETURNS INFORMATION ABOUT THE
C   INPUT DATA SET ITSELF (NAME, CENTER DATE, DUMP TIME).  VARIOUS
C   CONTINGENCIES ARE COVERED BY RETURN VALUE OF THE FUNCTION AND
C   PARAMETER 'IER' - FUNCTION AND IER HAVE SAME VALUE.  REPEATED
C   CALLS OF FUNCTION WILL RETURN A SEQUENCE OF UNPACKED REPORTS.
C   THE CALLING PROGRAM MAY SWITCH TO A NEW 'NUNIT' AT ANY TIME,
C   THAT DATASET WILL THEN BE READ IN SEQUENCE.  IF USER SWITCHES
C   BACK TO A PREVIOUS 'NUNIT', THAT DATA SET WILL BE READ FROM THE
C   BEGINNING, NOT FROM WHERE THE USER LEFT OFF (THIS IS A 'SOFTWARE
C   TOOL', NOT AN ENTIRE I/O SYSTEM).
C
C PROGRAM HISTORY LOG:
C 1998-02-17  D. A. KEYSER -- ORIGINAL AUTHOR (ADAPTED FROM W3LIB
C             ROUTINE IW3GAD)
C 1998-03-16  D. A. KEYSER -- BUMPED UP ARRAY SIZE FOR UNPACKED
C             REPORT IN "OBS" FROM 1608 TO 2500; SUBROUTINE NO
C             LONGER RETURNS WITH IER=999 IF UNPACKED ARRAY SIZE
C             IS EXCEEDED - RATHER IT SKIPS PROCESSING OF OFFENDING
C             REPORT AND PRINTS A DIAGNOSTIC; BUMP UP MAXIMUM
C             NUMBER OF LEVELS THAT CAN BE PROCESSED FOR A CATEGORY
C             FROM 150 TO 160
C 1998-06-15  D. A. KEYSER -- STREAMLINED THE PROCESSING OF FLIGHT-
C             LEVEL RECONNAISSANCE REPORTS; REDEFINED UNITS FOR
C             UNPACKED WORDS 1 (LATITUDE), 2 (LONGITUDE), 4 (OBS.
C             TIME) AND 11 (RECEIPT TIME) - ALL TO STREAMLINE
C             PROCESSING IN PREPDATA PROGRAM; IN ADPUPA PROCESSING,
C             REMOVED WRITING OF CAT. 8 C.F. 106, ADDED WRITING OF
C             NEW CAT. 8 C.F. 353 (SOLAR AND INFRARED RADIATION
C             CORRECTION INDICATOR) AND 354 (TRACKING TECHNIQUE/
C             STATUS OF SYSTEM USED INDICATOR)
C 1998-09-21  D. A. KEYSER -- SUBROUTINE NOW Y2K AND FORTRAN 90
C             COMPLIANT
C 1999-01-20 D. A. KEYSER -- INCORPORATED BOB KISTLER'S CHANGES NEEDED
C             TO PORT THE CODE TO THE IBM SP; ADDED REPORT SUBTYPE 1
C             FOR RECCOS UNDER REPORT TYPE 31 AND SUBTYPE 2 FOR DROPS
C             UNDER REPORT TYPE 31
C
C
C USAGE:    II = IW3UNPBF(NUNIT, OBS, STNID, CRES1, CRES2, 
C                         DSNAME, IDSDAT, IDSDMP, IER)
C   INPUT ARGUMENT LIST:
C     NUNIT    - FORTRAN UNIT NUMBER FOR SEQUENTIAL DATA SET CONTAINING
C              - PACKED JBUFR REPORTS
C
C   OUTPUT ARGUMENT LIST:
C     OBS      - ARRAY CONTAINING ONE REPORT IN SPECIFIED UNPACKED
C              - FORMAT.  FORMAT IS MIXED, USER MUST EQUIVALENCE
C              - INTEGER ARRAY TO THIS ARRAY (SEE REMARKS BELOW FOR
C              - UNPACKED FORMAT) THE LENGTH OF THE ARRAY SHOULD BE AT
C              - LEAST 2500 - (NOTE: DOES NOT INCLUDE STATION ID,
C              - CHARACTER RESERVE WORD 1, AND CHARACTER RESERVE WORD 2)
C     STNID    - CHARACTER*8 SINGLE REPORT STATION IDENTIFICATION (UP
C              - TO 8 CHARACTERS, LEFT-JUSTIFIED)
C     CRES1    - CHARACTER*8 SINGLE REPORT CHARACTER RESERVE WORD 1
C              - (SEE DOCUMENTATION/COMMENTS IN THIS PROGRAM (VARIES BY
C              - TYPE OF DATA)
C     CRES2    - CHARACTER*8 SINGLE REPORT CHARACTER RESERVE WORD 2
C              - (SEE DOCUMENTATION/COMMENTS IN THIS PROGRAM (VARIES BY
C              - TYPE OF DATA)
C     DSNAME   - CHARACTER*8 DATA SET NAME (SAME FOR ALL REPORTS IN
C              - A COMMON INPUT DATA SET - SEE REMARKS FOR IER=1)
C     IDSDAT   - INTEGER DATA SET CENTER DATE IN FORM YYYYMMDDHH (SAME
C              - FOR ALL REPORTS IN A COMMON INPUT DATA SET - SEE
C              - REMARKS FOR IER=1)
C     IDSDMP   - INTEGER DATA SET DUMP TIME IN FORM YYYYMMDDHHMM (SAME
C              - FOR ALL REPORTS IN A COMMON INPUT DATA SET - SEE
C              - REMARKS FOR IER=1)
C     IER      - RETURN FLAG (EQUAL TO FUNCTION VALUE) - SEE REMARKS
C
C   INPUT FILES:
C     UNIT AA  - SEQUENTIAL JBUFR DATA SET ("AA" IS UNIT NUMBER
C              - SPECIFIED BY INPUT ARGUMENT "NUNIT")
C
C   OUTPUT FILES:
C     UNIT 06  - PRINTOUT
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:    - R01UBF   R02UBF   R03UBF   R04UBF   R05UBF
C                - R06UBF   I02UBF   I04UBF   I05UBF   L01UBF
C                - L02UBF   L03UBF   S01UBF   S02UBF   S03UBF
C                - S04UBF   S05UBF   SE01UBF  EQSUBF   EQMUBF
C                - EWZUBF   ERTUBF   C01UBF
C     LIBRARY:
C       FORTRAN  - ORDERS (CRAY PLATFORM ONLY)
C       W3LIB    - ORDERS (SGI  PLATFORM ONLY)
C       BUFRLIB  - OPENBF   CLOSBF   DATELEN  STATUS   DUMPBF
C                - READMG   READNS   UFBINT
C
C REMARKS:
C     INPUT DATA SET SHOULD BE ASSIGNED IN THIS WAY:
C       Cray:
C        assign -a ADPUPA                       fort.XX
C       SGI:
C        assign -a ADPUPA     -F cos            fort.XX 
C
C     THE RETURN FLAGS IN IER (AND FUNCTION IW3UNPBF ITSELF) ARE:
C          =   0  OBSERVATION READ AND UNPACKED INTO LOCATIONS 'OBS',
C                   'STNID', 'CRES1', 'CRES2', 'DSNAME', 'IDSDAT', AND
C                   'IDSDMP'.  SEE REMARKS BELOW FOR CONTENTS. NEXT
C                   CALL TO IW3UNPBF WILL RETURN NEXT OBSERVATION IN
C                   DATA SET.
C          =   1  INFORMATION ABOUT THE BUFR DATASET IS RETURNED IN
C                   THE OUTPUT ARGUMENTS DSNAME, IDSDAT, IDSDMP (SEE
C                   OUTPUT ARGUMENT LIST ABOVE)
C
C                   THIS SHOULD ALWAYS OCCUR AFTER THE FIRST CALL TO
C                   THIS SUBROUTINE.  NO REPORT IS UNPACKED AT THIS
C                   POINT, AND ONLY DSNAME, IDSDAT, AND IDSDMP
C                   CONTAIN INFORMATION.  ALL SUBSEQUENT CALLS TO
C                   IW3UNPBF SHOULD RETURN THE OBSERVATIONS IN THIS
C                   DATA SET, SEQUENTIALLY, (IER=0) UNTIL THE END OF
C                   FILE IS ENCOUNTERED (IER=2 OR 3).  THE VALUES
C                   STORED IN DSNAME, IDSDAT, AND IDSDMP WILL
C                   CONTINUE TO BE RETURNED ALONG WITH EACH REPORT
C                   WHEN IER = 0.
C          =   2  THE PHYSICAL END OF FILE IS ENCOUNTERED (THIS IS
C                   FOR A VALID JBUFR FILE CONTAINING REPORTS).  ALL
C                   DONE.
C          =   3  THIS IS AN EMPTY (NULL) FILE.  ALL DONE.
C          = 999  ERROR: EITHER FILE IS NOT JBUFR, ERROR READING FILE,
C                        CENTER DATE DUMMY MESSAGE NOT FOUND AT TOP OF
C                        JBUFR FILE, OR SOME OTHER PROBLEM IN DECODING
C                        ONE OR MORE REPORTS IN A JBUFR FILE.
C                  NO USEFUL INFORMATION IN 'OBS', 'STNID', 'CRES1',
C                  'CRES2', 'DSNAME', 'IDSDAT', AND 'IDSDMP' ARRAYS.
C                  CALLING PROGRAM CAN CHOOSE TO STOP WITH NON-ZERO
C                  CONDITION CODE OR RESET 'NUNIT' TO POINT TO A NEW
C                  DATA SET (IN WHICH CASE NEXT CALL TO IW3UNPBF
C                  SHOULD RETURN WITH IER=1).
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C           CONTENTS OF AN UNPACKED REPORT IN THE "OBS" ARRAY
C                  (MISSING INTEGER DATA ARE SET TO 99999;
C                    MISSING REAL DATA ARE SET TO 99999.)
C
C     (NOTE: DOES NOT INCLUDE STATION IDENTIFICATION, CHARACTER RESERVE
C            WORD 1, AND CHARACTER RESERVE WORD 2 - SEE OUTPUT ARGUMENTS
C            "STNID", "CRES1", AND "CRES2" ABOVE)
C
C   ***************************************************************
C   WORD   CONTENT                   UNIT                 FORMAT
C   ----   ----------------------    -------------------  ---------
C     1    LATITUDE                  DEGREES (N+,S-)      REAL
C     2    LONGITUDE                 DEGREES EAST         REAL
C     3    NOT USED                  MISSING              REAL
C     4    OBSERVATION TIME          HOURS (UTC)          REAL
C     5    NOT USED                  MISSING              REAL
C     6    NOT USED                  MISSING              REAL
C     7    STATION ELEVATION         METERS               REAL
C     8    INSTRUMENT TYPE           ON29 TABLE R.2       INTEGER
C     9    REPORT TYPE               ON29 TABLE R.1 OR    INTEGER
C                                    ON124 TABLE S.3
C    10    REPORT SUBTYPE            (SEE %)              INTEGER
C    11    RECEIPT TIME              HOURS (UTC)          REAL
C    12    NOT USED                  MISSING              REAL
C
C    13    CATEGORY  1, NO. LEVELS   COUNT                INTEGER
C    14    CATEGORY  1, DATA INDEX   COUNT                INTEGER
C    15    CATEGORY  2, NO. LEVELS   COUNT                INTEGER
C    16    CATEGORY  2, DATA INDEX   COUNT                INTEGER
C    17    CATEGORY  3, NO. LEVELS   COUNT                INTEGER
C    18    CATEGORY  3, DATA INDEX   COUNT                INTEGER
C    19    CATEGORY  4, NO. LEVELS   COUNT                INTEGER
C    20    CATEGORY  4, DATA INDEX   COUNT                INTEGER
C    21    CATEGORY  5, NO. LEVELS   COUNT                INTEGER
C    22    CATEGORY  5, DATA INDEX   COUNT                INTEGER
C    23    CATEGORY  6, NO. LEVELS   COUNT                INTEGER
C    24    CATEGORY  6, DATA INDEX   COUNT                INTEGER
C    25    CATEGORY  7, NO. LEVELS   COUNT                INTEGER
C    26    CATEGORY  7, DATA INDEX   COUNT                INTEGER
C    27    CATEGORY  8, NO. LEVELS   COUNT                INTEGER
C    28    CATEGORY  8, DATA INDEX   COUNT                INTEGER
C    29    CATEGORY 51, NO. LEVELS   COUNT                INTEGER
C    30    CATEGORY 51, DATA INDEX   COUNT                INTEGER
C    31    CATEGORY 52, NO. LEVELS   COUNT                INTEGER
C    32    CATEGORY 52, DATA INDEX   COUNT                INTEGER
C    33    CATEGORY  9, NO. LEVELS   COUNT                INTEGER
C    34    CATEGORY  9, DATA INDEX   COUNT                INTEGER
C 35-42    ZEROED OUT - NOT USED                          INTEGER
C
C 43-END   UNPACKED DATA GROUPS      (SEE BELOW)          MIXED
C
C    % - WORD 10, REPORT SUBTYPE IS CURRENTLY SET TO MISSING FOR
C         ALL TYPES EXCEPT:
C             AIRCRAFT (REPORT TYPE 41), WHERE:
C                 = 1 - AIREP AIRCRAFT
C                 = 2 - PIREP AIRCRAFT
C                 = 3 - AMDAR/ASDAR AIRCRAFT
C                 = 4 - ACARS AIRCRAFT
C             RECONNAISSANCE/DROPWINSONDE (REPORT TYPE 31), WHERE:
C                 = 1 - RECONNAISSANCE AIRCRAFT
C                 = 2 - DROPWINSONDE
C
C
C
C   ***************************************************************
C
C
C     DATA ARE UNPACKED INTO FIXED LOCATIONS IN WORDS 1-12 AND INTO
C     INDEXED LOCATIONS IN WORD 43 AND FOLLOWING.  THE VERTICAL
C     SIGNIFICANCE DESCRIPTOR FOR EACH LEVEL IN THE BUFR REPORT IS
C     USED TO UNPACK THE LEVEL INTO THE "CATEGORIES" DESCRIBED BELOW.
C     EACH CATEGORY HAS A LAYOUT IN LOCATIONS IN ARRAY OBS THAT MAY
C     BE FOUND BY USING THE CORRESPONDING INDEX AMOUNT FROM WORDS 14,
C     16, ..., 34, IN ARRAY OBS.  FOR INSTANCE, IF A REPORT IS
C     UNPACKED INTO ONE OR MORE CATEGORY 3 DATA GROUPS (WIND DATA AT
C     VARIABLE PRESSURE LEVELS) THAT DATA WILL BE SPECIFIED IN THE
C     UNPACKED BINARY FORMAT AS DESCRIBED BELOW UNDER CATEGORY 3.
C     THE NUMBER OF LEVELS WILL BE STORED IN WORD 17 OF OBS AND THE
C     INDEX OF THE FIRST LEVEL OF UNPACKED DATA IN THE OUTPUT ARRAY
C     WILL BE STORED IN WORD 18.  THE SECOND LEVEL, IF ANY, WILL BE
C     STORED BEGINNING SIX WORDS FURTHER ON, AND SO FORTH UNTIL THE
C     COUNT IN WORD 17 IS EXHAUSTED.  THE FIELD LAYOUT IN EACH
C     CATEGORY IS GIVEN BELOW...
C
C     CATEGORY 1 - MANDATORY LEVEL DATA
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    PRESSURE             0.1 MILLIBARS       REAL
C       2    GEOPOTENTIAL         METERS              REAL
C       3    TEMPERATURE          0.1 DEGREES C       REAL
C       4    DEWPOINT DEPRESSION  0.1 DEGREES C       REAL
C       5    WIND DIRECTION       DEGREES             REAL
C       6    WIND SPEED           0.1 METERS/SEC      REAL
C       7    PRES. QUALITY MARKER (SEE $)             REAL
C       8    GEOP. QUALITY MARKER (SEE $)             REAL
C       9    TEMP. QUALITY MARKER (SEE $)             REAL
C      10    DDPR. QUALITY MARKER (SEE $)             REAL
C      11    WIND  QUALITY MARKER (SEE $)             REAL
C
C     CATEGORY 2 - TEMPERATURE AT VARIABLE PRESSURE
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    PRESSURE             0.1 MILLIBARS       REAL
C       2    TEMPERATURE          0.1 DEGREES C       REAL
C       3    DEWPOINT DEPRESSION  0.1 DEGREES C       REAL
C       4    PRES. QUALITY MARKER (SEE $)             REAL
C       5    TEMP. QUALITY MARKER (SEE $)             REAL
C       6    DDPR. QUALITY MARKER (SEE $)             REAL
C       7    SPECIAL INDICATOR    (SEE $$)            REAL
C
C     CATEGORY 3 - WINDS AT VARIABLE PRESSURE
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    PRESSURE             0.1 MILLIBARS       REAL
C       2    WIND DIRECTION       DEGREES             REAL
C       3    WIND SPEED           0.1 METERS/SEC      REAL
C       4    PRES. QUALITY MARKER (SEE $)             REAL
C       5    WIND  QUALITY MARKER (SEE $)             REAL
C       6    SPECIAL INDICATOR    (SEE $$$)           REAL
C
C     CATEGORY 4 - WINDS AT VARIABLE HEIGHTS
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    GEOPOTENTIAL         METERS              REAL
C       2    WIND DIRECTION       DEGREES             REAL
C       3    WIND SPEED           0.1 METERS/SEC      REAL
C       4    GEOP. QUALITY MARKER (SEE $)             REAL
C       5    WIND  QUALITY MARKER (SEE $)             REAL
C
C     CATEGORY 5 - TROPOPAUSE DATA
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    PRESSURE             0.1 MILLIBARS       REAL
C       2    TEMPERATURE          0.1 DEGREES C       REAL
C       3    DEWPOINT DEPRESSION  0.1 DEGREES C       REAL
C       4    WIND DIRECTION       DEGREES             REAL
C       5    WIND SPEED           0.1 METERS/SEC      REAL
C       6    PRES. QUALITY MARKER (SEE $)             REAL
C       7    TEMP. QUALITY MARKER (SEE $)             REAL
C       8    DDPR. QUALITY MARKER (SEE $)             REAL
C       9    WIND  QUALITY MARKER (SEE $)             REAL
C
C     CATEGORY 6 - CONSTANT-LEVEL DATA (AIRCRAFT, SAT. CLOUD-DRIFT)
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    PRESSURE             0.1 MILLIBARS       REAL
C       2    PRESSURE ALTITUDE    METERS              REAL
C       3    TEMPERATURE          0.1 DEGREES C       REAL
C       4    DEWPOINT DEPRESSION  0.1 DEGREES C       REAL
C       5    WIND DIRECTION       DEGREES             REAL
C       6    WIND SPEED           0.1 METERS/SEC      REAL
C       7    PRES. QUALITY MARKER (SEE $)             REAL
C       8    P-ALT QUALITY MARKER (SEE $)             REAL
C       9    TEMP. QUALITY MARKER (SEE $)             REAL
C      10    DDPR. QUALITY MARKER (SEE $)             REAL
C      11    WIND  QUALITY MARKER (SEE $)             REAL
C
C     CATEGORY 7 - CLOUD COVER
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    PRESSURE             0.1 MILLIBARS       REAL
C       2    AMOUNT OF CLOUDS     PER CENT            REAL
C       3    PRES. QUALITY MARKER (SEE $)             REAL
C       4    C-AMT QUALITY MARKER (SEE $)             REAL
C
C     CATEGORY 8 - ADDITIONAL DATA
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    SPECIFIED IN ON29    VARIABLE            REAL
C            TABLE 101.1 OR
C            ON124 TABLE SM.8A.1
C       2    FORM OF ADD'L DATA   CODE FIGURE FROM    REAL
C                                 ON29 TABLE 101 OR
C                                 ON124 TABLE SM.8A
C       3    INDICATOR 1          (SEE @)             REAL    
C       4    INDICATOR 2          (SEE @)             REAL    
C
C     CATEGORY 51 - SURFACE DATA
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    SEA-LEVEL PRESSURE   0.1 MILLIBARS       REAL
C       2    STATION PRESSURE     0.1 MILLIBARS       REAL
C       3    WIND DIRECTION       DEGREES             REAL
C       4    WIND SPEED           0.1 METERS/SEC      REAL
C       5    AIR TEMPERATURE      0.1 DEGREES C       REAL
C       6    DEWPOINT DEPRESSION  0.1 DEGREES C       REAL
C       7    MAXIMUM TEMPERATURE  0.1 DEGREES C       REAL
C       8    MINIMUM TEMPERATURE  0.1 DEGREES C       REAL
C       9    MSL-P QUALITY MARKER (SEE $)             REAL
C      10    STN-P QUALITY MARKER (SEE $)             REAL
C      11    WIND  QUALITY MARKER (SEE $)             REAL
C      12    ATEMP QUALITY MARKER (SEE $)             REAL
C      13    DDPR. QUALITY MARKER (SEE $)             REAL
C      14    HORIZ. VISIBILITY    WMO CODE TABLE 4300 INTEGER
C      15    PRESENT WEATHER      WMO CODE TABLE 4677 INTEGER
C      16    PAST WEATHER         WMO CODE TABLE 4561 INTEGER
C      17    PAST WEATHER 2       WMO CODE TABLE ???? INTEGER
C      18    TOTAL CLOUD COVER N  WMO CODE TABLE 2700 INTEGER
C      19    CLOUD COVER OF C/LN  WMO CODE TABLE 2700 INTEGER
C      20    CLOUD TYPE OF C/L    WMO CODE TABLE 0513 INTEGER
C      21    CLOUD HEIGHT OF C/L  WMO CODE TABLE 1600 INTEGER
C      22    CLOUD TYPE OF C/M    WMO CODE TABLE 0515 INTEGER
C      23    CLOUD TYPE OF C/H    WMO CODE TABLE 0509 INTEGER
C      24    CHARACTERISTIC OF    WMO CODE TABLE 0200 INTEGER
C            3-HR PRESS TENDENCY
C      25    AMT. PRESS TENDENCY  0.1 MILLIBARS       REAL
C            (50.0 WILL BE ADDED TO INDICATE 24-HR TENDENCY)
C
C     CATEGORY 52 - ADDITIONAL SURFACE DATA
C     WORD   PARAMETER            UNITS               FORMAT
C     ----   ---------            -----------------   -------------
C       1    6-HR PRECIPITATION   0.01 INCH           INTEGER
C       2    SNOW DEPTH           INCH                INTEGER
C       3    24-HR PRECIPITATION  0.01 INCH           INTEGER
C       4    DURATION OF PRECIP.  NO. 6-HR PERIODS    INTEGER
C       5    PERIOD OF WAVES      SECONDS             INTEGER
C       6    HEIGHT OF WAVES      0.5 METERS          INTEGER
C       7    SWELL DIRECTION      WMO CODE TABLE 0877 INTEGER
C       8    SWELL PERIOD         SECONDS             INTEGER
C       9    SWELL HEIGHT         0.5 METERS          INTEGER
C      10    SEA SFC TEMPERATURE  0.1 DEGREES C       INTEGER
C   **-11    SPECIAL PHEN, GEN'L                      INTEGER
C   **-12    SPECIAL PHEN, DET'L                      INTEGER
C      13    SHIP'S COURSE        WMO CODE TABLE 0700 INTEGER
C      14    SHIP'S AVERAGE SPEED WMO CODE TABLE 4451 INTEGER
C   **-15    WATER EQUIVALENT OF  0.01 INCH           INTEGER
C            SNOW AND/OR ICE
C
C     CATEGORY 9 - PLAIN LANGUAGE DATA (ALPHANUMERIC TEXT)
C       ==> CURRENTLY NOT POSSIBLE TO UNPACK IN THIS FUNCTION
C
C
C   **-CURRENTLY NOT POSSIBLE TO UNPACK THIS PARAMETER
C
C    $ - QUALITY MARKER CODE TABLE:
C             0. - MONITOR KEEP
C             1. - GOOD
C             2. - NEUTRAL OR NOT CHECKED (DEFAULT)
C             3. - SUSPECT
C             4. - GOOD - CORRECTED BY O.P.C. (SURFACE MARINE ONLY)
C            12. - REJECT LIST, DO NOT USE
C            13. - FAILED AUTOMATIC Q.C. TESTS, DO NOT USE
C            14. - MONITOR PURGE, DO NOT USE
C
C    $$  - CATEGORY 2 SPECIAL INDICATOR:
C             0. - NOTHING INDICATED (DEFAULT)
C             4. - LEVEL PRESSURE ESTIMATED
C
C    $$$ - CATEGORY 3 SPECIAL INDICATOR:
C             0. - NOTHING INDICATED (DEFAULT)
C             1. - TROPOPAUSE LEVEL
C             2. - MAXIMUM WIND LEVEL
C             3. - MAXIMUM WIND LEVEL AT TERMINATING LEVEL
C             4. - LEVEL PRESSURE ESTIMATED
C
C    @ - CATEGORY 8 INDICATORS 1 AND 2:
C             0. - NOTHING INDICATED (DEFAULT)
C           ELSE - SEE DOCUMENTATION IN THIS CODE FOR VARIOUS DATA
C                  DATA TYPE CATEGORY 8 PROCESSING
C
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM-SP, CRAY, SGI
C
C$$$
      FUNCTION IW3UNPBF(LUNIT,OBS,STNID,CRES1,CRES2,DSNAME,IDSDAT,
     $ IDSDMP,IER)
 
      COMMON/IUBFBB/KNDX,KSKACF(8),KSKUPA,KSKSFC,KSKSAT
      COMMON/IUBFCC/SUBSET
      COMMON/IUBFDD/HDR(12),RCATS(50,160,11),IKAT(11),MCAT(11),NCAT(11)
      COMMON/IUBFEE/ROBS(255,11)
      COMMON/IUBFFF/QMS(255,7)
      COMMON/IUBFGG/SFO(35)
      COMMON/IUBFHH/SFQ(5)
      COMMON/IUBFII/PWMIN
      COMMON/IUBFJJ/ISET,MANLIN(1001)
      COMMON/IUBFKK/KOUNT(8,9),KNTSAT(250:260),IFLSAT
      COMMON/IUBFLL/Q8(255,2)
      COMMON/IUBFMM/XIND(255)
      COMMON/IUBFNN/STNIDX,CRES1X,CRES2X
      COMMON/IUBFOO/DSNAMX,IDSDAX,IDSDMX
 
      DIMENSION    OBS(*),JWFILE(100)

      CHARACTER*8  STNID,STNIDX,CRES1,CRES1X,CRES2,CRES2X,DSNAME,DSNAMX,
     $ SUBSET

      INTEGER(8) IDSDMP,IDSDAX,IDSDMX

      SAVE

      DATA JWFILE/100*0/,LASTF/0/,ITIMES/0/

      IF(ITIMES.EQ.0)  THEN

C  THE FIRST TIME IN, INITIALIZE SOME DATA
C  (NOTE: FORTRAN 77/90 STANDARD DOES NOT ALLOW COMMON BLOCK VARIABLES
C         TO BE INITIALIZED VIA DATA STATEMENTS, AND, FOR SOME REASON,
C         THE BLOCK DATA DOES NOT INITIALIZE DATA IN THE W3LIB
C         A V O I D   B L O C K   D A T A   I N   W 3 L I B )
C  --------------------------------------------------------------------

         ITIMES = 1
         KNDX = 0
         KSKACF = 0
         KSKUPA = 0
         KSKSFC = 0
         KSKSAT = 0
         KOUNT =  0
         KNTSAT = 0
         IFLSAT = 0
         IKAT(1)  =  1
         IKAT(2)  =  2
         IKAT(3)  =  3
         IKAT(4)  =  4
         IKAT(5)  =  5
         IKAT(6)  =  6
         IKAT(7)  =  7
         IKAT(8)  =  8
         IKAT(9)  = 51
         IKAT(10) = 52
         IKAT(11) =  9

C  MCAT defines the number of parameters in a level for each category
C --> THIS NEEDS TO BE UPDATED WHEN ADDING MORE WORDS PER CAT LEVEL

         MCAT(1)  = 11
         MCAT(2)  =  7
         MCAT(3)  =  6
         MCAT(4)  =  5
         MCAT(5)  =  9
         MCAT(6)  = 11
         MCAT(7)  =  4
         MCAT(8)  =  4
         MCAT(9)  = 25
         MCAT(10) = 15
         MCAT(11) =  3

         ISET  = 0
      END IF

C  UNIT NUMBER OUT OF RANGE RETURNS A 999
C  --------------------------------------

      IF(LUNIT.LT.1 .OR. LUNIT.GT.100)  THEN
         PRINT *, '##IW3UNPBF - UNIT NUMBER ',LUNIT,' OUT OF RANGE -- ',
     $    'IER = 999'
         GO TO 9999
      END IF
      IF(LASTF.NE.LUNIT .AND. LASTF.GT.0) THEN
         CALL CLOSBF(LASTF)
         JWFILE(LASTF) = 0
      END IF
      LASTF = LUNIT
 
C  THE JWFILE INDICATOR: =0 IF UNOPENED; =2 IF OPENED AND JBUFR
C  ------------------------------------------------------------
 
      IF(JWFILE(LUNIT).EQ.0) THEN
         PRINT *,'===> IW3UNPBF - Y2K/F90 VERSION: 01-20-1999'
         IERR = I02UBF(LUNIT,OBS,IER)
         IF(IERR.EQ.1) THEN
            PRINT *, 'IW3UNPBF - OPENED A JBUFR FILE IN UNIT ',LUNIT
            JWFILE(LUNIT) = 2
            KNDX = 0
            KSKACF = 0
            KSKUPA = 0
            KSKSFC = 0
            KSKSAT = 0
            IFLSAT = 0
            IER = 1
            IW3UNPBF = 1
         ELSE  IF(IERR.NE.999) then
            IER = 3
            IW3UNPBF = 3                                
         ELSE
            IER = 999
            IW3UNPBF = 999
         END IF
      ELSEIF(JWFILE(LUNIT).EQ.2) THEN
         IF(I02UBF(LUNIT,OBS,IER).NE.0) JWFILE(LUNIT) = 0
         IF(IER.GT.0) CALL CLOSBF(LUNIT)
         IF(IER.EQ.2.OR.IER.EQ.3)  THEN
            IF(KSKACF(1).GT.0)  PRINT *, 'IW3UNPBF - NO. OF AIRCFT/',
     $      'AIRCAR REPORTS TOSSED DUE TO ZERO CAT. 6 LVLS = ',KSKACF(1)
            IF(KSKACF(2).GT.0)  PRINT *, 'IW3UNPBF - NO. OF AIRCFT ',
     $       'REPORTS TOSSED DUE TO BEING "LFPW" AMDAR = ',KSKACF(2)
            IF(KSKACF(8).GT.0)  PRINT *, 'IW3UNPBF - NO. OF AIRCFT ',
     $       'REPORTS TOSSED DUE TO BEING "PHWR" AIREP = ',KSKACF(8)
            IF(KSKACF(3).GT.0)  PRINT *, 'IW3UNPBF - NO. OF AIRCFT ',
     $       'REPORTS TOSSED DUE TO BEING CARSWELL AMDAR = ',KSKACF(3)
            IF(KSKACF(4).GT.0)  PRINT *, 'IW3UNPBF - NO. OF AIRCFT ',
     $       'REPORTS TOSSED DUE TO BEING CARSWELL ACARS = ',KSKACF(4)
            IF(KSKACF(5).GT.0)  PRINT *, 'IW3UNPBF - NO. OF AIRCFT/',
     $       'AIRCAR REPORTS TOSSED DUE TO HAVING MISSING WIND = ',
     $       KSKACF(5)
            IF(KSKACF(6).GT.0)  PRINT *, 'IW3UNPBF - NO. OF AIRCFT ',
     $       'REPORTS TOSSED DUE TO BEING AMDAR < 2286 M = ',KSKACF(6)
            IF(KSKACF(7).GT.0)  PRINT *, 'IW3UNPBF - NO. OF AIRCFT ',
     $       'REPORTS TOSSED DUE TO BEING AIREP <  100 M = ',KSKACF(7)
            IF(KSKACF(1)+KSKACF(2)+KSKACF(3)+KSKACF(4)+KSKACF(5)+
     $       KSKACF(6)+KSKACF(7)+KSKACF(8).GT.0)
     $       PRINT *, 'IW3UNPBF - TOTAL NO. OF AIRCFT/AIRCAR REPORTS ',
     $        'TOSSED = ',KSKACF(1)+KSKACF(2)+KSKACF(3)+KSKACF(4)+
     $        KSKACF(5)+KSKACF(6)+KSKACF(7)+KSKACF(8)
            IF(KSKUPA.GT.0)  PRINT *, 'IW3UNPBF - TOTAL NO. OF ADPUPA ',
     $       'REPORTS TOSSED = ',KSKUPA
            IF(KSKSFC.GT.0)  PRINT *, 'IW3UNPBF - TOTAL NO. OF ADPSFC/',
     $       'SFCSHP REPORTS TOSSED = ',KSKSFC
            IF(KSKSAT.GT.0)  PRINT *, 'IW3UNPBF - TOTAL NO. OF SATWND ',
     $       'REPORTS TOSSED = ',KSKSAT
            IF(IFLSAT.EQ.1)  THEN
               PRINT 8102
 8102 FORMAT(/' ---> IW3UNPBF: SUMMARY OF GOES REPORT COUNTS GROUPED',
     $ ' BY SATELLITE ID (PRIOR TO ANY FILTERING BY CALLING PROGRAM)'/)
               DO  IDSAT = 250,259
                  IF(KNTSAT(IDSAT).GT.0) PRINT 8103, IDSAT,KNTSAT(IDSAT)
               ENDDO
 8103 FORMAT(15X,'NUMBER FROM SAT. ID',I4,4X,':',I6)
               IF(KNTSAT(260).GT.0)  PRINT 8104, KNTSAT(260)
 8104 FORMAT(15X,'NUMBER FROM UNKNOWN SAT. ID:',I6)
               PRINT 8105
 8105          FORMAT(/)
            END IF
            KNDX = 0
            KSKACF = 0
            KSKUPA = 0
            KSKSFC = 0
            KSKSAT = 0
            IFLSAT = 0
         END IF
         IW3UNPBF = IER
      END IF
 
      STNID = STNIDX
      CRES1 = CRES1X
      CRES2 = CRES2X
      DSNAME = DSNAMX
      IDSDAT = IDSDAX
      IDSDMP = IDSDMX

      RETURN

 9999 CONTINUE
      IER = 999
      IW3UNPBF = 999
      RETURN

      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION I02UBF(LUNIT,OBS,IER)
 
      COMMON/IUBFCC/SUBSET
      COMMON/IUBFOO/DSNAMX,IDSDAX,IDSDMX
 
      CHARACTER*8  SUBSET,DSNAMX
      CHARACTER*6  C01UBF
      CHARACTER*4  CBUFR
      DIMENSION    OBS(2500),JDATE(5),JDUMP(5)
      INTEGER(8) IDSDAX,IDSDMX
 
      SAVE

      JDATE = -1
      JDUMP = -1
 
C  IF FILE IS CLOSED TRY TO OPEN IT AND RETURN INFORMATION ABOUT IT
C  IN THE FIRST THREE WORDS OF OBS
C  ----------------------------------------------------------------
 
      CALL STATUS(LUNIT,LUN,IL,IM)

      IF(IL.EQ.0) THEN
         IRET = -1
         I02UBF = 2
         REWIND LUNIT
         READ(LUNIT,END=11,ERR=12) CBUFR
         IF(CBUFR.NE.'BUFR') then
            PRINT *, '##IW3UNPBF/I02UBF - INPUT FILE ON UNIT ',LUNIT,
     $       ' IS NOT JBUFR -- IER = 999'
            I02UBF = 999
            GOTO 10
         END IF
         CALL DATELEN(10)
         CALL DUMPBF(LUNIT,JDATE,JDUMP)
cppppp
         print *, 'CENTER DATE (JDATE) = ',jdate
         print *, 'DUMP DATE (JDUMP) = ',jdump
cppppp
         IF(JDATE(1).LE.0)  then
            PRINT *, '##IW3UNPBF/I02UBF - CENTER DATE COULD NOT BE ',
     $       'OBTAINED FROM INPUT FILE ON UNIT ',LUNIT,' -- IER = 999'
            I02UBF = 999
            GO TO 10
         END IF
         IF(JDATE(1).LT.100)  THEN

C IF 2-DIGIT YEAR RETURNED IN JDATE(1), MUST USE "WINDOWING" TECHNIQUE
C  TO CREATE A 4-DIGIT YEAR

C IMPORTANT: IF DATELEN(10) IS CALLED, THE DATE HERE SHOULD ALWAYS
C            CONTAIN A 4-DIGIT YEAR, EVEN IF INPUT FILE IS NOT
C            Y2K COMPLIANT (BUFRLIB DOES THE WINDOWING HERE)

            PRINT *, '##IW3UNPBF/I02UBF - THE FOLLOWING SHOULD NEVER ',
     $       'HAPPEN!!!!!'
            PRINT *, '##IW3UNPBF/I02UBF - 2-DIGIT YEAR IN JDATE(1) ',
     $       'RETURNED FROM DUMPBF (JDATE IS: ',JDATE,') - USE ',
     $       'WINDWOING TECHNIQUE TO OBTAIN 4-DIGIT YEAR'
            IF(JDATE(1).GT.20)  THEN
               JDATE(1) = 1900 + JDATE(1)
            ELSE
               JDATE(1) = 2000 + JDATE(1)
            ENDIF
            PRINT *, '##IW3UNPBF/I02UBF - CORRECTED JDATE(1) WITH ',
     $       '4-DIGIT YEAR, JDATE NOW IS: ',JDATE
         ENDIF
!        IDSDAX = JDATE(1)*1E6+JDATE(2)*1E4+JDATE(3)*1E2+JDATE(4)
         IDSDAX = JDATE(1)*1000000+JDATE(2)*10000+JDATE(3)*100+JDATE(4)

         IF(JDUMP(1).LE.0)  THEN
            IDSDMX = 999999999999_8
         ELSE
            IF(JDUMP(1).LT.100)  THEN

C IF 2-DIGIT YEAR RETURNED IN JDUMP(1), MUST USE "WINDOWING" TECHNIQUE
C  TO CREATE A 4-DIGIT YEAR

C IMPORTANT: IF DATELEN(10) IS CALLED, THE DATE HERE SHOULD ALWAYS
C            CONTAIN A 4-DIGIT YEAR, EVEN IF INPUT FILE IS NOT
C            Y2K COMPLIANT (BUFRLIB DOES THE WINDOWING HERE)

               PRINT *, '##IW3UNPBF/I02UBF - THE FOLLOWING SHOULD ',
     $          'NEVER HAPPEN!!!!!'
               PRINT *, '##IW3UNPBF/I02UBF - 2-DIGIT YEAR IN JDUMP(1) ',
     $          'RETURNED FROM DUMPBF (JDUMP IS: ',JDUMP,') - USE ',
     $          'WINDWOING TECHNIQUE TO OBTAIN 4-DIGIT YEAR'
               IF(JDUMP(1).GT.20)  THEN
                  JDUMP(1) = 1900 + JDUMP(1)
               ELSE
                  JDUMP(1) = 2000 + JDUMP(1)
               ENDIF
               PRINT *, '##IW3UNPBF/I02UBF - CORRECTED JDUMP(1) WITH ',
     $          '4-DIGIT YEAR, JDUMP NOW IS: ',JDUMP
            END IF
!           IDSDMX = JDUMP(1)*1E8+JDUMP(2)*1E6+JDUMP(3)*1E4+
!    $       JDUMP(4)*1E2+JDUMP(5)
            IDSDMX = JDUMP(1)*100000000_8+JDUMP(2)*1000000_8+
     $       JDUMP(3)*10000_8+JDUMP(4)*100_8+JDUMP(5)
         ENDIF

         CALL OPENBF(LUNIT,'IN',LUNIT)

C This next call, I believe, is needed only because SUBSET is not
C  returned in DUMPBF ...
         call readmg(lunit,subset,idateb,iret)

         DSNAMX = C01UBF(SUBSET)//'  '

         I02UBF = 1
         GO TO 10
   11    CONTINUE
         PRINT *, '##IW3UNPBF/I02UBF - INPUT FILE ON UNIT ',LUNIT,
     $    ' IS EMPTY (NULL) -- ALL DONE WITH THIS FILE (IER = 3)'
         GO TO 10
   12    CONTINUE
         PRINT *, '##IW3UNPBF/I02UBF - ERROR READING INPUT FILE ON ',
     $    'UNIT ',LUNIT,' -- IER = 999'
         I02UBF = 999
   10    CONTINUE
         IER = I02UBF
         RETURN
      END IF
 
C  IF THE FILE IS ALREADY OPENED FOR INPUT TRY TO READ THE NEXT SUBSET
C  -------------------------------------------------------------------
 
      IF(IL.LT.0) THEN
 7822    CONTINUE
         CALL READNS(LUNIT,SUBSET,idateb,IRET)
         IF(IRET.EQ.0) I02UBF = R01UBF(SUBSET,LUNIT,OBS)
         IF(IRET.NE.0) I02UBF = 2
         IF(I02UBF.EQ.-9999)  GO TO 7822
         IER = I02UBF
         RETURN
      END IF
 
C  FILE MUST BE OPEN FOR INPUT!
C  ----------------------------
 
      PRINT *, '##IW3UNPBF/I02UBF - FILE ON UNIT ',LUNIT,' IS OPENED ',
     $ 'FOR OUTPUT -- IER = 999'
      I02UBF = 999
      IER = 999
      RETURN
 
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION C01UBF(SUBSET)
 
      CHARACTER*(*) SUBSET
      CHARACTER*6   C01UBF
 
      SAVE
 
      C01UBF = 'NONE'
 
      IF(SUBSET(1:5).EQ.'NC000')  THEN
         C01UBF = 'ADPSFC'
      ELSE  IF(SUBSET(1:5).EQ.'NC001')  THEN
         C01UBF = 'SFCSHP'
      ELSE  IF(SUBSET(1:5).EQ.'NC002')  THEN
         C01UBF = 'ADPUPA'
      ELSE  IF(SUBSET(1:5).EQ.'NC004')  THEN
         IF(SUBSET(6:8).EQ.'004')  THEN
            C01UBF = 'AIRCAR'
         ELSE  IF(SUBSET(6:8).EQ.'005')  THEN
            C01UBF = 'ADPUPA'
         ELSE
            C01UBF = 'AIRCFT'
         END IF
      ELSE  IF(SUBSET(1:5).EQ.'NC005')  THEN
         C01UBF = 'SATWND'
      END IF
 
      IF(C01UBF.EQ.'NONE') PRINT*,'##IW3UNPBF/C01UBF - UNKNOWN SUBSET ',
     $ '(=',SUBSET,') -- CONTINUE~~'
 
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION R01UBF(SUBSET,LUNIT,OBS)
 
      CHARACTER*(*) SUBSET
      CHARACTER*6   C01UBF,ADPSUB
      DIMENSION     OBS(*)

      SAVE
 
C  FIND SPECIFIED UNPACKED DATA TYPE AND CALL A TRANSLATOR
C  -------------------------------------------------------
 
      R01UBF = 4
      ADPSUB = C01UBF(SUBSET)
      IF(ADPSUB .EQ. 'ADPUPA')  THEN
         R01UBF = R03UBF(LUNIT,OBS)
      ELSE  IF(ADPSUB(1:3).EQ.'AIR')  THEN
         R01UBF = R05UBF(LUNIT,OBS)
      ELSE  IF(ADPSUB .EQ. 'SATWND')  THEN
         R01UBF = R06UBF(LUNIT,OBS)
      ELSE
         R01UBF = R04UBF(LUNIT,OBS)
      END IF

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE S01UBF(SID,XOB,YOB,RHR,RCTIM,RSV1,RSV2,ELV,ITP,RTP,
     $ RSTP)
C     ---> PROCESSES HEADER
 
      COMMON/IUBFDD/HDR(12),RCATS(50,160,11),IKAT(11),MCAT(11),NCAT(11)
      COMMON/IUBFNN/STNIDX,CRES1X,CRES2X
 
      CHARACTER*(*) RSV1,RSV2
      CHARACTER*8   SID,STNIDX,CRES1X,CRES2X
      DIMENSION     IHDR(12),RHDR(12)
      EQUIVALENCE   (IHDR(1),RHDR(1))

      SAVE
 
      DATA XMISS/99999./,IMISS/99999/,BMISS/10E10/
 
C  INITIALIZE THE UNPACK ARRAY TO MISSINGS
C  ---------------------------------------
 
C  NCAT will hold the number of unpacked levels in each category
      NCAT = 0
      RCATS = IMISS
 
C  STORE THE UNPACKED HEADER INFORMATION INTO UNP FORMAT
C  -----------------------------------------------------
 
      RHDR( 1) = MIN(YOB,XMISS)
cppppp
      IF(YOB.GE.BMISS)  print *, '~~IW3UNPBF/S01UBF: ID ',sid,' has a ',
     $ 'missing LATITUDE - unpked hdr, word 1 is set to ',RHDR(1)
cppppp

C Important: According to BUFR Manual, CLON (0-06-002) - represented
C  here by "XOB" - should be in units of Degrees West - and East +
C  (-180.0 to +180.0); however some BUFR data sets (e.g., PREPBUFR) are
C  known to encode CLON in units of Degree East (0.0 to 359.99) -- So
C  we use the following conversion to work in either case ...
      RHDR( 2) = XMISS
      IF(XOB.LT.BMISS)  RHDR( 2) = 360. - MOD(720.-XOB,360.)
      IF(RHDR(2).EQ.360.0)  RHDR(2) = 0.0
cppppp
      IF(XOB.GE.BMISS)  print *, '~~IW3UNPBF/S01UBF: ID ',sid,' has a ',
     $ 'missing LONGITUDE - unpked hdr, word 2 is set to ',RHDR(2)
cppppp
      RHDR( 3) = XMISS
      RHDR( 4) = MIN(RHR,XMISS)
cppppp
      IF(RHR.GE.BMISS)  print *, '~~IW3UNPBF/S01UBF: ID ',sid,' has a ',
     $ 'missing OB TIME - unpked hdr, word 4 is set to ',RHDR(4)
cppppp
      RHDR( 5) = XMISS
      RHDR( 6) = XMISS
      RHDR( 7) = NINT(ELV)
      IHDR( 8) = ITP
      IHDR( 9) = RTP
      IHDR(10) = RSTP
      RHDR(11) = XMISS
      IF(RCTIM.LT.24.01.AND.RCTIM.GT.-.01)  RHDR(11) = RCTIM
      RHDR(12) = XMISS
      STNIDX = SID
      CRES1X = RSV1
      CRES2X = RSV2
 
C  STORE THE HEADER INTO A HOLDING ARRAY
C  -------------------------------------
 
      HDR = RHDR
 
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE S02UBF(ICAT,N,*)
C     ---> PROCESSES CATEGORIES

C      Input argument N - level indicator (unless = 0, then signals
C                         subr. to write an empty cat. 2,3, or 4 level
 
      COMMON/IUBFDD/HDR(12),RCATS(50,160,11),IKAT(11),MCAT(11),NCAT(11)
      COMMON/IUBFEE/POB(255),QOB(255),TOB(255),ZOB(255),DOB(255),
     $               SOB(255),VSG(255),CLP(255),CLA(255),OB8(255),
     $               CF8(255)
      COMMON/IUBFFF/PQM(255),QQM(255),TQM(255),ZQM(255),WQM(255),
     $               QCP(255),QCA(255)
      COMMON/IUBFGG/PSL,STP,SDR,SSP,STM,DPD,TMX,TMI,HVZ,PRW,PW1,PW2,
     $               CCN,CHN,CTL,CTM,CTH,HCB,CPT,APT,PC6,SND,P24,
     $               DOP,POW,HOW,SWD,SWP,SWH,SST,SPG,SPD,SHC,SAS,WES
      COMMON/IUBFHH/PSQ,SPQ,SWQ,STQ,DDQ
      COMMON/IUBFII/PWMIN
      COMMON/IUBFLL/Q81(255),Q82(255)
      COMMON/IUBFMM/XIND(255)
      COMMON/IUBFNN/STNIDX,CRES1X,CRES2X

      CHARACTER*8 STNIDX,CRES1X,CRES2X
      DIMENSION   RCAT(50),JCAT(50)
      EQUIVALENCE (RCAT(1),JCAT(1))
      LOGICAL     SURF

      SAVE
 
      DATA IMISS/99999/,BMISS/10E10/
 
cppppp-ID
      iprint = 0

c     if(stnidx.eq.'68906   ')  iprint = 1
c     if(stnidx.eq.'59362   ')  iprint = 1
c     if(stnidx.eq.'57957   ')  iprint = 1
c     if(stnidx.eq.'74794   ')  iprint = 1
c     if(stnidx.eq.'74389   ')  iprint = 1
c     if(stnidx.eq.'96801A  ')  iprint = 1
cppppp-ID

      SURF = .FALSE.
      GOTO 1

C  ENTRY POINT SE01UBF FORCES DATA INTO THE SURFACE (FIRST) LEVEL
C  --------------------------------------------------------------

      ENTRY SE01UBF(ICAT,N)
C     ---> PROCESSES DATA INTO SURFACE LEVEL
      SURF = .TRUE.
 
C  CHECK THE PARAMETERS COMING IN
C  ------------------------------
 
1     KCAT = 0
      DO I = 1,11
         IF(ICAT.EQ.IKAT(I))  THEN
            KCAT = I
            GO TO 991
         END IF
      ENDDO

  991 CONTINUE

C  PARAMETER ICAT (UNPACKED CATEGORY) OUT OF BOUNDS RETURNS A 999
C  --------------------------------------------------------------

      IF(KCAT.EQ.0)  THEN
         PRINT *, '##IW3UNPBF/S02UBF - UNPACKED CATEGORY ',ICAT,' OUT ',
     $    'OF BOUNDS -- IER = 999'
         RETURN 1
      END IF

C  PARAMETER N (LEVEL INDEX) OUT OF BOUNDS RETURNS A 999
C  -----------------------------------------------------

      IF(N.GT.255)  THEN
         PRINT *, '##IW3UNPBF/S02UBF - LEVEL INDEX ',N,' EXCEEDS 255 ',
     $    '-- IER = 999'
         RETURN 1
      END IF
 
C  MAKE A MISSING LEVEL AND RETURN WHEN N=0 (NOT ALLOWED FOR CAT 01)
C   (NOTE: QUALITY MARKERS ARE SET TO 2 AND SPECIAL LEVEL INDICATORS
C          ARE SET TO 0)
C  -----------------------------------------------------------------
 
      IF(N.EQ.0) THEN
         IF(KCAT.EQ.1) RETURN
         NCAT(KCAT) = MIN(159,NCAT(KCAT)+1)
cppppp
         if(iprint.eq.1)  then
            print *, 'To prepare for sfc. data, write all missings on ',
     $       'lvl ',ncat(kcat),' for cat ',kcat
            print *, '  also write default =2 for all q. markers'
           print *, '  also write default =0 for special lvl indicators'
         end if
cppppp
         IF(KCAT.EQ.2)  THEN
            RCATS(4:6,NCAT(KCAT),2) = 2.
            RCATS(7,NCAT(KCAT),2)   = 0
         ELSE  IF(KCAT.EQ.3)  THEN
            RCATS(4:5,NCAT(KCAT),3) = 2.
            RCATS(6,NCAT(KCAT),3)   = 0
         ELSE  IF(KCAT.EQ.4)  THEN
            RCATS(4:5,NCAT(KCAT),4) = 2.
         END IF
         RETURN
      END IF
 
C  FIGURE OUT WHICH LEVEL TO UPDATE AND RESET THE LEVEL COUNTER
C  ------------------------------------------------------------
 
      IF(KCAT.EQ.1) THEN
         L = I04UBF(POB(N)*.1)

C  MANDATORY LEVEL WITH NON-MANDATORY PRESSURE RETURNS A 999
C  ---------------------------------------------------------

         IF(L.LE.0)  THEN
            PRINT *, '##IW3UNPBF/S02UBF - MANDATORY LEVEL WITH NON-',
     $       'MANDATORY PRESSURE (P = ',POB(N)*.1,'mb) -- IER = 999'
            RETURN 1
         END IF
         NCAT(KCAT) = MAX(NCAT(KCAT),L)
cppppp
         if(iprint.eq.1)
     $    print *, 'Will write cat. 1 data on lvl ',L,' for cat ',kcat,
     $    ', - total no. cat. 1 lvls processed so far = ',ncat(kcat)
cppppp
      ELSEIF(SURF) THEN
         L = 1
         NCAT(KCAT) = MAX(NCAT(KCAT),1)
cppppp
         if(iprint.eq.1)
     $    print *, 'Will write cat. ',kcat,' SURFACE data on lvl ',L,
     $    ', - total no. cat. ',kcat,' lvls processed so far = ',
     $    ncat(kcat)
cppppp
      ELSE
         L = MIN(159,NCAT(KCAT)+1)
         IF(L.EQ.159) THEN
cppppp
            print *, '~~IW3UNPBF/S02UBF: ID ',stnidx,
     $ ' - This cat. ',kcat,', level cannot be processed because ',
     $ 'the limit has already been reached'
cppppp
            RETURN
         END IF
         NCAT(KCAT) = L
cppppp
         if(iprint.eq.1)
     $    print *, 'Will write cat. ',kcat,' NON-SFC data on lvl ',L,
     $    ', - total no. cat. ',kcat,' lvls processed so far = ',
     $    ncat(kcat)
cppppp
      END IF
 
C  EACH CATEGORY NEEDS A SPECIFIC DATA ARRANGEMENT
C  -----------------------------------------------
 
      IF(ICAT.EQ.1) THEN
         RCAT(1)  = MIN(NINT(POB(N)),NINT(RCATS( 1,L,KCAT)))
         RCAT(2)  = MIN(NINT(ZOB(N)),NINT(RCATS( 2,L,KCAT)))
         RCAT(3)  = MIN(NINT(TOB(N)),NINT(RCATS( 3,L,KCAT)))
         RCAT(4)  = MIN(NINT(QOB(N)),NINT(RCATS( 4,L,KCAT)))
         RCAT(5)  = MIN(NINT(DOB(N)),NINT(RCATS( 5,L,KCAT)))
         RCAT(6)  = MIN(NINT(SOB(N)),NINT(RCATS( 6,L,KCAT)))
         IF(RCATS(7,L,KCAT).LT.IMISS)  THEN
            RCAT(7)  = MAX(NINT(PQM(N)),NINT(RCATS( 7,L,KCAT)))
         ELSE
            RCAT(7)  = NINT(PQM(N))
         END IF
         IF(RCATS(8,L,KCAT).LT.IMISS)  THEN
            RCAT(8)  = MAX(NINT(ZQM(N)),NINT(RCATS( 8,L,KCAT)))
         ELSE
            RCAT(8)  = NINT(ZQM(N))
         END IF
         IF(RCATS(9,L,KCAT).LT.IMISS)  THEN
            RCAT(9)  = MAX(NINT(TQM(N)),NINT(RCATS( 9,L,KCAT)))
         ELSE
            RCAT(9)  = NINT(TQM(N))
         END IF
         IF(RCATS(10,L,KCAT).LT.IMISS)  THEN
            RCAT(10) = MAX(NINT(QQM(N)),NINT(RCATS(10,L,KCAT)))
         ELSE
            RCAT(10) = NINT(QQM(N))
         END IF
         IF(RCATS(11,L,KCAT).LT.IMISS)  THEN
            RCAT(11) = MAX(NINT(WQM(N)),NINT(RCATS(11,L,KCAT)))
         ELSE
            RCAT(11) = NINT(WQM(N))
         END IF
      ELSEIF(ICAT.EQ.2) THEN
         RCAT(1) = MIN(NINT(POB(N)),IMISS)
         RCAT(2) = MIN(NINT(TOB(N)),IMISS)
         RCAT(3) = MIN(NINT(QOB(N)),IMISS)
         RCAT(4) = NINT(PQM(N))
         RCAT(5) = NINT(TQM(N))
         RCAT(6) = NINT(QQM(N))
         RCAT(7) = NINT(XIND(N))
      ELSEIF(ICAT.EQ.3) THEN
         RCAT(1) = MIN(NINT(POB(N)),IMISS)
         RCAT(2) = MIN(NINT(DOB(N)),IMISS)
         RCAT(3) = MIN(NINT(SOB(N)),IMISS)
         RCAT(4) = NINT(PQM(N))
         RCAT(5) = NINT(WQM(N))
         IF(NINT(VSG(N)).EQ.16)  THEN

C  MARK THE TROPOPAUSE LEVEL IN CAT. 3

            XIND(N) = 1
         ELSE  IF(NINT(VSG(N)).EQ. 8)  THEN

C  MARK THE MAXIMUM WIND LEVEL IN CAT. 3

            XIND(N) = 2
            IF(NINT(POB(N)).EQ.NINT(PWMIN))  XIND(N) = 3
         END IF
         RCAT(6) = NINT(XIND(N))
      ELSEIF(ICAT.EQ.4) THEN
         RCAT(1) = MIN(NINT(ZOB(N)),IMISS)
         RCAT(2) = MIN(NINT(DOB(N)),IMISS)
         RCAT(3) = MIN(NINT(SOB(N)),IMISS)
         RCAT(4) = NINT(ZQM(N))
         RCAT(5) = NINT(WQM(N))
      ELSEIF(ICAT.EQ.5) THEN
         RCAT(1) = MIN(NINT(POB(N)),IMISS)
         RCAT(2) = MIN(NINT(TOB(N)),IMISS)
         RCAT(3) = MIN(NINT(QOB(N)),IMISS)
         RCAT(4) = MIN(NINT(DOB(N)),IMISS)
         RCAT(5) = MIN(NINT(SOB(N)),IMISS)
         RCAT(6) = NINT(PQM(N))
         RCAT(7) = NINT(TQM(N))
         RCAT(8) = NINT(QQM(N))
         RCAT(9) = NINT(WQM(N))
      ELSEIF(ICAT.EQ.6) THEN
         RCAT(1)  = MIN(NINT(POB(N)),IMISS)
         RCAT(2)  = MIN(NINT(ZOB(N)),IMISS)
         RCAT(3)  = MIN(NINT(TOB(N)),IMISS)
         RCAT(4)  = MIN(NINT(QOB(N)),IMISS)
         RCAT(5)  = MIN(NINT(DOB(N)),IMISS)
         RCAT(6)  = MIN(NINT(SOB(N)),IMISS)
         RCAT(7)  = NINT(PQM(N))
         RCAT(8)  = NINT(ZQM(N))
         RCAT(9)  = NINT(TQM(N))
         RCAT(10) = NINT(QQM(N))
         RCAT(11) = NINT(WQM(N))
      ELSEIF(ICAT.EQ.7) THEN
         RCAT(1) = MIN(NINT(CLP(N)),IMISS)
         RCAT(2) = MIN(NINT(CLA(N)),IMISS)
         RCAT(3) = NINT(QCP(N))
         RCAT(4) = NINT(QCA(N))
      ELSEIF(ICAT.EQ.8) THEN
         RCAT(1) = MIN(NINT(OB8(N)),IMISS)
         RCAT(2) = MIN(NINT(CF8(N)),IMISS)
         RCAT(3) = NINT(Q81(N))
         RCAT(4) = NINT(Q82(N))
      ELSEIF(ICAT.EQ.51) THEN
         RCAT( 1) = MIN(NINT(PSL),IMISS)
         RCAT( 2) = MIN(NINT(STP),IMISS)
         RCAT( 3) = MIN(NINT(SDR),IMISS)
         RCAT( 4) = MIN(NINT(SSP),IMISS)
         RCAT( 5) = MIN(NINT(STM),IMISS)
         RCAT( 6) = MIN(NINT(DPD),IMISS)
         RCAT( 7) = MIN(NINT(TMX),IMISS)
         RCAT( 8) = MIN(NINT(TMI),IMISS)
         RCAT( 9) = NINT(PSQ)
         RCAT(10) = NINT(SPQ)
         RCAT(11) = NINT(SWQ)
         RCAT(12) = NINT(STQ)
         RCAT(13) = NINT(DDQ)
         JCAT(14) = MIN(NINT(HVZ),IMISS)
         JCAT(15) = MIN(NINT(PRW),IMISS)
         JCAT(16) = MIN(NINT(PW1),IMISS)
         JCAT(17) = MIN(NINT(PW2),IMISS)
         JCAT(18) = MIN(NINT(CCN),IMISS)
         JCAT(19) = MIN(NINT(CHN),IMISS)
         JCAT(20) = MIN(NINT(CTL),IMISS)
         JCAT(21) = MIN(NINT(HCB),IMISS)
         JCAT(22) = MIN(NINT(CTM),IMISS)
         JCAT(23) = MIN(NINT(CTH),IMISS)
         JCAT(24) = MIN(NINT(CPT),IMISS)
         RCAT(25) = MIN(IABS(NINT(APT)),IMISS)
         IF(CPT.GE.BMISS.AND.APT.LT.0.)
     $    RCAT(25) = MIN(IABS(NINT(APT))+500,IMISS)
      ELSEIF(ICAT.EQ.52) THEN
         JCAT( 1) = MIN(NINT(PC6),IMISS)
         JCAT( 2) = MIN(NINT(SND),IMISS)
         JCAT( 3) = MIN(NINT(P24),IMISS)
         JCAT( 4) = MIN(NINT(DOP),IMISS)
         JCAT( 5) = MIN(NINT(POW),IMISS)
         JCAT( 6) = MIN(NINT(HOW),IMISS)
         JCAT( 7) = MIN(NINT(SWD),IMISS)
         JCAT( 8) = MIN(NINT(SWP),IMISS)
         JCAT( 9) = MIN(NINT(SWH),IMISS)
         JCAT(10) = MIN(NINT(SST),IMISS)
         JCAT(11) = MIN(NINT(SPG),IMISS)
         JCAT(12) = MIN(NINT(SPD),IMISS)
         JCAT(13) = MIN(NINT(SHC),IMISS)
         JCAT(14) = MIN(NINT(SAS),IMISS)
         JCAT(15) = MIN(NINT(WES),IMISS)
      ELSE

C  UNSUPPORTED CATEGORY RETURNS A 999
C  ----------------------------------

         PRINT *, '##IW3UNPBF/S02UBF - CATEGORY ',ICAT,' NOT SUPPORTED',
     $    ' -- IER = 999'
         RETURN 1
      END IF
 
C  TRANSFER THE LEVEL DATA INTO THE HOLDING ARRAY AND EXIT
C  -------------------------------------------------------
 
      RCATS(1:MCAT(KCAT),L,KCAT) = RCAT(1:MCAT(KCAT))
 
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE S03UBF(UNP,*,*,*)
C     ---> PACKS DATA INTO UNP ARRAY
 
      COMMON/IUBFDD/HDR(12),RCATS(50,160,11),IKAT(11),MCAT(11),NCAT(11)
      COMMON/IUBFNN/STNIDX,CRES1X,CRES2X
 
      CHARACTER*8 STNIDX,CRES1X,CRES2X
      DIMENSION   RCAT(50),JCAT(50),UNP(*)
      EQUIVALENCE (RCAT(1),JCAT(1))

      SAVE

C  CALL TO SORT CATEGORIES 02, 03, 04, AND 08 LEVELS
C  -------------------------------------------------
 
      CALL S04UBF
 
C  TRANSFER DATA FROM ALL CATEGORIES INTO UNP ARRAY & SET POINTERS
C  ---------------------------------------------------------------
 
      INDX = 43
      JCAT = 0
      NLEVTO = 0
 
      DO K = 1,11
         JCAT(2*K+11) = NCAT(K)
         IF(K.NE.7.AND.K.NE.8.AND.K.NE.11)  NLEVTO = NLEVTO + NCAT(K)
         IF(NCAT(K).GT.0) THEN
            JCAT(2*K+12) = INDX
         ELSE
            JCAT(2*K+12) = 0
         END IF
         DO J = 1,NCAT(K)
            DO I = 1,MCAT(K)

C  UNPACKED REPORT CONTAINS MORE THAN 2500 WORDS - RETURNS A 999
C  -------------------------------------------------------------

               IF(INDX.GT.2500)  THEN
                  PRINT *, '~~IW3UNPBF/S03UBF: RPT with ID= ',STNIDX,
     $             ' TOSSED - CONTAINS ',INDX,' WORDS, > LIMIT OF 2500'
                  RETURN 3
               END IF
               UNP(INDX) = RCATS(I,J,K)
               INDX = INDX+1
            ENDDO
         ENDDO
      ENDDO

C  RETURN WITHOUT PROCESSING THIS REPORT IF NO DATA IN CAT. 1-6, 51, 52
C  --------------------------------------------------------------------

      IF(NLEVTO.EQ.0)  RETURN 2

C  TRANSFER THE HEADER AND POINTER ARRAYS INTO UNP
C  -----------------------------------------------
 
      UNP(1:12)  = HDR
      UNP(13:42) = RCAT(13:42)
 
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE S04UBF
C     ---> SORTS CATEGORIES
 
      COMMON/IUBFDD/HDR(12),RCATS(50,160,11),IKAT(11),MCAT(11),NCAT(11)
      COMMON/IUBFNN/STNIDX,CRES1X,CRES2X

      CHARACTER*8 STNIDX,CRES1X,CRES2X
 
      DIMENSION RCAT(50,160),IORD(160),IWORK(65536),SCAT(50,160),RCTL(3)
      DIMENSION PMAND(21)

      SAVE

      DATA OMISS/99999/,PMAND/10000.,9250.,8500.,7000.,5000.,4000.,
     $ 3000.,2500.,2000.,1500.,1000.,700.,500.,300.,200.,100.,70.,50.,
     $ 30.,20.,10./
 
C  INSERT CATEGORY 1 PRESSURE & DEF. Q.M.'S INTO THOSE LVLS WHERE MSSING
C  ---------------------------------------------------------------------

      DO I=1,NCAT(1)
         IF(RCATS(1,I,1).GE.OMISS)  THEN
            RCATS(1,I,1) = PMAND(I)
            RCATS(7:11,I,1) = 2.0
         END IF
      ENDDO

C  SORT CATEGORIES 2, 3, AND 4 - LEAVE THE FIRST LEVEL IN EACH INTACT
C  ------------------------------------------------------------------
 
      DO K=2,4
         IF(NCAT(K).GT.1) THEN
            DO J=1,NCAT(K)-1
               DO I=1,MCAT(K)
                  SCAT(I,J) = RCATS(I,J+1,K)
               ENDDO
            ENDDO
            CALL ORDERS(2,IWORK,SCAT(1,1),IORD,NCAT(K)-1,50,8,2)
            RCTL = 10E9
            DO J=1,NCAT(K)-1
               IF(K.LT.4) JJ = IORD((NCAT(K)-1)-J+1)
               IF(K.EQ.4) JJ = IORD(J)
               DO I=1,MCAT(K)
                  RCAT(I,J) = SCAT(I,JJ)
               ENDDO
               IDUP = 0
               IF(NINT(RCAT(1,J)).EQ.NINT(RCTL(1)))  THEN
                  IF(NINT(RCAT(2,J)).EQ.NINT(RCTL(2)).AND.
     $               NINT(RCAT(3,J)).EQ.NINT(RCTL(3)))  THEN
cppppp
                     if(k.ne.4)  then
                        print *,'~~@@IW3UNPBF/S04UBF: ID ',stnidx,
     $ ' has a dupl. cat. ',k,' lvl (all data) at ',rcat(1,j)*.1,
     $ ' mb -- lvl will be excluded from processing'
                     else
                     print *,'~~@@IW3UNPBF/S04UBF: ID ',stnidx,' has ',
     $ 'a dupl. cat. ',k,' lvl (all data) at ',rcat(1,j),' m -- lvl',
     $ ' will be excluded from processing'
                     end if
cppppp
                     IDUP = 1
                  ELSE
cppppp
                     if(k.ne.4)  then
                        print *,'~~@@#IW3UNPBF/S04UBF: ID ',stnidx,
     $ ' has a dupl. cat. ',k,' press. lvl (data differ) at ',
     $ rcat(1,j)*.1,' mb -- lvl will NOT be excluded'
                     else
                        print *,'~~@@#IW3UNPBF/S04UBF: ID ',stnidx,
     $ ' has a dupl. cat. ',k,' height lvl (data differ) at ',rcat(1,j),
     $ ' m -- lvl will NOT be excluded'
                     end if
cppppp
                  END IF
               END IF
               RCTL = RCAT(1:3,J)
               IF(IDUP.EQ.1)  RCAT(1,J) = 10E8
            ENDDO
            JJJ = 1
            DO J=2,NCAT(K)
               IF(RCAT(1,J-1).GE.10E8)  GO TO 887
               JJJ = JJJ + 1
               DO I=1,MCAT(K)
                  RCATS(I,JJJ,K) = RCAT(I,J-1)
               ENDDO
  887          CONTINUE
            ENDDO
cppppp
            if(jjj.ne.NCAT(K))
     $       print *,'~~@@IW3UNPBF/S04UBF: ID ',stnidx,' has had ',
     $       NCAT(K)-jjj,' lvls removed due to their being duplicates'
cppppp
            ncat(k) = jjj
         end if
         IF(NCAT(K).EQ.1)  THEN
            IF(MIN(RCATS(1,1,K),RCATS(2,1,K),RCATS(3,1,K)).GT.99998.8)
     $       NCAT(K) = 0
         END IF
      ENDDO
 
C  SORT CATEGORY 08 BY CODE FIGURE
C  -------------------------------
 
      DO K=8,8
         IF(NCAT(K).GT.1) THEN
            CALL ORDERS(2,IWORK,RCATS(2,1,K),IORD,NCAT(K),50,8,2)
            DO J=1,NCAT(K)
               DO I=1,MCAT(K)
                  RCAT(I,J) = RCATS(I,IORD(J),K)
               ENDDO
            ENDDO
            DO J=1,NCAT(K)
               DO I=1,MCAT(K)
                  RCATS(I,J,K) = RCAT(I,J)
               ENDDO
            ENDDO
         END IF
      ENDDO
 
C  NORMAL EXIT
C  -----------
 
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE S05UBF
C     ---> INITIALIZES INPUT ARRAYS
 
      CHARACTER*8  STNIDX,CRES1X,CRES2X

      COMMON/IUBFEE/OBS(255,11)
      COMMON/IUBFFF/QMS(255,7)
      COMMON/IUBFGG/SFO(35)
      COMMON/IUBFHH/SFQ(5)
      COMMON/IUBFLL/Q8(255,2)
      COMMON/IUBFMM/XIND(255)
      COMMON/IUBFNN/STNIDX,CRES1X,CRES2X
 
      SAVE
 
      DATA BMISS/10E10/
 
C  SET THE INPUT OBS DATA ARRAYS TO MISSING, INPUT Q.M. DATA ARRAYS
C  TO 2, INPUT CAT. 8 INDICATORS TO 99999 AND INPUT SPECIAL LEVEL
C  INDICATORS TO 0
C  ----------------------------------------------------------------
 
      OBS  = BMISS
      QMS  = 2
      SFO  = BMISS
      SFQ  = 2
      XIND = 0
      Q8   = 99999
      STNIDX = '        '
      CRES1X = '        '
      CRES2X = '        '
 
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION I04UBF(P)
C     ---> ASSIGNS NUMBER TO MANDATORY PRESSURE LEVEL
 
      COMMON/IUBFJJ/ISET,MANLIN(1001)

      SAVE
 
      IF(ISET.EQ.0) THEN
         MANLIN = 0

         MANLIN(1000) =  1
         MANLIN(925)  =  2
         MANLIN(850)  =  3
         MANLIN(700)  =  4
         MANLIN(500)  =  5
         MANLIN(400)  =  6
         MANLIN(300)  =  7
         MANLIN(250)  =  8
         MANLIN(200)  =  9
         MANLIN(150)  = 10
         MANLIN(100)  = 11
         MANLIN(70)   = 12
         MANLIN(50)   = 13
         MANLIN(30)   = 14
         MANLIN(20)   = 15
         MANLIN(10)   = 16
         MANLIN(7)    = 17
         MANLIN(5)    = 18
         MANLIN(3)    = 19
         MANLIN(2)    = 20
         MANLIN(1)    = 21

         ISET = 1
      END IF
 
      IP = NINT(P*10.)
 
      IF(IP.GT.10000 .OR. IP.LT.10 .OR. MOD(IP,10).NE.0) THEN
         I04UBF = 0
      ELSE
         I04UBF = MANLIN(IP/10)
      END IF
 
      RETURN

      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION R02UBF()
 
      CHARACTER*8 SUBSET,RPID
      LOGICAL     L02UBF,L03UBF

      SAVE
 
      DATA BMISS/10E10/
 
      R02UBF = 0
 
      RETURN

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENTRY ERTUBF(SUBSET,RPID)
         ERTUBF = BMISS
         IF(SUBSET(1:5).EQ.'NC000')  THEN
            IF(L02UBF(RPID))  ERTUBF = 511
            IF(L03UBF(RPID))  ERTUBF = 512
         ELSE  IF(SUBSET(1:5).EQ.'NC001')  THEN
            IF(SUBSET(6:8).EQ.'001')  THEN
               IF(RPID.NE.'SHIP')  THEN
                  ERTUBF = 522
               ELSE
                  ERTUBF = 523
               END IF
            ELSE  IF(SUBSET(6:8).EQ.'002') THEN
               ERTUBF = 562
            ELSE  IF(SUBSET(6:8).EQ.'003') THEN
               ERTUBF = 561
            ELSE  IF(SUBSET(6:8).EQ.'004') THEN
               ERTUBF = 531
            ELSE  IF(SUBSET(6:8).EQ.'005') THEN   ! add tide gauge station
               ERTUBF = 532
            END IF
         ELSE  IF(SUBSET(1:5).EQ.'NC002')  THEN
            IF(SUBSET(6:8).EQ.'001')  THEN

C  LAND RADIOSONDE - FIXED
C  -----------------------

               ERTUBF = 011
               IF(L03UBF(RPID)) ERTUBF = 012
               IF(RPID(1:4).EQ.'CLAS') ERTUBF = 013
            ELSE  IF(SUBSET(6:8).EQ.'002') THEN

C  LAND RADIOSONDE - MOBILE
C  ------------------------

               ERTUBF = 013
            ELSE  IF(SUBSET(6:8).EQ.'003') THEN

C  SHIP RADIOSONDE
C  ---------------

               ERTUBF = 022
               IF(RPID(1:4).EQ.'SHIP') ERTUBF = 023
            ELSE  IF(SUBSET(6:8).EQ.'004') THEN

C  DROPWINSONDE
C  -------------

               ERTUBF = 031
            ELSE  IF(SUBSET(6:8).EQ.'005') THEN

C  PIBAL
C  -----

               ERTUBF = 011
               IF(L03UBF(RPID)) ERTUBF = 012
            END IF

         ELSE  IF(SUBSET.EQ.'NC004005')  THEN

C  RECCOS/DROPS
C  ------------

            ERTUBF = 031
         END IF
         RETURN

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENTRY EWZUBF(HGT,Z100)

C  ALL WINDS-BY-HEIGHT HEIGHTS ARE TRUNCATED DOWN TO THE NEXT
C   10 METER LEVEL IF PART DD (ABOVE 100 MB LEVEL) (ON29 CONVENTION)
C  -----------------------------------------------------------------

         IF(HGT.GT.Z100)  THEN
            IF(MOD(NINT(HGT),10).NE.0)  HGT = INT(HGT/10.) * 10
            EWZUBF = NINT(HGT)
         ELSE
            IF(MOD(NINT(HGT/1.016),1500).EQ.0) HGT = NINT(HGT - 1.0)
            EWZUBF = INT(HGT)
         END IF
         RETURN

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENTRY EQMUBF(QMK)
         IF((NINT(QMK).GE.0 .AND. NINT(QMK).LT.4) .OR.
     $      (NINT(QMK).GT.9 .AND. NINT(QMK).LT.15))  THEN
            EQMUBF = NINT(QMK)
         ELSE  IF(NINT(QMK).EQ.6)  THEN
C  To get around the 3-bit limit to ON29 pressure q.m. mnemonic "QMPR",
C   a purge or reject flag on pressure in changed from 12 or 14 to 6 to
C   fit in 3-bits)
            EQMUBF = 14
         ELSE
            EQMUBF = 2
         END IF
         RETURN

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENTRY EQSUBF(QMK)
         IF((NINT(QMK).GE.0 .AND. NINT(QMK).LT.5) .OR.
     $      (NINT(QMK).GT.9 .AND. NINT(QMK).LT.15))  THEN
            EQSUBF = NINT(QMK)
         ELSE  IF(NINT(QMK).EQ.6)  THEN
C  To get around the 3-bit limit to ON29 pressure q.m. mnemonic "QMPR",
C   a purge or reject flag on pressure in changed from 12 or 14 to 6 to
C   fit in 3-bits)
            EQSUBF = 14
         ELSE
            EQSUBF = 2
         END IF
         RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION L01UBF()
      CHARACTER*8 RPID
      LOGICAL L01UBF,L02UBF,L03UBF

      SAVE

      L01UBF = .TRUE.

      RETURN

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENTRY L02UBF(RPID)
         L02UBF = .FALSE.
         READ(RPID,'(I5)',ERR=1) IBKS
         L02UBF = .TRUE.
1        RETURN

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ENTRY L03UBF(RPID)
         L03UBF = .TRUE.
         READ(RPID,'(I5)',ERR=2) IBKS
         L03UBF = .FALSE.
2        RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION R03UBF(LUNIT,OBS)
C     ---> PROCESSES ADPUPA DATA
 
      COMMON/IUBFBB/KNDX,KSKACF(8),KSKUPA,KSKSFC,KSKSAT
      COMMON/IUBFCC/SUBSET
      COMMON/IUBFDD/HDR(12),RCATS(50,160,11),IKAT(11),MCAT(11),NCAT(11)
      COMMON/IUBFEE/POB(255),QOB(255),TOB(255),ZOB(255),DOB(255),
     $               SOB(255),VSG(255),CLP(255),CLA(255),OB8(255),
     $               CF8(255)
      COMMON/IUBFFF/PQM(255),QQM(255),TQM(255),ZQM(255),WQM(255),
     $               QCP(255),QCA(255)
      COMMON/IUBFII/PWMIN
      COMMON/IUBFLL/Q81(255),Q82(255)
      COMMON/IUBFMM/XIND(255)
 
      CHARACTER*80 HDSTR,LVSTR,QMSTR,CLSTR,RCSTR
      CHARACTER*8  SUBSET,SID,RSV1,RSV2
      INTEGER    ICLA(0:8)
      REAL(8) RID,HDR_8(12),VSG_8(255)
      REAL(8)  RCT_8(5,255),ARR_8(10,255)
      REAL(8)  RAT_8(255),RMORE_8(4),RGP10_8(255),PRGP10_8(255)
      REAL(8)  RPMSL_8,RPSAL_8
      DIMENSION  OBS(*),RCT(5,255),ARR(10,255)
      DIMENSION  RAT(255),RMORE(4),RGP10(255),PRGP10(255)
      DIMENSION  P2(255),P8(255),P16(255)

      EQUIVALENCE  (RID,SID)
      LOGICAL      L02UBF

      SAVE
 
      DATA HDSTR/'RPID CLON CLAT HOUR MINU SELV               '/
      DATA LVSTR/'PRLC TMDP TMDB GP07 GP10 WDIR WSPD          '/
      DATA QMSTR/'QMPR QMAT QMDD QMGP QMWN                    '/
      DATA CLSTR/'HOCB CLAM QMCA                              '/
      DATA RCSTR/'RCHR RCMI RCTS                              '/
 
      DATA ICLA/0,5,25,40,50,60,75,95,100/
      DATA BMISS/10E10/,GRAV/9.8/
 
      PRS1(Z) = 1013.25 * (((288.15 - (.0065 * Z))/288.15)**5.256)
      PRS2(Z) = 226.3 * EXP(1.576106E-4 * (11000. - Z))

      R03UBF = 0
 
      CALL S05UBF
 
C  VERTICAL SIGNIFICANCE DESCRIPTOR TO ASSIGN UNPACKED LEVEL CATEGORY
C  ------------------------------------------------------------------

C NOTE: MNEMONIC "VSIG" 008001 IS DEFINED AS VERTICAL SOUNDING
C       SIGNIFICANCE -- CODE TABLE FOLLOWS:
C        64   Surface 
C              processed as unpacked category 2 and/or 3 and/or 4
C        32   Standard (mandatory) level
C              processed as unpacked category 1
C        16   Tropopause level
C              processed as unpacked category 5
C         8   Maximum wind level
C              processed as unpacked category 3 or 4
C         4   Significant level, temperature
C              processed as unpacked category 2
C         2   Significant level, wind
C              processed as unpacked category 3 or 4
C
C NOTE: THIS SUBR. ASSIGNS VSIG=1 TO LEVELS THAT SHOULD BE
C       PROCESSED AS UNPACKED CATEGORY 6 (ONLY APPLIES TO
C       RECCOS)
C
C  anything else - the level is not processed
 
      CALL UFBINT(LUNIT,VSG_8,1,255,NLEV,'VSIG');VSG=VSG_8
 
C  PUT THE HEADER INFORMATION INTO UNPACKED FORMAT
C  -----------------------------------------------
 
      CALL UFBINT(LUNIT,HDR_8,12,  1,IRET,HDSTR);HDR(2:)=HDR_8(2:)
      IF(HDR(5).GE.BMISS) HDR(5) = 0
      RID = HDR_8(1)
cppppp-ID
      iprint = 0
c     if(sid.eq.'68906   ')  iprint = 1
c     if(sid.eq.'61094   ')  iprint = 1
c     if(sid.eq.'62414   ')  iprint = 1
c     if(sid.eq.'59362   ')  iprint = 1
c     if(sid.eq.'57957   ')  iprint = 1
c     if(sid.eq.'74794   ')  iprint = 1
c     if(sid.eq.'74389   ')  iprint = 1
c     if(sid.eq.'96801A  ')  iprint = 1
      if(iprint.eq.1)
     $ print *, '@@@ START DIAGNOSTIC PRINTOUT FOR ID ',sid
cppppp-ID

      RGP10  = BMISS
      PRGP10 = BMISS
      CALL UFBINT(LUNIT,RPMSL_8,1,  1,IRET,'PMSL');RPMSL=RPMSL_8
      IF(SUBSET.EQ.'NC004005')  THEN
         CALL UFBINT(LUNIT,RGP10_8,1,255,NLEV,'GP10');RGP10=RGP10_8
         CALL UFBINT(LUNIT,RPSAL_8,1,1,IRET,'PSAL');RPSAL=RPSAL_8
         IF(NINT(VSG(1)).EQ.32.AND.RPMSL.GE.BMISS.AND.
     $    MAX(RGP10(1),RPSAL).LT.BMISS)  THEN
cppppp
cdak        print *, '~~IW3UNPBF/R03UBF: ID ',sid,' is a Cat. 6 type ',
cdak $       'RECCO with a mand-lvl GEOPOT'
cppppp
            VSG(1) = 1
            CALL UFBINT(LUNIT,PRGP10_8,1,255,NLEV,'PRLC')
            PRGP10=PRGP10_8
            HDR(6) = RPSAL + SIGN(0.0000001,RPSAL)
         ELSE  IF(MIN(VSG(1),RPMSL,RGP10(1)).GE.BMISS.AND.RPSAL.LT.
     $    BMISS)  THEN
cppppp
cdak        print *, '~~IW3UNPBF/R03UBF: ID ',sid,' is a Cat. 6 type ',
cdak $       'RECCO'
cppppp
            VSG(1) = 1
            HDR(6) = RPSAL + SIGN(0.0000001,RPSAL)
         ELSE  IF(MIN(VSG(1),RGP10(1)).GE.BMISS.AND.MAX(RPMSL,RPSAL)
     $    .LT.BMISS)  THEN
cppppp
cdak        print *, '~~IW3UNPBF/R03UBF: ID ',sid,' is a Cat. 6 type ',
cdak $       'RECCO with a valid PMSL'
cppppp
            VSG(1) = 1
            HDR(6) = RPSAL + SIGN(0.0000001,RPSAL)
         ELSE
cppppp
            print *, '~~IW3UNPBF/R03UBF: ID ',sid,' is currently an ',
     $       'unknown type of RECCO - VSIG =',VSG(1),'; PMSL =',RPMSL,
     $       '; GP10 =',RGP10(1),' -- SKIP IT for now'
cppppp
            R03UBF = -9999
            KSKUPA =KSKUPA + 1
            RETURN
         END IF
      END IF

      XOB = HDR(2)
      YOB = HDR(3)

      RHR = BMISS
      IF(HDR(4).LT.BMISS)  RHR = NINT(HDR(4))+NINT(HDR(5))/60.

      RCTIM = BMISS
      RSV1 = '        '
      RSV2 = '        '
      ELV = HDR(6)
     
      CALL UFBINT(LUNIT,RAT_8, 1,255,NLEV,'RATP');RAT=RAT_8
      ITP = MIN(99,NINT(RAT(1)))
      RTP = ERTUBF(SUBSET,SID)
      IF(ELV.GE.BMISS)  THEN
         IF((RTP.GT.20.AND.RTP.LT.24).OR.SUBSET.EQ.'NC002004')  THEN
cppppp
            print *, 'IW3UNPBF/R03UBF: ID ',sid,' has a missing elev, ',
     $       'so elevation set to ZERO'
cppppp
            ELV = 0
         ELSE
cppppp
            print *, '~~IW3UNPBF/R03UBF: ID ',sid,' has a missing elev'
cppppp
         END IF
      END IF
cdak  if(sid(5:5).eq.' ') print*,sid
      IF(RID.EQ.BMISS) SID = 'MISSING'
      IF(L02UBF(SID).AND.SID(5:5).EQ.' ') SID = '0'//SID
      RSTP = 99999
      IF(RTP.EQ.31)  THEN
         IF(SUBSET.EQ.'NC004005')  THEN
            RSTP = 1
         ELSE
            RSTP = 2
         END IF
      END IF
      IF (SUBSET.EQ.'NC002005') THEN
         RSTP = 3
      END IF
      CALL S01UBF(SID,XOB,YOB,RHR,RCTIM,RSV1,RSV2,ELV,ITP,RTP,RSTP)
 
C  PUT THE LEVEL DATA INTO SPECIFIED UNPACKED FORMAT
C  -------------------------------------------------
 
      CALL UFBINT(LUNIT,ARR_8,10,255,NLEV,LVSTR);ARR=ARR_8

      PWMIN = 999999.

      DO L=1,NLEV

         POB(L) = BMISS
         IF(ARR(1,L).LT.BMISS)  POB(L) = NINT(ARR(1,L)*.1)
         
         IF(NINT(ARR(1,L)).LE.0) THEN
            POB(L) =  BMISS
cppppp
            print *,'~~@@IW3UNPBF/R03UBF: ID ',sid,' has a ZERO or ',
     $       'negative reported pressure that is reset to missing'
cppppp
         END IF

         QOB(L) = BMISS
         IF(ARR(2,L).LT.BMISS .AND. ARR(3,L).LT.BMISS)
     $    QOB(L) = (ARR(3,L)-ARR(2,L))*10.

         TOB(L) = BMISS
         ITMP = NINT(ARR(3,L)*100.)
         IF(ARR(3,L).LT.BMISS)  TOB(L) = NINT((ITMP-27315)*0.1)

         XXX = BMISS
         IF(ARR(4,L).LT.BMISS) XXX = (ARR(4,L)/GRAV)
         YYY = BMISS
         IF(ARR(5,L).LT.BMISS) YYY = (ARR(5,L)/GRAV)
         ZOB(L) = MIN(XXX,YYY)

         DOB(L) = ARR(6,L)
         SOB(L) = MIN(ARR(7,L)*10.,BMISS)
         IF(NINT(DOB(L)).EQ.0.AND.NINT(SOB(L)).GT.0)  THEN
            DOB(L) = 360
         ELSE  IF(NINT(DOB(L)).EQ.360.AND.NINT(SOB(L)).EQ.0)  THEN
            DOB(L) = 0
         END IF
cppppp
      if(iprint.eq.1)  then
         print *, 'At lvl=',L,'; VSG=',vsg(L),'; POB = ',pob(L),
     $    '; QOB = ',qob(L),'; TOB = ',tob(L),'; ZOB = ',zob(L),
     $    '; DOB = ',dob(L),'; SOB = ',sob(L)
      end if
cppppp
         IF(MAX(POB(L),DOB(L),SOB(L)).LT.BMISS) PWMIN =MIN(PWMIN,POB(L))
      ENDDO

      CALL UFBINT(LUNIT,ARR_8,10,255,NLEV,QMSTR);ARR=ARR_8

      DO L=1,NLEV
         PQM(L) = EQMUBF(ARR(1,L))
         TQM(L) = EQMUBF(ARR(2,L))
         QQM(L) = EQMUBF(ARR(3,L))
         ZQM(L) = EQMUBF(ARR(4,L))
         WQM(L) = EQMUBF(ARR(5,L))
      ENDDO

C  SURFACE DATA MUST GO FIRST
C  --------------------------
 
      CALL S02UBF(2,0,*9999)
      CALL S02UBF(3,0,*9999)
      CALL S02UBF(4,0,*9999)

      INDX2  = 0
      INDX8  = 0
      INDX16 = 0
      P2  = BMISS
      P8  = BMISS
      P16 = BMISS

      DO L=1,NLEV
      IF(NINT(VSG(L)).EQ.64) THEN
cppppp
      if(iprint.eq.1)  then
         print *, 'Lvl=',L,' is a surface level'
      end if
      if(iprint.eq.1.and.POB(L).LT.BMISS.AND.TOB(L).LT.BMISS)  then
         print *, ' --> valid cat. 2 sfc. lvl '
      end if
cppppp
         IF(POB(L).LT.BMISS.AND.TOB(L).LT.BMISS)  CALL SE01UBF(2,L)
cppppp
      if(iprint.eq.1.and.POB(L).LT.BMISS.AND.DOB(L).LT.BMISS)  then
         print *, ' --> valid cat. 3 sfc. lvl '
      end if
cppppp
         IF(POB(L).LT.BMISS.AND.DOB(L).LT.BMISS)  CALL SE01UBF(3,L)
         IF(MAX(ZOB(L),DOB(L)).LT.BMISS) THEN
cppppp
            if(iprint.eq.1)  print *, ' --> valid cat. 4 sfc. lvl '
cppppp

C  CAT. 4 HEIGHT DOES NOT PASS ON A KEEP, PURGE, OR REJECT LIST Q.M.

            ZQM(L) = 2
            CALL SE01UBF(4,L)
         END IF
         VSG(L) = 0
      ELSE  IF(NINT(VSG(L)).EQ.2)  THEN
         P2(L) = POB(L)
         INDX2 = L
         IF(INDX8.GT.0)  THEN
            DO II = 1,INDX8
               IF(NINT(POB(L)).EQ.NINT(P8(II)).AND.POB(L).LT.BMISS) THEN
cppppp
                  if(iprint.eq.1)  then
                     print *, ' ## This cat. 3 level, on lvl ',L,
     $                ' will have already been processed as a cat. 3 ',
     $                'MAX wind lvl (on lvl ',II,') - skip this Cat. ',
     $                '3 lvl'
                  end if
cppppp
                  IF(MAX(SOB(II),DOB(II)).GE.BMISS)  THEN
                     SOB(II) = SOB(L)
                     DOB(II) = DOB(L)
cppppp
                     if(iprint.eq.1)  then
                        print *, ' ...... also on lvl ',L,' - transfer',
     $                  ' wind data to dupl. MAX wind lvl because its ',
     $                  'missing there'
                     end if
cppppp
                  END IF
                  VSG(L) = 0
                  GO TO 7732
               END IF
            ENDDO
         END IF
      ELSE  IF(NINT(VSG(L)).EQ.8)  THEN
         P8(L) = POB(L)
         INDX8 = L
         IF(INDX2.GT.0)  THEN
            DO II = 1,INDX2
               IF(NINT(POB(L)).EQ.NINT(P2(II)).AND.POB(L).LT.BMISS)THEN
cppppp
                  if(iprint.eq.1)  then
                     print *, ' ## This MAX wind level, on lvl ',L,
     $                ' will have already been processed as a cat. 3 ',
     $                'lvl (on lvl ',II,') - skip this MAX wind lvl ',
     $                'but set'
                     print *, '     cat. 3 lvl XIND to "2"'
                  end if
cppppp
                  XIND(II) = 2
                  IF(NINT(POB(L)).EQ.NINT(PWMIN))  XIND(II) = 3
                  IF(MAX(SOB(II),DOB(II)).GE.BMISS)  THEN
                     SOB(II) = SOB(L)
                     DOB(II) = DOB(L)
cppppp
                     if(iprint.eq.1)  then
                        print *, ' ...... also on lvl ',L,' - transfer',
     $                  ' wind data to dupl. cat. 3 lvl because its ',
     $                  'missing there'
                     end if
cppppp
                  END IF
                  VSG(L) = 0
                  GO TO 7732
               END IF
            ENDDO
         END IF
         IF(INDX8-1.GT.0)  THEN
            DO II = 1,INDX8-1
               IF(NINT(POB(L)).EQ.NINT(P8(II)).AND.POB(L).LT.BMISS)THEN
cppppp
                  if(iprint.eq.1)  then
                     print *, ' ## This cat. 3 MAX wind lvl, on lvl ',L,
     $                ' will have already been processed as a cat. 3 ',
     $                'MAX wind lvl (on lvl ',II,') - skip this Cat. ',
     $                '3 MAX wind lvl'
                  end if
cppppp
                  IF(MAX(SOB(II),DOB(II)).GE.BMISS)  THEN
                     SOB(II) = SOB(L)
                     DOB(II) = DOB(L)
cppppp
                     if(iprint.eq.1)  then
                        print *, ' ...... also on lvl ',L,' - transfer',
     $                  ' wind data to dupl. MAX wind lvl because its ',
     $                  'missing there'
                     end if
cppppp
                  END IF
                  VSG(L) = 0
                  GO TO 7732
               END IF
            ENDDO
         END IF
      ELSE  IF(NINT(VSG(L)).EQ.16)  THEN
         INDX16 = INDX16 + 1
         P16(INDX16) = POB(L)
      END IF
 7732 CONTINUE
      ENDDO
 
 
C  REST OF THE DATA
C  ----------------
 
      Z100 = 16000
      DO L=1,NLEV
      IF(NINT(VSG(L)).EQ.32) THEN
         IF(MIN(DOB(L),ZOB(L),TOB(L)).GE.BMISS)  THEN
cppppp
            if(iprint.eq.1)  then
               print *,' ==> For lvl ',L,'; VSG=32 & DOB,ZOB,TOB all ',
     $          'missing --> this level not processed'
            end if
            VSG(L) = 0
         ELSE  IF(MIN(ZOB(L),TOB(L)).LT.BMISS) THEN
cppppp
            if(iprint.eq.1)  then
               print *,' ==> For lvl ',L,'; VSG=32 & one or both of ',
     $          'ZOB,TOB non-missing --> valid cat. 1 lvl'
            end if
cppppp
            CALL S02UBF(1,L,*9999)
            IF(NINT(POB(L)).EQ.1000.AND.ZOB(L).LT.BMISS)  Z100 = ZOB(L)
            VSG(L) = 0
         END IF
      END IF
      ENDDO
      DO L=1,NLEV
      IF(NINT(VSG(L)).EQ.32) THEN
         IF(DOB(L).LT.BMISS.AND.MIN(ZOB(L),TOB(L)).GE.BMISS) THEN
            LL = I04UBF(POB(L)*.1)
            IF(LL.LE.0)  THEN
cppppp
               print *, '~~IW3UNPBF/R03UBF: ID ',sid,' has VSG=32 for ',
     $          'lvl ',L,' but pressure (=',POB(L)*.1,'mb) not mand.!!',
     $          ' --> this level not processed'
cppppp
            ELSE  IF(MIN(RCATS(2,LL,1),RCATS(3,LL,1)).LT.99999.)  THEN
               IF(RCATS(5,LL,1).GE.99998.)  THEN
cppppp
                  if(iprint.eq.1)  then
                     print *,' ==> For lvl ',L,'; VSG=32 & ZOB,TOB ',
     $                'both missing while DOB non-missing BUT one or ',
     $                'both of Z, T non-missing while wind missing in'
                     print *,'      earlier cat. 1 processing of this ',
     $                POB(L)*.1,'mb level --> valid cat. 1 lvl'
                  end if
cppppp
                  CALL S02UBF(1,L,*9999)
               ELSE
cppppp
                  if(iprint.eq.1)  then
                     print *,' ==> For lvl ',L,'; VSG=32 & ZOB,TOB ',
     $                'both missing while DOB non-missing BUT one or ',
     $                'both of Z, T non-missing while wind non-missing',
     $                ' in'
                     print *,'      earlier cat. 1 processing of this ',
     $                POB(L)*.1,'mb level --> valid cat. 3 lvl'
                  end if
cppppp
                  CALL S02UBF(3,L,*9999)
               END IF
            ELSE
cppppp
               if(iprint.eq.1)  then
                  print *,' ==> For lvl ',L,'; VSG=32 & ZOB,TOB both ',
     $             'missing while DOB non-missing AND both Z, T ',
     $             'missing on'
                  print *,'      this ',POB(L)*.1,'mb level in cat. 1 ',
     $             ' --> valid cat. 3 lvl'
               end if
cppppp
               CALL S02UBF(3,L,*9999)
            END IF
         ELSE
cppppp
            print *, '~~IW3UNPBF/R03UBF: ID ',sid,' has VSG=32 for ',
     $       'lvl ',L,' & should never come here!! - by default output',
     $       ' as cat. 1 lvl'
cppppp
            CALL S02UBF(1,L,*9999)
         END IF
         VSG(L) = 0
      END IF
      ENDDO

      DO L=1,NLEV
      IF(NINT(VSG(L)).EQ. 4) THEN
cppppp
         if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG= 4 --> valid cat. 2 lvl'
         end if
cppppp
         IF(INDX16.GT.0)  THEN
            DO II = 1,INDX16
               IF(NINT(POB(L)).EQ.NINT(P16(II)).AND.POB(L).LT.BMISS)THEN
cppppp
                  if(iprint.eq.1)  then
                     print *, ' ## This cat. 2 level, on lvl ',L,' is',
     $                ' also the tropopause level, as its pressure ',
     $                'matches that of trop. lvl no. ',II,' - ',
     $                'set this cat. 2'
                     print *, '    lvl XIND to "1"'
                  end if
cppppp
                  XIND(L) = 1
                  GO TO 7738
               END IF
            ENDDO
         END IF
 7738    CONTINUE
         CALL S02UBF(2,L,*9999)
         VSG(L) = 0
      ELSEIF(NINT(VSG(L)).EQ.16) THEN
cppppp
         if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG=16 --> valid cat. 3/5 lvl'
         end if
cppppp
         IF(MIN(SOB(L),DOB(L)).LT.BMISS)  CALL S02UBF(3,L,*9999)
         CALL S02UBF(5,L,*9999)
         VSG(L) = 0
      ELSEIF(NINT(VSG(L)).EQ. 1) THEN
cppppp
         if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG=1 --> valid cat. 6 lvl',
     $      ' - can only be a Recco'
         end if
cppppp
         POB(L) = BMISS
         ZOB(L) = ELV
         CALL S02UBF(6,L,*9999)
         VSG(L) = 0
      ELSEIF(NINT(VSG(L)).EQ. 2)  THEN
         IF(POB(L).LT.BMISS) THEN
            IF(MAX(SOB(L),DOB(L)).LT.BMISS)  THEN
cppppp
            if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG= 2 & POB .ne. missing ',
     $       '--> valid cat. 3 lvl (expect that ZOB is missing)'
            end if
cppppp
            CALL S02UBF(3,L,*9999)
            ELSE
cppppp
            if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG= 2 & POB .ne. missing ',
     $       '--> Cat. 3 level not processed - wind is missing'
            end if
cppppp
            END IF
            VSG(L) = 0
         ELSE IF(ZOB(L).LT.BMISS) THEN
            IF(MAX(SOB(L),DOB(L)).LT.BMISS)  THEN

C  CERTAIN U.S. WINDS-BY-HEIGHT ARE CORRECTED TO ON29 CONVENTION
C  -------------------------------------------------------------

            IF((SID(1:2).GE.'70'.AND.SID(1:2).LE.'72').OR.SID(1:2).EQ.
     $       '74')  ZOB(L) = EWZUBF(ZOB(L),Z100)
cppppp
            if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG= 2 & ZOB .ne. missing ',
     $       '--> valid cat. 4 lvl (POB must always be missing)'
            if(sid(1:2).eq.'70'.or.sid(1:2).eq.'71'.or.sid(1:2).eq.'72'
     $       .or.sid(1:2).eq.'74')  print *, '   .... ZOB at this ',
     $       'U.S. site adjusted to ',zob(L)
            end if
cppppp

C  CAT. 4 HEIGHT DOES NOT PASS ON A KEEP, PURGE, OR REJECT LIST Q.M.
C  -----------------------------------------------------------------

            ZQM(L) = 2

            CALL S02UBF(4,L,*9999)
            ELSE
cppppp
            if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG= 2 & ZOB .ne. missing ',
     $       '--> Cat. 4 level not processed - wind is missing'
            end if
cppppp
            END IF
            VSG(L) = 0
         END IF
      ELSEIF(NINT(VSG(L)).EQ. 8)  THEN
         IF(POB(L).LT.BMISS) THEN
cppppp
            if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG= 8 & POB .ne. missing ',
     $       '--> valid cat. 3 lvl (expect that ZOB is missing)'
            end if
cppppp
            CALL S02UBF(3,L,*9999)
            VSG(L) = 0
         ELSE  IF(ZOB(L).LT.BMISS) THEN
            IF(MAX(SOB(L),DOB(L)).LT.BMISS)  THEN

C  CERTAIN U.S. WINDS-BY-HEIGHT ARE CORRECTED TO ON29 CONVENTION
C  -------------------------------------------------------------

               IF((SID(1:2).GE.'70'.AND.SID(1:2).LE.'72').OR.SID(1:2)
     $          .EQ.'74')  ZOB(L) = EWZUBF(ZOB(L),Z100)
cppppp
               if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG= 8 & ZOB .ne. missing ',
     $       '--> valid cat. 4 lvl (POB must always be missing)'
            if(sid(1:2).eq.'70'.or.sid(1:2).eq.'71'.or.sid(1:2).eq.'72'
     $       .or.sid(1:2).eq.'74')  print *, '   .... ZOB at this ',
     $       'U.S. site adjusted to ',zob(L)
               end if
cppppp

C  CAT. 4 HEIGHT DOES NOT PASS ON A KEEP, PURGE, OR REJECT LIST Q.M.
C  -----------------------------------------------------------------

               ZQM(L) = 2

               CALL S02UBF(4,L,*9999)
            ELSE
cppppp
               if(iprint.eq.1)  then
            print *, ' ==> For lvl ',L,'; VSG= 8 & ZOB .ne. missing ',
     $       '--> Cat. 4 level not processed - wind is missing'
               end if
cppppp
            END IF
            VSG(L) = 0
         END IF
      END IF
      ENDDO
 
C  CHECK FOR LEVELS WHICH GOT LEFT OUT
C  -----------------------------------
 
      DO L=1,NLEV
      IF(NINT(VSG(L)).GT.0)  THEN
         PRINT 887, L,SID,NINT(VSG(L))
  887 FORMAT(' ##IW3UNPBF/R03UBF - ~~ON LVL',I4,' OF ID ',A8,', A ',
     $    'VERTICAL SIGNIFICANCE OF',I3,' WAS NOT SUPPORTED - LEAVE ',
     $    'THIS LEVEL OUT OF THE PROCESSING')
         print *, ' ..... at lvl=',L,'; POB = ',pob(L),'; QOB = ',
     $    qob(L),'; TOB = ',tob(L),'; ZOB = ',zob(L),'; DOB = ',dob(L),
     $    ';'
         print *, '                  SOB = ',sob(L)
      END IF
      ENDDO
 
C  CLOUD DATA GOES INTO CATEGORY 07
C  --------------------------------
 
      CALL UFBINT(LUNIT,ARR_8,10,255,NLEV,CLSTR);ARR=ARR_8
      DO L=1,NLEV

         IF(ELV+ARR(1,L).GE.BMISS)  THEN
            CLP(L) =  BMISS
         ELSE  IF(ELV+ARR(1,L).LE.11000)  THEN
            CLP(L) = (PRS1(ELV+ARR(1,L))*10.)
         ELSE
            CLP(L) = (PRS2(ELV+ARR(1,L))*10.)
         END IF

         IF(NINT(ARR(2,L)).GT.-1.AND.NINT(ARR(2,L)).LT.9)  THEN
            CLA(L) = ICLA(NINT(ARR(2,L)))
         ELSE
            CLA(L) = BMISS
         END IF

         QCP(L) = 2
         QCA(L) = EQMUBF(ARR(3,L))
         IF(MIN(CLP(L),CLA(L)).LT.BMISS)  CALL S02UBF(7,L,*9999)
      ENDDO
 
C  -----------------------------------------------------
C  MISC DATA GOES INTO CATEGORY 08
C  -----------------------------------------------------
C  CODE FIGURE 104 - RELEASE TIME IN .01*HR
C  CODE FIGURE 105 - RECEIPT TIME IN .01*HR
C  CODE FIGURE 351 - GEOPOTENTIAL HEIGHT IN METERS FOR
C                    PRESSURE LEVEL DEFINED IN INDICATOR 1 (MB)
C                    (QUALITY MARKER FOR GEOPOTENTIAL STORED
C                     IN INDICATOR 2)
C  CODE FIGURE 352 - MEAN-SEA LEVEL PRESSURE IN .1*MB
C  CODE FIGURE 353 - SOLAR AND INFRARED RADIATION CORRECTION
C                    INDICATOR
C  CODE FIGURE 354 - TRACKING TECHNIQUE/STATUS OF SYSTEM USED
C                    INDICATOR
C  -----------------------------------------------------------
 
      CALL UFBINT(LUNIT,RCT_8, 5,255,NRCT,RCSTR);RCT=RCT_8

C NOTE: MNEMONIC "RCTS" 008202 IS A LOCAL DESCRIPTOR DEFINED AS
C       RECEIPT TIME SIGNIFICANCE - THIS IS STORED IN INDICATOR 1
C       FOR CATEGORY 8 CODE FIGURE 105 -- CODE TABLE FOLLOWS:
C         0   General decoder receipt time
C         1   NCEP receipt time
C         2   OSO  receipt time
C         3   ARINC ground station receipt time
C         4   Radiosonde TEMP AA part receipt time
C         5   Radiosonde TEMP BB part receipt time
C         6   Radiosonde TEMP CC part receipt time
C         7   Radiosonde TEMP DD part receipt time
C         8   Radiosonde PILOT AA part receipt time
C         9   Radiosonde PILOT BB part receipt time
C        10   Radiosonde PILOT CC part receipt time
C        11   Radiosonde PILOT DD part receipt time
C      12-62  Reserved for future use
C        63   Missing

      DO L=1,NRCT
         CF8(L) = 105
         OB8(L) = NINT((NINT(RCT(1,L))+NINT(RCT(2,L))/60.) * 100.)
         Q81(L) = 99999
         IF(RCT(3,L).LT.BMISS)  Q81(L) = NINT(RCT(3,L))
         Q82(L) = 99999
         CALL S02UBF(8,L,*9999)
      ENDDO

      CALL UFBINT(LUNIT,RMORE_8,4,1,NRMORE,'SIRC TTSS UALNHR UALNMN')
      RMORE=RMORE_8
      IF(MAX(RMORE(3),RMORE(4)).LT.BMISS)  THEN
         CF8(1) = 104
         OB8(1) = NINT((RMORE(3)+RMORE(4)/60.) * 100.)
         Q81(1) = 99999
         Q82(1) = 99999
         CALL S02UBF(8,1,*9999)
      END IF
      IF(SUBSET.EQ.'NC004005')  THEN
         IF(MAX(PRGP10(1),RGP10(1)).LT.BMISS)  THEN
            RGP10(1) = (NINT(RGP10(1))/GRAV)
cppppp
            if(iprint.eq.1)  print *, ' orig. RGP10 = ',rgp10(1)
cppppp
            IF(MOD(NINT(RGP10(1)),10).NE.0)  RGP10(1) =
     $       INT(RGP10(1)/10.) * 10
            CF8(1) = 351
            OB8(1) = NINT(RGP10(1))
            Q81(1) = NINT(PRGP10(1)/100.)
            Q82(1) = 3
            CALL S02UBF(8,1,*9999)
         END IF
         IF(RPMSL.LT.BMISS)  THEN
            CF8(1) = 352
            OB8(1) = NINT(RPMSL*.1)
            Q81(1) = 99999
            Q82(1) = 99999
            CALL S02UBF(8,1,*9999)
         END IF
      END IF
      IF(NINT(RAT(1)).LT.100)  THEN
         IF(RMORE(1).LT.BMISS)  THEN
            CF8(1) = 353
            OB8(1) = NINT(RMORE(1))
            Q81(1) = 99999
            Q82(1) = 99999
            CALL S02UBF(8,1,*9999)
         END IF
         IF(RMORE(2).LT.BMISS)  THEN
            CF8(1) = 354
            OB8(1) = NINT(RMORE(2))
            Q81(1) = 99999
            Q82(1) = 99999
            CALL S02UBF(8,1,*9999)
         END IF
      END IF
 
C  PUT THE UNPACKED REPORT INTO OBS
C  --------------------------------
 
      CALL S03UBF(OBS,*9999,*9998,*9997)
 
      RETURN
 9999 CONTINUE
      R03UBF = 999
      RETURN
 9998 CONTINUE
      print *,'IW3UNPBF/R03UBF: RPT with ID= ',SID,' TOSSED - ZERO ',
     $ 'CAT.1-6,51,52 LVLS'
 9997 CONTINUE
      R03UBF = -9999
      KSKUPA =KSKUPA + 1
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION R04UBF(LUNIT,OBS)
C     ---> PROCESSES SURFACE DATA
 
      COMMON/IUBFBB/KNDX,KSKACF(8),KSKUPA,KSKSFC,KSKSAT
      COMMON/IUBFCC/SUBSET
      COMMON/IUBFEE/POB(255),QOB(255),TOB(255),ZOB(255),DOB(255),
     $               SOB(255),VSG(255),CLP(255),CLA(255),OB8(255),
     $               CF8(255)
      COMMON/IUBFGG/PSL,STP,SDR,SSP,STM,DPD,TMX,TMI,HVZ,PRW,PW1,PW2,
     $               CCN,CHN,CTL,CTM,CTH,HCB,CPT,APT,PC6,SND,P24,
     $               DOP,POW,HOW,SWD,SWP,SWH,SST,SPG,SPD,SHC,SAS,WES
      COMMON/IUBFHH/PSQ,SPQ,SWQ,STQ,DDQ
      COMMON/IUBFLL/Q81(255),Q82(255)
 
      CHARACTER*80 HDSTR,RCSTR
      CHARACTER*8  SUBSET,SID,RSV1,RSV2
      INTEGER KKK(0:99),KKKK(49),ICCN(0:100)
      REAL(8) RID,UFBINT_8
      REAL(8) HDR_8(20),RCT_8(5,255),RRSV_8(3),CLDS_8(4,255),
     $ TMXMNM_8(4,255)
      DIMENSION  OBS(*),HDR(20),RCT(5,255),RRSV(3),CLDS(4,255),JTH(0:9),
     $ JTL(0:9),LTL(0:9),TMXMNM(4,255)
      EQUIVALENCE  (RID,SID)

      SAVE
 
      DATA HDSTR/'RPID CLON CLAT HOUR MINU SELV AUTO          '/
      DATA RCSTR/'RCHR RCMI RCTS                              '/
 
      DATA BMISS /10E10  /
      
      DATA JTH/0,1,2,3,4,5,6,8,7,9/,JTL/0,1,5,8,7,2,3,4,6,9/
      DATA LTL/0,1,5,6,7,2,8,4,3,9/
      DATA KKK /5*90,16*91,30*92,49*93/
      DATA KKKK/94,2*95,6*96,10*97,30*98/
      DATA ICCN/0,14*1,20*2,10*3,10*4,10*5,20*6,15*7,8/

      R04UBF = 0
 
      CALL S05UBF
 
C  PUT THE HEADER INFORMATION INTO UNPACKED FORMAT
C  -----------------------------------------------
 
      CALL UFBINT(LUNIT,HDR_8,20,  1,IRET,HDSTR);HDR(2:)=HDR_8(2:)
      CALL UFBINT(LUNIT,RCT_8, 5,255,NRCT,RCSTR);RCT=RCT_8
      IF(HDR(5).GE.BMISS) HDR(5) = 0
      RCTIM = NINT(RCT(1,1))+NINT(RCT(2,1))/60.
      RID = HDR_8(1)
      XOB = HDR(2)
      YOB = HDR(3)
      RHR = BMISS
      IF(HDR(4).LT.BMISS)  RHR = NINT(HDR(4))+NINT(HDR(5))/60.
      ELV = HDR(6)

C  I1 DEFINES SYNOPTIC FORMAT FLAG (SUBSET NC000001, NC000009)
C  I1 DEFINES AUTOMATED STATION TYPE (SUBSET NC000003-NC000008,NC000010)
C  I2 DEFINES CONVERTED HOURLY FLAG (SUBSET NC000xxx)
C  I2 DEFINES SHIP LOCATION FLAG (SUBSET NC001xxx)

      I1  = 9
      I2  = 9
      IF(SUBSET(1:5).EQ.'NC000')  THEN
         IF(SUBSET(6:8).EQ.'001'.OR.SUBSET(6:8).EQ.'009')  THEN
            I1 = 1
            IF(SUBSET(6:8).EQ.'009')  I2 = 1
         ELSE  IF(SUBSET(6:8).NE.'002')  THEN
            IF(HDR(7).LT.15)  THEN
               IF(HDR(7).GT.0.AND.HDR(7).LT.5) THEN
                  I1 = 2
               ELSE  IF(HDR(7).EQ.8) THEN
                  I1 = 3
               ELSE
                  I1 = 4
               END IF
            END IF
         END IF
      END IF
      ITP = (10 * I1) + I2
      RTP = ERTUBF(SUBSET,SID)

C  INDICATOR FOR PRECIPITATION (INCLUDED/EXCLUDED STORED IN BYTE 1 OF
C   HEADER RESERVE CHARACTER WORD 5
C  INDICATOR FOR WIND SPEED (SOURCE/UNITS) STORED IN BYTE 3 OF HEADER
C   RESERVE CHARACTER WORD 5
C  INDICATOR FOR STATION OPERATION/PAST WEATHER DATA STORED IN BYTE 5 OF
C   HEADER RESERVE CHARACTER WORD 5

      CALL UFBINT(LUNIT,RRSV_8,3,1,NRSV,'INPC SUWS ITSO');RRSV=RRSV_8
      RSV1 = '        '
      RSV2 = '        '
      J = -1
      DO  I=1,3
         J = J + 2
         IF(RRSV(I).LT.BMISS)  WRITE(RSV1(J:J),'(I1)') NINT(RRSV(I))
      ENDDO

C  READ THE CATEGORY 51 SURFACE DATA FROM BUFR
C  -------------------------------------------
 
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'PMSL');PSL=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'PRES');STP=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'WDIR');SDR=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'WSPD');SSP=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TMDB');STM=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TMDP');DPD=UFBINT_8
      IF(SUBSET.NE.'NC000007')  THEN
         CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'MXTM');TMX=UFBINT_8
         CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'MITM');TMI=UFBINT_8
      ELSE
         TMX = BMISS
         TMI = BMISS
      END IF
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'QMPR');QSL=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'QMPR');QSP=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'QMWN');QMW=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'QMAT');QMT=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'QMDD');QMD=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'HOVI');HVZ=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'PRWE');PRW=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'PSW1');PW1=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'PSW2');PW2=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TOCC');CCN=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'CHPT');CPT=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'3HPC');APT=UFBINT_8
      IF(MAX(APT,CPT).GE.BMISS) THEN
         APT = BMISS
         CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'24PC');APT24=UFBINT_8
         IF(APT24.LT.BMISS)  THEN
            APT = APT24
            CPT = BMISS
         END IF
      END IF
      
 
C  READ THE CATEGORY 52 SURFACE DATA FROM BUFR
C  -------------------------------------------
 
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TP06');PC6=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TOSD');SND=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TP24');P24=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TOPC');PTO=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'.DTHTOPC');DOP=UFBINT_8
      IF(PTO.LT.BMISS)  THEN
         IF(PC6.GE.BMISS.AND.NINT(DOP).EQ. 6)  PC6 = PTO
cppppp
         IF(PC6.GE.BMISS.AND.NINT(DOP).EQ. 6)
     $    print *, '~~IW3UNPBF/R04UBF: PTO used for PC6 since latter ',
     $    'missing &  6-hr DOP'
cppppp
         IF(P24.GE.BMISS.AND.NINT(DOP).EQ.24)  P24 = PTO
cppppp
         IF(P24.GE.BMISS.AND.NINT(DOP).EQ.24)
     $    print *, '~~IW3UNPBF/R04UBF: PTO used for P24 since latter ',
     $    'missing & 24-hr DOP'
cppppp
      END IF
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'POWW');POW=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'HOWW');HOW=UFBINT_8
      IF(SUBSET(1:5).EQ.'NC001')  THEN
         IF(MIN(POW,HOW).GE.BMISS)  THEN
            CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'POWV');POW=UFBINT_8
            CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'HOWV');HOW=UFBINT_8
         END IF
      END IF
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'DOSW');SWD=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'POSW');SWP=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'HOSW');SWH=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'SST1');SST=UFBINT_8
      IF(SST.GE.BMISS) THEN
          CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'STMP');SST=UFBINT_8
      ENDIF
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'????');SPG=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'????');SPD=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TDMP');SHC=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'ASMP');SAS=UFBINT_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'????');WES=UFBINT_8
      I52FLG = 0
      IF(MIN(SND,P24,POW,HOW,SWD,SWP,SWH,SST,SPG,SPD,SHC,SAS,WES)
     $ .GE.BMISS.AND.(PC6.EQ.0..OR.PC6.GE.BMISS))  I52FLG= 1
 
C  SOME CLOUD DATA IS NEEDED FOR LOW, MIDDLE, AND HIGH CLOUDS IN CAT. 51
C  ---------------------------------------------------------------------
 
      CALL UFBINT(LUNIT,CLDS_8,4,255,NCLD,'VSSO CLAM CLTP HOCB')
      CLDS=CLDS_8
      CTH = -9999.
      CTM = -9999.
      CTL = -9999.
      CHH = BMISS
      CHM = BMISS
      CHL = BMISS
      IF(NCLD.EQ.0)  THEN
         CCM = BMISS
         CCL = BMISS
      ELSE
         CCM = 0.
         CCL = 0.
         DO L=1,NCLD
            VSS = CLDS(1,L)
            CAM = CLDS(2,L)
            CTP = CLDS(3,L)
            CHT = CLDS(4,L)
            IF(CHT.LT.BMISS)  CHT = CHT * 3.2808
            IF(NINT(VSS).EQ.0)  THEN
               IF(NINT(CTP).GT.9.AND.NINT(CTP).LT.20)  THEN
                  ITH = MOD(NINT(CTP),10)
                  KTH = JTH(ITH)
                  CTH = MAX(KTH,NINT(CTH))
                  CHH = MIN(NINT(CHT),NINT(CHH))
               ELSE  IF(NINT(CTP).LT.30)  THEN
                  ITM = MOD(NINT(CTP),10)
                  CTM = MAX(ITM,NINT(CTM))
                  IF(ITM.EQ.0)  CAM = 0.
                  CCM = MAX(NINT(CAM),NINT(CCM))
                  CHM = MIN(NINT(CHT),NINT(CHM))
               ELSE  IF(NINT(CTP).LT.40)  THEN
                  ITL = MOD(NINT(CTP),10)
                  KTL = JTL(ITL)
                  CTL = MAX(KTL,NINT(CTL))
                  IF(ITL.EQ.0)  CAM = 0
                  CCL = MAX(NINT(CAM),NINT(CCL))
                  CHL = MIN(NINT(CHT),NINT(CHL))
               ELSE  IF(NINT(CTP).EQ.59)  THEN
                  CTH = 10
                  CTM = 10
                  IF(CCM.EQ.0.)  CCM = 15
                  CTL = 10
                  IF(CCL.EQ.0.)  CCL = 15
               ELSE  IF(NINT(CTP).EQ.60)  THEN
                  CTH = 10
               ELSE  IF(NINT(CTP).EQ.61)  THEN
                  CTM = 10
                  IF(CCM.EQ.0.)  CCM = 15
               ELSE  IF(NINT(CTP).EQ.62)  THEN
                  CTL = 10
                  IF(CCL.EQ.0.)  CCL = 15
               END IF
            END IF
         ENDDO
      END IF
      IF(NINT(CTH).GT.-1.AND.NINT(CTH).LT.10)  THEN
         CTH = JTH(NINT(CTH))
      ELSE  IF(NINT(CTH).NE.10)  THEN
         CTH = BMISS
      END IF
      IF(NINT(CTM).LT.0.OR.NINT(CTM).GT.10)  THEN
         CTM = BMISS
         CCM = BMISS
      END IF
      IF(NINT(CTL).GT.-1.AND.NINT(CTL).LT.10)  THEN
         CTL = LTL(NINT(CTL))
      ELSE  IF(NINT(CTL).NE.10)  THEN
         CTL = BMISS
         CCL = BMISS
      END IF
 
C  CALL FUNCTIONS TO TRANSFORM TO SPECIFIED UNPACKED FORMAT
C  --------------------------------------------------------

      PSQ = EQSUBF(QSL)
      SPQ = EQSUBF(QSP)
      SWQ = EQSUBF(QMW)
      STQ = EQSUBF(QMT)
      DDQ = EQSUBF(QMD)
 
      IF(PSL.LT.BMISS)  THEN
         PSL = NINT(PSL*.1)
      ELSE
         PSQ = 2
       END IF
      IF(SUBSET(1:5).EQ.'NC001'.AND.NINT(QSL).EQ.4)  STP = BMISS
      IF(STP.LT.BMISS)  THEN
         STP = NINT(STP*.1)
      ELSE
         SPQ = 2
      END IF
      SSP = MIN(SSP*10.,BMISS)
      IF(NINT(SDR).EQ.0)   SDR = 360
      IF(SDR.GE.BMISS.AND.NINT(SSP).EQ.0)   SDR = 360
      IF(MAX(SDR,SSP).GE.BMISS)  SWQ = 2
      IF(MAX(DPD,STM).LT.BMISS)  THEN
         DPD = (STM-DPD)*10.
      ELSE
         DPD = BMISS
C - note this was not set before  & could make a bufr chg in checking
         DDQ = 2
      END IF
      IF(STM.LT.BMISS)  THEN
         ISTM = NINT(STM*100.)
         STM = NINT((ISTM-27315)*0.1)
      ELSE
         STQ = 2
      END IF
      ITMX = NINT(TMX*100.)
      IF(TMX.LT.BMISS)  TMX = NINT((ITMX-27315)*0.1)
      ITMI = NINT(TMI*100.)
      IF(TMI.LT.BMISS)  TMI = NINT((ITMI-27315)*0.1)

      IF(SUBSET(1:5).EQ.'NC000'.OR.SUBSET.EQ.'NC001004')  THEN
         IF(HVZ.GE.BMISS.OR.HVZ.LT.0.) THEN
            HVZ = BMISS
         ELSE  IF(NINT(HVZ).LT.6000)  THEN
            HVZ = MIN(INT(NINT(HVZ)/100),50)
         ELSE  IF(NINT(HVZ).LT.30000)  THEN
            HVZ = INT(NINT(HVZ)/1000) + 50
         ELSE  IF(NINT(HVZ).LE.70000)  THEN
            HVZ = INT(NINT(HVZ)/5000) + 74
         ELSE
            HVZ = 89
         END IF
      ELSE
         IF(HVZ.GE.BMISS.OR.HVZ.LT.0.) THEN
            HVZ = BMISS
         ELSE  IF(NINT(HVZ).LT.1000)  THEN
            KK = MIN(INT(NINT(HVZ)/10),99)
            HVZ = KKK(KK)
         ELSE  IF(NINT(HVZ).LT.50000)  THEN
            KK = MIN(INT(NINT(HVZ)/1000),49)
            HVZ = KKKK(KK)
         ELSE
            HVZ = 99
         END IF
      END IF

      IF(PRW.LT.BMISS)  PRW = NINT(MOD(PRW,100.))
      IF(PW1.LT.BMISS)  PW1 = NINT(MOD(PW1,10.))
      IF(PW2.LT.BMISS)  PW2 = NINT(MOD(PW2,10.))

      IF(NINT(CCN).GT.-1.AND.NINT(CCN).LT.101)  THEN
         CCN = ICCN(NINT(CCN))
      ELSE
         CCN = BMISS
      END IF

      CHN = CCL
      IF(NINT(CHN).EQ.0)  CHN = CCM
      IF(NINT(CHN).GT.9)  THEN
         IF(NINT(CHN).EQ.10)  THEN
            CHN = 9.
         ELSE  IF(NINT(CHN).EQ.15)  THEN
            CHN = 10.
         ELSE
            CHN = BMISS
         END IF
      END IF

      IF(NINT(MAX(CTL,CTM,CTH)).EQ.0)  THEN
         HCB = 9
      ELSE
         HCB = BMISS
         IF(CHH.LT.BMISS) HCB = CHH
         IF(CHM.LT.BMISS) HCB = CHM
         IF(CHL.LT.BMISS) HCB = CHL
         IF(HCB.GE.0.AND.HCB.LT.BMISS)  THEN
            IF(HCB.LT. 150)  THEN
               HCB = 0
            ELSE  IF(HCB.LT. 350)  THEN
               HCB = 1
            ELSE  IF(HCB.LT. 650)  THEN
               HCB = 2
            ELSE  IF(HCB.LT. 950)  THEN
               HCB = 3
            ELSE  IF(HCB.LT.1950)  THEN
               HCB = 4
            ELSE  IF(HCB.LT.3250)  THEN
               HCB = 5
            ELSE  IF(HCB.LT.4950)  THEN
               HCB = 6
            ELSE  IF(HCB.LT.6750)  THEN
               HCB = 7
            ELSE  IF(HCB.LT.8250)  THEN
               HCB = 8
            ELSE
               HCB = 9
            END IF
         END IF
      END IF

      IF(NINT(CPT).LT.0.OR.NINT(CPT).GT.8)  CPT = BMISS

      IF(APT.LT.BMISS)  APT = NINT(APT*.1)
 
      IF(PC6.LT.0.) THEN
         PC6 = 9998
      ELSE  IF(PC6.LT.BMISS) THEN
         PC6 = NINT(PC6*3.937)
      END IF

      IF(SND.LT.0.)  THEN
         SND = 998
      ELSE  IF(SND.LT.BMISS)  THEN
         SND = NINT(SND*39.37)
      END IF
 
      IF(P24.LT.0.) THEN
         P24 = 9998
      ELSE  IF(P24.LT.BMISS) THEN
         P24 = NINT(P24*3.937)
      END IF

      DOP = BMISS
      IF(PC6.LT.BMISS)  DOP = 1

      POW = NINT(POW)
      HOW = NINT(MIN(2.*HOW,BMISS))

      IF(NINT(SWD).EQ.0)  THEN
         SWD = 0
      ELSE  IF(SWD.LT.5)  THEN
         SWD = 36
      ELSE  IF(SWD.LT.BMISS)  THEN
         SWD = NINT((SWD+.001)*.1)
      END IF

      SWP = NINT(SWP)

      SWH = NINT(MIN(2.*SWH,BMISS))

      ISST = NINT(SST*100.)
      IF(SST.LT.BMISS)  SST = NINT((ISST-27315)*0.1)

      SHC = NINT(SHC)
      IF(NINT(SHC).LT.0.OR.NINT(SHC).GT.8)  SHC = BMISS

      SAS = NINT(SAS)
      IF(NINT(SAS).LT.0.OR.NINT(SAS).GT.9)  SAS = BMISS

C  MAKE THE UNPACKED REPORT INTO OBS
C  ---------------------------------
 
      RSTP = 99999
      CALL S01UBF(SID,XOB,YOB,RHR,RCTIM,RSV1,RSV2,ELV,ITP,RTP,RSTP)
      CALL S02UBF(51,1,*9999)
      IF(I52FLG.EQ.0)  CALL S02UBF(52,1,*9999)

C  ------------------------------------------------------------------
C  MISC DATA GOES INTO CATEGORY 08
C  ------------------------------------------------------------------
C  CODE FIGURE 020 - ALTIMETER SETTING IN 0.1*MB
C  CODE FIGURE 081 - CALENDAR DAY MAXIMUM TEMPERATURE
C  CODE FIGURE 082 - CALENDAR DAY MINIMUM TEMPERATURE
C  CODE FIGURE 083 - SIX HOUR MAXIMUM TEMPERATURE
C  CODE FIGURE 084 - SIX HOUR MINIMUM TEMPERATURE
C  CODE FIGURE 085 - PRECIPITATION OVER PAST HOUR IN 0.01*INCHES
C  CODE FIGURE 098 - DURATION OF SUNSHINE FOR CALENDAR DAY IN MINUTES
C  ------------------------------------------------------------------

      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'ALSE');ALS=UFBINT_8
      IF(ALS.LT.BMISS) THEN
         OB8(1) =  NINT(ALS*.1)
         CF8(1) = 20 
         Q81(1) = 99999
         Q82(1) = 99999
         CALL S02UBF(8,1,*9999)
      END IF
      IF(SUBSET.EQ.'NC000007')  THEN
         CALL UFBINT(LUNIT,TMXMNM_8,4,255,NTXM,
     $    '.DTHMXTM MXTM .DTHMITM MITM');TMXMNM=TMXMNM_8
         IF(NTXM.GT.0)  THEN
            DO I = 1,NTXM
               DO J = 1,3,2
                  IF(NINT(TMXMNM(J,I)).EQ.24) THEN
                     IF(TMXMNM(J+1,I).LT.BMISS)  THEN
                        ITMX = NINT(TMXMNM(J+1,I)*100.)
                        TMX = NINT((ITMX-27315)*0.1)
                        IF(TMX.LT.0)  THEN
                           OB8(1) = 1000 + ABS(NINT(TMX))
                        ELSE 
                           OB8(1) = NINT(TMX)
                        END IF
                        CF8(1) = 81 + INT(J/2)
                        Q81(1) = 99999
                        Q82(1) = 99999
                        CALL S02UBF(8,1,*9999)
                     END IF
                  ELSE  IF(NINT(TMXMNM(J,I)).EQ.6) THEN
                     IF(TMXMNM(J+1,I).LT.BMISS)  THEN
                        ITMX = NINT(TMXMNM(J+1,I)*100.)
                        TMX = NINT((ITMX-27315)*0.1)
                        IF(TMX.LT.0)  THEN
                           OB8(1) = 1000 + ABS(NINT(TMX))
                        ELSE 
                           OB8(1) = NINT(TMX)
                        END IF
                        CF8(1) = 83 + INT(J/2)
                        Q81(1) = 99999
                        Q82(1) = 99999
                        CALL S02UBF(8,1,*9999)
                     END IF
                  END IF
               ENDDO
            ENDDO
         END IF
      END IF
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TP01');PC1=UFBINT_8
      IF(PC1.LT.10000) THEN
         IF(PC1.GE.0.) THEN
            OB8(1) = NINT(PC1*3.937)
         ELSE
            OB8(1) = 9998
         END IF
         CF8(1) = 85
         Q81(1) = 99999
         Q82(1) = 99999
         CALL S02UBF(8,1,*9999)
      END IF
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'TOSS');DUS=UFBINT_8
      IF(NINT(DUS).LT.1000) THEN
         OB8(1) = NINT(98000. + DUS)
         CF8(1) = 98 
         Q81(1) = 99999
         Q82(1) = 99999
         CALL S02UBF(8,1,*9999)
      END IF
 
      CALL S03UBF(OBS,*9999,*9998,*9997)
 
      RETURN

 9999 CONTINUE
      R04UBF = 999
      RETURN

 9998 CONTINUE
      print *,'IW3UNPBF/R04UBF: RPT with ID= ',SID,' TOSSED - ZERO ',
     $ 'CAT.1-6,51,52 LVLS'
 9997 CONTINUE
      R04UBF = -9999
      KSKSFC =KSKSFC + 1
      RETURN

      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION R05UBF(LUNIT,OBS)
C     ---> PROCESSES AIRCRAFT DATA
 
      COMMON/IUBFBB/KNDX,KSKACF(8),KSKUPA,KSKSFC,KSKSAT
      COMMON/IUBFCC/SUBSET
      COMMON/IUBFEE/POB(255),QOB(255),TOB(255),ZOB(255),DOB(255),
     $               SOB(255),VSG(255),CLP(255),CLA(255),OB8(255),
     $               CF8(255)
      COMMON/IUBFFF/PQM(255),QQM(255),TQM(255),ZQM(255),WQM(255),
     $               QCP(255),QCA(255)
      COMMON/IUBFLL/Q81(255),Q82(255)
 
      CHARACTER*80 HDSTR,LVSTR,QMSTR,RCSTR,CRAWR
      CHARACTER*8  SUBSET,SID,RSV1,RSV2,CCL,CRAW(1,255)
      REAL(8) RID,RCL,UFBINT_8
      REAL(8) HDR_8(20),RCT_8(5,255),ARR_8(10,255),RAW_8(1,255)
      DIMENSION    OBS(*),HDR(20),RCT(5,255),ARR(10,255),RAW(1,255)
      EQUIVALENCE  (RID,SID),(RCL,CCL),(RAW_8,CRAW)

      SAVE
 
      DATA HDSTR/'RPID CLON CLAT HOUR MINU SECO               '/
      DATA LVSTR/'PRLC TMDP TMDB WDIR WSPD                    '/
      DATA QMSTR/'QMPR QMAT QMDD QMGP QMWN                    '/
      DATA RCSTR/'RCHR RCMI RCTS                              '/
 
      DATA BMISS /10E10  /
 
      HGTF(P) = (1.-(P/1013.25)**(1./5.256))*(288.15/.0065)

      R05UBF = 0
 
      CALL S05UBF
 
C  PUT THE HEADER INFORMATION INTO UNPACKED FORMAT
C  -----------------------------------------------
 
      CALL UFBINT(LUNIT,HDR_8,20,  1,IRET,HDSTR);HDR(2:)=HDR_8(2:)
      IF(IRET.EQ.0)  SID = '        '
      CALL UFBINT(LUNIT,RCT_8, 5,255,NRCT,RCSTR);RCT=RCT_8
      IF(HDR(5).GE.BMISS) HDR(5) = 0
      IF(HDR(6).GE.BMISS) HDR(6) = 0
      RCTIM = NINT(RCT(1,1))+NINT(RCT(2,1))/60.
      RID = HDR_8(1)
      XOB = HDR(2)
      YOB = HDR(3)
      RHR = BMISS
      IF(HDR(4).LT.BMISS) RHR = (NINT(HDR(4)) + ((NINT(HDR(5)) * 60.) +
     $  NINT(HDR(6)))/3600.) + 0.0000000001
 
C  TRY TO FIND FIND THE FLIGHT LEVEL HEIGHT
C  ----------------------------------------
 
      CALL UFBINT(LUNIT,HDR_8,20,1,IRET,'PSAL FLVL IALT HMSL PRLC')
      HDR=HDR_8
      IF(HDR(4).LT.BMISS)  THEN
         ELEV = HDR(4)
      ELSE  IF(HDR(5).LT.BMISS)  THEN
         ELEV = HGTF(HDR(5)*.01)
      ELSE
         ELEV = BMISS
      END IF

C FOR MDCARS ACARS DATA ONLY:
C  UNCOMMENTING THE 2 LINES BELOW WILL SET P-ALT TO RPTD "IALT" VALUE --
C    IN THIS CASE, PREPDATA WILL LATER GET PRESS. VIA STD. ATMOS. FCN.
C  COMMENTING THE 2 LINES BELOW WILL USE RPTD PRESSURE "PRLC" TO GET
C    P-ALT VIA INVERSE STD. ATMOS. FCN. -- IN THIS CASE, PREPDATA WILL
C    LATER RETURN THIS SAME PRESS. VIA STD. ATMOS. FCN.
C  (this may not be correct now that we are storing reported pressure
C    which is later read by prepdata program - ???)

      IF(HDR(1).LT.BMISS) THEN
         ELEV = HDR(1) + SIGN(0.0000001,HDR(1))
      ELSE  IF(HDR(2).LT.BMISS)  THEN
         ELEV = HDR(2) + SIGN(0.0000001,HDR(2))
cdak  ELSE  IF(HDR(3).LT.BMISS)  THEN
cdak     ELEV = HDR(3)
      END IF

      ELV = ELEV

C  ACFT NAVIGATION SYSTEM STORED IN INSTR. TYPE LOCATION (AS WITH ON29)
C  --------------------------------------------------------------------

      ITP = 99
      CALL UFBINT(LUNIT,RNS_8,1,1,IRET,'ACNS');RNS=RNS_8
      IF(RNS.LT.BMISS)  THEN
         IF(NINT(RNS).EQ.0)  THEN
            ITP = 97
         ELSE  IF(NINT(RNS).EQ.1)  THEN
            ITP = 98
         END IF
      END IF

      RTP = 41
      READ(SUBSET(8:8),'(F1.0)') RSTP

      RSV1 = '        '
      RSV2 = '        '

C  ICAO LOCATION ID STORED IN HEADER RESERVE CHARACTER WORD 5
C  ----------------------------------------------------------

      CALL UFBINT(LUNIT,RCL_8,1,1,IRET,'ICLI');RCL=RCL_8
      IF(IRET.NE.0)  RSV1 = CCL

      POF = BMISS
      PCT = BMISS

      IF(SUBSET.EQ.'NC004003')  THEN
 
C  ------------------------------------
C  ASDAR/AMDAR AIRCRAFT TYPE COME HERE
C  ------------------------------------

C  GET PHASE OF FLIGHT FOR LATER STORAGE INTO CAT. 8
C  -------------------------------------------------

         CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'POAF');POF=UFBINT_8

C  GET TEMPERATURE PRECISION FOR LATER STORAGE INTO CAT. 8
C  -------------------------------------------------------

         CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'PCAT');PCT=UFBINT_8

C  CARSWELL (NEVER HAPPENS) INDICATOR STORED IN BYTE 1 OF HEADER
C  RESERVE CHARACTER WORD 6
C   (NOTE: NAS9000 ONLY ASSIGNED HEADER "KAWN" AS CARSWELL, ALTHOUGH
C          "PHWR" AND "EGWR" ARE ALSO APPARENTLY ALSO CARSWELL)
C  ------------------------------------------------------------------

         IF(RSV1(1:4).EQ.'KAWN')  RSV2(1:1) = 'C'

      ELSE IF(SUBSET.EQ.'NC004004') THEN
 
C  ------------------------------
C  ACARS AIRCRAFT TYPE COME HERE
C  ------------------------------
 
         CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'ACRN');RID=UFBINT_8
         IF(IRET.EQ.0)  SID = 'ACARS   '
         KNDX = KNDX + 1

      ELSE IF(SUBSET.EQ.'NC004001'.OR.SUBSET.EQ.'NC004002') THEN
 
C  -----------------------------------------
C  AIREP AND PIREP AIRCRAFT TYPES COME HERE
C  -----------------------------------------

C  MAY POSSIBLY NEED TO MODIFY THE RPID HERE
C  (for later tests in prepacqc program)
C  -----------------------------------------

cdak no more!IF(SID(6:6).EQ.'Z')  SID(6:6) = 'X'
         IF(SID.EQ.'A       '.OR.SID.EQ.'        '.OR.SID(1:3).EQ.'ARP'
     $    .OR.SID(1:3).EQ.'ARS')  SID = 'AIRCFT  '

cvvvvv temporary?
C  Determined that Hickum AFB reports are much like Carswell - they have
C   problems!  They also are usually duplicates of either Carswell or
C   non-Carswell reports.  Apparently the front-end processing filters
C   them out (according to B. Ballish).  So, to make things match,
C   we will do the same here.
C   ACTUALLY, JEFF ATOR HAS REMOVED THESE FROM THE DECODER, SO WE
C    SHOULD NEVER EVEN SEE THEM IN THE DATABASE, but it won't hurt
C    anything to keep this in here.
C   (NOTE: These all have headers of "PHWR")

         if(RSV1(1:4).eq.'PHWR')  then
cppppp
cdak  print *, 'IW3UNPBF/R05UBF: TOSS "PHWR" AIREP with ID = ',SID,
cdak $ '; CCL = ',RSV1(1:4)
cppppp
            R05UBF = -9999
            kskacf(8) = kskacf(8) + 1
            return
         end if
caaaaa temporary?

cvvvvv temporary?
C        1) Carswell/Tinker AMDARS are processed as AIREP subtypes.
C      Nearly all of them are duplicated as true non-Carswell AMDARS in
C      the AMDAR subtype.  The earlier version of the aircraft dup-
C      checker could not remove such duplicates; the new verison now
C      in operations can remove these. SO, WE HAVE COMMENTED THIS OUT.
C
C      The Carswell AMDARS can be identified by the string " Sxyz" in
C      the raw report (beyond byte 40), where y is 0,1, or 2.
C      (NOTE: Apparently Carswell here applies to more headers than
C             just "KAWN", so report header is not even checked.)

C        2) Carswell/Tinker ACARS are processed as AIREP subtypes.
C      These MAY duplicate true non-Carswell ACARS in the ACARS
C      subtype.  The NAS9000 decoder always excluded this type (no
C      dup-checking was done).  All of these will be removed here.
C      The Carswell ACARS can be identified by the string " Sxyz" in
C      the raw report (beyond byte 40), where y is 3 or greater.
C      (NOTE: Apparently Carswell here applies to more headers than
C             just "KAWN", so report header is not even checked.)

         call ufbint(lunit,raw_8,1,255,nlev,'RRSTG');raw=raw_8
         if(nlev.gt.5)  then
            ni = -7
            do mm = 6,nlev
               ni = ni + 8
               crawr(ni:ni+7) = craw(1,mm)
               if(ni+8.gt.80)  go to 556
            enddo
  556       continue
            do mm = 1,ni+7
               if(crawr(mm:mm+1).eq.' S')  then
                  if((crawr(mm+2:mm+2).ge.'0'.and.crawr(mm+2:mm+2).le.
     $             '9').or.crawr(mm+2:mm+2).eq.'/')  then
                     if((crawr(mm+3:mm+3).ge.'0'.and.crawr(mm+3:mm+3)
     $                .le.'9').or.crawr(mm+3:mm+3).eq.'/')  then
                        if((crawr(mm+4:mm+4).ge.'0'.and.
     $                   crawr(mm+4:mm+4).le.'9').or.crawr(mm+4:mm+4)
     $                   .eq.'/')  then
cppppp
cdak  print *, 'IW3UNPBF/R05UBF: For ',SID,', raw_8(',ni+7,') = ',
cdak $ crawr(1:ni+7)
cppppp
                           if(crawr(mm+3:mm+3).lt.'3')  then

C  THIS IS A CARSWELL/TINKER AMDAR REPORT --> THROW OUT
C   (NOT ANYMORE, DUP-CHECKER IS HANDLING THESE OKAY NOW)
C  ----------------------------------------------------

cppppp
cdak  print *, 'IW3UNPBF/R05UBF: Found a Carswell AMDAR for ',SID,
cdak $ '; CCL = ',RSV1(1:4)
cppppp
cdak  R05UBF = -9999
cdak  KSKACF(3) = KSKACF(3) + 1
cdak  RETURN
                           else

C  THIS IS A CARSWELL/TINKER ACARS REPORT --> THROW OUT
C  ----------------------------------------------------

cppppp
cdak  print *, 'IW3UNPBF/R05UBF: Found a Carswell ACARS for ',SID,
cdak $ '; CCL = ',RSV1(1:4)
cppppp
      R05UBF = -9999
      KSKACF(4) = KSKACF(4) + 1
      RETURN

                           end if
                        end if
                     end if
                  end if
               end if
               if(mm+5.gt.ni+7)  go to 557
            enddo
  557       continue
         END IF
caaaaa temporary?

C  CARSWELL INDICATOR STORED IN BYTE 1 OF HEADER RESERVE CHARACTER WRD 6
C   (NOTE: NAS9000 ONLY ASSIGNED HEADER "KAWN" AS CARSWELL, ALTHOUGH
C          "PHWR" AND "EGWR" ARE ALSO APPARENTLY ALSO CARSWELL)
C  ------------------------------------------------------------------

         IF(RSV1(1:4).EQ.'KAWN')  RSV2(1:1) = 'C'

      END IF
 
C  -----------------------------
C  ALL AIRCRAFT TYPES COME HERE
C  -----------------------------

      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'DGOT');DGT=UFBINT_8

C  PUT THE LEVEL DATA INTO SPECIFIED UNPACKED FORMAT
C  -------------------------------------------------
 
      CALL UFBINT(LUNIT,ARR_8,10,255,NLEV,LVSTR);ARR=ARR_8

      DO L=1,NLEV

         POB(L) = BMISS
         IF(ARR(1,L).LT.BMISS)  POB(L) =  NINT(ARR(1,L)*.1)

         QOB(L) = BMISS
         IF(MAX(ARR(2,L),ARR(3,L)).LT.BMISS)
     $    QOB(L) = (ARR(3,L)-ARR(2,L))*10.

         TOB(L) = BMISS
         ITMP = NINT(ARR(3,L)*100.)
         IF(ARR(3,L).LT.BMISS)  TOB(L) = NINT((ITMP-27315)*0.1)

         ZOB(L) = ELEV
         DOB(L) = ARR(4,L)
         SOB(L) = MIN(ARR(5,L)*10.,BMISS)
      ENDDO

      CALL UFBINT(LUNIT,ARR_8,10,255,NLEV,QMSTR);ARR=ARR_8

      IF(SUBSET.EQ.'NC004004') THEN
 
C  ---------------------------------------------------------
C  ACARS AIRCRAFT TYPE COME HERE FOR QUALITY MARK ASSIGNMENT
C  ---------------------------------------------------------
 
         DO L=1,NLEV
            PQM(L) = EQMUBF(ARR(1,L))
            TQM(L) = EQMUBF(ARR(2,L))
            QQM(L) = EQMUBF(ARR(3,L))
            ZQM(L) = EQMUBF(ARR(4,L))
            WQM(L) = EQMUBF(ARR(5,L))
         ENDDO

C  DEFAULT Q.MARK FOR WIND: "1"
C  (note: This should be moved to prepdata)
C  ----------------------------------------

         IF(NLEV.EQ.0.OR.ARR(5,1).GE.BMISS)  WQM(1) = 1

      ELSE

C  --------------------------------------------------------------
C  ALL OTHER AIRCRAFT TYPES COME HERE FOR QUALITY MARK ASSIGNMENT
C   NOTE: This is out of date and needs to be fixed to look like
C         ACARS above (except no default q.m. on wind)!!!
C  (note: This should be moved to prepdata)
C  --------------------------------------------------------------
 
         DO L=1,NLEV

C  DEFAULT Q.MARK FOR ALL VARIABLES IS 2
C  -------------------------------------

            PQM(L) =  2
            TQM(L) =  2
            QQM(L) =  2
            ZQM(L) =  2
            WQM(L) =  2

C  IF KEEP  FLAG ON WIND ( 0), ALL Q. MARKS GET KEEP  FLAG ( 0)
C   -- unless....
C  IF PURGE FLAG ON WIND (14), ALL Q. MARKS GET PURGE FLAG (14)
C  IF PURGE FLAG ON TEMP (14), ALL Q. MARKS GET PURGE FLAG (14)
C  -----------------------------------------------------------------

            IF(ARR(5,L).EQ.0.AND.(ARR(2,L).LT.10.OR.ARR(2,L).GT.15))THEN
               PQM(L) =  0
               TQM(L) =  0
               QQM(L) =  0
               ZQM(L) =  0
               WQM(L) =  0
            ELSE  IF(ARR(5,L).EQ.14.OR.ARR(2,L).EQ.14)  THEN
               PQM(L) = 14
               TQM(L) = 14
               QQM(L) = 14
               ZQM(L) = 14
               WQM(L) = 14
            END IF
         ENDDO
      END IF
 
C  PUT THE UNPACKED REPORT INTO OBS
C  --------------------------------
 
      CALL S01UBF(SID,XOB,YOB,RHR,RCTIM,RSV1,RSV2,ELV,ITP,RTP,RSTP)
      CALL S02UBF(6,1,*9999)
 
C  ------------------------------------------------------------------
C  MISC DATA GOES INTO CATEGORY 08
C  ------------------------------------------------------------------
C  CODE FIGURE 021 - REPORT SEQUENCE NUMBER
C  CODE FIGURE 914 - PHASE OF FLIGHT
C                     (CURRENTLY ONLY FOR AMDAR/ASDAR)
C  CODE FIGURE 915 - TEMPERATURE PRECISION
C                     (CURRENTLY ONLY FOR AMDAR/ASDAR)
C  CODE FIGURE 916 - DEGREE OF TURBULENCE INDICATOR
C  ------------------------------------------------------------------

      IF(SUBSET.EQ.'NC004004') THEN
         OB8(1) = KNDX
         CF8(1) = 21
         Q81(1) = 99999
         Q82(1) = 99999
         CALL S02UBF(8,1,*9999)
      ELSE IF(SUBSET.EQ.'NC004003') THEN
         IF(POF.LT.BMISS) THEN
            OB8(1) = NINT(POF)
            CF8(1) = 914
            Q81(1) = 99999
            Q82(1) = 99999
            CALL S02UBF(8,1,*9999)
         END IF
         IF(PCT.LT.BMISS) THEN
            OB8(1) = NINT(PCT)
            CF8(1) = 915
            Q81(1) = 99999
            Q82(1) = 99999
            CALL S02UBF(8,1,*9999)
         END IF
      END IF
      IF(NINT(DGT).LT.16.) THEN
         OB8(3) = NINT(DGT)
         CF8(3) = 916
         Q81(3) = 99999
         Q82(3) = 99999
         CALL S02UBF(8,3,*9999)
      END IF

      CALL S03UBF(OBS,*9999,*9998,*9997)
 
      RETURN

 9999 CONTINUE
      R05UBF = 999
      RETURN

 9998 CONTINUE
      print *,'IW3UNPBF/R05UBF: RPT with ID= ',SID,' TOSSED - ZERO ',
     $ 'CAT.1-6,51,52 LVLS'
 9997 CONTINUE
      R05UBF = -9999
      KSKACF(1) = KSKACF(1) + 1
      RETURN

      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      FUNCTION R06UBF(LUNIT,OBS)
C     ---> PROCESSES SATWIND DATA
 
      COMMON/IUBFBB/KNDX,KSKACF(8),KSKUPA,KSKSFC,KSKSAT
      COMMON/IUBFCC/SUBSET
      COMMON/IUBFEE/POB(255),QOB(255),TOB(255),ZOB(255),DOB(255),
     $               SOB(255),VSG(255),CLP(255),CLA(255),OB8(255),
     $               CF8(255)
      COMMON/IUBFFF/PQM(255),QQM(255),TQM(255),ZQM(255),WQM(255),
     $               QCP(255),QCA(255)
      COMMON/IUBFKK/KOUNT(8,9),KNTSAT(250:260),IFLSAT
 
      CHARACTER*80 HDSTR,LVSTR,QMSTR,RCSTR
      CHARACTER*8  SUBSET,SID,RSV1,RSV2
      CHARACTER*3  CINDX3
      CHARACTER*1  CSAT(8),CPRD(9),CINDX7,C7(26)
      REAL(8) RID,UFBINT_8
      REAL(8) HDR_8(20),RCT_8(5,255),ARR_8(10,255)
      DIMENSION    OBS(*),HDR(20),RCT(5,255),ARR(10,255)
      EQUIVALENCE  (RID,SID)

      integer inum
      integer lityp(69)
      real*8  AQCM

      SAVE
 
      DATA HDSTR/'RPID CLON CLAT HOUR MINU SAID               '/
      DATA LVSTR/'PRLC TMDP TMDB WDIR WSPD                    '/
      DATA QMSTR/'QMPR QMAT QMDD QMGP SWQM MAQC               '/
      DATA RCSTR/'RCHR RCMI RCTS                              '/
 
      DATA BMISS /10E10  /
      DATA CSAT  /'A','B','C','D','X','Z','V','Y'/
      DATA CPRD  /'C','V','I','W','P','T','L','Z','G'/
      DATA C7    /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     $            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/


!  Add translation table for 'instrument type' based on last two
!  digits of subset name.  Undefined instrument type -> 99
!
      DATA lityp /1, 6, 4, 3, 16, 18, 99, 17, 3, 16, 18, 17, 3,
     &      99, 6*99, 1, 6, 4, 3, 16*99, 1, 6, 4, 3, 5*99, 99, 99,
     &      9*99, 1, 6, 4, 64, 65, 66, 3*99/

 
      HGTF(P) = (1.-(P/1013.25)**(1./5.256))*(288.15/.0065)

      R06UBF = 0
      IFLSAT = 1
 
      CALL S05UBF
 
C  TRY TO FIND FIND THE HEIGHT ASSIGNMENT
C  --------------------------------------
 
      CALL UFBINT(LUNIT,HDR_8,20,1,IRET,'HGHT PRLC');HDR=HDR_8
      IF(HDR(1).LT.BMISS)  THEN
         ELEV = HDR(1)
      ELSE  IF(HDR(2).LT.BMISS)  THEN
         ELEV = HGTF(HDR(2)*.01)
      ELSE
         ELEV = BMISS
      END IF
 
C  PUT THE HEADER INFORMATION INTO UNPACKED FORMAT
C  -----------------------------------------------
 
      CALL UFBINT(LUNIT,HDR_8,20,  1,IRET,HDSTR);HDR(2:)=HDR_8(2:)
      CALL UFBINT(LUNIT,RCT_8, 5,255,NRCT,RCSTR);RCT=RCT_8
      IF(HDR(5).GE.BMISS) HDR(5) = 0
      RCTIM = NINT(RCT(1,1))+NINT(RCT(2,1))/60.
      RID = HDR_8(1)
      XOB = HDR(2)
      YOB = HDR(3)
      RHR = BMISS
      IF(HDR(4).LT.BMISS)  RHR = NINT(HDR(4))+NINT(HDR(5))/60.
      RSV1 = '        '
      RSV2 = '        '

C  SATWIND PRODUCER INDICATOR STORED IN BYTE 1 OF HEADER RESERVE
C   CHARACTER WORD 5
C  ------------------------------------------------------------------

C  CLOUD MASK/DEEP LAYER INDICATOR STORED IN BYTE 3 OF HEADER RESERVE
C   CHARACTER WORD 5
C   {=2 - CLOUD TOP (NORMAL CLOUD DRIFT), =1 - DEEP LAYER,
C    =9 - INDICATOR MISSING, THUS REVERTS TO DEFAULT CLOUD TOP}
C   (=9 FOR ALL BUT U.S. HIGH-DENSITY SATWND TYPES)
C  --------------------------------------------------------------------

      IF(SID(1:1).GE.'A'.AND.SID(1:1).LE.'D')  THEN
         CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'SWPR');SWPR=UFBINT_8
         IF(NINT(SWPR).GT.0.AND.NINT(SWPR).LT.10)
     $    WRITE(RSV1(1:1),'(I1)') NINT(SWPR)

C  FOR U.S. SATWNDS: REPROCESS THE STN. ID TO CONFORM WITH OLD ON29
C  APPEARANCE, EXCEPT THAT THE ID WILL BE 8-CHAR RATHER THAN 6-CHAR
C  ----------------------------------------------------------------

C    REPROCESSED CHAR 1 -----> JBUFR CHAR 1
C    REPROCESSED CHAR 2 -----> RETURNED VALUE IN BUFR FOR 'SWPR'
C                               (PRODUCER)
C    REPROCESSED CHAR 3-5 ---> SEQUENTIAL SERIAL INDEX (001 - 999)
C                              (UNIQUE FOR EACH JBUFR CHAR 1/6 COMB.)
C    REPROCESSED CHAR 6 -----> JBUFR CHAR 6
C    REPROCESSED CHAR 7 -----> GROUP NUMBER FOR SERIAL INDEX IN
C                              REPROCESSED CHAR 3-5 (0 - 9, A - Z)
C    REPROCESSED CHAR 8 -----> ALWAYS BLANK (' ') FOR NOW

         DO I = 1,8
            IF(SID(1:1).EQ.CSAT(I))  THEN
               DO J = 1,9
                  IF(SID(6:6).EQ.CPRD(J))  THEN
                     KOUNT(I,J) = KOUNT(I,J) + 1
                     KOUNT(I,J) = MIN(KOUNT(I,J),35999)
                     KOUNT3 = MOD(KOUNT(I,J),1000)
                     KOUNT7 = INT(KOUNT(I,J)/1000)
                     WRITE(CINDX3,'(I3.3)')  KOUNT3
                     IF(KOUNT7.LT.10)  THEN
                        WRITE(CINDX7,'(I1.1)')  KOUNT7
                     ELSE
                        CINDX7 = C7(KOUNT7-9)
                     END IF
                     GO TO 9945
                  END IF
               ENDDO
               GO TO 9944
            END IF
         ENDDO
 9944    CONTINUE
         CINDX3 = '???'
         CINDX7 = '?'
 9945    CONTINUE
         SID = SID(1:1)//RSV1(1:1)//CINDX3//SID(6:6)//CINDX7//' '

         CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'SWDL');SWDL=UFBINT_8
         IF(NINT(SWDL).GT.-1.AND.NINT(SWDL).LT.10)
     $    WRITE(RSV1(3:3),'(I1)') NINT(SWDL)

      ELSE  IF(SID(1:1).EQ.'V')  THEN    ! India 
         RSV1(1:1) = 'E'
      ELSE  IF(SID(1:1).EQ.'Z')  THEN    ! ESA/EUMETSAT
         RSV1(1:1) = 'C'
      ELSE  IF(SID(1:1).EQ.'Y')  THEN    ! Japanese
         RSV1(1:1) = 'D'
      END IF
!
!
!  Note: the above code should be modified so that RSV1 is set according
!  to inum value below (SUBSET(7:8))
!     inum   0 - 19   US  (wisc or Nesdis)            '1'
!     inum   21- 24   India                           'E'
!     inum   41-44    Japan                           'D'
!     inum   61-66    European                        'C'
!     {inum in 50's will be NESDIS processing of GMS data, possibly}
!
!   then we can use RSV1(1:1) to determine producer, independent of
!   station IDs
!
      read(subset(7:8),'(I2)') inum
      if (inum .lt. 20) then
         RSV1(1:1) = '1'
      else if (inum .gt. 20 .and. inum .lt. 25) then
         RSV1(1:1) = 'E'
      else if (inum .gt. 40 .and. inum .lt. 45) then
         RSV1(1:1) = 'D'
      else if (inum .gt. 60 .and. inum .lt. 67) then
         RSV1(1:1) = 'C'
      else
         RSV1(1:1) = ' '
      endif 

C  THE INSTRUMENT TYPE INDICATES THE PRODUCT TYPE
C  ----------------------------------------------

      ITP = 99

      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'SWTP');SWTP=UFBINT_8
      IF(SWTP.LT.BMISS)  THEN
         ITP = NINT(SWTP)
      ELSE

!  Modifications for changes in satwind July 2000 & beyond 
!  Set instrument type based on last 2 digits of subset name

!        read(subset(7:8),'(I2)') inum
         if (inum .gt. 0 .and. inum .lt. 70) then
           itp = lityp(inum)
         endif
      endif

!        IF(SUBSET(8:8).EQ.'1') THEN

C  IR AUTOMATED (CONVENTIONAL)
C  ---------------------------

!           ITP = 1
!        ELSE  IF(SUBSET(8:8).EQ.'2') THEN

C  VISIBLE AUTOMATED (CONVENTIONAL)
C  --------------------------------

!           ITP = 6
!        ELSE  IF(SUBSET(8:8).EQ.'3') THEN

C  WATER VAPOR AUTOMATED (CONVENTIONAL)
C  ------------------------------------

!           ITP = 4
!        ELSE  IF(SUBSET(8:8).EQ.'4') THEN

C  PICTURE TRIPLET (OLD PROCESSING/FORMAT)
C  ---------------------------------------

!           ITP = 3
!        END IF

!     END IF

      ELV = ELEV
      RTP = 63
 
C  PUT THE LEVEL DATA INTO SPECIFIED UNPACKED FORMAT
C  -------------------------------------------------
 
      CALL UFBINT(LUNIT,ARR_8,10,255,NLEV,LVSTR);ARR=ARR_8
      DO L=1,NLEV

         POB(L) = BMISS
         IF(ARR(1,L).LT.BMISS)  POB(L) =  NINT(ARR(1,L)*.1)

C  GROSS CHECK ON PRESSURE
C  -----------------------

         IF(NINT(POB(L)).EQ.0)  THEN
            print *,'~~IW3UNPBF/R06UBF: RPT with ID= ',SID,' TOSSED - ',
     $       'PRES. IS ZERO MB'
            R06UBF = -9999
            KSKSAT = KSKSAT + 1
            RETURN
         END IF

         QOB(L) = BMISS
         IF(MAX(ARR(2,L),ARR(3,L)).LT.BMISS)
     $    QOB(L) = (ARR(3,L)-ARR(2,L))*10.

         TOB(L) = BMISS
         ITMP = NINT(ARR(3,L)*100.)
         IF(ARR(3,L).LT.BMISS)  TOB(L) = NINT((ITMP-27315)*0.1)

         ZOB(L) = ELEV
         DOB(L) = ARR(4,L)
         SOB(L) = MIN(ARR(5,L)*10.,BMISS)
      ENDDO

C  DETERMINE QUALITY MARKERS
C  (later should be moved to prepdata, or let prepbufr just store RFFL)
C  --------------------------------------------------------------------
 
      CALL UFBINT(LUNIT,ARR_8,10,255,NLEV,QMSTR);ARR=ARR_8
      CALL UFBINT(LUNIT,UFBINT_8,1,1,IRET,'RFFL');RFFL=UFBINT_8
      IF(RFFL.LT.BMISS.AND.(NINT(ARR(5,1)).EQ.2.OR.NINT(ARR(5,1)).GE.
     $ BMISS))  THEN
         IF(NINT(RFFL).GT.84)  THEN
            ARR(5,1) = 1
         ELSE  IF(NINT(RFFL).GT.55)  THEN
            ARR(5,1) = 2
         ELSE  IF(NINT(RFFL).GT.49)  THEN
            ARR(5,1) = 3
         ELSE
            ARR(5,1) = 13
         END IF
      END IF

!  New quality code for ELW winds  as of 25 Apr. 2001 (from iw3gad changes)
!
!  Note that the setting of WQM,PQM, etc. below is performed for NLEV
!   levels.  NLEV is set by reading QMSTR above.  If none of the
!   mnemonics listed in QMSTR are found, NLEV = 0.  So when new
!   QM mnemonics are added NLEV may need to be set manually (or
!   add the new mnemonic to QMSTR if desired).

      IF (ITP .ge. 64 .and. ITP .le. 66) then   ! new ELW winds
         CALL UFBINT(LUNIT,AQCM,1,1,IRET,'MAQC')
         if (NINT(AQCM) .EQ. 3) THEN
             ARR(5,1) = 14         ! decided to purge these
         else if (NINT(AQCM) .EQ. 0) THEN
             ARR(5,1) =  2         ! these get QM = ' '
         end if
      END IF


      DO L=1,NLEV
         WQM(L) = EQMUBF(ARR(5,L))

         IF(NINT(WQM(L)).GT.11.AND.NINT(WQM(L)).LT.15)  THEN

C  A REJECT, PURGE, OR FAIL FLAG ON WIND IS TRANSFERRED TO ALL VARIABLES
C  (later move to prepdata?)
C  ---------------------------------------------------------------------

            PQM(L) = WQM(L)
            TQM(L) = WQM(L)
            QQM(L) = WQM(L)
            ZQM(L) = WQM(L)

         ELSE

            PQM(L) = EQMUBF(ARR(1,L))
            TQM(L) = EQMUBF(ARR(2,L))
            QQM(L) = EQMUBF(ARR(3,L))
            ZQM(L) = EQMUBF(ARR(4,L))

         END IF

      ENDDO
 
C  PUT THE UNPACKED REPORT INTO OBS
C  --------------------------------
 
      RSTP = 99999
      CALL S01UBF(SID,XOB,YOB,RHR,RCTIM,RSV1,RSV2,ELV,ITP,RTP,RSTP)
      CALL S02UBF(6,1,*9999)
 
C  ---------------------------------------------------------------------
C  MISC DATA GOES INTO CATEGORY 08
C  ---------------------------------------------------------------------
C  CURRENTLY NONE
C  ---------------------------------------------------------------------

      CALL S03UBF(OBS,*9999,*9998,*9997)

      IF(SUBSET(7:8).LT.'21') THEN

C  KEEP TRACK OF GOES SATELLITE WIND COUNTS BY SATELLITE ID
C  --------------------------------------------------------

         IF(HDR(6).LT.BMISS)  THEN
            IF(NINT(HDR(6)).GT.249.AND.NINT(HDR(6)).LT.260)  THEN
               KNTSAT(NINT(HDR(6))) = KNTSAT(NINT(HDR(6))) + 1
            ELSE
               KNTSAT(260) = KNTSAT(260) + 1
            END IF
         END IF
      END IF
 
      RETURN

 9999 CONTINUE
      R06UBF = 999
      RETURN

 9998 CONTINUE
      print *,'IW3UNPBF/R06UBF: RPT with ID= ',SID,' TOSSED - ZERO ',
     $ 'CAT.1-6,51,52 LVLS'
 9997 CONTINUE
      R06UBF = -9999
      KSKSAT =KSKSAT + 1
      RETURN

      END
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    I05UBF      FINDS LOCATION OF NEXT NUMERIC
C   PRGMMR: RAY CRAYTON      ORG: W/NMC41    DATE: 1989-07-07
C
C ABSTRACT: FINDS THE LOCATION OF THE NEXT NUMERIC CHARACTER
C           IN A STRING OF CHARACTERS.
C
C PROGRAM HISTORY LOG:
C 1989-07-07  RAY CRAYTON
C
C USAGE:    LOC=I05UBF(STRING,NUM,CHAR)
C   INPUT ARGUMENT LIST:
C     STRING   - CHARACTER ARRAY.
C     NUM      - NUMBER OF CHARACTERS TO SEARCH IN STRING.
C
C   OUTPUT ARGUMENT LIST:
C     I05UBF   - INTEGER*4 LOCATION OF ALPHANUMERIC CHARACTER.
C                = 0 IF NOT FOUND.
C     CHAR     - CHARACTER FOUND.
C
C REMARKS: NONE
C
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM-SP, CRAY, SGI
C
C$$$
      FUNCTION I05UBF(STRING,NUM,CHAR)
      CHARACTER*1  STRING(1),CHAR

      SAVE

      DO I = 1,NUM
         IF(STRING(I).GE.'0'.AND.STRING(I).LE.'9')  THEN
            I05UBF = I
            CHAR = STRING(I)
            GO TO 200
         END IF
      ENDDO
      I05UBF = 0
      CHAR = '?'
 200  CONTINUE
      RETURN
      END
