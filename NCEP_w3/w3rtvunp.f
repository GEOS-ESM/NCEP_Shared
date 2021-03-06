      SUBROUTINE W3RTVUNP(IUNIT,IBDATE,PP,TT,QQ,CLAL,CLAM,NLEV,
     $                    IRTCHN,RTRAD,NCHN,STNID,ISATOB,RSATOB,
     $                    DSNAME,IDSDAT,IDSDMP,IERR)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    W3RTVUNP    DECODES RTOVS RETRIEVAL DATA FROM BUFR
C   PRGMMR: KEYSER           ORG: NP22        DATE: 1999-01-20
C
C ABSTRACT: READS AND UNPACKS ONE REPORT FROM INPUT RTOVS RETRIEVAL
C   JBUFR FILE INTO SPECIFIED FORMAT. RETURNS INFORMATION IN THE
C   FORMAT DESCRIBED BELOW.  INCLUDES SOUNDING DATA ON 40 LEVELS AND
C   RTOVS RADIANCES FROM 32 CHANNELS (20 HIRS, 4 MSU, 3 SSU, AND 5
C   AVHRR CHANNELS).  ALSO, INFORMATION ABOUT THE INPUT DATA SET
C   ITSELF (NAME, CENTER DATE, DUMP TIME) IS RETURNED.
C
C PROGRAM HISTORY LOG:
C 1997-11-21  BERT B. KATZ/GSC -- ORIGINAL AUTHOR
C 1998-05-05  D.A. KEYSER/NP22 -- STREAMLINED, CORRECTED ERRORS
C 1998-06-15  D.A. KEYSER/NP22 -- MODIFIED TO READ IN LONGITUDE
C     (MNEMONIC "CLON") CORRECTLY WHETHER STORED IN DEG. E, OR IN
C     DEG. W (-) AND DEG. E (+)
C 1998-09-21  D.A. KEYSER/NP22 -- SUBROUTINE NOW Y2K AND FORTRAN 90
C      COMPLIANT
C 1999-01-20 D. A. KEYSER -- INCORPORATED BOB KISTLER'S CHANGES NEEDED
C      TO PORT THE CODE TO THE IBM SP
C
C USAGE:    CALL W3RTVUNP(IUNIT,IBDATE,PP,TT,QQ,CLAL,CLAM,NLEV,
C                         IRTCHN,RTRAD,NCHN,STNID,ISATOB,RSATOB,
C                         DSNAME,IDSDAT,IDSDMP,IERR)
C   INPUT ARGUMENT LIST:
C     IUNIT    - UNIT NUMBER OF INPUT FILE CONTAINING RETRIEVAL DATA
C              - IN BUFR FORMAT.
C
C   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
C     IBDATE   - DATE IN SECTION 1 OF BUFR MESSAGE CONTAINING THIS
C              - RETRIEVAL (IN FORM YYYYMMDDHH)
C     PP       - PRESSURE IN MB AT RETRIEVAL LEVELS (MUST BE
C              - DIMENSIONED TO AT LEAST 40 WORDS BY CALLING PROGRAM
C              - - ONLY 40 WORDS WILL BE FILLED HERE)
C     TT       - VIRTUAL TEMPERATURE IN DEG. K AT RETRIEVAL LEVELS
C              - (MUST BE DIMENSIONED TO AT LEAST 40 WORDS BY CALLING
C              - PROGRAM - ONLY 40 WORDS WILL BE FILLED HERE)
C     QQ       - SPECIFIC HUMIDITY IN G/G AT RETRIEVAL LEVELS (MUST
C              - BE DIMENSIONED TO AT LEAST 40 WORDS BY CALLING
C              - PROGRAM - ONLY 40 WORDS WILL BE FILLED HERE)
C     CLAL     - CLOUD TOP ALBEDO IN PERCENT AT RETRIEVAL LEVELS (MUST
C              - BE DIMENSIONED TO AT LEAST 40 WORDS BY CALLING
C              - PROGRAM - ONLY 40 WORDS WILL BE FILLED HERE)
C     CLAM     - CLOUD AMOUNT AT RETRIEVAL LEVELS (C.T. 0-20-011)
C              - (MUST BE DIMENSIONED TO AT LEAST 40 WORDS BY CALLING
C              - PROGRAM - ONLY 40 WORDS WILL BE FILLED HERE)
C     NLEV     - NUMBER OF RETRIEVAL LEVELS (SHOULD ALWAYS BE 40)
C     IRTCHN   - RTOVS RADIANCE CHANNELS (NUMERIC) (MUST BE
C              - DIMENSIONED TO AT LEAST 32 WORDS BY CALLING PROGRAM -
C              - ONLY 32 WORDS WILL BE FILLED HERE)
C     RTRAD    - RTOVS RADIANCES IN DEG. K (MUST BE DIMENSIONED TO AT
C              - LEAST 32 WORDS BY CALLING PROGRAM - ONLY 32 WORDS
C              - WILL BE FILLED HERE)
C     NCHN     - NUMBER OF RADIANCE CHANNELS (SHOULD ALWAYS BE 32)
C     STNID    - STATION ID (CHAR*8)
C     ISATOB   - 11-WORD INTEGER ARRAY CONTAINING RETURNED DATA
C              - (SEE REMARKS FOR CONTENTS)
C     RSATOB   - 17-WORD REAL ARRAY CONTAINING RETURNED DATA
C              - (SEE REMARKS FOR CONTENTS)
C     DSNAME   - CHARACTER*8 DATA SET NAME (SAME FOR ALL REPORTS IN
C              - A COMMON INPUT DATA SET - SEE REMARKS FOR IERR=1)
C     IDSDAT   - INTEGER DATA SET CENTER DATE IN FORM YYYYMMDDHH (SAME
C              - FOR ALL REPORTS IN A COMMON input data set - see
C              - REMARKS FOR IERR=1)
C     IDSDMP   - INTEGER DATA SET DUMP TIME IN FORM YYYYMMDDHHMM (SAME
C              - FOR ALL REPORTS IN A COMMON INPUT DATA SET - SEE
C              - REMARKS FOR IERR=1)
C     IERR     - ERROR RETURN CODE
C                 =  0 OBSERVATION READ AND UNPACKED INTO OUTPUT
C                        ARGUMENT LOCATIONS (SEE ABOVE).  SEE REMARKS
C                        FOR CONTENTS. NEXT CALL TO W3RTVUNP WILL
C                        RETURN NEXT OBSERVATION IN DATA SET.
C                 =  1 INFORMATION ABOUT THE BUFR DATASET IS RETURNED
C                        IN THE OUTPUT ARGUMENTS DSNAME, IDSDAT, IDSDMP
C                        (SEE OUTPUT ARGUMENT LIST ABOVE)
C
C                        THIS SHOULD ALWAYS OCCUR AFTER THE FIRST CALL
C                        TO THIS SUBROUTINE.  NO REPORT IS UNPACKED AT
C                        THIS POINT, AND ONLY DSNAME, IDSDAT, AND
C                        IDSDMP CONTAIN INFORMATION.  ALL SUBSEQUENT
C                        CALLS TO W3RTVUNP SHOULD RETURN THE
C                        OBSERVATIONS IN THIS DATA SET, SEQUENTIALLY,
C                        (IERR=0) UNTIL THE END OF FILE IS ENCOUNTERED
C                        (IERR=2).  THE VALUES STORED IN DSNAME,
C                        IDSDAT, AND IDSDMP WILL CONTINUE TO BE
C                        RETURNED ALONG WITH EACH REPORT WHEN IERR = 0.
C                 =  2 FOR NORMAL END-OF-FILE ENCOUNTERED.
C                 = -1 FOR END-OF-FILE ON FIRST READ -- EMPTY FILE
C                 = -2 FOR INPUT BUFR FILE NOT Y2K COMPLIANT -- NO
C                        PROCESSING DONE
C                 = -3 CENTER DATE COULD NOT BE OBTAINED FROM INPUT
C                        FILE -- NO PROCESSING DONE
C
C   INPUT FILES:
C     UNIT "IUNIT" - FILE CONTAINING BUFR RTOVS DATA
C
C   OUTPUT FILES:
C     UNIT 06 - STANDARD OUTPUT PRINT
C
C REMARKS:
C
C     CONTENTS OF OUTPUT ARGUMENT ISATOB (11 INTEGER WORDS) FOR EACH
C      RETRIEVAL:
C        WORD  1 - NESDIS SATELLITE ID NUMBER (=1 FOR NOAA-11)
C                                             (=2 FOR NOAA-12)
C                                             (=3 FOR NOAA-14)
C        WORD  2 - PREPBUFR REPORT TYPE
C                     (=161 FOR CLEAR PATH OVER LAND)
C                     (=163 FOR CLOUDY PATH OVER LAND)
C                     (=171 FOR CLEAR PATH OVER OCEAN)
C                     (=173 FOR CLOUDY PATH OVER OCEAN)
C                     {NOTE: R. TYPE 162 (NSTAR PATH OVER LAND) AND
C                            R. TYPE 172 (NSTAR PATH OVER OCEAN)
C                            ARE NOT VALID FOR RTOVS}
C        WORD  3 - SWATH LOCATION :  1 INDICATES LEFT  LIMB OF SWATH
C                                   56 INDICATES RIGHT LIMB OF SWATH
C                           USED FOR DETERMINING DISTANCE FROM NADIR
C        WORD  4 - ORBIT NUMBER
C        WORD  5 - NESDIS LAND-SEA FLAG
C        WORD  6 - SATELLITE DATA PROCESSING TECHNIQUE
C                     (=12 FOR CLOUDY PATH)
C                     (=36 FOR CLEAR PATH)
C                     (= 4 FOR UNKNOWN PATH????)
C        WORD  7 - TOVS FILTER FLAG (=0 FOR GOOD)
C                                   (=1 FOR REDUNDANT)
C                   IMPORTANT: REDUNDANT RETRIEVALS ARE AT FULL 80 KM
C                              RESOLUTION; GOOD RETRIEVALS ARE AT
C                              250 KM RESOLUTION
C        WORD  8 - SUPERADIABATIC FLAG (=0 FOR NOT SUPERADIABATIC)
C                                      (.NE. 0 FOR SUPERADIABATIC)
C        WORD  9 - DAY/NIGHT QUALIFIER (=0 FOR NIGHT)
C                                      (=1 FOR DAY)
C        WORD 10 - VERTICAL SIGNIFICANCE (SAT OBSERVATION)
C                                      (BUFR C.T. 0-08-003)
C        WORD 11 - SNOW/ICE FLAGS (TERRAIN TYPE)
C                                      (BUFR C.T. 0-13-039)
C
C     CONTENTS OF OUTPUT ARGUMENT RSATOB (17 REAL WORDS) FOR EACH
C      RETRIEVAL:
C        WORD  1 - OBSERVATION TIME IN HOURS GOOD TO 0.01 HOUR.
C        WORD  2 - LATITUDE IN DEGREES (N+,S-)
C        WORD  3 - LONGITUDE IN DEGREES EAST (0.0 to 359.99)
C        WORD  4 - SURFACE HEIGHT IN METERS
C        WORD  5 - FIRST ABOVE-GROUND TOVS PRESSURE LEVEL IN MB
C        WORD  6 - SKIN TEMPERATURE IN DEGREES K
C        WORD  7 - TOTAL PRECIPITABLE WATER IN MM
C        WORD  8 - TOTAL OZONE IN DOBSONS
C        WORD  9 - SEA SURFACE TEMPERATURE IN DEGREES K
C        WORD 10 - CLOUD-TOP PRESSURE IN MB
C        WORD 11 - SOLAR ZENITH ANGLE IN DEGREES
C        WORD 12 - SATELLITE ELEVATION IN DEGREES
C        WORD 13 - NESDIS AVERAGE NSTAR (CLOUDINESS PARAMETER)
C                   (NOTE: THIS IS NOT OF ANY USE SINCE NSTAR
C                         PATH NOT VALID FOR RTOVS
C                         = -1 FOR CLOUDY PATH
C                         = MISSING FOR CLEAR PATH)
C        WORD 14 - HIRS CHANNEL UTILIZATION FLAGS
C        WORD 15 - MSU CHANNEL UTILIZATION FLAGS
C        WORD 16 - SSU CHANNEL UTILIZATION FLAGS
C        WORD 17 - AVHRR CHANNEL UTILIZATION FLAGS
C
C FOR ALL DATA: MISSING VALUES ARE RETURNED AS 99999.0 (REAL) OR
C               99999 (INTEGER)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM-SP, CRAY, SGI
C
C$$$
 
      DATA IFIRST/0/,KOUNTR/0/

      REAL  PP(*),TT(*),QQ(*),CLAL(*),CLAM(*),RTRAD(*)

      INTEGER  IRTCHN(*),JDATE(5),JDUMP(5)

      INTEGER(8) IDSDMP

      CHARACTER*8 SUBSET,DSNAME

      CHARACTER*50 IDTSTR
      DATA IDTSTR/'SAID CLAT CLON TOFF HOUR MINU SECO FOVN ORBN SADF '/
      CHARACTER*5 IDTST2
      DATA IDTST2/'DNQL '/

      CHARACTER*50 INFSTR
      DATA INFSTR/'SOEL ELEV NSAV OZON SDPT HIRC MSUC SSUC AVHC TPWT '/
      CHARACTER*10 INFST1
      DATA INFST1/'CDTP VSAT '/

      CHARACTER*30 SRFSTR
      DATA SRFSTR/'PRES SELV LSQL TMSK SST1 TERR '/

      CHARACTER*25 RETSTR
      DATA RETSTR/'PRLC TMDB MIXR CLAL CLAM '/

      CHARACTER*10 BRTSTR
      DATA BRTSTR/'CHNM TMBR '/

      REAL(8) XIDENT(11),XINFO(12),SRFDAT(6),RETDAT(5,40),BRTDAT(2,32)

      DIMENSION  KOUNT(0:3,4),ISATOB(11),RSATOB(17)
      DATA  KOUNT/16*0/,XMSG/99999./,IMSG/99999/,BMISS/10.E10/

      CHARACTER*8   STNID

      CHARACTER*1   NMCHAR(0:35)
      DATA          NMCHAR/'0','1','2','3','4','5','6','7','8','9',
     $                     'A','B','C','D','E','F','G','H','I','J',
     $                     'K','L','M','N','O','P','Q','R','S','T',
     $                     'U','V','W','X','Y','Z'/

      IF(IFIRST.EQ.0) THEN

C  FIRST TIME IN, SET DATELEN AND GET CENTER AND DUMP TIME FOR FILE
C  ----------------------------------------------------------------

         PRINT *, ' ==> W3RTVUNP -- Y2K/F90 VERSION 01-20-1999'
         IFIRST = 1
         CALL DATELEN(10)
         CALL DUMPBF(IUNIT,JDATE,JDUMP)
cppppp
         print *, 'CENTER DATE (JDATE) = ',jdate
         print *, 'DUMP DATE (JDUMP) = ',jdump
         print *, ' '
cppppp
         IF(JDATE(1).LE.0)  then
            PRINT *, '##W3RTVUNP - CENTER DATE COULD NOT BE ',
     $       'OBTAINED FROM INPUT FILE ON UNIT ',IUNIT,' -- IERR = -3'
            IERR = -3
            RETURN
         END IF
         IF(JDATE(1).LT.100)  THEN

C IF 2-DIGIT YEAR RETURNED IN JDATE(1), MUST USE "WINDOWING" TECHNIQUE
C  TO CREATE A 4-DIGIT YEAR

C IMPORTANT: IF DATELEN(10) IS CALLED, THE DATE HERE SHOULD ALWAYS
C            CONTAIN A 4-DIGIT YEAR, EVEN IF INPUT FILE IS NOT
C            Y2K COMPLIANT (BUFRLIB DOES THE WINDOWING HERE)

            PRINT *, '##W3RTVUNP - THE FOLLOWING SHOULD NEVER ',
     $       'HAPPEN!!!!!'
            PRINT *, '##W3RTVUNP - 2-DIGIT YEAR IN JDATE(1) RETURNED ',
     $       'FROM DUMPBF (JDATE IS: ',JDATE,') - USE WINDOWING ',
     $       'TECHNIQUE TO OBTAIN 4-DIGIT YEAR'
            IF(JDATE(1).GT.20)  THEN
               JDATE(1) = 1900 + JDATE(1)
            ELSE
               JDATE(1) = 2000 + JDATE(1)
            ENDIF
            PRINT *, '##W3RTVUNP - CORRECTED JDATE(1) WITH 4-DIGIT ',
     $       'YEAR, JDATE NOW IS: ',JDATE
         ENDIF
         IDSDAT = JDATE(1)*1E6+JDATE(2)*1E4+JDATE(3)*1E2+JDATE(4)
         IF(JDUMP(1).LE.0)  THEN
            IDSDMP = 999999999999_8
         ELSE
            IF(JDUMP(1).LT.100)  THEN

C IF 2-DIGIT YEAR RETURNED IN JDUMP(1), MUST USE "WINDOWING" TECHNIQUE
C  TO CREATE A 4-DIGIT YEAR

C IMPORTANT: IF DATELEN(10) IS CALLED, THE DATE HERE SHOULD ALWAYS
C            CONTAIN A 4-DIGIT YEAR, EVEN IF INPUT FILE IS NOT
C            Y2K COMPLIANT (BUFRLIB DOES THE WINDOWING HERE)

               PRINT *, '##W3RTVUNP - THE FOLLOWING SHOULD NEVER ',
     $          'HAPPEN!!!!!'
               PRINT *, '##W3RTVUNP - 2-DIGIT YEAR IN JDUMP(1) ',
     $          'RETURNED FROM DUMPBF (JDUMP IS: ',JDUMP,') - USE ',
     $          'WINDWOING TECHNIQUE TO OBTAIN 4-DIGIT YEAR'
               IF(JDUMP(1).GT.20)  THEN
                  JDUMP(1) = 1900 + JDUMP(1)
               ELSE
                  JDUMP(1) = 2000 + JDUMP(1)
               ENDIF
               PRINT *, '##W3RTVUNP - CORRECTED JDUMP(1) WITH 4-DIGIT ',
     $          'YEAR, JDUMP NOW IS: ',JDUMP
            END IF
            IDSDMP = JDUMP(1)*1E8+JDUMP(2)*1E6+JDUMP(3)*1E4+
     $       JDUMP(4)*1E2+JDUMP(5)
         ENDIF
         DSNAME = 'RTOVS   '
         IERR = 1
         RETURN
      ELSE  IF(IFIRST.EQ.1) THEN

C  SECOND TIME IN, OPEN BUFR DATASET FOR INPUT AND DECODE FIRST MESSAGE
C  --------------------------------------------------------------------

         IFIRST = 2
         CALL OPENBF(IUNIT,'IN',IUNIT)
         CALL READMG(IUNIT,SUBSET,IBDATE,IRET)
         IF(IRET.NE.0) THEN
            WRITE(6,1009) IUNIT
 1009       FORMAT('##W3RTVUNP ERROR: EMPTY FILE IN UNIT ',I5)
            IERR = -1
            RETURN
         ENDIF
         IF(IBDATE.LT.100000000)  THEN
C IF INPUT BUFR FILE DOES NOT RETURN MESSAGES WITH A 4-DIGIT YEAR,
C  SOMETHING IS WRONG (EVEN NON-COMPLIANT BUFR MESSAGES SHOULD
C  CONSTRUCT A 4-DIGIT YEAR AS LONG AS DATELEN(10) HAS BEEN CALLED
            WRITE(6,1209) IUNIT
 1209       FORMAT('##W3RTVUNP ERROR: A 10-DIGIT SECT. 1 BUFR MESSAGE ',
     $       'DATE WAS NOT RETURNED IN UNIT',I5,' - PROBLEM WITH BUFR ',
     $       'FILE')
            IERR = -2
            RETURN
         END IF
      ENDIF

 1000 CONTINUE                                                          

C  EACH CALL TO READSB INCREASES "KOUNTR" BY 1 (REGARDLESS OF RESULT)
C  ------------------------------------------------------------------

      KOUNTR = KOUNTR + 1

      CALL READSB(IUNIT,IRET)
      IF(IRET.NE.0) THEN
        CALL READMG(IUNIT,SUBSET,IBDATE,IRET)
        IF(IRET.NE.0) THEN

C  ALL BUFR MESSAGES HAVE BEEN READ AND DECODED -- ALL DONE
C  --------------------------------------------------------

           WRITE(6,1001) IUNIT
 1001 FORMAT(//' ==> W3RTVUNP: END OF FILE ON UNIT',I3,' -- ALL DONE'/)
           IERR = 2
           RETURN
        ENDIF
         IF(IBDATE.LT.100000000)  THEN
C IF INPUT BUFR FILE DOES NOT RETURN MESSAGES WITH A 4-DIGIT YEAR,
C  SOMETHING IS WRONG (EVEN NON-COMPLIANT BUFR MESSAGES SHOULD
c  CONSTRUCT A 4-DIGIT YEAR AS LONG AS DATELEN(10) HAS BEEN CALLED
            WRITE(6,1209) IUNIT
            IERR = -2
            RETURN
         END IF
        GO TO 1000
      ENDIF

      XIDENT = BMISS
      CALL UFBINT(IUNIT,XIDENT,11, 1,IRET,IDTSTR//IDTST2)

      XINFO = BMISS
      CALL UFBINT(IUNIT,XINFO ,12, 1,IRET,INFSTR//INFST1)

      SRFDAT = BMISS
      CALL UFBINT(IUNIT,SRFDAT, 6, 1,IRET,SRFSTR)

      RETDAT = BMISS
      CALL UFBINT(IUNIT,RETDAT, 5,40,NLEV,RETSTR)

      BRTDAT = BMISS
      CALL UFBINT(IUNIT,BRTDAT, 2,32,NCHN,BRTSTR)

      ISATOB = IMSG
      RSATOB = XMSG
      PP(1:40) = XMSG
      TT(1:40) = XMSG
      QQ(1:40) = XMSG
      CLAL(1:40) = XMSG
      CLAM(1:40) = XMSG
      IRTCHN(1:32) = IMSG
      RTRAD(1:32) = XMSG

C  SATELLITE ID
C  ------------

      IF(XIDENT(1).LT.BMISS)  THEN
         ISATOB(1) = NINT(XIDENT(1))
         IF(ISATOB(1).EQ.203) THEN
            ISATOB(1)   = 1
         ELSE IF(ISATOB(1).EQ.204) THEN
            ISATOB(1)   = 2
         ELSE IF(ISATOB(1).EQ.205) THEN
            ISATOB(1)   = 3
         ELSE
            WRITE(6,2071) KOUNTR,NINT(XIDENT(1))
 2071 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') HAS ',
     $ 'AN UNRECOGNIZED SATELLITE ID (=',I5,') - SKIP IT')
            GO TO 1000
         ENDIF
      ELSE
         WRITE(6,2072) KOUNTR
 2072 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') HAS A',
     $ ' MISSING SATELLITE ID - SKIP IT')
         GO TO 1000
      ENDIF

C  LATITUDE
C  --------

      IF(XIDENT(2).LT.BMISS)  THEN
         RSATOB(2) = XIDENT(2)
      ELSE
         WRITE(6,2073) KOUNTR
 2073 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') HAS A',
     $ ' MISSING LATITUDE - SKIP IT')
         GO TO 1000
      ENDIF

C  LONGITUDE
C  ---------

      IF(XIDENT(3).LT.BMISS)  THEN

C Important: According to BUFR Manual, CLON (0-06-002) - represented
C  here by "XIDENT(3)" - should be in units of Degrees West - and East +
C  (-180.0 to +180.0); however some BUFR data sets (e.g., PREPBUFR) are
C  known to encode 0-06-002 in units of Degrees East (0.0 to 359.99) --
C  So we use the following conversion to work in either case ...
         RSATOB(3) = 360._8 - MOD(360._8-XIDENT(3),360._8)
         IF(RSATOB(3).EQ.360.0)  RSATOB(3) = 0.0
      ELSE
         WRITE(6,2074) KOUNTR
 2074 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') HAS A',
     $ ' MISSING LONGITUDE - SKIP IT')
         GO TO 1000
      ENDIF

C  FILTER FLAG
C  -----------

      IF(XIDENT(4).LT.BMISS)  THEN
         ISATOB(7) = NINT(XIDENT(4))
      ELSE
         WRITE(6,2075) KOUNTR
 2075 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') HAS A',
     $ ' MISSING FILTER FLAG - SKIP IT')
         GO TO 1000
       ENDIF

C  OBSERVATION TIME
C  ----------------

      IF(MAX(XIDENT(5),XIDENT(6),XIDENT(7)).LE.BMISS) THEN
         HRFRAC = (60.0 * XIDENT(6) + XIDENT(7)) / 3600.0
         RSATOB(1) = XIDENT(5) + HRFRAC
         RSATOB(1) = 0.01 * AINT(100.0 * RSATOB(1) + 0.5)
      ELSE
         WRITE(6,2076) KOUNTR
 2076 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') HAS ',
     $ 'ONE OR MORE MISSING TIME UNITS - SKIP IT')
         GO TO 1000
      ENDIF

C  SWATH LOCATION
C  --------------

      IF(XIDENT(8).LT.BMISS)  THEN
         ISATOB(3) = NINT(XIDENT(8))
      ELSE
         WRITE(6,2077) KOUNTR
 2077 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') HAS A',
     $ ' MISSING SWATH LOCATION - SKIP IT')
         GO TO 1000
      ENDIF

C  ORBIT NUMBER
C  ------------

      IF(XIDENT(9).LT.BMISS)  ISATOB(4) = NINT(XIDENT(9))

C  SUPERADIABATIC FLAG
C  -------------------

      IF(XIDENT(10).LT.BMISS)  ISATOB(8) = NINT(XIDENT(10))
      IF(ISATOB(8).NE.0)  THEN
         WRITE(6,2078) KOUNTR,ISATOB(8)
 2078 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') IS ',
     $ 'SUPERADIABATIC (FLAG=',I5,') - SKIP IT')
         GO TO 1000
      ENDIF

C  DAY-NIGHT INDICATOR
C  -------------------

      IF(XIDENT(11).LT.BMISS)  ISATOB(9) = NINT(XIDENT(11))

C  SOLAR ZENITH ANGLE
C  ------------------

      IF(XINFO(1).LT.BMISS)  RSATOB(11) = XINFO(1)

C  SOLAR ELEVATION
C  ---------------

      IF(XINFO(2).LT.BMISS)  RSATOB(12) = XINFO(2)

C  AVERAGE NSTAR VALUE
C  -------------------

      IF(XINFO(3).LT.BMISS)  RSATOB(13) = XINFO(3)

C  OZONE
C  -----

      IF(XINFO(4).LT.BMISS)  RSATOB(8) = XINFO(4)

C  SATELLITE DATA PROCESSING TECHNIQUE
C  -----------------------------------

      IF(XINFO(5).LT.BMISS)  THEN
         ISATOB(6) = NINT(XINFO(5))
      ELSE
         WRITE(6,2079) KOUNTR
 2079 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') HAS A',
     $ ' MISSING SATELLITE DATA PROCESSING TECHNIQUE - SKIP IT')
         GO TO 1000
      ENDIF

C  HIRS CHANNEL UTILIZATION FLAGS
C  ------------------------------

      IF(XINFO(6).LT.BMISS)  RSATOB(14) = XINFO(6)

C  MSU CHANNEL UTILIZATION FLAGS
C  ------------------------------

      IF(XINFO(7).LT.BMISS)  RSATOB(15) = XINFO(7)

C  SSU CHANNEL UTILIZATION FLAGS
C  -----------------------------

      IF(XINFO(8).LT.BMISS)  RSATOB(16) = XINFO(8)

C  AVHRR CHANNEL UTILIZATION FLAGS
C  -------------------------------

      IF(XINFO(9).LT.BMISS)  RSATOB(17) = XINFO(9)

C  TOTAL PRECIPITABLE WATER
C  ------------------------

      IF(XINFO(10).LT.BMISS)  RSATOB(7) = XINFO(10) * 10.

C  CLOUD-TOP PRESSURE
C  ------------------

      IF(XINFO(11).LT.BMISS)  RSATOB(10) = XINFO(11) * 0.01

C  VERTICAL SIGNIFICANCE INDICATOR (SATELLITE SOUNDINGS)
C  -----------------------------------------------------

      IF(XINFO(12).LT.BMISS)  ISATOB(10) = NINT(XINFO(12))

C  STATION ID
C  ----------

      STNID = '      R '
      MODSAT = MOD(ISATOB(1),4) + 1
      IF(MODSAT.EQ.4) THEN
         NF = 13
      ELSE IF(MODSAT.EQ.2) THEN
         NF = 17
      ELSE IF(MODSAT.EQ.3) THEN
         NF = 31
      ELSE IF(MODSAT.EQ.1) THEN
         NF = 35
      ENDIF
      INST = IMSG
      IF(ISATOB(6).GE.32) THEN

C  CLEAR PATH
C  ----------

         INST = 1
         STNID(8:8) = NMCHAR(NF-3)
      ELSE IF(ISATOB(6).GE.16) THEN

C  NSTAR PATH ===> NOT VALID FOR RTOVS!!!! - SHOULD NEVER SEE THIS
C  ---------------------------------------------------------------

         INST = 2
         STNID(8:8) = NMCHAR(NF-2)
      ELSE IF(ISATOB(6).GE.8) THEN

C  CLOUDY PATH
C  -----------

         INST = 3
         STNID(8:8) = NMCHAR(NF-1)
      ELSE

C  UNKNOWN PATH
C  ------------

         INST = 0
         STNID(8:8) = NMCHAR(NF)
         WRITE(6,2080) KOUNTR,STNID,ISATOB(6)
 2080 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,',ID=',A8,
     $ ') HAS AN UNKNOWN SATELLITE DATA PROCESSING TECH. (=',I5,') - ',
     $ 'SKIP IT')
         GO TO 1000
      ENDIF

      KOUNT(INST,MODSAT) = MIN(99999,KOUNT(INST,MODSAT) + 1)
      WRITE(STNID(2:6),'(I5.5)') KOUNT(INST,MODSAT)

      NADIR = INT(ABS(FLOAT(ISATOB(3)) - 28.5)) / 3
      IF(NADIR.LT.9) NADIR = NADIR + 1
      IF(NADIR.GT.0.AND.NADIR.LT.10)  THEN
         STNID(1:1) = NMCHAR(NADIR)
      ELSE
         WRITE(6,2081) KOUNTR,STNID,NADIR
 2081 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,',ID=',A8,
     $ ') HAS AN INVALID NADIR PROX. IND. (=',I5,') - SKIP IT')
         GO TO 1000
      ENDIF

      IF(NLEV.EQ.40)  THEN

C  FILL IN THE LEVEL SOUNDING DATA - MUST HAVE 40 TEMP LEVELS TO ACCEPT
C  --------------------------------------------------------------------

         DO L = 1,NLEV
            LREV = NLEV + 1 - L
            PP(LREV) = RETDAT(1,L) * 0.01
            IF(RETDAT(2,L).LT.BMISS) THEN
               IF(RETDAT(3,L).GE.BMISS) THEN
                  TT(LREV) = RETDAT(2,L)
               ELSE
                  QQ(LREV) = RETDAT(3,L) / (1.0 + RETDAT(3,L))
                  TT(LREV) = RETDAT(2,L) * (1.0 + 0.61 * QQ(LREV))
               ENDIF
            ELSE
               WRITE(6,3082) KOUNTR
 3082 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,') ',
     $ 'CONTAINS 1 OR MORE MISSING TEMERATURES IN ITS SOUNDING - SKIP ',
     $ 'IT')
               GO TO 1000
            ENDIF
            IF(RETDAT(4,L).LT.BMISS)  CLAL(LREV) = RETDAT(4,L)
            IF(RETDAT(5,L).LT.BMISS)  CLAM(LREV) = RETDAT(5,L)
         ENDDO
      ELSE
         WRITE(6,2082) STNID,KOUNTR,NLEV
 2082 FORMAT(/'##W3RTVUNP WARNING: ID ',A8,' (#',I5,') DOES NOT ',
     $ 'CONTAIN 40 SOUNDING LEVELS (HAS ',I3,' LEVELS) -- NO SOUNDING ',
     $ 'LEVEL DATA PROCESSED')
      ENDIF

      IF(NCHN.EQ.32)  THEN

C  FILL IN THE RADIANCE DATA
C  -------------------------

         DO L = 1,NCHN
            IF(BRTDAT(1,L).LT.BMISS)  IRTCHN(L) = NINT(BRTDAT(1,L))
            IF(BRTDAT(2,L).LT.BMISS)  RTRAD(L)  = BRTDAT(2,L)
         ENDDO
      ELSE
         WRITE(6,2083) STNID,KOUNTR,NCHN
 2083 FORMAT(/'##W3RTVUNP WARNING: ID ',A8,' (#',I5,') DOES NOT ',
     $ 'CONTAIN 32 RADIANCE CHANNELS (HAS ',I3,' CHANNELS) -- NO ',
     $ 'RADIANCE DATA PROCESSED')
      ENDIF

C  BOTTOM LEVEL PRESSURE
C  ---------------------

      IF(SRFDAT(1).LT.BMISS)  THEN
         RSATOB(5) = 0.01 * SRFDAT(1)
      ELSE
         WRITE(6,2084) KOUNTR,STNID
 2084 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,',ID=',A8,
     $ ') HAS A MISSING BOTTOM PRESSURE LEVEL - SKIP IT')
         GO TO 1000
      ENDIF

C  SURFACE HEIGHT
C  --------------

      IF(SRFDAT(2).LT.BMISS)  THEN
         RSATOB(4) = SRFDAT(2)
      ELSE
         WRITE(6,2085) KOUNTR,STNID
 2085 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,',ID=',A8,
     $ ') HAS A MISSING SURFACE HEIGHT - SKIP IT')
         GO TO 1000
      ENDIF

C  REPORT TYPE AND LAND/SEA TAG
C  ----------------------------

      IF(SRFDAT(3).LT.BMISS)  THEN
         ISATOB(5) = NINT(SRFDAT(3))
         IF(MOD(ISATOB(5),2).EQ.0)  THEN
            IF(INST.LT.IMSG)  ISATOB(2) = 170 + INST
            ISATOB(5) = 0
         ELSE 
            IF(INST.LT.IMSG)  ISATOB(2) = 160 + INST
            ISATOB(5) = 1
         ENDIF
      ELSE
         WRITE(6,2086) KOUNTR,STNID
 2086 FORMAT(/'##W3RTVUNP WARNING: A DECODED RETRIEVAL (#',I5,',ID=',A8,
     $ ') HAS A MISSING LAND/SEA INDICATOR - SKIP IT')
            GO TO 1000
      ENDIF

C  SKIN TEMPERATURE
C  ----------------

      IF(SRFDAT(4).LT.BMISS)  RSATOB(6) = SRFDAT(4)

C  SEA-SURFACE  TEMPERATURE
C  ------------------------

      IF(SRFDAT(5).LT.BMISS)  RSATOB(9) = SRFDAT(5)

C  TERRAIN INDICATOR
C  -----------------

      IF(SRFDAT(6).LT.BMISS)  ISATOB(11) = NINT(SRFDAT(6))

C  RETURN WITH DECODED RTOVS RETRIEVAL
C  -----------------------------------

      IERR = 0

      RETURN
      END
