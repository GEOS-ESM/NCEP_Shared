!------------------------------------------------------------------------------
!M+
! NAME:
!       TauCoeff_SA_TraceGas_Binary_IO
!
! PURPOSE:
!       Module containing routines to read and write Binary format
!       TauCoeff_SARTA_TraceGas files.
!       
! CATEGORY:
!       Optical Depth : Coefficients
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       USE TauCoeff_SA_TraceGas_Binary_IO
!
! MODULES:
!       Type_Kinds:            Module containing definitions for kinds
!                              of variable types.
!
!       File_Utility:          Module containing generic file utility routines
!
!       Message_Handler:         Module to define simple error codes and
!                              handle error conditions
!                              USEs: FILE_UTILITY module
!
!       Binary_File_Utility:   Module containing utility routines for "Binary" 
!                              format datafiles.
!                              USEs: TYPE_KINDS module
!                                    FILE_UTILITY module
!                                    Message_Handler module
!
!       TauCoeff_SARTA_TraceGas_Define:Module defining the TauCoeff_SARTA_TraceGas data structure and
!                              containing routines to manipulate it.
!                              USEs: TYPE_KINDS module
!                                    FILE_UTILITY module
!                                    Message_Handler module
!
! CONTAINS:
!       Inquire_TauCoeff_TraceGas_Binary: Function to inquire a Binary format
!                                       TauCoeff_SARTA_TraceGas file.
!
!       Read_TauCoeff_TraceGas_Binary:    Function to read a Binary format
!                                       TauCoeff_SARTA_TraceGas file.
!
!       Write_TauCoeff_TraceGas_Binary:   Function to write a Binary format
!                                       TauCoeff_SARTA_TraceGas file.
!
! INCLUDE FILES:
!       None.
!
! EXTERNALS:
!       None.
!
! COMMON BLOCKS:
!       None.
!
! FILES ACCESSED:
!       User specified Binary format TauCoeff_SARTA_TraceGas data files for both
!       input and output.
!
! CREATION HISTORY:
!       Written by:     Yong Chen, CSU/CIRA 03-May-2006
!                       Yong.Chen@noaa.gov
!
!  Copyright (C) 2006 Yong Chen
!
!M-
!------------------------------------------------------------------------------

MODULE TauCoeff_SA_TraceGas_Binary_IO


  ! ----------
  ! Module use
  ! ----------

  USE Type_Kinds
  USE File_Utility
  USE Message_Handler
  USE Binary_File_Utility

  USE TauCoeff_SARTA_TraceGas_Define



  ! -----------------------
  ! Disable implicit typing
  ! -----------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE
  
  PUBLIC :: Inquire_TauCoeff_TraceGas_Bin
  PUBLIC :: Read_TauCoeff_TraceGas_Binary
  PUBLIC :: Write_TauCoeff_TraceGas_Binary

  ! ---------------------
  ! Procedure overloading
  ! ---------------------

  INTERFACE Read_TauCoeff_TraceGas_Binary 
    MODULE PROCEDURE Read_TauCoeff_TraceGas_Scalar
    MODULE PROCEDURE Read_TauCoeff_TraceGas_Rank1
  END INTERFACE Read_TauCoeff_TraceGas_Binary 

  INTERFACE Write_TauCoeff_TraceGas_Binary 
    MODULE PROCEDURE Write_TauCoeff_TraceGas_Scalar
    MODULE PROCEDURE Write_TauCoeff_TraceGas_Rank1
  END INTERFACE Write_TauCoeff_TraceGas_Binary
  
  ! -----------------
  ! Module parameters
  ! -----------------

  ! -- Module RCS Id string
  CHARACTER( * ), PRIVATE, PARAMETER :: MODULE_RCS_ID = &

  ! -- Keyword set value
  INTEGER, PRIVATE, PARAMETER :: UNSET = 0
  INTEGER, PRIVATE, PARAMETER :: SET   = 1

  INTEGER, PRIVATE, PARAMETER :: NO    = UNSET
  INTEGER, PRIVATE, PARAMETER :: YES   = SET


CONTAINS



  FUNCTION Read_TauCoeff_TraceGas_Record( FileID,           &  ! Input
                                     	TauCoeff_TraceGas,  &  ! Output
                                     	No_Allocate,      &  ! Optional input
                                     	Message_Log )     &  ! Error messaging
                                      RESULT ( Error_Status )



    !#--------------------------------------------------------------------------#
    !#                          -- TYPE DECLARATIONS --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    INTEGER,                       INTENT( IN )     :: FileID

    ! -- Outut
    TYPE( TauCoeff_SARTA_TraceGas_type ),  INTENT( IN OUT ) :: TauCoeff_TraceGas

    ! -- Optional input
    INTEGER,        OPTIONAL, INTENT( IN )     :: No_Allocate

    ! -- Error handler Message log
    CHARACTER( * ), OPTIONAL, INTENT( IN )     :: Message_Log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: Error_Status


    ! -------------------
    ! Function parameters
    ! -------------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'Read_TauCoeff_SARTA_TraceGas_Binary(Record)'


    ! ------------------
    ! Function variables
    ! ------------------

    CHARACTER( 256 ) :: Message
    LOGICAL :: Yes_Allocate
    INTEGER :: IO_Status
    INTEGER :: Destroy_Status
    INTEGER( Long ) :: n_Layers
    INTEGER( Long ) :: n_Predictors  
    INTEGER( Long ) :: n_Channels         
  

    !#--------------------------------------------------------------------------#
    !#                     -- INITIALISE THE ERROR STATUS --                    #
    !#--------------------------------------------------------------------------#

    Error_Status = SUCCESS



    !#--------------------------------------------------------------------------#
    !#                            -- CHECK INPUT --                             #
    !#--------------------------------------------------------------------------#

    ! -------------------------------------------
    ! Process the optional No_Allocate argument
    ! -------------------------------------------

    ! -- Default action is to allocate the structure....
    Yes_Allocate = .TRUE.

    ! -- ...unless the No_Allocate optional argument is set.
    IF ( PRESENT( No_Allocate ) ) THEN
      IF ( No_Allocate == SET ) Yes_Allocate = .FALSE.
    END IF


    !#--------------------------------------------------------------------------#
    !#                     -- READ THE TauCoeff_SARTA_TraceGas DIMENSIONS --                      #
    !#--------------------------------------------------------------------------#

    READ( FileID, IOSTAT = IO_Status ) n_Layers, &                    
                                       n_Predictors, &
                                       n_Channels

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error reading TauCoeff_SARTA_TraceGas data dimensions. IOSTAT = ", i5 )' ) &
                      IO_Status
      GOTO 1000  ! Clean up
    END IF



    !#--------------------------------------------------------------------------#
    !#              -- ALLOCATE THE TauCoeff_SARTA_TraceGas STRUCTURE IF REQUIRED --              #
    !#--------------------------------------------------------------------------#

    IF ( Yes_Allocate ) THEN


      ! ----------------------
      ! Perform the allocation
      ! ----------------------

      Error_Status = Allocate_TauCoeff_TraceGas( n_Layers, &
                                               n_Predictors, &
                                               n_Channels, &
                                               TauCoeff_TraceGas, &
                                          Message_Log = Message_Log )

      IF ( Error_Status /= SUCCESS ) THEN
        Message = 'Error allocating TauCoeff_SARTA_TraceGas structure.'
        GOTO 1000
      END IF


    ELSE


      ! ---------------------------------------------------------
      ! Structure already allocated. Check the association status
      ! ---------------------------------------------------------

      IF ( .NOT. Associated_TauCoeff_TraceGas( TauCoeff_TraceGas ) ) THEN
        Message = 'TauCoeff_SARTA_TraceGas structure components are not associated.'
        GOTO 1000  ! Clean up
      END IF


      ! --------------------------
      ! Check the dimension values
      ! --------------------------

      IF ( n_Layers /= TauCoeff_TraceGas%n_Layers ) THEN
        WRITE( Message, '( "TauCoeff_SARTA_TraceGas data dimensions, ", i5, &
                          &" are inconsistent with structure definition, ", i5, "." )' ) &
                        n_Layers, TauCoeff_TraceGas%n_Layers
        GOTO 1000  ! Clean up
      END IF

    END IF



    !#--------------------------------------------------------------------------#
    !#                          -- READ THE TauCoeff_SARTA_TraceGas DATA --             #
    !#--------------------------------------------------------------------------#

    READ( FileID, IOSTAT = IO_Status ) TauCoeff_TraceGas%Absorber_ID
 
    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error reading TauCoeff_SARTA_TraceGas Absorber_ID data. IOSTAT = ", i5 )' ) &
                      IO_Status
      GOTO 1000  ! Clean up
    END IF


    !#--------------------------------------------------------------------------#
    !#                       -- READ THE CHANNEL INDEX DATA --                  #
    !#--------------------------------------------------------------------------#

    READ( FileID, IOSTAT = IO_Status ) TauCoeff_TraceGas%Channel_Index

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error reading TauCoeff_SARTA_TraceGas Channel_Index data. IOSTAT = ", i5 )' ) &
                      IO_Status
      GOTO 1000  ! Clean up
    END IF

    !#--------------------------------------------------------------------------#
    !#                -- READ THE GAS ABSORPTION COEFFICIENTS --                #
    !#--------------------------------------------------------------------------#
 
    READ( FileID, IOSTAT = IO_Status ) TauCoeff_TraceGas%Trace_Coeff

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error reading TauCoeff_SARTA_TraceGas absorption coefficients data. IOSTAT = ", i5 )' ) &
                      IO_Status
      GOTO 1000  ! Clean up
    END IF


    RETURN
        
    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    !#                      -= CLEAN UP AFTER AN ERROR -=                       #
    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    1000 CONTINUE
    Error_Status = FAILURE
    CALL Display_Message( ROUTINE_NAME, &
                          TRIM( Message ), &
                          Error_Status, &
                          Message_Log = Message_Log )

    IF ( Associated_TauCoeff_TraceGas( TauCoeff_TraceGas ) ) THEN
      Destroy_Status = Destroy_TauCoeff_TraceGas( TauCoeff_TraceGas, &
                                           Message_Log = Message_Log )
    END IF

    CLOSE( FileID, IOSTAT = IO_Status )

  END FUNCTION Read_TauCoeff_TraceGas_Record


  FUNCTION Write_TauCoeff_TraceGas_Record( FileID,           &  ! Input
                                     	 TauCoeff_TraceGas,  &  ! Input
                                     	 Message_Log )     &  ! Error messaging
                                       RESULT ( Error_Status )



    !#--------------------------------------------------------------------------#
    !#                          -- TYPE DECLARATIONS --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    INTEGER,                       INTENT( IN )  :: FileID
    TYPE( TauCoeff_SARTA_TraceGas_type ),  INTENT( IN )  :: TauCoeff_TraceGas

    ! -- Error handler Message log
    CHARACTER( * ), OPTIONAL, INTENT( IN )  :: Message_Log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: Error_Status


    ! -------------------
    ! Function parameters
    ! -------------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'Write_TauCoeff_SARTA_TraceGas_Binary(Record)'
    CHARACTER( * ), PARAMETER :: FILE_STATUS_ON_ERROR = 'DELETE'


    ! ------------------
    ! Function variables
    ! ------------------

    CHARACTER( 256 ) :: Message

    INTEGER :: IO_Status
 


    !#--------------------------------------------------------------------------#
    !#                     -- INITIALISE THE ERROR STATUS --                    #
    !#--------------------------------------------------------------------------#

    Error_Status = SUCCESS



    !#--------------------------------------------------------------------------#
    !#             -- CHECK STRUCTURE POINTER ASSOCIATION STATUS --             #
    !#--------------------------------------------------------------------------#

    IF ( .NOT. Associated_TauCoeff_TraceGas( TauCoeff_TraceGas ) ) THEN
      Message = 'Some or all INPUT TauCoeff_SARTA_TraceGas pointer members are NOT associated.'
      GOTO 1000
    END IF


    !#--------------------------------------------------------------------------#
    !#                    -- WRITE THE TauCoeff_SARTA_TraceGas DIMENSIONS --            #
    !#--------------------------------------------------------------------------#

    WRITE( FileID, IOSTAT = IO_Status ) TauCoeff_TraceGas%n_Layers, &
                                        TauCoeff_TraceGas%n_Predictors, &
                                        TauCoeff_TraceGas%n_Channels

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error writing TauCoeff_SARTA_TraceGas data dimensions. IOSTAT = ", i5 )' ) &
                      IO_Status
      GOTO 1000  ! Clean up
    END IF



    !#--------------------------------------------------------------------------#
    !#                         -- WRITE THE TauCoeff_SARTA_TraceGas DATA --             #
    !#--------------------------------------------------------------------------#
 
    WRITE( FileID, IOSTAT = IO_Status ) TauCoeff_TraceGas%Absorber_ID
 
    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error writing TauCoeff_SARTA_TraceGas Absorber_ID data. IOSTAT = ", i5 )' ) &
                      IO_Status
      GOTO 1000  ! Clean up
    END IF


    !#--------------------------------------------------------------------------#
    !#                       -- WRITE THE CHANNEL INDEX DATA --                 #
    !#--------------------------------------------------------------------------#


    WRITE( FileID, IOSTAT = IO_Status )  TauCoeff_TraceGas%Channel_Index 
    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error writing TauCoeff_SARTA_TraceGas Channel_Index data. IOSTAT = ", i5 )' ) &
                      IO_Status
      GOTO 1000  ! Clean up
    END IF
 
    !#--------------------------------------------------------------------------#
    !#                -- WRITE THE GAS ABSORPTION COEFFICIENTS --               #
    !#--------------------------------------------------------------------------#
 
    WRITE( FileID, IOSTAT = IO_Status ) TauCoeff_TraceGas%Trace_Coeff

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error writing TauCoeff_SARTA_TraceGas absorption coefficients data. IOSTAT = ", i5 )' ) &
                      IO_Status
      GOTO 1000  ! Clean up
    END IF


    RETURN


    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    !#                      -= CLEAN UP AFTER AN ERROR -=                       #
    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    1000 CONTINUE
    Error_Status = FAILURE
    CALL Display_Message( ROUTINE_NAME, &
                          TRIM( Message ), &
                          Error_Status, &
                          Message_Log = Message_Log )
    CLOSE( FileID, STATUS = FILE_STATUS_ON_ERROR, IOSTAT = IO_Status )

  END FUNCTION Write_TauCoeff_TraceGas_Record





  FUNCTION Inquire_TauCoeff_TraceGas_Bin( Filename,     &  ! Input
                                           n_Tracegases,   &  ! Optional output
                                           RCS_Id,	 &  ! Revision control
                                           Message_Log ) &  ! Error messaging
                                         RESULT ( Error_Status )



    !#--------------------------------------------------------------------------#
    !#                          -- TYPE DECLARATIONS --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    CHARACTER( * ),           INTENT( IN )  :: Filename

    ! -- Optional output
    INTEGER,        OPTIONAL, INTENT( OUT ) :: n_Tracegases
 
    ! -- Revision control
    CHARACTER( * ), OPTIONAL, INTENT( OUT ) :: RCS_Id

    ! -- Error Message log file
    CHARACTER( * ), OPTIONAL, INTENT( IN )  :: Message_Log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: Error_Status


    ! -------------------
    ! Function parameters
    ! -------------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'Inquire_TauCoeff_SARTA_TraceGas_Binary'


    ! ------------------
    ! Function variables
    ! ------------------

    CHARACTER( 256 ) :: Message
    INTEGER :: IO_Status
    INTEGER :: FileID
    INTEGER( Long ) :: File_n_Tracegases
 



    !#--------------------------------------------------------------------------#
    !#                  -- DEFINE A SUCCESSFUL EXIT STATUS --                   #
    !#--------------------------------------------------------------------------#

    Error_Status = SUCCESS


    ! -----------------------------------
    ! Set the RCS Id argument if supplied
    ! -----------------------------------

    IF ( PRESENT( RCS_Id ) ) THEN
      RCS_Id = ' '
      RCS_Id = MODULE_RCS_ID
    END IF

    !#--------------------------------------------------------------------------#
    !#                            -- CHECK INPUT --                             #
    !#--------------------------------------------------------------------------#

    ! --------------------------
    ! Check that the file exists
    ! --------------------------

    IF ( .NOT. File_Exists( TRIM( Filename ) ) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'File '//TRIM( Filename )//' not found.', &
                            Error_Status, &
                            Message_Log = Message_Log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#             -- OPEN THE BINARY FORMAT TauCoeff_SARTA_TraceGas DATA FILE --              #
    !#--------------------------------------------------------------------------#

    Error_Status = Open_Binary_File( TRIM( Filename ), &
                                     FileID,   &
                                     Message_Log = Message_Log )

    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error opening TauCoeff_SARTA_TraceGas file '//TRIM( Filename ), &
                            Error_Status, &
                            Message_Log = Message_Log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#                      -- READ THE n_Tracegases of TauCoeff --                #
    !#--------------------------------------------------------------------------#

    ! ------------------------------------
    ! Read the Release/Version information
    ! ------------------------------------

    READ( FileID, IOSTAT = IO_Status ) File_n_Tracegases 

    IF ( IO_Status /= 0 ) THEN
      Error_Status = FAILURE
      WRITE( Message, '( "Error reading File_n_Tracegases data dimension from ", a, &
                        &". IOSTAT = ", i5 )' ) &
                      TRIM( Filename ), IO_Status
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM( Message ), &
                            Error_Status, &
                            Message_Log = Message_Log )
      CLOSE( FileID, IOSTAT = IO_Status )
      RETURN
    END IF

    !#--------------------------------------------------------------------------#
    !#                    -- SAVE THE NUMBER OF TRACE GASES  --                 #
    !#--------------------------------------------------------------------------#

    IF ( PRESENT( n_Tracegases ) ) n_Tracegases = File_n_Tracegases

    !#--------------------------------------------------------------------------#
    !#                          -- CLOSE THE FILE --                            #
    !#--------------------------------------------------------------------------#

    CLOSE( FileID, IOSTAT = IO_Status )

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error closing ", a, ". IOSTAT = ", i5 )' ) &
                      TRIM( Filename ), IO_Status
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM( Message ), &
                            WARNING, &
                            Message_Log = Message_Log )
    END IF


 
  END FUNCTION Inquire_TauCoeff_TraceGas_Bin





  FUNCTION Read_TauCoeff_TraceGas_Scalar( Filename,   &  ! Input
                                 TauCoeff_TraceGas,   &  ! Output
                                 No_File_Close,     &  ! Optional input
                                 No_Allocate,       &  ! Optional input
                                 n_Tracegases,         &  ! Optional output
                                 RCS_Id,            &  ! Revision control
                                 Message_Log )      &  ! Error messaging
                               RESULT ( Error_Status )



    !#--------------------------------------------------------------------------#
    !#                          -- TYPE DECLARATIONS --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    CHARACTER( * ),           INTENT( IN )     :: Filename

    ! -- Output
    TYPE( TauCoeff_SARTA_TraceGas_type ),    INTENT( IN OUT ) :: TauCoeff_TraceGas

    ! -- Optional input
    INTEGER,        OPTIONAL, INTENT( IN )     :: No_File_Close 
    INTEGER,        OPTIONAL, INTENT( IN )     :: No_Allocate

    ! -- Optional output
    INTEGER,        OPTIONAL, INTENT( OUT )    :: n_Tracegases

    ! -- Revision control
    CHARACTER( * ), OPTIONAL, INTENT( OUT )    :: RCS_Id

    ! -- Error message log file
    CHARACTER( * ), OPTIONAL, INTENT( IN )     :: Message_Log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: Error_Status


    ! -------------------
    ! Function parameters
    ! -------------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'Read_TauCoeff_SARTA_TraceGas_Binary(Scalar)'


    ! ------------------
    ! Function variables
    ! ------------------

    CHARACTER( 256 ) :: Message
    LOGICAL :: Yes_File_Close
    INTEGER :: IO_Status
    INTEGER :: Destroy_Status
    INTEGER :: FileID
    INTEGER( Long ):: n_Input_Traces
 


    !#--------------------------------------------------------------------------#
    !#                     -- INITIALISE THE ERROR STATUS --                    #
    !#--------------------------------------------------------------------------#

    Error_Status = SUCCESS


    !#--------------------------------------------------------------------------#
    !#                 -- SET THE RCS ID ARGUMENT IF SUPPLIED --                #
    !#--------------------------------------------------------------------------#
 
    IF ( PRESENT( RCS_Id ) ) THEN
      RCS_Id = ' '
      RCS_Id = MODULE_RCS_ID
    END IF


    !#--------------------------------------------------------------------------#
    !#                            -- CHECK INPUT --                             #
    !#--------------------------------------------------------------------------#
    ! ---------------------------------------------
    ! Check that the file is open. If not, open it.
    ! Otherwise get its file id
    ! ---------------------------------------------

    IF ( .NOT. File_Open( FileName ) ) THEN

      ! -- Check that the file exists
      IF ( .NOT. File_Exists( TRIM( Filename ) ) ) THEN
        Message = 'File '//TRIM( Filename )//' not found.'
        GOTO 2000  ! Clean up
      END IF 

      ! -- Open the file
      Error_Status = Open_Binary_File( TRIM( Filename ), &
                                     FileID,   &
                                     Message_Log = Message_Log )
 
      IF ( Error_Status /= SUCCESS ) THEN
       Message ='Error opening '//TRIM( Filename )  
        GOTO 2000  ! Clean up
      END IF

    ELSE

      ! -- Inquire for the logical unit number
      INQUIRE( FILE = Filename, NUMBER = FileID )

      ! -- Ensure it's valid
      IF ( FileID == -1 ) THEN
        Message = 'Error inquiring '//TRIM( Filename )//' for its FileID'
        GOTO 2000  ! Clean up
      END IF
 
    END IF


    ! -------------------------------------------
    ! Process the optional No_File_Close argument
    ! -------------------------------------------

    ! -- Default action is to close the file on exit....
    Yes_File_Close = .TRUE.

    ! -- ...unless the No_File_Close optional argument is set.
    IF ( PRESENT( No_File_Close ) ) THEN
      IF ( No_File_Close == SET ) Yes_File_Close = .FALSE.
    END IF


    !#--------------------------------------------------------------------------#
    !#                       -- READ THE NUMBER OF TRACE GASES --                   #
    !#--------------------------------------------------------------------------#

    READ( FileID, IOSTAT = IO_Status ) n_Input_Traces

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error reading n_Input_Traces data dimension from ", a, &
                        &". IOSTAT = ", i5 )' ) &
                      TRIM( Filename ), IO_Status
      GOTO 1000  ! Clean up
    END IF

    ! ---------------------
    ! Check if n_Tracegases > 1
    ! ---------------------

    IF ( n_Input_Traces > 1 ) THEN
      WRITE( Message, '( "Number of trace gases > 1 (",i5,") and output TauCoeff_SARTA_TraceGas structure ", &
                        &"is scalar." )' ) n_Input_Traces
      GOTO 1000  ! Clean up
    END IF

 
    !#--------------------------------------------------------------------------#
    !#                        -- READ THE STRUCTURE DATA --                     #
    !#--------------------------------------------------------------------------#
    Error_Status = Read_TauCoeff_TraceGas_Record( FileID, &
                                      TauCoeff_TraceGas, &
                                      No_Allocate = No_Allocate, &
                                      Message_Log = Message_Log )

    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error reading Cloud record from '//TRIM( Filename )
      GOTO 1000  ! Clean up  
    END IF

    !#--------------------------------------------------------------------------#
    !#                   -- SAVE THE NUMBER OF TRACE GASES READ --                  #
    !#--------------------------------------------------------------------------#

    IF ( PRESENT( n_Tracegases ) ) n_Tracegases = 1

    !#--------------------------------------------------------------------------#
    !#                          -- CLOSE THE FILE --                            #
    !#--------------------------------------------------------------------------#

    IF ( Yes_File_Close ) THEN

      CLOSE( FileID, IOSTAT = IO_Status )

      IF ( IO_Status /= 0 ) THEN
        WRITE( Message, '( "Error closing ", a, ". IOSTAT = ", i5 )' ) &
                        TRIM( Filename ), IO_Status
        CALL Display_Message( ROUTINE_NAME, &
                              TRIM( Message ), &
                              WARNING, &
                              Message_Log = Message_Log )
      END IF

    END IF

    RETURN
 

    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    !#                      -= CLEAN UP AFTER AN ERROR -=                       #
    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    1000 CONTINUE
    CLOSE( FileID, IOSTAT = IO_Status )

    2000 CONTINUE
    Error_Status = FAILURE
    CALL Display_Message( ROUTINE_NAME, &
                          TRIM( Message ), &
                          Error_Status, &
                          Message_Log = Message_Log )
    IF ( Associated_TauCoeff_TraceGas( TauCoeff_TraceGas ) ) THEN

        Destroy_Status = Destroy_TauCoeff_TraceGas( TauCoeff_TraceGas, &
                                           Message_Log = Message_Log )
     
    END IF
    
  END FUNCTION Read_TauCoeff_TraceGas_Scalar


  FUNCTION Read_TauCoeff_TraceGas_Rank1( Filename,    &  ! Input
                                 TauCoeff_TraceGas,   &  ! Output
                                 No_File_Close,     &  ! Optional input
                                 No_Allocate,       &  ! Optional input
                                 n_Tracegases,         &  ! Optional output
                                 RCS_Id,            &  ! Revision control
                                 Message_Log )      &  ! Error messaging
                               RESULT ( Error_Status )



    !#--------------------------------------------------------------------------#
    !#                          -- TYPE DECLARATIONS --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    CHARACTER( * ),           INTENT( IN )     :: Filename

    ! -- Output
    TYPE( TauCoeff_SARTA_TraceGas_type ), DIMENSION( : ), INTENT( IN OUT ) :: TauCoeff_TraceGas

    ! -- Optional input
    INTEGER,        OPTIONAL, INTENT( IN )     :: No_File_Close
    INTEGER,        OPTIONAL, INTENT( IN )     :: No_Allocate

    ! -- Optional output
    INTEGER,        OPTIONAL, INTENT( OUT )    :: n_Tracegases

    ! -- Revision control
    CHARACTER( * ), OPTIONAL, INTENT( OUT )    :: RCS_Id

    ! -- Error message log file
    CHARACTER( * ), OPTIONAL, INTENT( IN )     :: Message_Log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: Error_Status


    ! -------------------
    ! Function parameters
    ! -------------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'Read_TauCoeff_SARTA_TraceGas_Binary(Rank-1)'


    ! ------------------
    ! Function variables
    ! ------------------

    CHARACTER( 256 ) :: Message
    LOGICAL :: Yes_File_Close 
    INTEGER :: IO_Status
    INTEGER :: Destroy_Status
    INTEGER :: FileID
    INTEGER( Long )::m, n_Input_Traces
 


    !#--------------------------------------------------------------------------#
    !#                     -- INITIALISE THE ERROR STATUS --                    #
    !#--------------------------------------------------------------------------#

    Error_Status = SUCCESS


    !#--------------------------------------------------------------------------#
    !#                 -- SET THE RCS ID ARGUMENT IF SUPPLIED --                #
    !#--------------------------------------------------------------------------#

    IF ( PRESENT( RCS_Id ) ) THEN
      RCS_Id = ' '
      RCS_Id = MODULE_RCS_ID
    END IF


 
    !#--------------------------------------------------------------------------#
    !#                            -- CHECK INPUT --                             #
    !#--------------------------------------------------------------------------#
    ! ---------------------------------------------
    ! Check that the file is open. If not, open it.
    ! Otherwise get its file id
    ! ---------------------------------------------

    IF ( .NOT. File_Open( FileName ) ) THEN

      ! -- Check that the file exists
      IF ( .NOT. File_Exists( TRIM( Filename ) ) ) THEN
        Message = 'File '//TRIM( Filename )//' not found.'
        GOTO 2000  ! Clean up
      END IF 
  
      ! -- Open the file
      Error_Status = Open_Binary_File( TRIM( Filename ), &
                                     FileID,   &
                                     Message_Log = Message_Log )

      IF ( Error_Status /= SUCCESS ) THEN
       Message ='Error opening '//TRIM( Filename )  
        GOTO 2000  ! Clean up
      END IF

    ELSE

      ! -- Inquire for the logical unit number
      INQUIRE( FILE = Filename, NUMBER = FileID )

      ! -- Ensure it's valid
      IF ( FileID == -1 ) THEN
        Message = 'Error inquiring '//TRIM( Filename )//' for its FileID'
        GOTO 2000  ! Clean up
      END IF
 
    END IF

    ! -------------------------------------------
    ! Process the optional No_File_Close argument
    ! -------------------------------------------

    ! -- Default action is to close the file on exit....
    Yes_File_Close = .TRUE.

    ! -- ...unless the No_File_Close optional argument is set.
    IF ( PRESENT( No_File_Close ) ) THEN
      IF ( No_File_Close == SET ) Yes_File_Close = .FALSE.
    END IF


    !#--------------------------------------------------------------------------#
    !#                       -- READ THE NUMBER OF TRACE GASES --                   #
    !#--------------------------------------------------------------------------#

    READ( FileID, IOSTAT = IO_Status ) n_Input_Traces

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error reading n_Input_Traces data dimension from ", a, &
                        &". IOSTAT = ", i5 )' ) &
                      TRIM( Filename ), IO_Status
      GOTO 1000  ! Clean up
    END IF

    ! -----------------------------------------
    ! Check if n_Tracegases > size of output array
    ! -----------------------------------------

    IF ( n_Input_Traces > size( TauCoeff_TraceGas ) ) THEN
      WRITE( Message, '( "Number of trace gases ",i5," > size of the output ", &
                         &"structure array, ", i5, "." )' ) &
                      n_Input_Traces, SIZE( TauCoeff_TraceGas ) 
      GOTO 1000  ! Clean up
    END IF

    !#--------------------------------------------------------------------------#
    !#                           -- LOOP OVER SUMSETS --                        #
    !#--------------------------------------------------------------------------#

    Trace_Loop: DO m = 1, n_Input_Traces

      Error_Status = Read_TauCoeff_TraceGas_Record( FileID, &
                                        TauCoeff_TraceGas(m), &
                                        No_Allocate = No_Allocate, &
                                        Message_Log = Message_Log )!

      IF ( Error_Status /= SUCCESS ) THEN
        WRITE( Message, '( "Error reading Cloud element #", i5, " from ", a )' ) &
                        m, TRIM( Filename )
        GOTO 1000  ! Clean up
      END IF
 
    END DO Trace_Loop

 
    !#--------------------------------------------------------------------------#
    !#                   -- SAVE THE NUMBER OF TRACE GASES READ --                  #
    !#--------------------------------------------------------------------------#

    IF ( PRESENT( n_Tracegases ) ) n_Tracegases = n_Input_Traces

 
    !#--------------------------------------------------------------------------#
    !#                          -- CLOSE THE FILE --                            #
    !#--------------------------------------------------------------------------#

    IF ( Yes_File_Close ) THEN

      CLOSE( FileID, IOSTAT = IO_Status )

      IF ( IO_Status /= 0 ) THEN
        WRITE( Message, '( "Error closing ", a, ". IOSTAT = ", i5 )' ) &
                        TRIM( Filename ), IO_Status
        CALL Display_Message( ROUTINE_NAME, &
                              TRIM( Message ), &
                              WARNING, &
                              Message_Log = Message_Log )
      END IF

    END IF

    RETURN


    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    !#                      -= CLEAN UP AFTER AN ERROR -=                       #
    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    1000 CONTINUE
    CLOSE( FileID, IOSTAT = IO_Status )

    2000 CONTINUE
    Error_Status = FAILURE
    CALL Display_Message( ROUTINE_NAME, &
                          TRIM( Message ), &
                          Error_Status, &
                          Message_Log = Message_Log )
    Destroy_Status = Destroy_TauCoeff_TraceGas( TauCoeff_TraceGas, &
                                           Message_Log = Message_Log )
    
  END FUNCTION Read_TauCoeff_TraceGas_Rank1





  FUNCTION Write_TauCoeff_TraceGas_Scalar( Filename,     &  ! Input
                                  TauCoeff_TraceGas,     &  ! Input
                                  No_File_Close,       &  ! Optional input
                                  RCS_Id,              &  ! Revision control
                                  Message_Log )        &  ! Error messaging
                                RESULT ( Error_Status )



    !#--------------------------------------------------------------------------#
    !#                          -- TYPE DECLARATIONS --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    CHARACTER( * ),           INTENT( IN )  :: Filename
    TYPE( TauCoeff_SARTA_TraceGas_type ),    INTENT( IN )  :: TauCoeff_TraceGas

    ! -- Optional input
    INTEGER,        OPTIONAL, INTENT( IN )  :: No_File_Close

    ! -- Revision control
    CHARACTER( * ), OPTIONAL, INTENT( OUT ) :: RCS_Id

    ! -- Error Message log file
    CHARACTER( * ), OPTIONAL, INTENT( IN )  :: Message_Log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: Error_Status


    ! -------------------
    ! Function parameters
    ! -------------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'Write_TauCoeff_SARTA_TraceGas_Binary(Scalar)'
    CHARACTER( * ), PARAMETER :: FILE_STATUS_ON_ERROR = 'DELETE'


    ! ------------------
    ! Function variables
    ! ------------------

    CHARACTER( 256 ) :: Message
    LOGICAL :: Yes_File_Close
    INTEGER :: IO_Status
    INTEGER :: FileID



    !#--------------------------------------------------------------------------#
    !#                     -- INITIALISE THE ERROR STATUS --                    #
    !#--------------------------------------------------------------------------#

    Error_Status = SUCCESS



    !#--------------------------------------------------------------------------#
    !#                 -- SET THE RCS ID ARGUMENT IF SUPPLIED --                #
    !#--------------------------------------------------------------------------#

    IF ( PRESENT( RCS_Id ) ) THEN
      RCS_Id = ' '
      RCS_Id = MODULE_RCS_ID
    END IF



    !#--------------------------------------------------------------------------#
    !#                            -- CHECK INPUT --                             #
    !#--------------------------------------------------------------------------#

    ! ---------------------------------------------
    ! Check that the file is open. If not, open it.
    ! Otherwise get the file ID.
    ! ---------------------------------------------

    IF ( .NOT. File_Open( FileName ) ) THEN

      Error_Status = Open_Binary_File( TRIM( Filename ), &
                                       FileID, &
                                       For_Output  = SET, &
                                       Message_Log = Message_Log )

      IF ( Error_Status /= SUCCESS ) THEN
        Message = 'Error opening '//TRIM( Filename )
        GOTO 2000  ! Clean up
      END IF


    ELSE

      ! -- Inquire for the logical unit number
      INQUIRE( FILE = Filename, NUMBER = FileID )

      ! -- Ensure it's valid
      IF ( FileID == -1 ) THEN
        Message = 'Error inquiring '//TRIM( Filename )//' for its FileID'
        GOTO 1000  ! Clean up
      END IF

    END IF


    !#--------------------------------------------------------------------------#
    !#             -- CHECK STRUCTURE POINTER ASSOCIATION STATUS --             #
    !#                                                                          #
    !#                ALL structure pointers must be associated                 #
    !#--------------------------------------------------------------------------#

    IF ( .NOT. Associated_TauCoeff_TraceGas( TauCoeff_TraceGas ) ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'Some or all INPUT TauCoeff_SARTA_TraceGas pointer '//&
                            'members are NOT associated.', &
                            Error_Status,    &
                            Message_Log = Message_Log )
      RETURN
    END IF


    ! ---------------------------------------
    ! Check the TauCoeff_SARTA_TraceGas structure dimensions
    ! ---------------------------------------

    IF ( TauCoeff_TraceGas%n_Layers           < 1 .OR. &
         TauCoeff_TraceGas%n_Predictors < 1 .OR. &
         TauCoeff_TraceGas%n_Channels         < 1      ) THEN
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            'One or more dimensions of TauCoeff_SARTA_TraceGas structure are < or = 0.', &
                            Error_Status, &
                            Message_Log = Message_Log )
      RETURN
    END IF

    ! -------------------------------------------
    ! Process the optional No_File_Close argument
    ! -------------------------------------------

    ! -- Default action is to close the file on exit....
    Yes_File_Close = .TRUE.

    ! -- ...unless the No_File_Close optional argument is set.
    IF ( PRESENT( No_File_Close ) ) THEN
      IF ( No_File_Close == SET ) Yes_File_Close = .FALSE.
    END IF


    !#--------------------------------------------------------------------------#
    !#                     -- WRITE THE NUMBER OF TRACE GASES --                #
    !#--------------------------------------------------------------------------#

    WRITE( FileID, IOSTAT = IO_Status ) 1

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error writing n_Clouds data dimension to ", a, &
                        &". IOSTAT = ", i5 )' ) &
                      TRIM( Filename ), IO_Status
      CLOSE( FileID, STATUS = FILE_STATUS_ON_ERROR )
      GOTO 1000
    END IF


    !#--------------------------------------------------------------------------#
    !#                        -- WRITE THE STRUCTURE DATA --                    #
    !#--------------------------------------------------------------------------#

    Error_Status = Write_TauCoeff_TraceGas_Record( FileID, &
                                       TauCoeff_TraceGas, &
                                       Message_Log = Message_Log )

    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error writing TauCoeff_SARTA_TraceGas record to '//TRIM( Filename )
      GOTO 1000
    END IF

    !#--------------------------------------------------------------------------#
    !#                 -- CLOSE THE OUTPUT FILE IF REQUIRED --                  #
    !#--------------------------------------------------------------------------#

    IF ( Yes_File_Close ) THEN

      CLOSE( FileID, IOSTAT = IO_Status )

      IF ( IO_Status /= 0 ) THEN
        WRITE( Message, '( "Error closing ", a, ". IOSTAT = ", i5 )' ) &
                        TRIM( Filename ), IO_Status
        CALL Display_Message( ROUTINE_NAME, &
                              TRIM( Message ), &
                              WARNING, &
                              Message_Log = Message_Log )
      END IF

    END IF

    RETURN


    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    !#                      -= CLEAN UP AFTER AN ERROR -=                       #
    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    1000 CONTINUE
    CLOSE( FileID, IOSTAT = IO_Status )

    2000 CONTINUE
    Error_Status = FAILURE
    CALL Display_Message( ROUTINE_NAME, &
                          TRIM( Message ), &
                          Error_Status, &
                          Message_Log = Message_Log )

  END FUNCTION Write_TauCoeff_TraceGas_Scalar
  

  FUNCTION Write_TauCoeff_TraceGas_Rank1( Filename,     &  ! Input
                                  TauCoeff_TraceGas,    &  ! Input
                                  No_File_Close,      &  ! Optional input
                                  RCS_Id,             &  ! Revision control
                                  Message_Log )       &  ! Error messaging
                                RESULT ( Error_Status )



    !#--------------------------------------------------------------------------#
    !#                          -- TYPE DECLARATIONS --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    CHARACTER( * ),           INTENT( IN )  :: Filename
    TYPE( TauCoeff_SARTA_TraceGas_type ), DIMENSION( : ), INTENT( IN )  :: TauCoeff_TraceGas

    ! -- Optional input
    INTEGER,        OPTIONAL, INTENT( IN )  :: No_File_Close 

    ! -- Revision control
    CHARACTER( * ), OPTIONAL, INTENT( OUT ) :: RCS_Id

    ! -- Error Message log file
    CHARACTER( * ), OPTIONAL, INTENT( IN )  :: Message_Log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: Error_Status


    ! -------------------
    ! Function parameters
    ! -------------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'Write_TauCoeff_SARTA_TraceGas_Binary(Scalar)'
    CHARACTER( * ), PARAMETER :: FILE_STATUS_ON_ERROR = 'DELETE'


    ! ------------------
    ! Function variables
    ! ------------------

    CHARACTER( 256 ) :: Message
    LOGICAL :: Yes_File_Close 
    INTEGER :: IO_Status
    INTEGER :: FileID
    INTEGER :: m



    !#--------------------------------------------------------------------------#
    !#                     -- INITIALISE THE ERROR STATUS --                    #
    !#--------------------------------------------------------------------------#

    Error_Status = SUCCESS


    !#--------------------------------------------------------------------------#
    !#                 -- SET THE RCS ID ARGUMENT IF SUPPLIED --                #
    !#--------------------------------------------------------------------------#
 
    IF ( PRESENT( RCS_Id ) ) THEN
      RCS_Id = ' '
      RCS_Id = MODULE_RCS_ID
    END IF


    !#--------------------------------------------------------------------------#
    !#                            -- CHECK INPUT --                             #
    !#--------------------------------------------------------------------------#

    ! ---------------------------------------------
    ! Check that the file is open. If not, open it.
    ! Otherwise get the file ID.
    ! ---------------------------------------------

    IF ( .NOT. File_Open( FileName ) ) THEN

      Error_Status = Open_Binary_File( TRIM( Filename ), &
                                       FileID, &
                                       For_Output  = SET, &
                                       Message_Log = Message_Log )

      IF ( Error_Status /= SUCCESS ) THEN
        Message = 'Error opening '//TRIM( Filename )
        GOTO 2000  ! Clean up
      END IF


    ELSE

      ! -- Inquire for the logical unit number
      INQUIRE( FILE = Filename, NUMBER = FileID )

      ! -- Ensure it's valid
      IF ( FileID == -1 ) THEN
        Message = 'Error inquiring '//TRIM( Filename )//' for its FileID'
        GOTO 1000  ! Clean up
      END IF

    END IF

    ! ---------------------------------------
    ! Check the TauCoeff_SARTA_TraceGas structure dimensions
    ! ---------------------------------------

    IF ( ANY( TauCoeff_TraceGas%n_Layers           < 1 ) .OR. &
         ANY( TauCoeff_TraceGas%n_Predictors < 1 ) .OR. &
         ANY( TauCoeff_TraceGas%n_Channels         < 1 )     ) THEN
      Message = 'One or more dimensions of TauCoeff_SARTA_TraceGas structure are < or = 0.'
      GOTO 1000
    END IF

    ! -------------------------------------------
    ! Process the optional No_File_Close argument
    ! -------------------------------------------

    ! -- Default action is to close the file on exit....
    Yes_File_Close = .TRUE.

    ! -- ...unless the No_File_Close optional argument is set.
    IF ( PRESENT( No_File_Close ) ) THEN
      IF ( No_File_Close == SET ) Yes_File_Close = .FALSE.
    END IF


    !#--------------------------------------------------------------------------#
    !#                     -- WRITE THE NUMBER OF TRACE GASES --                    #
    !#--------------------------------------------------------------------------#

    WRITE( FileID, IOSTAT = IO_Status ) SIZE(TauCoeff_TraceGas)

    IF ( IO_Status /= 0 ) THEN
      WRITE( Message, '( "Error writing n_Clouds data dimension to ", a, &
                        &". IOSTAT = ", i5 )' ) &
                      TRIM( Filename ), IO_Status
      CLOSE( FileID, STATUS = FILE_STATUS_ON_ERROR )
      GOTO 2000
    END IF


    !#--------------------------------------------------------------------------#
    !#                        -- WRITE THE STRUCTURE DATA --                    #
    !#--------------------------------------------------------------------------#
    
    TauCoeff_TraceGas_Loop: DO m = 1, SIZE( TauCoeff_TraceGas )

      Error_Status = Write_TauCoeff_TraceGas_Record( FileID, &			    
        				 TauCoeff_TraceGas(m), &			    
        				 Message_Log = Message_Log )		    

      IF ( Error_Status /= SUCCESS ) THEN					    
        WRITE( Message, '( "Error writing Cloud element #", i5, " to ", a )' ) &
                      m, TRIM( Filename )
        GOTO 1000								    
      END IF									    
     
    END DO TauCoeff_TraceGas_Loop
  

    !#--------------------------------------------------------------------------#
    !#                 -- CLOSE THE OUTPUT FILE IF REQUIRED --                  #
    !#--------------------------------------------------------------------------#

    IF ( Yes_File_Close ) THEN

      CLOSE( FileID, IOSTAT = IO_Status )

      IF ( IO_Status /= 0 ) THEN
        WRITE( Message, '( "Error closing ", a, ". IOSTAT = ", i5 )' ) &
                        TRIM( Filename ), IO_Status
        CALL Display_Message( ROUTINE_NAME, &
                              TRIM( Message ), &
                              WARNING, &
                              Message_Log = Message_Log )
      END IF

    END IF

    RETURN



    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
    !#                      -= CLEAN UP AFTER AN ERROR -=                       #
    !#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

    1000 CONTINUE
    CLOSE( FileID, IOSTAT = IO_Status )

    2000 CONTINUE
    Error_Status = FAILURE
    CALL Display_Message( ROUTINE_NAME, &
                          TRIM( Message ), &
                          Error_Status, &
                          Message_Log = Message_Log )

  END FUNCTION Write_TauCoeff_TraceGas_Rank1

END MODULE TauCoeff_SA_TraceGas_Binary_IO



 
