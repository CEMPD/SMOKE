
        SUBROUTINE WRINVEMIS( FNAME )

C***********************************************************************
C  subroutine body starts at line 111
C
C  DESCRIPTION:
C      This subroutine writes the average inventory emissions to the inventory
C      files
C
C  PRECONDITIONS REQUIRED:
C      Logical name of output file defined
C      Emission arrays populated
C      MODINFO values assigned
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
C
C*************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  EXTERNAL FUNCTIONS
        CHARACTER*2   CRLF
        LOGICAL       ENVYN
        INTEGER       INDEX1

        EXTERNAL      CRLF, ENVYN, INDEX1

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME   ! logical file name of output file

C...........   LOCAL PARAMETERS
        CHARACTER*16, PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'

C.........  Inventory temporay arrays
        INTEGER, ALLOCATABLE:: IPPTR ( : ) ! position in POLVAL sparse array
        INTEGER, ALLOCATABLE:: IPMAX ( : ) ! max IPPTR by source
        REAL   , ALLOCATABLE:: SRCPOL( :,: )  ! data-spec values by source
        REAL   , ALLOCATABLE:: COMPUTED( :,: )! computed data-spec values by src

C...........   Other local variables
        INTEGER         S, L, L2, VA, VB     ! counters and indices

        INTEGER         IOS       ! i/o status
        INTEGER         LEQU      ! position of '=' in formula
        INTEGER         LDIV      ! position of '-' or '+' in formula
        INTEGER         LMNS      ! position of '-' in formula
        INTEGER         LPLS      ! position of '+' in formula

        INTEGER         RIMISS3          ! real value of integer missing

        LOGICAL      :: CHKPLUS  = .FALSE. ! true: formula uses a + sign
        LOGICAL      :: CHKMINUS = .FALSE. ! true: formula uses a - sign
        LOGICAL      :: EFLAG    = .FALSE. ! true: error found
        LOGICAL      :: FFLAG    = .FALSE. ! true: formula in use
        LOGICAL      :: SFLAG    = .FALSE. ! true: error on missing species

        CHARACTER*60    VAR_FORMULA      ! formula
        CHARACTER*100   BUFFER           ! message buffer
        CHARACTER*300   MESG             ! message buffer

        CHARACTER(LEN=IOVLEN3 ) VIN_A
        CHARACTER(LEN=IOVLEN3 ) VIN_B
        CHARACTER(LEN=IOVLEN3 ) VNAME

        CHARACTER*16 :: PROGNAME = 'WRINVEMIS' !  program name

C***********************************************************************
C   begin body of program WRINVEMIS

C.........  Get environment variables that control the program
        CALL ENVSTR( FORMEVNM, MESG, ' ', VAR_FORMULA, IOS )

        SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                 'Flag to check for missing species-records',
     &                 .FALSE., IOS )

C.........  Compute real value of integer missing
        RIMISS3 = REAL( IMISS3 )

C.........  Allocate memory for local source-specific arrays used for output   
C.........  Allocate memory for indices IPPTR & IPMAX for pointing to position
C           in sparsely stored data array.  
        ALLOCATE( IPPTR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPPTR', PROGNAME )
        ALLOCATE( IPMAX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPMAX', PROGNAME )
        ALLOCATE( SRCPOL( NSRC, NPTPPOL3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCPOL', PROGNAME )
       
C.........  Initialize global index based on count of number of sources per 
C           pollutant
C.........  Since the emission values are already sorted in the output order of
C           the pollutants/activities, can use the IPPTR and IPMAX indices 
C           to keep track of what pollutant/activity we are on for each source
C           (sources can have different numbers and types of output pollutants
C           and activities).
        IPPTR( 1 ) = 1
        IPMAX( 1 ) = NPCNT( 1 )
        DO S = 2, NSRC
            IPPTR( S ) = IPPTR( S-1 ) + NPCNT( S-1 )
            IPMAX( S ) = IPPTR( S )   + NPCNT( S ) - 1
        ENDDO

C.........  If there is a computed output variable, get set up for that
        IF( VAR_FORMULA .NE. ' ' ) THEN

            FFLAG = .TRUE.

C.............  Make sure formula makes sense
            LEQU = INDEX( VAR_FORMULA, '=' )
            LPLS = INDEX( VAR_FORMULA, '+' )
            LMNS = INDEX( VAR_FORMULA, '-' )

            CHKPLUS  = ( LPLS .GT. 0 )
            CHKMINUS = ( LMNS .GT. 0 )

            LDIV = LPLS
            IF( CHKMINUS ) LDIV = LMNS

            IF( LEQU .LE. 0 .OR. 
     &        ( .NOT. CHKPLUS .AND. .NOT. CHKMINUS ) ) THEN

                L = LEN_TRIM( FORMEVNM )
                MESG = 'Could not interpret formula for extra ' //
     &                 'pollutant from environment variable ' //
     &                 CRLF() // BLANK10 // '"' // FORMEVNM( 1:L ) //
     &                 '": ' // VAR_FORMULA
                EFLAG = .TRUE.
            END IF

C.............  Extract formula variable names
            L     = LEN_TRIM( VAR_FORMULA )
            VNAME = ADJUSTL ( VAR_FORMULA(      1:LEQU-1 ) )
            VIN_A = ADJUSTL ( VAR_FORMULA( LEQU+1:LDIV-1 ) )
            VIN_B = ADJUSTL ( VAR_FORMULA( LDIV+1:L      ) )

C.............  Find formula inputs in existing variable list
            VA = INDEX1( VIN_A, NIPPA, EANAM )
            VB = INDEX1( VIN_B, NIPPA, EANAM )

            IF( VA .LE. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( VIN_A )
                MESG = 'Variable "'// VIN_A( 1:L ) // 
     &                 '" from formula was not found in inventory.'
                CALL M3MSG2( MESG )
            END IF

            IF( VB .LE. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( VIN_B )
                MESG = 'Variable "'// VIN_B( 1:L ) // 
     &                 '" from formula was not found in inventory.'
                CALL M3MSG2( MESG )
            END IF

            IF( EFLAG ) THEN
                MESG = 'Problem processing formula.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Allocate memory for computed variable
            ALLOCATE( COMPUTED( NSRC, NPTPPOL3 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'COMPUTED', PROGNAME )
            COMPUTED = 0.   ! array

        END IF

C.........  Loop through pollutants, store, and write to inventory file
        IF( NIPOL .GT. 0 ) 
     &      CALL LOOP_FOR_OUTPUT( NIPOL, NPPOL, EIIDX, EINAM )

C.........  Loop through activity data, store, and write to inventory file
        IF( NIACT .GT. 0 ) 
     &      CALL LOOP_FOR_OUTPUT( NIACT, NPACT, AVIDX, ACTVTY )

        IF( EFLAG ) THEN
            MESG = 'Missing data for some sources is not allowed ' //
     &             CRLF() // BLANK5 //
     &             'because the environment variable ' //
     &             'RAW_SRC_CHECK was set to "N".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If needed, check for negative values and output computed variable
        IF( FFLAG ) THEN

C.............  Check for negative values
            DO S = 1, NSRC

                IF( COMPUTED( S,1 ) .LT. 0 ) THEN
                    L = LEN_TRIM( VNAME )
                    CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )
                    WRITE( MESG,94020 ) 
     &                     'WARNING: Resetting negative value of "' //
     &                     VNAME( 1:L ) // '" from', COMPUTED( S,1 ),
     &                     'to 0. for source:'// CRLF() // BLANK5 // 
     &                     BUFFER( 1:L2 )
                    CALL M3MESG( MESG )

                    COMPUTED( S,1 ) = 0.

                END IF

            END DO

C.............  Output data
            CALL WRINVPOL( FNAME, CATEGORY, NSRC, 1, NPPOL, 
     &                     VNAME, COMPUTED )
        END IF

C.........  Deallocate local arrays
        IF( ALLOCATED( IPPTR ) )    DEALLOCATE( IPPTR )
        IF( ALLOCATED( IPMAX ) )    DEALLOCATE( IPMAX )
        IF( ALLOCATED( SRCPOL ) )   DEALLOCATE( SRCPOL )
        IF( ALLOCATED( COMPUTED ) ) DEALLOCATE( COMPUTED )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( 10( A, :, E10.2, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram is for writing out the inventory
C               data, whether it is the pollutant data or the activity data.
C.............  Most variables are defined from the main subroutine
            SUBROUTINE LOOP_FOR_OUTPUT( NOUT, NPVAR, INDX, NAMES )

C.............  Subroutine arguments 
            INTEGER     , INTENT (IN) :: NOUT          ! no. pols/act for output
            INTEGER     , INTENT (IN) :: NPVAR         ! no. vars per data
            INTEGER     , INTENT (IN) :: INDX ( NOUT ) ! index to master list
            CHARACTER(*), INTENT (IN) :: NAMES( NOUT ) ! names of pols/act

C.............  Local variables
            INTEGER    I, J, K, L, S

C----------------------------------------------------------------------

            DO I = 1, NOUT

C.................  Initialize SRCPOL pollutant/activity-specific data array
C                   as missing
                SRCPOL        = BADVAL3 ! array
                SRCPOL( :,1 ) = 0.      ! array
                IF( NOZ .GT. 0 ) SRCPOL( :,NOZ ) = 0.      ! array
                IF( NC1 .GT. 0 ) SRCPOL( :,NC1 ) = RIMISS3 ! array
                IF( NC2 .GT. 0 ) SRCPOL( :,NC2 ) = RIMISS3 ! array
                IF( NCE .GT. 0 ) SRCPOL( :,NCE ) = 0.      ! array
                IF( NRE .GT. 0 ) SRCPOL( :,NRE ) = 100.    ! array
                IF( NRP .GT. 0 ) SRCPOL( :,NRP ) = 100.    ! array

C.................  Transfer emissions or activity data to output SRCPOL array
                K = 0
                DO S = 1, NSRC

C.....................  Set position in sparse array by comparing pointer to max
                    K = MIN( IPPTR( S ), IPMAX( S ) ) 

C.....................  Retrieve emissions from the pollutant array if the 
C                       source pollutant ID equals the pollutant ID of the 
C                       pollutant of iteration I.
                    IF( IPOSCOD( K ) .EQ. INDX( I ) ) THEN

                        IPPTR( S ) = K + 1  ! pointer for source S

                        DO J = 1, NPVAR     ! rearrange pollutant-specific info
                            SRCPOL( S,J ) = POLVAL( K,J )
                        END DO

C.........................  If missing emissions or VMT data records are fatal, 
C                           then write message and set error flag
                        IF( SFLAG .AND. SRCPOL( S,1 ) .LT. 0 ) THEN
                            EFLAG = .TRUE.
                            L = LEN_TRIM( NAMES( I ) )
                            CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )
                            MESG = 'ERROR: Missing data for "' // 
     &                             NAMES( I )( 1:L )// '" for source:'//
     &                             CRLF() // BLANK5 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                        END IF

                    END IF

                END DO  ! end of loop through sources

C.................  Write data to inventory file
                CALL WRINVPOL( FNAME, CATEGORY, NSRC, 1, NPVAR, 
     &                         NAMES( I ), SRCPOL )

C.................  If current data variable is the first variable in the
C                   formula, then store data in formula arrays
                IF( NAMES( I ) .EQ. VIN_A ) THEN
                    COMPUTED( :,1 ) = COMPUTED( :,1 ) + SRCPOL( :,1 )
                    IF ( NPVAR .GT. 1 ) THEN
                        COMPUTED( :,2 )= COMPUTED( :,2 ) + SRCPOL( :,2 )
                    END IF

                    IF ( NPVAR .GT. 2 ) THEN
                        COMPUTED( :,3:NPVAR ) = SRCPOL( :,3:NPVAR )
                    END IF
                END IF

C.................  If current data variable is the second variable in the
C                   formula, then use data in formula to compute output value
                IF( NAMES( I ) .EQ. VIN_B ) THEN

                    IF( CHKPLUS ) THEN
                        COMPUTED( :,1 )= COMPUTED( :,1 ) + SRCPOL( :,1 )
                        IF ( NPVAR .GT. 1 ) THEN
                            COMPUTED(:,2)= COMPUTED(:,2) + SRCPOL(:,2)
                        END IF

                    ELSE IF( CHKMINUS ) THEN
                        COMPUTED( :,1 )= COMPUTED( :,1 ) - SRCPOL( :,1 )
                        IF ( NPVAR .GT. 1 ) THEN
                             COMPUTED(:,2)= COMPUTED(:,2) - SRCPOL(:,2)
                        END IF

                   END IF

                END IF

            END DO  ! end loop through actual output data

            END SUBROUTINE LOOP_FOR_OUTPUT
  
        END SUBROUTINE WRINVEMIS
