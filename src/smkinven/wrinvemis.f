
        SUBROUTINE WRINVEMIS( FNAME )

C***********************************************************************
C  subroutine body starts at line 
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
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

        EXTERNAL      CRLF, ENVYN

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME   ! logical file name of output file

C.........  Inventory temporay arrays
        INTEGER, ALLOCATABLE:: IPPTR ( : ) ! position in POLVAL sparse array
        INTEGER, ALLOCATABLE:: IPMAX ( : ) ! max IPPTR by source
        REAL   , ALLOCATABLE:: SRCPOL( :,: )!  pollutant-spec values by source

C...........   Other local variables
        INTEGER         S, L2            ! counters and indices
        INTEGER         IOS              ! i/o status

        INTEGER         RIMISS3          ! real value of integer missing

        LOGICAL      :: EFLAG = .FALSE.  ! TRUE iff ERROR
        LOGICAL         SFLAG            ! true: error on missing species

        CHARACTER*100   BUFFER           !  message buffer
        CHARACTER*300   MESG             !  message buffer

        CHARACTER*16 :: PROGNAME = 'WRINVEMIS' !  program name

C***********************************************************************
C   begin body of program WRINVEMIS

C.........  Get environment variables that control the program
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

C.........  Loop through pollutants, store, and write to inventory file
        CALL LOOP_FOR_OUTPUT( NIPOL, EIIDX, EINAM )

C.........  Loop through activity data, store, and write to inventory file
        CALL LOOP_FOR_OUTPUT( NIACT, AVIDX, ACTVTY )

        IF( EFLAG ) THEN
            MESG = 'Missing data for some sources is not allowed ' //
     &             CRLF() // BLANK5 //
     &             'because the environment variable ' //
     &             'RAW_SRC_CHECK was set to "N".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Deallocate local arrays
        DEALLOCATE( IPPTR, IPMAX, SRCPOL )

C.........  Deallocate global emissions arrays
        CALL SRCMEM( CATEGORY, 'SORTED', .FALSE., .TRUE., 1, 1, 1 )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram is for writing out the inventory
C               data, whether it is the pollutant data or the activity data.
C.............  Most variables are defined from the main subroutine
            SUBROUTINE LOOP_FOR_OUTPUT( NOUT, INDX, NAMES )

C.............  Subroutine arguments 
            INTEGER     , INTENT (IN) :: NOUT          ! no. pols/act for output
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

                        DO J = 1, NPPOL     ! rearrange pollutant-specific info
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

                CALL WRINVPOL( FNAME, CATEGORY, NSRC, 1, NPPOL, 
     &                         NAMES( I ), SRCPOL )

            END DO  ! end loop through actual output data

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

            END SUBROUTINE LOOP_FOR_OUTPUT
 
        END SUBROUTINE WRINVEMIS
