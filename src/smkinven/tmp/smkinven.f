
         PROGRAM RAWPOINT

C***********************************************************************
C  program body starts at line 151
C
C  DESCRIPTION:
C    The rawpoint program reads the point source inventory in one of four
C    formats: EMS-95, EPS, IDA, and SMOKE list format.  It permits a flexible
C    definition of a point source, which depends on the inventory input
C    formats.  It allows any number of inventory pollutants, within the limit
C    of I/O API (this works out to only 15 pollutants per file). 
C
C  PRECONDITIONS REQUIRED:
C    Set environment variables:
C      PROMPTFLAG:    If N, default inputs are used
C      RAW_SRC_CHECK: If Y, missing species disallowed
C      RAW_DUP_CHECK: If Y, duplicate sources disallowed
C      VELOC_RECALC:  If Y, recalculates velocity based on flow
C      WEST_HSPHERE:  If N, does not reset positive longitudes to negative
C    Input files:  
C      PTINV: ASCII point sources inventory
C      PSTK: Replacement stack parameters file
C      ZONES: Time zones files
C      SIPOLS: Master list of pollutant codes and names (in output order)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C    Subroutines: I/O API subroutines, INITEM, CHECKMEM, RDTZONE, RDSIPOLS, 
C       RDPTINV, FMTCSRC, FIXSTK, WRPTSCC, WRPTREF, OPENPNTS, WPNTSCHR, 
C       WPNTSPOL
C    Functions: I/O API functions, GETFLINE, GETTZONE 
C
C  REVISION  HISTORY:
C    copied by: mhouyoux 10/98
C    origin: emspoint.F 4.3
C
C****************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2             CRLF
        LOGICAL                 ENVYN
        INTEGER                 GETFLINE
        INTEGER                 GETTZONE
        INTEGER                 PROMPTFFILE

        EXTERNAL CRLF, ENVYN, GETFLINE, GETTZONE, PROMPTFFILE

C...........  LOCAL PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C.........  LOCAL VARIABLES and their descriptions:
C.........  Time zone tables:  FIP-independent; state-only; state-county

        INTEGER                TZONE0
        INTEGER                NZS         !  no of state-specific time zones
        INTEGER                NZF         !  no of FIP-specific time zones
        INTEGER, ALLOCATABLE:: TZONST( : ) !  state-specific zones
        INTEGER, ALLOCATABLE:: TFIPST( : ) !  state FIPS codes (2 digit)
        INTEGER, ALLOCATABLE:: TZONEF( : ) !  fip-specific zones
        INTEGER, ALLOCATABLE:: TFIPEF( : ) !  state/county FIPS codes (5 digit)

C.........  Full list of inventory pollutants (in output order)
        INTEGER                     MXIPOL       ! Max no of inv pollutants
        INTEGER     , ALLOCATABLE:: INVPCOD( : ) ! 5-digit pollutant code
        INTEGER     , ALLOCATABLE:: INVSTAT( : ) ! Status (0=not in inventory)
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: INVPNAM( : ) ! Name of pollutant

C.........  Inventory pollutants actually in the inventory
        INTEGER                               NIPOL      ! Actual no of inv pols
        INTEGER               , ALLOCATABLE:: EIIDX( : ) ! pos in full inven arr
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: EINAM( : ) ! Name of actual pols

C.........  Inventory temporay arrays
        INTEGER, ALLOCATABLE:: IPPTR ( : ) ! position in POLVAL sparse array
        INTEGER, ALLOCATABLE:: IPMAX ( : ) ! max IPPTR by source
        INTEGER, ALLOCATABLE:: IDXSCC( : ) ! sorting index for output SCC file
        REAL   , ALLOCATABLE:: SRCPOL( :,: )!  pollutant-spec values by source

C.........  File units and logical/physical names

        INTEGER         ADEV    !  Unit number for output actual SCCs
        INTEGER         IDEV    !  Inventory file (various formats)
        INTEGER         LDEV    !  log-device
        INTEGER         PDEV    !  Unit number for pollutants codes/names file
        INTEGER         RDEV    !  Unit number for default stack parameters
        INTEGER         TDEV    !  Unit number for output temporal profile #s
        INTEGER         SDEV    !  for ASCII output inventory file
        INTEGER         ZDEV    !  for time zone file

        CHARACTER(LEN=NAMLEN3) ENAME !  emissions output inventory logical name

C...........   Other local variables
                                
        INTEGER         S, I, J, K, L, LK, LS, V !  counters and indices
        INTEGER         L1, L2           !  counters and indices

        INTEGER         FIP     !  Temporary FIPS code
        INTEGER         INVFMT  !  Inventory format code
        INTEGER         IOS     !  I/O status
        INTEGER         MAXK    !  test for maximum value of K in output loop
        INTEGER         NFIPLIN !  number of lines in ZDEV
        INTEGER         NPSRC   !  actual source count
        INTEGER         TZONE   !  tmp time zone

        LOGICAL      :: EFLAG = .FALSE.  !  TRUE iff ERROR
        LOGICAL         SFLAG   !  input verification:  report missing species
        LOGICAL         TFLAG   !  TRUE if temporal x-ref output

        CHARACTER*300   BUFFER  !  input line from POINT file
        CHARACTER*300   MESG    !  text for M3EXIT()

        CHARACTER*16  :: PROGNAME = 'RAWPOINT'   !  program name

C***********************************************************************
C   begin body of program RAWPOINT

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get environment variables that control the program
        SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                 'Flag to check for missing species-records',
     &                 .FALSE., IOS )

C.........  Get file name for opening input raw point source file
        IDEV = PROMPTFFILE( 
     &         'Enter the name of the RAW POINT INVENTORY file',
     &          .TRUE., .TRUE., 'PTINV', PROGNAME )

C.........  Get file name for input replacement stack parameters file
        RDEV = PROMPTFFILE(
     &          'Enter REPLACEMENT STACK PARAMETERS file',
     &          .TRUE., .TRUE., 'PSTK', PROGNAME )

C.........  Get file name for time zones files
        ZDEV = PROMPTFFILE( 
     &         'Enter logical name for TIME ZONE file',
     &         .TRUE., .TRUE., 'ZONES', PROGNAME )

C.........  Get file name for inventory pollutants codes/names
        PDEV = PROMPTFFILE( 
     &         'Enter logical name for POLLUTANT CODES & NAMES file',
     &         .TRUE., .TRUE., 'SIPOLS', PROGNAME )

C.........  Allocate memory for time zone tables based on no. of lines in file
        NFIPLIN = GETFLINE( ZDEV, 'Time zones file')

C.........  Get number of lines of pollutant codes/names file 
        MXIPOL = GETFLINE( PDEV, 'Pollutant codes and names file')

        CALL M3MSG2( 'Setting up to read inventory data...' )

C.........  Allocate memory for time zones tables 
        ALLOCATE( TZONST( NFIPLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONST', PROGNAME )
        ALLOCATE( TFIPST( NFIPLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TFIPST', PROGNAME )
        ALLOCATE( TZONEF( NFIPLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONEF', PROGNAME )
        ALLOCATE( TFIPEF( NFIPLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TFIPEF', PROGNAME )

C.........  Allocate memory for storing contents of pollutants file
        ALLOCATE( INVPCOD( MXIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVPCOD', PROGNAME )
        ALLOCATE( INVPNAM( MXIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVPNAM', PROGNAME )
        ALLOCATE( INVSTAT( MXIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTAT', PROGNAME )

C.........  Read time zone file, separate into categories, and sort tables
        CALL RDTZONE( ZDEV, NFIPLIN, NZS, NZF, TZONE0,
     &                TZONST, TFIPST, TZONEF, TFIPEF )

C.........  Read and sort pollutant codes/names file
        CALL RDSIPOLS( PDEV, MXIPOL, INVPCOD, INVPNAM )

C.........  Initialize pollutant status (present in inventory or not)
        INVSTAT = 0  ! array

        CALL M3MSG2( 'Reading raw inventory data...' )

C.........  Read the raw inventory data, and store in sorted order
C.........  The inventory arrays that are populated by this subroutine call
C           are contained in the module MODSOURC

        CALL RDPTINV( IDEV, MXIPOL, INVPCOD, INVPNAM, INVFMT, 
     &                NPSRC, TFLAG, EFLAG, INVSTAT )

        IF( EFLAG ) THEN
           MESG = 'Error reading raw inventory file(s)'         
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        CALL M3MSG2( 'Processing inventory data...' )

C.........  Determine memory needed for actual pollutants list and allocate it
        NIPOL = 0
        DO I = 1, MXIPOL
            IF( INVSTAT( I ) .NE. 0 ) NIPOL = NIPOL + 1
        ENDDO

        ALLOCATE( EIIDX( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EIIDX', PROGNAME )
        ALLOCATE( EINAM( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EINAM', PROGNAME )

C.........  Create list of actual pollutants and an index to the master list.
C           The order in EINAM will be the output order. The index is for
C           accessing INVPCOD and INVSTAT if needed.
        J = 0
        DO I = 1, MXIPOL
           IF( INVSTAT( I ) .NE. 0 ) THEN
               J = J + 1
               EIIDX( J ) = I
               EINAM( J ) = INVPNAM( I )
           ENDIF
        ENDDO

C.........   Fix stack parameters. 
C.........   Some of these arguments are variables that are defined in the
C            module MODSOURC
        CALL FIXSTK( RDEV, NPSRC )

C.........  Set time zones based on state and county FIPS code. Note that a
C           few counties in the Western U.S. are divided by a time zone, so this
C           is not perfectly accurate for all counties.
        ALLOCATE( TZONES( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONES', PROGNAME )

        DO S = 1, NPSRC
            FIP   = IFIP( S )
            TZONE = GETTZONE( FIP, NZS, NZF, TZONE0,
     &                        TZONST, TFIPST, TZONEF, TFIPEF )

            TZONES( S ) = TZONE

        ENDDO

C.........  Allocate memory for local source-specific arrays used for output      
C.........  Allocate memory for indices IPPTR & IPMAX for pointing to position
C           in sparsely stored emissions array.  
        ALLOCATE( IPPTR( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPPTR', PROGNAME )
        ALLOCATE( IPMAX( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPMAX', PROGNAME )
        ALLOCATE( SRCPOL( NPSRC, NPTPPOL3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCPOL', PROGNAME )
       
C.........  Get unique-SCC output file
        ADEV = PROMPTFFILE( 
     &          'Enter the name of the ACTUAL SCC output file',
     &          .FALSE., .TRUE., 'PSCC', PROGNAME )

C.........  Write out SCCs list (this subroutine expect the data structure
C           that is being provided.

        CALL M3MSG2( 'Writing out ACTUAL SCC file...' )

        CALL WRCHRSCC( ADEV, NPSRC, CSCC )

C.........  Write out temporal x-ref file. (TFLAG is true for EMS-95 format)

        IF( TFLAG ) THEN

           TDEV = PROMPTFFILE( 
     &            'Enter the name of the TEMPORAL X-REF output file',
     &            .FALSE., .TRUE., 'PTREF', PROGNAME )

           CALL M3MSG2( 'Writing out TEMPORAL CROSS-REFERENCE file...' )

           CALL WRPTREF( TDEV, NPSRC, IDIU, IWEK, IWEK ) ! no monthly 

        ENDIF

C.........  Open output I/O API and ASCII files for PNTS
        CALL OPENPNTS( NPSRC, NIPOL, EINAM, ENAME, SDEV )

        CALL M3MSG2( 'Writing SMOKE POINT SOURCE INVENTORY file...' )

C.........  Write source characteristics to PNTS file (I/O API and ASCII)
        CALL WPNTSCHR( ENAME, SDEV, NPSRC )

C.........  Deallocate memory to potentially speed the rest of the program up
        DEALLOCATE( IFIP, CSCC, ISIC, IORIS, TZONES, TPFLAG, 
     &              INVYR, XLOCA, YLOCA, STKHT, STKDM, STKTK, STKVE, 
     &              CBLRID, CPDESC )

C.........  Initialize global index based on count of number of sources per 
C           pollutant
C.........  Since the emission values are already sorted in the output order of
C           the pollutants, can use the IPPTR and IPMAX indices to keep track 
C           of what pollutant we are on for each source (sources can have 
C           different numbers and types of output pollutants).
        IPPTR( 1 ) = 1
        IPMAX( 1 ) = NPCNT( 1 )
        DO S = 2, NPSRC
            IPPTR( S ) = IPPTR( S-1 ) + NPCNT( S-1 )
            IPMAX( S ) = IPPTR( S )   + NPCNT( S ) - 1
        ENDDO

C.........  Loop through inventory pollutants, store, and write to PNTS file
        DO I = 1, NIPOL

C.............  Initialize SRCPOL pollutant-specific data array
            SRCPOL = 0. ! array

C.............  Transfer emissions data to output SRCPOL array
            K = 0
            DO S = 1, NPSRC

C.................  Set position in sparse array by compring pointer to max
                K = MIN( IPPTR( S ), IPMAX( S ) ) 

C.................  Retrieve emissions from the pollutant array
C                   if the source pollutant ID equals the pollutant ID of the 
C                   pollutant of iteration I.
                IF( IPOSCOD( K ) .EQ. EIIDX( I ) ) THEN

                    IPPTR( S ) = K + 1  ! increase pointer for this source by 1

                    DO J = 1, NPTPPOL3    ! store pollutant-specific info
                        SRCPOL( S,J ) = POLVAL( K,J )
                    ENDDO

C.....................  If emissions data records are fatal, then set error flag
                    IF( SFLAG .AND. SRCPOL( S,1 ) .LT. 0 ) THEN
                        EFLAG = .TRUE.
                        CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )
                        MESG = 'ERROR: Missing emissions for source:' //
     &                         CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                    END IF

                END IF

            ENDDO  ! end of loop through sources x pollutants

            CALL WPNTSPOL( ENAME, NPSRC, 1, EINAM( I ), SRCPOL )

        ENDDO  ! end loop through actual output pollutants

        IF( EFLAG ) THEN
            MESG = 'Missing species for some sources but this is' //
     &             CRLF() // BLANK5 //
     &             'not allowed because the environment variable ' //
     &             'RAW_SRC_CHECK was set to "N".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  End program successfully
        MESG = ' '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C93041   FORMAT( I5, X, I5, X, I3, X, I8, X, I5.5, 3( X, I3 ) )

93060   FORMAT( 10( A, :, E10.3, :, 1X ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,
     &          A, I6, 1X, A, I6, 1X, A, :, I6 )

94040   FORMAT( A, I2.2 )

94060   FORMAT( 10( A, :, E10.3, :, 1X ) )

94080   FORMAT( '************  ', A, I7, ' ,  ' , A, I12 )
 
        END PROGRAM RAWPOINT
