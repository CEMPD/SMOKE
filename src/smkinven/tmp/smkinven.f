
         PROGRAM SMKINVEN

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C    The smkinven program reads the any source inventory in one of four
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
C    started 10/98 by M Houyoux as rawpoint.f from emspoint.F 4.3
C    smkinven changes started 4/98
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

C.........  This module contains the information about the source category
        USE MODINFO

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
        INTEGER               , ALLOCATABLE:: EIIDX( : ) ! pos in full inven arr

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

        CHARACTER(LEN=NAMLEN3) ANAME !  inventory ASCII output logical name
        CHARACTER(LEN=NAMLEN3) ENAME !  inventory I/O API output logical name
        CHARACTER(LEN=NAMLEN3) INAME !  inventory input logical name

C...........   Other local variables
                                
        INTEGER         S, I, J, K, L, L2, V !  counters and indices

        INTEGER         FILFMT  !  input file(s) format code
        INTEGER         FIP     !  Temporary FIPS code
        INTEGER         IOS     !  I/O status
        INTEGER         MAXK    !  test for maximum value of K in output loop
        INTEGER         NFIPLIN !  number of lines in ZDEV
        INTEGER         NRAWBP  !  number of sources x pollutants
        INTEGER         TZONE   !  tmp time zone

        REAL            PRATIO  !  position ratio

        LOGICAL      :: EFLAG = .FALSE.  ! TRUE iff ERROR
        LOGICAL         SFLAG            ! input check:  report missing species
        LOGICAL      :: TFLAG = .FALSE.  ! TRUE if temporal x-ref output

        CHARACTER*300   BUFFER  !  input line from POINT file
        CHARACTER*300   MESG    !  text for M3EXIT()

        CHARACTER*16  :: PROGNAME = 'SMKINVEN'   !  program name

        integer bdev
C***********************************************************************
C   begin body of program SMKINVEN

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get environment variables that control the program
        SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                 'Flag to check for missing species-records',
     &                 .FALSE., IOS )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get names of input files
        CALL OPENINVIN( CATEGORY, IDEV, RDEV, PDEV, ZDEV, INAME )

        CALL M3MSG2( 'Setting up to read inventory data...' )

C.........  Get no. lines in time zone file for allocating memory
        NFIPLIN = GETFLINE( ZDEV, 'Time zones file')

C.........  Allocate memory for time zones tables 
        ALLOCATE( TZONST( NFIPLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONST', PROGNAME )
        ALLOCATE( TFIPST( NFIPLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TFIPST', PROGNAME )
        ALLOCATE( TZONEF( NFIPLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONEF', PROGNAME )
        ALLOCATE( TFIPEF( NFIPLIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TFIPEF', PROGNAME )

C.........  Get no. lines in pollutant codes file for allocating memory
        MXIPOL = GETFLINE( PDEV, 'Pollutant codes and names file')

C.........  Allocate memory for storing contents of pollutants file.
C.........  Increase MXIPOL by +1 to append 'VMT' to the list, if needed
        ALLOCATE( INVPCOD( MXIPOL+1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVPCOD', PROGNAME )
        ALLOCATE( INVPNAM( MXIPOL+1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVPNAM', PROGNAME )
        ALLOCATE( INVSTAT( MXIPOL+1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTAT', PROGNAME )

C.........  Read time zone file, separate into categories, and sort tables
        CALL RDTZONE( ZDEV, NFIPLIN, NZS, NZF, TZONE0,
     &                TZONST, TFIPST, TZONEF, TFIPEF   )

C.........  Read and sort pollutant codes/names file
        CALL RDSIPOLS( PDEV, MXIPOL, INVPCOD, INVPNAM )

C.........  Insert VMT as a pollutant name for use by mobile sources
        IF( CATEGORY .EQ. 'MOBILE' ) THEN
            MXIPOL = MXIPOL + 1
            INVPNAM( MXIPOL ) = 'VMT'
            INVPCOD( MXIPOL ) = 0
        END IF

C.........  Initialize pollutant status (present in inventory or not)
        INVSTAT = 0  ! array

        CALL M3MSG2( 'Reading raw inventory data...' )

C.........  Read the raw inventory data, and store in unsorted order
C.........  The arrays that are populated by this subroutine call
C           are contained in the module MODSOURC

        CALL RDINVEN( IDEV, INAME, MXIPOL, INVPCOD, INVPNAM, FILFMT,
     &                NRAWBP, PRATIO, TFLAG )

        CALL M3MSG2( 'Sorting raw inventory data...' )

C.........  Sort inventory and pollutants (sources x pollutants). Note that
C           sources are sorted based on character string definition of the 
C           source so that source definition can be consistent with that of
C           the input format.

        CALL SORTIC( NRAWBP, INDEXA, CSOURCA )

        CALL M3MSG2( 'Processing inventory data...' )

C.........  Processing inventory records and store in sorted order

        CALL PROCINVEN( NRAWBP, MXIPOL, FILFMT, PRATIO, INVSTAT )

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

C.........   Fix stack parameters for point sources
C.........   Some of these arguments are variables that are defined in the
C            module MODSOURC
        IF( CATEGORY .EQ. 'POINT' ) CALL FIXSTK( RDEV, NSRC )

C.........  Set time zones based on country/state/county code. Note that a
C           few counties in the Western U.S. are divided by a time zone, so this
C           is not perfectly accurate for all counties.
        ALLOCATE( TZONES( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONES', PROGNAME )

        DO S = 1, NSRC
            FIP   = IFIP( S )
            TZONE = GETTZONE( FIP, NZS, NZF, TZONE0,
     &                        TZONST, TFIPST, TZONEF, TFIPEF )

            TZONES( S ) = TZONE

        ENDDO

C.........  Get unique-SCC output file
        ADEV = PROMPTFFILE( 
     &          'Enter the name of the ACTUAL SCC output file',
     &          .FALSE., .TRUE., CRL // 'SCC', PROGNAME )

C.........  Write out SCCs list (this subroutine expect the data structure
C           that is being provided.

        CALL M3MSG2( 'Writing out ACTUAL SCC file...' )

        CALL WRCHRSCC( ADEV, NSRC, CSCC )

C.........  Write out temporal x-ref file. (TFLAG is true for EMS-95 format for
C           point sources only)
        IF( TFLAG ) THEN

           TDEV = PROMPTFFILE( 
     &            'Enter the name of the TEMPORAL X-REF output file',
     &            .FALSE., .TRUE., CRL // 'TREF', PROGNAME )

           CALL M3MSG2( 'Writing out TEMPORAL CROSS-REFERENCE file...' )

           CALL WRPTREF( TDEV, NSRC, IDIU, IWEK, IWEK ) ! no monthly 

        ENDIF

C.........  Get output inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Generate message to use just before writing out inventory files
C.........  Open output I/O API and ASCII files, then write source 
C           characteristics to inventory files (I/O API and ASCII)
        CALL OPENINVOUT( ENAME, ANAME, SDEV )

        MESG = 'Writing SMOKE ' // CATEGORY( 1:CATLEN ) // 
     &         ' SOURCE INVENTORY file...'

        CALL M3MSG2( MESG )

        CALL WRINVCHR( ENAME, SDEV )

C.........  Deallocate sorted inventory info arrays
        CALL SRCMEM( CATEGORY, 'SORTED', .FALSE., .FALSE., 1, 1, 1 )

C.........  Allocate memory for local source-specific arrays used for output      
C.........  Allocate memory for indices IPPTR & IPMAX for pointing to position
C           in sparsely stored emissions array.  
        ALLOCATE( IPPTR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPPTR', PROGNAME )
        ALLOCATE( IPMAX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPMAX', PROGNAME )
        ALLOCATE( SRCPOL( NSRC, NPTPPOL3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCPOL', PROGNAME )
       
C.........  Initialize global index based on count of number of sources per 
C           pollutant
C.........  Since the emission values are already sorted in the output order of
C           the pollutants, can use the IPPTR and IPMAX indices to keep track 
C           of what pollutant we are on for each source (sources can have 
C           different numbers and types of output pollutants).
        IPPTR( 1 ) = 1
        IPMAX( 1 ) = NPCNT( 1 )
        DO S = 2, NSRC
            IPPTR( S ) = IPPTR( S-1 ) + NPCNT( S-1 )
            IPMAX( S ) = IPPTR( S )   + NPCNT( S ) - 1
        ENDDO

C.........  Loop through inventory pollutants, store, and write to PNTS file
        DO I = 1, NIPOL

C.............  Initialize SRCPOL pollutant-specific data array
            SRCPOL = 0. ! array

C.............  Transfer emissions data to output SRCPOL array
            K = 0
            DO S = 1, NSRC

C.................  Set position in sparse array by compring pointer to max
                K = MIN( IPPTR( S ), IPMAX( S ) ) 

C.................  Retrieve emissions from the pollutant array
C                   if the source pollutant ID equals the pollutant ID of the 
C                   pollutant of iteration I.
                IF( IPOSCOD( K ) .EQ. EIIDX( I ) ) THEN

                    IPPTR( S ) = K + 1  ! increase pointer for this source by 1

                    DO J = 1, NPPOL     ! rearrange pollutant-specific info
                        SRCPOL( S,J ) = POLVAL( K,J )
                    END DO

C.....................  If missing emissions or VMT data records are fatal, 
C                       then write message and set error flag
                    IF( SFLAG .AND. SRCPOL( S,1 ) .LT. 0 ) THEN
                        EFLAG = .TRUE.
                        L = LEN_TRIM( EINAM( I ) )
                        CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )
                        MESG = 'ERROR: Missing data for "' // 
     &                         EINAM( I )( 1:L ) // '" for source:' //
     &                         CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                    END IF

                END IF

            END DO  ! end of loop through sources

            CALL WRINVPOL( ENAME, CATEGORY, NSRC, 1, NPPOL, 
     &                     EINAM( I ), SRCPOL )

        END DO  ! end loop through actual output pollutants

        IF( EFLAG ) THEN
            MESG = 'Missing data for some sources is not allowed ' //
     &             CRLF() // BLANK5 //
     &             'because the environment variable ' //
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
 
        END PROGRAM SMKINVEN
