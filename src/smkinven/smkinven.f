
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

C.........  Full list of inventory pollutants/activities (in output order)
        INTEGER                  :: MXIDAT = 0   ! Max no of inv pols & acvtys
        INTEGER     , ALLOCATABLE:: INVDCOD( : ) ! 5-digit pollutant/actvty code
        INTEGER     , ALLOCATABLE:: INVSTAT( : ) ! Status (<0 activity; >0 pol)
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE:: INVDNAM( : ) ! Name of pollutant

C.........  Index to master list for pollutant and activity names
        INTEGER, ALLOCATABLE:: AVIDX( : ) ! for activity names
        INTEGER, ALLOCATABLE:: EIIDX( : ) ! for pol names

C.........  Inventory temporay arrays
        INTEGER, ALLOCATABLE:: IPPTR ( : ) ! position in POLVAL sparse array
        INTEGER, ALLOCATABLE:: IPMAX ( : ) ! max IPPTR by source
        INTEGER, ALLOCATABLE:: IDXSCC( : ) ! sorting index for output SCC file
        REAL   , ALLOCATABLE:: SRCPOL( :,: )!  pollutant-spec values by source

C.........  File units and logical/physical names

        INTEGER         ADEV    !  unit no. for output actual SCCs
        INTEGER         IDEV    !  inventory file (various formats)
        INTEGER         LDEV    !  log-device
        INTEGER         PDEV    !  unit number for pollutants codes/names file
        INTEGER         RDEV    !  unit no. for def stack pars or mobile codes
        INTEGER         TDEV    !  unit no. for output temporal profile no.s
        INTEGER         SDEV    !  unit no. for ASCII output inventory file
        INTEGER         VDEV    !  unit no. for activity codes/names file
        INTEGER         ZDEV    !  unit no. for time zone file

        CHARACTER(LEN=NAMLEN3) ANAME !  inventory ASCII output logical name
        CHARACTER(LEN=NAMLEN3) ENAME !  inventory I/O API output logical name
        CHARACTER(LEN=NAMLEN3) INAME !  inventory input logical name

C...........   Other local variables
                                
        INTEGER         S, I, J, J1, J2, K, L, L2, V !  counters and indices

        INTEGER         FILFMT  !  input file(s) format code
        INTEGER         FIP     !  Temporary FIPS code
        INTEGER         IOS     !  I/O status
        INTEGER         MAXK    !  test for maximum value of K in output loop
        INTEGER      :: NDAT = 0!  tmp no. actual pols & activities
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
        CALL OPENINVIN( CATEGORY, IDEV, RDEV, PDEV, VDEV, ZDEV, INAME )

        CALL M3MSG2( 'Setting up to read inventory data...' )

C.........  Get no. lines in pollutant codes & activities files for allocating
C           memory     
        IF( PDEV .GT. 0 ) THEN
            MXIDAT = GETFLINE( PDEV, 'Pollutant codes and names file' )
        END IF
        IF( VDEV .GT. 0 ) THEN
            I = GETFLINE( VDEV, 'Activity names file' )
            MXIDAT = MXIDAT + I
        END IF
        
C.........  Allocate memory for storing contents of pollutants & activities
C           files
        ALLOCATE( INVDCOD( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDCOD', PROGNAME )
        ALLOCATE( INVDNAM( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVDNAM', PROGNAME )
        ALLOCATE( INVSTAT( MXIDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVSTAT', PROGNAME )

C.........  Read country, state, and county file for time zones
        CALL RDSTCY( ZDEV, 1, I )   !  "I" used as a dummy

C.........  Initialize inventory data status.  PROCINVEN will expect this
C           type of initialization.
        INVSTAT = 1   ! array

C.........  Read, sort, and store pollutant codes/names file
        IF( PDEV .GT. 0 ) THEN
            CALL RDCODNAM( PDEV, MXIDAT, NDAT, INVDCOD, 
     &                     INVDNAM )
        END IF

C.........  Read, sort, and store activity codes/names file
        IF( VDEV .GT. 0 ) THEN

            INVSTAT( NDAT+1:MXIDAT ) = -1

            CALL RDCODNAM( VDEV, MXIDAT, NDAT, INVDCOD, 
     &                     INVDNAM )
        END IF

        MXIDAT = NDAT

C.........  Fill tables for translating mobile road classes and vehicle types
C.........  The tables are passed through MODINFO
        IF( CATEGORY .EQ. 'MOBILE' ) THEN

            CALL RDMVINFO( RDEV )

        END IF

        CALL M3MSG2( 'Reading raw inventory data...' )

C.........  Read the raw inventory data, and store in unsorted order
C.........  The arrays that are populated by this subroutine call
C           are contained in the module MODSOURC

        CALL RDINVEN( IDEV, INAME, MXIDAT, INVDCOD, INVDNAM, FILFMT,
     &                NRAWBP, PRATIO, TFLAG )

        CALL M3MSG2( 'Sorting raw inventory data...' )

C.........  Sort inventory and pollutants (sources x pollutants). Note that
C           sources are sorted based on character string definition of the 
C           source so that source definition can be consistent with that of
C           the input format.

        CALL SORTIC( NRAWBP, INDEXA, CSOURCA )

        CALL M3MSG2( 'Processing inventory data...' )

C.........  Processing inventory records and store in sorted order

        CALL PROCINVEN( NRAWBP, MXIDAT, FILFMT, PRATIO, INVSTAT )

C.........  Determine memory needed for actual pollutants list and actual
C           activities list and allocate them. Invstat has been updated
C           to be +/- 2 depending on whether the pollutant or activity was
C           present in the inventory.
        NIPOL = 0
        NIACT = 0
        DO I = 1, MXIDAT
            IF( INVSTAT( I ) .GT.  1 ) NIPOL = NIPOL + 1
            IF( INVSTAT( I ) .LT. -1 ) NIACT = NIACT + 1
        ENDDO

        ALLOCATE( EIIDX( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EIIDX', PROGNAME )
        ALLOCATE( EINAM( NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EINAM', PROGNAME )
        ALLOCATE( AVIDX( NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVIDX', PROGNAME )
        ALLOCATE( ACTVTY( NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ACTVTY', PROGNAME )

C.........  Create list of actual pollutants and activities and indexes to the 
C           master list. The order in EINAM and ACTVTY will be the output 
C           order. The indexes are for accessing INVDCOD, if needed.
        J1 = 0
        J2 = 0
        DO I = 1, MXIDAT

           IF( INVSTAT( I ) .GT. 0 ) THEN
               J1 = J1 + 1
               EIIDX( J1 ) = I
               EINAM( J1 ) = INVDNAM( I )
           ENDIF

           IF( INVSTAT( I ) .LT. 0 ) THEN
               J2 = J2 + 1
               AVIDX ( J2 ) = I
               ACTVTY( J2 ) = INVDNAM( I )
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
            TZONE = GETTZONE( FIP )

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
     &            .FALSE., .TRUE., CRL // 'TREF_ALT', PROGNAME )

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

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram is for writing out the inventory
C               data, whether it is the pollutant data or the activity data.
C.............  Most variables are defined through host association
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
                SRCPOL = 0. ! array

C.................  Transfer emissions or activity data to output SRCPOL array
                K = 0
                DO S = 1, NSRC

C.....................  Set position in sparse array by compring pointer to max
                    K = MIN( IPPTR( S ), IPMAX( S ) ) 

C.....................  Retrieve emissions from the pollutant array if the 
C                       source pollutant ID equals the pollutant ID of the 
C                       pollutant of iteration I.
                    IF( IPOSCOD( K ) .EQ. INDX( I ) ) THEN

                        IPPTR( S ) = K + 1  ! pointer fo source S

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

                CALL WRINVPOL( ENAME, CATEGORY, NSRC, 1, NPPOL, 
     &                         NAMES( I ), SRCPOL )

            END DO  ! end loop through actual output data

C------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

            END SUBROUTINE LOOP_FOR_OUTPUT
 
        END PROGRAM SMKINVEN
