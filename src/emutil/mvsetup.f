
        PROGRAM MVSETUP

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C    Program MVSETUP uses the SMOKE mobile source file to create an condensed
C    list of all sources in the inventory by FIPS code, roadtype, vehicle type,
C    and including the speed from the inventory (if any).
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C    Original by M. Houyoux 5/2001
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                
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
C smoke@ncsc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************
 
C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE
 
C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'
        INCLUDE 'PARMS3.EXT'
        INCLUDE 'IODECL3.EXT'
        INCLUDE 'FDESC3.EXT'
      
C...........   EXTERNAL FUNCTIONS 
        CHARACTER*2   CRLF
        REAL          ENVREAL
        LOGICAL       ENVYN  
        INTEGER       FIND1
        INTEGER       FIND2
        INTEGER       GETFLINE
        INTEGER       JUNIT
        INTEGER       PROMPTFFILE

        EXTERNAL      CRLF, ENVREAL, ENVYN, FIND1, FIND2, JUNIT,
     &                GETFLINE, PROMPTFFILE

C...........   PARAMETERS and their descriptions:
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C.........  Array that contains the names of the inventory variables needed 
C           for this program
        CHARACTER(LEN=IOVLEN3) IVARNAMS( MXINVARR )

C...........   Local allocatable arrays...
C.........  Condensed source list
        INTEGER, ALLOCATABLE :: INDX     ( : )  ! Sorting index
        INTEGER, ALLOCATABLE :: SRCFIP   ( : )  ! FIPS code
        INTEGER, ALLOCATABLE :: SRCROAD  ( : )  ! Roadtype code
        INTEGER, ALLOCATABLE :: SRCVTC   ( : )  ! Vehicle codes
        REAL   , ALLOCATABLE :: SRCSPD   ( : )  ! Speeds
        REAL   , ALLOCATABLE :: SRCREFFIP( : )  ! Ref FIPS for MOBILE inputs
        INTEGER, ALLOCATABLE :: SRCFILIDX( : )  ! Index to MOBILE input file
        INTEGER, ALLOCATABLE :: SRCPSI   ( : )  ! parameter scheme index

C.........  Speeds file
        INTEGER, ALLOCATABLE :: SPDFIP   ( : )  ! FIPS code
        INTEGER, ALLOCATABLE :: SPDROAD  ( : )  ! Roadtype code
        REAL   , ALLOCATABLE :: SPDSPD   ( : )  ! Speeds
      
C.........  Alias file
        INTEGER, ALLOCATABLE :: MVAFIP   ( : )  ! FIPS code
        INTEGER, ALLOCATABLE :: MVAREFFIP( : )  ! Ref FIPS for MOBILE inputs

C.........  MOBILE input file list
        INTEGER, ALLOCATABLE      :: LISTFIP ( : ) ! FIPS code
        CHARACTER*256, ALLOCATABLE:: LISTNAMS( : ) ! Ref FIPS for MOBILE inputs

C.........  FIPS/road combinations
        INTEGER, ALLOCATABLE :: FRDFIP   ( : )  ! FIPS code
        INTEGER, ALLOCATABLE :: FRDROAD  ( : )  ! Roadtype code
        LOGICAL, ALLOCATABLE :: FRDSTAT  ( : )  ! true: 1 PSIs per FIPS/road

C...........   File units and logical/physical names
        INTEGER      :: ADEV = 0 !  alias file unit no.
        INTEGER      :: CDEV = 0 !  condensed sources file unit no.
        INTEGER         LDEV     !  log-device
        INTEGER      :: MDEV = 0 !  MPLIST output file unit no.
        INTEGER      :: NDEV = 0 !  MOBILE input file list unit no.
        INTEGER      :: PDEV = 0 !  pre-MPREF output file unit no.
        INTEGER      :: SDEV = 0 !  speeds file unit no.
        INTEGER      :: TDEV = 0 !  tmp MOBILE input file unit no.

C...........   LOCAL VARIABLES and their descriptions:
        INTEGER         I, J, JJ, K, KK, L, N      ! counters and indices

        INTEGER         CREF            ! tmp reference county ID
        INTEGER         CYID            ! tmp county ID
        INTEGER         IOS             ! i/o status
        INTEGER         IREC            ! record counter
        INTEGER         FIP             ! tmp country/state/county code
        INTEGER         LFIP            ! previous country/state/county code
        INTEGER         LM5P            ! MOBILE5 path name length
        INTEGER         LPSI            ! previous PSI code
        INTEGER         LREFFIP         ! previous reference FIPS code
        INTEGER         LROAD           ! previous roadtype code
        INTEGER         LSTA            ! previous state ID
        INTEGER         LVTC            ! previous vehicle type code
        INTEGER         NCSRC           ! number of condensed sources
        INTEGER         NFRD            ! no. FIPS/road class combos
        INTEGER         NLIST           ! no. recs in MOBILE inputs list file
        INTEGER         NMVA            ! number of recs in mobile aliases file
        INTEGER         NSPD            ! number of speeds entries
        INTEGER         PSI             ! tmp parameter scheme index
        INTEGER         REFFIP          ! tmp reference FIPS code
        INTEGER         ROAD, ROAD2     ! tmp roadtype code
        INTEGER         SREF            ! tmp reference state ID
        INTEGER         STID            ! tmp state ID
        INTEGER         VTC             ! tmp vehicle type code

        REAL            DFSPD           ! default speed from e.v.
        REAL            LSPD            ! previous speed
        REAL            SPD             ! tmp speed

        LOGICAL :: EFLAG = .FALSE.      ! true: error found
        LOGICAL :: MFLAG = .FALSE.      ! true: message should be printed
        LOGICAL :: SFLAG = .FALSE.      ! true: override speeds w/ speeds file

        CHARACTER*300           MESG    ! message field
        CHARACTER*256           FNAME   ! tmp MOBILE file name
        CHARACTER*512           M5PATH  ! tmp MOBILE file path
        CHARACTER*786           INFILE  ! tmp MOBILE path and file

        CHARACTER*16, PARAMETER :: PROGNAME = 'MVSETUP'

C***********************************************************************
C   begin body of program MVSETUP
 
        LDEV = INIT3()
 
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )
 
        MESG = 'Override inventory speeds with separate speeds file'
        SFLAG = ENVYN( 'OVERRIDE_INV_SPEEDS', MESG, .FALSE., IOS )

        MESG = 'Default average speed'
        DFSPD = ENVREAL( 'DEFAULT_SPEED', MESG, 30., IOS )

        MESG = 'Path for MOBILE input files'
        CALL ENVSTR( 'M5PATH', MESG, './', M5PATH, IOS )
        LM5P = LEN_TRIM( M5PATH )

C.........  Limit source category to mobile sources
        CATEGORY = 'MOBILE'
        CATDESC  = 'Mobile'
        CATLEN   = LEN_TRIM( CATEGORY )

C.........  Get unit number to use later for checking MOBILE input files
        TDEV = JUNIT()

C.........  Open input files
        
        CDEV = PROMPTFFILE( 
     &         'Enter logical name for CONDENSED INVENTORY FILE',
     &         .TRUE., .TRUE., 'MVCONDS', PROGNAME )

        SDEV = PROMPTFFILE( 
     &         'Enter logical name for SPEEDS FILE',
     &         .TRUE., .TRUE., 'MVSPEED', PROGNAME )

        ADEV = PROMPTFFILE( 
     &         'Enter logical name for ALIAS FILE',
     &         .TRUE., .TRUE., 'MVALIAS', PROGNAME )

        NDEV = PROMPTFFILE( 
     &         'Enter logical name for MOBILE INPUTS LIST FILE',
     &         .TRUE., .TRUE., 'MVINLIST', PROGNAME )

C.........  Get size of files
        NCSRC = GETFLINE( CDEV, 'Condensed inventory file' )
        NSPD  = GETFLINE( SDEV, 'Speeds file' )
        NMVA  = GETFLINE( ADEV, 'Alias file' )
        NLIST = GETFLINE( NDEV, 'MOBILE inputs list' )

C.........  Allocate memory for input files
        
        ALLOCATE( INDX( NCSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX', PROGNAME )
        ALLOCATE( SRCFIP( NCSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCFIP', PROGNAME )
        ALLOCATE( SRCROAD( NCSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCROAD', PROGNAME )
        ALLOCATE( SRCVTC( NCSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCVTC', PROGNAME )
        ALLOCATE( SRCSPD( NCSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCSPD', PROGNAME )
        ALLOCATE( SRCREFFIP( NCSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCREFFIP', PROGNAME )
        ALLOCATE( SRCFILIDX( NCSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCFILIDX', PROGNAME )
        ALLOCATE( SRCPSI( NCSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCPSI', PROGNAME )
        INDX      = 0   ! array
        SRCFIP    = 0   ! array
        SRCROAD   = 0   ! array
        SRCVTC    = 0   ! array
        SRCSPD    = 0.  ! array
        SRCREFFIP = 0.  ! array
        SRCFILIDX = 0   ! array
        SRCPSI    = 0   ! array

        ALLOCATE( SPDFIP( NSPD ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPDFIP', PROGNAME )
        ALLOCATE( SPDROAD( NSPD ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPDROAD', PROGNAME )
        ALLOCATE( SPDSPD( NSPD ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPDSPD', PROGNAME )
        SPDFIP  = 0   ! array
        SPDROAD = 0   ! array
        SPDSPD  = 0.  ! array

        ALLOCATE( MVAFIP( NMVA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MVAFIP', PROGNAME )
        ALLOCATE( MVAREFFIP( NMVA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MVAREFFIP', PROGNAME )
        MVAFIP    = 0   ! array
        MVAREFFIP = 0   ! array

        ALLOCATE( LISTFIP( NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LISTFIP', PROGNAME )
        ALLOCATE( LISTNAMS( NLIST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LISTNAMS', PROGNAME )
        LISTFIP  = 0     ! array
        LISTNAMS = ' '   ! array

C.........  Read input files
C.........  Assume list formatting and sorted files for this version.

        CALL M3MSG2( 'Reading condensed sources file...' )

C.........  Read condensed source list
        IREC  = 0
        I     = 0
        NFRD  = 0
        LFIP  = -9
        LROAD = -9
        LVTC  = -9
        DO N = 1, NCSRC

            READ( CDEV, *, END=9999, IOSTAT=IOS ) FIP, ROAD, VTC, SPD
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading condended sources file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( FIP .LT. LFIP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: FIPS code not in sorted ' //
     &                 'order at line', IREC
                CALL M3MESG( MESG )
            END IF

            IF( FIP  .EQ. LFIP .AND. 
     &          ROAD .LT. LROAD      ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: road code not in sorted ' //
     &                 'order at line', IREC
                CALL M3MESG( MESG )
            END IF

            IF( FIP  .EQ. LFIP  .AND. 
     &          ROAD .EQ. LROAD .AND. 
     &          VTC .LE. LVTC        ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: vehicle code not in ' //
     &                 'sorted order at line', IREC, CRLF()// BLANK10 //
     &                 'or records are not unique.'
                CALL M3MESG( MESG )
            END IF

            IF( .NOT. EFLAG ) THEN
                I = I + 1
                SRCFIP ( I ) = FIP
                SRCROAD( I ) = ROAD
                SRCVTC ( I ) = VTC
                SRCSPD ( I ) = SPD

C.................  Count the number of unique FIPS/road combinations
                IF ( FIP .NE. LFIP .OR. ROAD .NE. LROAD ) 
     &               NFRD = NFRD + 1

            END IF

            LFIP  = FIP
            LROAD = ROAD
            LVTC  = VTC

        END DO
        NCSRC = I

C.........  Read speeds file
        CALL M3MSG2( 'Reading speeds file...' )

        IREC = 0
        I    = 0
        LFIP  = -9
        LROAD = -9
        DO N = 1, NSPD

            READ( SDEV, *, END=9999, IOSTAT=IOS ) FIP, ROAD, SPD
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading speeds file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( FIP .LT. LFIP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: FIPS code not in sorted ' //
     &                 'order at line', IREC
                CALL M3MESG( MESG )
            END IF

            IF( FIP  .EQ. LFIP .AND. 
     &          ROAD .LT. LROAD      ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: road code not in sorted ' //
     &                 'order at line', IREC
                CALL M3MESG( MESG )
            END IF

            IF( .NOT. EFLAG ) THEN
                I = I + 1
                SPDFIP ( I ) = FIP
                SPDROAD( I ) = ROAD
                SPDSPD ( I ) = SPD

            END IF

            LFIP  = FIP
            LROAD = ROAD

        END DO
        NSPD = I

        CALL M3MSG2( 'Reading alias file...' )

C.........  Read mvalias file
        IREC = 0
        I    = 0
        LFIP  = -9
        DO N = 1, NMVA

            READ( ADEV, *, END=9999, IOSTAT=IOS ) STID, CYID, SREF, CREF
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading alias file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            FIP    = STID * 1000 + CYID
            REFFIP = SREF * 1000 + CREF

            IF( FIP .LT. LFIP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: FIPS code not in sorted ' //
     &                 'order at line', IREC
                CALL M3MESG( MESG )
            END IF

            IF( .NOT. EFLAG ) THEN
                I = I + 1
                MVAFIP   ( I ) = FIP
                MVAREFFIP( I ) = REFFIP

            END IF

            LFIP  = FIP

        END DO
        NMVA = I

        CALL M3MSG2( 'Reading MOBILE input list file...' )

C.........  Read file list of MOBILE inputs
        IREC = 0
        I    = 0
        LFIP = -9
        DO N = 1, NLIST

            READ( NDEV, *, END=9999, IOSTAT=IOS ) FIP, FNAME
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading MOBILE input list file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( FIP .LT. LFIP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: FIPS code not in sorted ' //
     &                 'order at line', IREC
                CALL M3MESG( MESG )
            END IF

            IF( .NOT. EFLAG ) THEN
                I = I + 1
                LISTFIP ( I ) = FIP
                LISTNAMS( I ) = FNAME

            END IF

            LFIP  = FIP

        END DO
        NLIST = I

C.........  Create list of FIPS/road combinations for later use
        ALLOCATE( FRDFIP( NFRD ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FRDFIP', PROGNAME )
        ALLOCATE( FRDROAD( NFRD ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FRDROAD', PROGNAME )
        ALLOCATE( FRDSTAT( NFRD ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FRDSTAT', PROGNAME )
        FRDFIP  = 0        ! array
        FRDROAD = 0        ! array
        FRDSTAT = .TRUE.   ! array

        LFIP = -9
        LROAD = -9
        I = 0
        DO N = 1, NCSRC

            FIP  = SRCFIP ( N )
            ROAD = SRCROAD( N )

            IF ( FIP .NE. LFIP .OR. ROAD .NE. LROAD ) THEN
                I = I + 1
                FRDFIP ( I ) = FIP
                FRDROAD( I ) = ROAD
            END IF

            LFIP = FIP
            LROAD = ROAD

        END DO

        CALL M3MSG2( 'Checking inputs for states...' )

C.........  Check to ensure sources have enough information in other inputs
        LSTA = -9
        LFIP = -9
        DO N = 1, NCSRC

            FIP  = SRCFIP ( N )
            STID = FIP / 1000
            ROAD = SRCROAD( N )
            VTC  = SRCVTC ( N )
            MFLAG = .FALSE.

C.............  Write out state code to screen and log file
            IF ( STID .NE. LSTA ) THEN
                WRITE( MESG,94010 ) BLANK10, STID
                CALL M3MSG2( MESG )
            END IF

C.............  Set sorting index
            INDX( N ) = N

C.............  Check that all sources have valid speeds, and replace speeds 
C               if requested by user...

C.............  If default is to try to use speeds file
            IF ( SFLAG ) THEN

C.................  Try to find source in speeds file...
                J = FIND2( FIP, ROAD, NSPD, SPDFIP, SPDROAD )

C.................  If speeds not found, check original value and write 
C                   message
                IF ( J .LE. 0 ) THEN

                    IF( SRCSPD( N ) .LE. 0. ) THEN
                        MFLAG = .TRUE.

                    ELSE
                        WRITE( MESG,94010 ) 'WARNING: Using '//
     &                     'fallback speed from inventory because' //
     &                     CRLF() // BLANK10 // 
     &                     'speeds file not available for:' //
     &                     CRLF() // BLANK10 // 'FIPS=', FIP, 
     &                     'Roadtype=', ROAD, 'Vehicle code=', VTC
                        CALL M3MESG( MESG )

                    END IF
                     
C.................  Otherwise, use speed from speeds file
                ELSE
                    SRCSPD( N ) = SPDSPD( J )

                END IF

C.............  If default is to try to use inventory speeds, then check
C               for valid speed value...
            ELSE IF( SRCSPD( N ) .LE. 0. ) THEN

C.................  Try to find speed in speeds file
                J = FIND2( FIP, ROAD, NSPD, SPDFIP, SPDROAD )

C.................  If speeds not found, write error
                IF ( J .LE. 0 ) THEN
                    MFLAG = .TRUE.

C.................  Otherwise, use speed from speeds file
                ELSE
                    WRITE( MESG,94010 ) 'WARNING: Speed '//
     &                 'from speeds file used for:' //
     &                 CRLF() // BLANK10 // 'FIPS=', FIP, 
     &                 'Roadtype=', ROAD, 'Vehicle code=', VTC
                    CALL M3MESG( MESG )
                    
                    SRCSPD( N ) = SPDSPD( J )

                END IF

            END IF

C.............  If message flag set on this iteration, write message
            IF ( MFLAG ) THEN
                WRITE( MESG,94013 ) 'WARNING: No valid speed '//
     &                 'available in source or speeds file for:' //
     &                 CRLF() // BLANK10 // 'FIPS=', FIP, 
     &                 'Roadtype=', ROAD, 'Vehicle code=', VTC, 
     &                 CRLF() // BLANK10 // 'Using default ' //
     &                 'speed of ', DFSPD
                CALL M3MESG( MESG )
                SRCSPD( N ) = DFSPD
            END IF

            IF ( FIP .NE. LFIP ) THEN

C.................  Check that all sources have valid entries in the alias 
C                   file, and that the MOBILE files being referenced are 
C                   available
                JJ = FIND1( FIP, NMVA, MVAFIP )

C.................  If there isn't an entry, try to fallback to county 001
                IF( JJ .LE. 0 ) THEN

                    JJ = FIND1( STID*1000 + 1, NMVA, MVAFIP )

C.....................  If still no luck, then write error
                    IF( JJ .LE. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: No valid alias '//
     &                      'available for:' //
     &                  CRLF() // BLANK10 // 'FIPS=', FIP
                        CALL M3MESG( MESG )

                    ELSE

                        WRITE( MESG,94010 ) 'WARNING: No valid alias '//
     &                      'available, so assuming ', MVAREFFIP( JJ ),
     &                      'for:' //
     &                      CRLF() // BLANK10 // 'FIPS=', FIP
                        CALL M3MESG( MESG )

                        KK = FIND1( STID*1000 + 1, NLIST, LISTFIP )

                    END IF

C.................  Otherwise, check to see that reference FIPS code has a
C                   MOBILE file associated with it
                ELSE
                    KK = FIND1( MVAREFFIP( JJ ), NLIST, LISTFIP )

                END IF

C.................  Write error if FIPS code is not available in input list
                IF( KK .LE. 0 ) THEN

                    WRITE( MESG,94010 ) 'ERROR: No alias code ' //
     &                 'available in MOBILE file list for:' //
     &                 CRLF() // BLANK10 // 'FIPS=', FIP
                    CALL M3MESG( MESG )

C.................  Otherwise, check to see if the file exists
                ELSE

                    INFILE = M5PATH( 1:LM5P ) // '/' // LISTNAMS( KK )
                    OPEN( TDEV, ERR=1001, FILE=INFILE, STATUS='OLD' )

C.....................  If file exists, close it 
                    CLOSE( TDEV )

                END IF
                
            END IF

            LFIP = FIP
            LSTA = STID

            SRCREFFIP( N ) = REAL( MVAREFFIP( JJ ) )
            SRCFILIDX( N ) = KK
            CYCLE

1001        L = LEN_TRIM( FNAME )
            MESG = 'ERROR: Could not find file:'// CRLF() //
     &                     BLANK10 // '"' // FNAME( 1:L ) // '"'
            CLOSE( TDEV )
            LFIP = FIP
            LSTA = STID
                
        END DO

        IF ( EFLAG ) THEN
            MESG = 'Problem with input files (see above).'
            CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Sort sources by speed and MOBILE reference
        CALL SORTR2( NCSRC, INDX, SRCREFFIP, SRCSPD )

C.........  Assign PSIs to sources
        PSI = 0
        LREFFIP = -9
        LSPD    = -9.
        DO N = 1, NCSRC

            J = INDX( N )
            REFFIP = INT( SRCREFFIP( J ) )
            SPD    = SRCSPD   ( J )

            IF ( REFFIP .NE. LREFFIP .OR.
     &           SPD    .NE. LSPD          ) THEN
    
                PSI = PSI + 1
                LREFFIP = REFFIP
                LSPD = SPD

            END IF

            SRCPSI( J ) = PSI ! set PSI for source 

        END DO

C.........  Flag FIPS/road combos if they do not have different PSIs across
C           multiple records with same FIPS/road
        LFIP = -9
        LROAD = -9
        LPSI = -9
        DO N = 1, NCSRC
            FIP    = SRCFIP ( N )
            ROAD   = SRCROAD( N )
            PSI    = SRCPSI ( N )

            IF( FIP .EQ. LFIP .AND. ROAD .EQ. LROAD ) THEN

                FRDSTAT( J ) = ( FRDSTAT( J ) .AND. PSI .EQ. LPSI )

            ELSE
                J = FIND2( FIP, ROAD, NFRD, FRDFIP, FRDROAD )

            END IF

            LFIP  = FIP
            LROAD = ROAD
            LPSI  = PSI

        END DO        

C.........  Open output files
        MDEV = PROMPTFFILE( 
     &         'Enter logical name for MPLIST OUTPUT FILE',
     &         .FALSE., .TRUE., 'MPLIST', PROGNAME )

        PDEV = PROMPTFFILE( 
     &         'Enter logical name for PRE-MPREF OUTPUT FILE',
     &         .FALSE., .TRUE., 'PREMPREF', PROGNAME )

C NOTE: This needs to be updated to write out fewer MPLIST records (when all
c    n: vehicle types in a county/road have the same PSI, don't write out).
C.........  Write out MPLIST and pre-MPREF file
        LPSI = -9
        DO N = 1, NCSRC

C.............  MPLIST file...
            FIP    = SRCFIP ( N )
            ROAD   = SRCROAD( N )
            VTC    = SRCVTC ( N )
            PSI    = SRCPSI ( N )

C.............  Find FIPS/road in list
            K = FIND2( FIP, ROAD, NFRD, FRDFIP, FRDROAD )

C.............  If same PSI for all entries with same FIP/road, then write
C               out only once, otherwise, write out all
            IF( FRDSTAT( K ) ) THEN

                IF( FIP .NE. LFIP .OR. ROAD .NE. LROAD ) THEN
                    WRITE( MDEV, 93888 ) 
     &                     FIP, ROAD, 0, 0, 'VMT 24*', PSI
                END IF

            ELSE 

                WRITE( MDEV, 93888 ) FIP, ROAD, 0, VTC, 'VMT 24*', PSI

            END IF

            LFIP  = FIP
            LROAD = ROAD

C.............  pre-MPREF file...
            J = INDX( N )
            ROAD2  = SRCROAD  ( J )
            REFFIP = INT( SRCREFFIP( J ) )
            SPD    = SRCSPD   ( J )
            PSI    = SRCPSI   ( J )

            IF( PSI .NE. LPSI ) THEN
                WRITE( PDEV, 93890 ) REFFIP, SPD, PSI, ROAD2
            END IF

            LPSI = PSI

        END DO

C.........  Normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

9999    WRITE( MESG, 94010) 'End of file reached unexpectedly at line',
     &         IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93888   FORMAT( I6.6, 1X, I8.2, 1X, I1, 1X, I8.4, 1X, A, I4.4 )

93890   FORMAT( I5.5, 1X, F5.1, 1X, I4.4, 1X, I8.2 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )
 
94013   FORMAT( 3( A, :, I7, :, 1X ), A, 1X, F5.1 )
 
	END PROGRAM MVSETUP
