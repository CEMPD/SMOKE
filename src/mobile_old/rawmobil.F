
         PROGRAM RAWMOBIL

C***********************************************************************
C  program body starts at line 164
C
C  DESCRIPTION:
C       This program reads the county and link VMT data and creates
C       the SMOKE mobile inventory file MOBL.
C
C  PRECONDITIONS REQUIRED:
C       Input files exist and correct format
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C	PROMPTFFILE, PROMPTMFILE, STR2INT, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  12/95 by MRH
C
C***********************************************************************
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
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
       
        INCLUDE 'MBDIMS3.EXT'
        INCLUDE 'TMDIMS3.EXT'
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         FIND1, FIND2, FIND3
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         JUNIT
        INTEGER         LEN2
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        INTEGER         TRIMLEN
        
        EXTERNAL  CRLF, ENVYN, FIND1, FIND2, FIND3, GETYN, INDEX1, 
     &            JUNIT, LEN2, PROMPTFFILE, PROMPTMFILE, STR2INT, 
     &            TRIMLEN

C...........   LOCAL PARAMETERS and their descriptions:

        CHARACTER*5     BLANK5
        PARAMETER(      BLANK5 = ' ' )

C...........   LOCAL VARIABLES and their descriptions:

C.........  Sorted Source variables
        INTEGER     IFIP   ( NMSRC )  ! Source FIP
        INTEGER     IRCLAS ( NMSRC )  ! Source Road Class
        INTEGER     TZONES ( NMSRC )  ! Source time zone number
        INTEGER     TPFLAG ( NMSRC )  ! Temporal profile type flag
        INTEGER     INVYR  ( NMSRC )  ! Source inventory year
        INTEGER     ILINK  ( NMSRC )  ! Source Link number
        REAL        XLOC1  ( NMSRC )  ! UTM X coord of link start [m]   
        REAL        XLOC2  ( NMSRC )  ! UTM X coord of link end [m]
        REAL        YLOC1  ( NMSRC )  ! UTM Y coord of link start [m]
        REAL        YLOC2  ( NMSRC )  ! UTM Y coord of link end [m]
        REAL        VMT    ( NMSRC )  ! Source VMT [mi/typical day]

        COMMON  / EMISREC / IFIP, IRCLAS, TZONES, TPFLAG, INVYR, ILINK,
     &                      XLOC1, YLOC1, XLOC2, YLOC2, VMT

C.........  Unsorted non-link and link source variables
        INTEGER     INDEXA ( MXNSRC + MXLSRC )  ! Non-link and link sorting index
        INTEGER     IFIPA  ( MXNSRC + MXLSRC )
        INTEGER     IRCLSA ( MXNSRC + MXLSRC )
        INTEGER     TPFLGA ( MXNSRC + MXLSRC )
        INTEGER     INVYRA ( MXNSRC + MXLSRC )
        INTEGER     ILINKA ( MXNSRC + MXLSRC )
        REAL        XLOC1A ( MXNSRC + MXLSRC )
        REAL        XLOC2A ( MXNSRC + MXLSRC )
        REAL        YLOC1A ( MXNSRC + MXLSRC )
        REAL        YLOC2A ( MXNSRC + MXLSRC )
        REAL        VMTA   ( MXNSRC + MXLSRC )

C.........  File units and logical names

        INTEGER       ADEV              ! Unit for non-link data
        INTEGER       KDEV              ! unit for link data
        INTEGER       IDEV              ! Link source raw datafile list
        INTEGER       LDEV              ! IO/API initialization unit
        INTEGER       MDEV              ! Non-link source raw datafile list
        INTEGER       ZDEV              ! Zones file unit
        CHARACTER*16  ENAME             ! Mobile source output file name

C.........  Other Local variables

        INTEGER     IREC, I, J, N, S  ! Counters and pointers
        INTEGER     L1, L2                    ! Counters and pointers

        INTEGER     FIP               ! Temporary read in FIP
        INTEGER     IOS               ! Error status on read
        INTEGER     INY               ! Temporary inventory year
        INTEGER     LFIP              ! Previous FIP in loop
        INTEGER     LINY              ! Previous inventory year
        INTEGER     LLNK              ! Previous link ID in loop
        INTEGER     LNK               ! Temporary link number
        INTEGER     LRCL              ! Previous road class in loop
        INTEGER     NLFIP             ! Fip for non-link source
        INTEGER     NLRCL             ! Road class for non-link source
        INTEGER     MLNUM             ! Actual number link sources
        INTEGER     MNNUM             ! Actual number non-link sources
        INTEGER     MTOT              ! Preliminary number of total sources
        INTEGER     NSRC              ! Actual number of total sources
        INTEGER     RCL               ! Temporary road class
        INTEGER     ROAD              ! Temporary read-in road code
        INTEGER     TZONE             ! Input TZONE for output data file
        INTEGER     UZONE             ! UTM zone to put into file header
        INTEGER     ZONE              ! Temporary UTM zone for link coords

        REAL        TVMT              ! Temporary read-in VMT
        REAL        X1, X2, Y1, Y2    ! Temporary link start and end coord [m]
        REAL        XX, YY            ! Temporary LAT/LON coords

        LOGICAL     DFLAG             ! True if link FIP/RDT can = non-link
        LOGICAL     EFLAG             ! Error flag
        LOGICAL     LFLAG             ! True if link sources exist
        LOGICAL     NFLAG             ! True if non-link sources exist
        LOGICAL     OVERFLAG          ! True if okay to replace non-link w/link 

        DATA    EFLAG  / .FALSE. /

        CHARACTER*2   RR              ! Temporary area/facility type
        CHARACTER*256 FNAME           ! Temporary logical raw datafile name
        CHARACTER*256 LINE            ! Line-read buffer
        CHARACTER*256 MESG
       
C***********************************************************************
C   begin body of program RAWMOBIL
        
        LDEV = INIT3()

        CALL INITEM( LDEV )
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program RAWMOBIL ,',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  'Optional checking that all species are reported for each ',
     &  'source may be turned on via "setenv RAW_SRC_CHECK Y".',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'RAWMOBIL', 0, 0, 'Ending program.', 2 )
        END IF

C.........  Get environment variables
        DFLAG = ENVYN( 'ADD_LINK_NONLINK',
     &                 'RAWMOBIL permits link/non-link in same FIP/RDT',
     &                 .FALSE., IOS )

C.........  Get file name; open input time zone file
        ZDEV = PROMPTFFILE( 
     &           'Enter logical name for TIME ZONE file',
     &           .TRUE., .TRUE., 'ZONES', 'RAWMOBIL' )

C.........  Assign the UTM zone to put into file header
C.........  (This should go away after update to store all coordinates as LAT-LON)

        UZONE = 17

C.......   Get file name for opening input raw mobile link-source files
        IDEV = PROMPTFFILE( 
     &        'Enter the name of the RAW LINK FILENAMES list or "NONE"',
     &         .TRUE., .TRUE., 'MLKLST', 'RAWMOBIL' )

C.............  End list if NONE entered
        IF( IDEV .LE. 0 ) GO TO 51

        CALL M3MSG2( 'Reading RAW MOBIL LINK files...' )

C.........  Loop over link input files

        INY    = IMISS3
        LFLAG  = .FALSE.
        S      = 0
   22   CONTINUE    ! head of MLKLST read loop

C.............  Read file names
C.............  Exit to line 51 if read is EOF
            READ( IDEV, 93000, END=51 ) LINE

            L1 = INDEX( LINE, 'INVYEAR' )
            L2 = TRIMLEN( LINE )

            IF( L1 .GT. 0 ) THEN
 
                INY = STR2INT( LINE( L1+7:L2 ) )
 
                IF( INY .LE. 0 ) THEN
 
                    CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &               'Must set inventory year using INVYEAR packet ' //
     &               'in MLKLST file', 2 )
 
                ELSEIF( INY .LT. 1950 .OR. INY .GT. 2030 ) THEN
 
                    CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &               'INVYEAR packet set 4-digit year outside range ' //
     &               '1950 to 2030 in MLKLST file', 2 )

                ENDIF

                GO TO 22                        ! To head of MLKLST read loop

            ELSE

                FNAME = LINE( 1:L2 )   ! Set inventory file name

            ENDIF

C.............  Get file name; open input link sources file
            KDEV = JUNIT()
            OPEN( KDEV, ERR=1006, FILE=FNAME, STATUS='OLD' )
 
            WRITE( MESG,94010 )
     &             'Successful OPEN using year', INY,
     &             'for inventory file:' // CRLF() // BLANK5 //
     &             FNAME( 1:TRIMLEN( FNAME ) )
            CALL M3MSG2( MESG )

C.............  If have gotten this far, then we have link data for sure
            LFLAG = .TRUE.
             
C.............  Read file for sources with link specific VMT (i.e. LINK sources)

            IREC = 0

   44       CONTINUE  ! Head of loop for reading link specific VMT data

                READ( KDEV, *, END=50, IOSTAT=IOS ) 
     &                FIP, ROAD, LNK, X1, Y1, X2, Y2, ZONE, TVMT

                IREC = IREC + 1

                IF( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 'Error ', IOS, 
     &                   'reading ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                   ' at line', IREC
                    CALL M3EXIT( 'RAWMOBIL', 0, 0, MESG, 2 )

C.................  Determine IRCLAS based on Urban/Rural and Facility type
C.................  and check for bad area/facility combination
                ELSEIF ( ROAD .LT. 100 ) THEN
     
                    WRITE( RR, '( I2.2 )' ) ROAD 
                    J = INDEX1( RR, NRCLAS, MARFAC3 )

                    IF( J .GT. 0 ) THEN
                        RCL = MROADS3( J )

                    ELSE      ! Error, road class not found

                        WRITE( MESG, 94010 )
     &                     'Skipping record', IREC,
     &                     ': Bad area/facility combination ' // RR
                        CALL M3MESG( MESG )

                        GO TO 44  ! To head of current KDEV read loop

                    ENDIF

C.................  Check for bad road class number
                ELSE

                    J = FIND1( ROAD, NRCLAS, MROADS3 )

                    IF( J .LE. 0 ) THEN
                        WRITE( MESG, 94010 )
     &                         'Skipping record', IREC,
     &                         ': Bad road class', ROAD
                        CALL M3MESG( MESG )

                        GO TO 44  ! To head of current KDEV read loop
                    ELSE
                        RCL = ROAD

                    ENDIF
                ENDIF

C.................  If link zone is not file zone, convert coords to zone
                IF( X1 .EQ. X2 .AND. Y1 .EQ. Y2 ) THEN
               
                    WRITE( MESG, 94010 ) 'Bad coordinates for Link ', 
     &                    LNK, 'record', IREC, '. Adjusting...'
                    CALL M3MESG( MESG )
                    X2 = X1 + 1.
                    Y2 = Y1 + 1.

                ELSEIF( ZONE .NE. UZONE ) THEN

                    CALL UTM2LL( X1, Y1,  ZONE, XX, YY )  ! X1, Y1 to LAT/LON
                    CALL LL2UTM( XX, YY, UZONE, X1, Y1 )  ! X1, Y1 to UZONE
                    CALL UTM2LL( X2, Y2,  ZONE, XX, YY )  ! X2, Y2 to LAT/LON
                    CALL LL2UTM( XX, YY, UZONE, X2, Y2 )  ! X2, Y2 to UZONE

                ENDIF

                S = S + 1

C.................  Ensure no overflow
                IF( S .GT. MXLSRC ) THEN  ! Error

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                     'Maximum number of LINK SOURCES', MXLSRC,
     &                     'exceeded at line', IREC, ', file ' // 
     &                     FNAME( 1:TRIMLEN( FNAME ) )
                    CALL M3MESG( MESG )

                ELSE

                    INDEXA( S ) = S       ! Index for sorting

                    IFIPA ( S ) = FIP
                    IRCLSA( S ) = RCL
                    ILINKA( S ) = LNK
                    TPFLGA( S ) = MTPRFAC*WTPRFAC
                    INVYRA( S ) = INY
                    XLOC1A( S ) = X1
                    XLOC2A( S ) = X2
                    YLOC1A( S ) = Y1
                    YLOC2A( S ) = Y2
                    VMTA  ( S ) = TVMT

                ENDIF

            GO TO 44  ! To head of current KDEV read loop

   50       CLOSE( KDEV )

        GO TO 22      ! To head of loop on LNKRAW files

   51   CONTINUE

        MLNUM = S

        IF( EFLAG ) THEN   ! Internal dimensions exceeded

            WRITE( MESG, 94010 ) 
     &             'Records count for LINK sources       :', MLNUM ,
     &             CRLF() // BLANK5 //
     &             'Max dimensioned record-count (MXLSRC):', MXLSRC,
     &             CRLF() // BLANK5 //
     &             'Max record-count exceeded.'

            CALL M3MSG2( MESG )
            CALL M3MESG( 'Max record-count exceeded in RAW VMT files')

        ENDIF

C.......   Get file name for opening input raw mobile non-link-source
C.......   (i.e., county-based VMT) files

        MDEV = PROMPTFFILE(
     &        'Enter the name of the RAW NON-LINK FILENAMES list ' //
     &        'or "NONE"',
     &         .TRUE., .TRUE., 'MNLST', 'RAWMOBIL' )

C.........  Loop over non-link input files

        INY   = IMISS3
        NFLAG = .FALSE.
        S     = MLNUM
        N     = 0
   55   CONTINUE

C.............  Read file names
C.............  Exit to line 111 if read is EOF
            READ( MDEV, 93000, END=111 ) LINE
 
            L1 = INDEX( LINE, 'INVYEAR' )
            L2 = TRIMLEN( LINE )
 
            IF( L1 .GT. 0 ) THEN
 
                INY = STR2INT( LINE( L1+7:L2 ) )

                IF( INY .LE. 0 ) THEN
 
                    CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &               'Must set inventory year using INVYEAR packet ' //
     &               'in MNLST file', 2 )
 
                ELSEIF( INY .LT. 1950 .OR. INY .GT. 2030 ) THEN
 
                    CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &               'INVYEAR packet set 4-digit year outside range ' //
     &               '1950 to 2030 in MNLST file', 2 )
 
                ENDIF

                GO TO 55                        ! To head of MNLST read loop
 
            ELSE
 
                FNAME = LINE( 1:L2 )   ! Set inventory file name
 
            ENDIF

C.............  Get file name; open input non-link mobile source file
            ADEV = JUNIT()
            OPEN( ADEV, ERR=1006, FILE=FNAME, STATUS='OLD' )

            WRITE( MESG,94010 )
     &             'Successful OPEN using year', INY,
     &             'for inventory file:' // CRLF() // BLANK5 //
     &             FNAME( 1:TRIMLEN( FNAME ) )
            CALL M3MSG2( MESG )

C.............  If have gotten this far, then we have non-link data for sure
            NFLAG = .TRUE.

C.............  Read county-VMT file(s) 

	    IREC   = 0     !  Counter for MRAW records

   66       CONTINUE       !  Head of the ADEV-read loop

		READ( ADEV, *, END=110, IOSTAT=IOS ) FIP, ROAD, TVMT

		IREC = IREC + 1

C.................  Abort if error reading file
		IF( IOS .NE. 0 ) THEN
		    WRITE( MESG, 94010 ) 'Error ', IOS, 
     &                   'reading ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                   ' at line', IREC
		    CALL M3EXIT( 'RAWMOBIL', 0, 0, MESG, 2 )

C.................  Determine IRCLAS based on Urban/Rural and Facility type
C.................  and check for bad area/facility combination
                ELSEIF ( ROAD .LT. 100 ) THEN
     
                    WRITE( RR, '( I2.2 )' ) ROAD 
                    J = INDEX1( RR, NRCLAS, MARFAC3 )

                    IF( J .GT. 0 ) THEN
                        RCL = MROADS3( J )
    
                    ELSE      ! Error, road class not found

                        WRITE( MESG, 94010 )
     &                         'Skipping record', IREC,
     &                         ': Bad area/facility combination ' // RR
                        CALL M3MESG( MESG )

                        GO TO 66  ! To head of current ADEV read loop

                    ENDIF

C.................  Check for bad road class number
                ELSE

                    J = FIND1( ROAD, NRCLAS, MROADS3 )

                    IF( J .LE. 0 ) THEN
                        WRITE( MESG, 94010 )
     &                         'Skipping record', IREC,
     &                         ': Bad road class', ROAD
                        CALL M3MESG( MESG )

                        GO TO 66  ! To head of current ADEV read loop
                    ELSE
                        RCL = ROAD

                    ENDIF
                ENDIF

                S = S + 1
                N = N + 1   

                IF( N .GT. MXNSRC ) THEN  ! Error

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                     'Maximum number of NON-LINK SOURCES', MXNSRC,
     &                     'exceeded at line', IREC, ', file ' // 
     &                      FNAME( 1:TRIMLEN( FNAME ) )
                    CALL M3MESG( MESG )

                ELSE

                    INDEXA( S ) = S

                    IFIPA ( S ) = FIP
                    IRCLSA( S ) = RCL
                    ILINKA( S ) = 0
                    TPFLGA( S ) = MTPRFAC*WTPRFAC
                    INVYRA( S ) = INY
                    XLOC1A( S ) = BADVAL3
                    XLOC2A( S ) = BADVAL3
                    YLOC1A( S ) = BADVAL3
                    YLOC2A( S ) = BADVAL3
                    VMTA  ( S ) = TVMT

                ENDIF !  overflow or not

            GO TO 66  ! To head of current ADEV-read loop

  110       CLOSE( ADEV )

        GO TO 55      ! To head of MNLST read loop

  111   CONTINUE

        MNNUM = N
        MTOT  = S 

C.........  Check for overflow and drop sources if necessary
        IF ( MNNUM .GT. MXNSRC ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'Records count for NON-LINK sources   ', MNNUM ,
     &              CRLF() // BLANK5 //
     &             'Max dimensioned record-count (MXNSRC)', MXNSRC,
     &              CRLF() // BLANK5 //
     &             'Max record-count exceeded.'
            CALL M3MSG2( MESG )

            CALL M3MESG( 'Max record-count exceeded in RAW VMT files')

        ELSEIF( MTOT .EQ. 0 ) THEN

            EFLAG = .TRUE.
            MESG = 'Blank or missing VMT input data - no sources found'
            CALL M3MSG2( MESG )

        ENDIF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &                   'Error reading RAW VMT file(s)', 2 )
        END IF

        CALL M3MSG2( 'Processing the RAW MOBILE SOURCE data...' )

C.........  Sort raw VMT data (both link and non-link) and check for cases
C           with link VMT for FIP/RDT using county data

        CALL SORTI3( MTOT, INDEXA, IFIPA, IRCLSA, ILINKA )

C.........  Merge and store sorted link and non-link sources

        LFIP  = -9      ! Previous FIP
        LRCL  = -9      ! Previous road class
        LLNK  = -9      ! Previous link ID
        LINY  = -9      ! Previous inventory year
        NLFIP = -9      ! Previous non-link FIP
        NLRCL = -9      ! Previous non-link road class
        S     = 0
        DO 122 I = 1, MTOT

            J   = INDEXA( I )

            FIP = IFIPA ( J )
            RCL = IRCLSA( J )
            LNK = ILINKA( J )
            INY = INVYRA( J )

C.............  Same source as previous iteration
            IF( FIP .EQ. LFIP .AND.
     &          RCL .EQ. LRCL .AND.
     &          LNK .EQ. LLNK       ) THEN  ! Error

C.................  Same inventory year
                IF( INY .EQ. LINY ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                'ERROR: Duplicate sources for'// CRLF()// BLANK5//
     &                BLANK5// 'FIP:', FIP, 'RCL:', RCL,
     &                'LNK:', LNK, 'INY:', INY
                    CALL M3MESG( MESG )
                    GO TO 122  ! To end of loop

C.................  Different inventory year
                ELSE

                    WRITE( MESG,94010 ) 
     &                'WARNING: Not recommended to use same source ' //
     &                'for multiple years...' // CRLF() // BLANK5 // 
     &                BLANK5 // 'FIP:', FIP, 'RCL:', RCL,
     &                'LNK:', LNK, 'YR(1):', INY, 'YR(2):', LINY
                    CALL M3MESG( MESG )

                ENDIF

C.............  Same FIP/RDT for link as previous non-link is okay, with warning
            ELSEIF( DFLAG           .AND.
     &              OVERFLAG        .AND. ! use here to prevent excessive output
     &              FIP  .EQ. NLFIP .AND.
     &              RCL  .EQ. NLRCL .AND.
     &              LNK  .NE. 0           ) THEN  ! Warning

                OVERFLAG = .FALSE. ! Turn off reporting

                WRITE( MESG,94010 )
     &               'WARNING: Link VMT added to non-link VMT for' 
     &                 // CRLF() // BLANK5 //
     &                 'FIP:', FIP, 'RCL:', RCL
                CALL M3MESG( MESG )

C.............  Same FIP/RDT for link as previous non-link, then non-link unused
            ELSEIF( .NOT. DFLAG      .AND.   ! link superceeds non-link
     &              OVERFLAG         .AND.   ! iff okay to overwrite (1st link)
     &              FIP   .EQ. NLFIP .AND.   ! same as previous non-link FIP
     &              RCL   .EQ. NLRCL .AND.   ! same as previous non-link RDT
     &              LNK   .NE. 0           ) THEN  ! Warning

                S = S - 1          ! Back up one so link will overwrite non-link
                OVERFLAG = .FALSE. ! Turn off one-back-step per FIP/RDT

                WRITE( MESG,94010 )
     &               'WARNING: Link VMT overwriting non-link VMT for'
     &               // CRLF() // BLANK5 //
     &               'FIP:', FIP, 'RCL:', RCL
                CALL M3MESG( MESG )

            ELSE

C.................  If new FIP/RDT combo, store FIP/RDT for last encountered
C                   non-link source
                IF( ( FIP .NE. LFIP .OR.
     &                RCL .NE. LRCL     ) .AND.
     &                LNK .EQ. 0                ) THEN

                    NLFIP = FIP
                    NLRCL = RCL
                    OVERFLAG = .TRUE.  ! Reset back-step flag

                ENDIF

            ENDIF

C.............  Update previous iteration stored values with this iteration
            LFIP = FIP
            LRCL = RCL
            LLNK = LNK
            LINY = INY

            S = S + 1

C.............  Store source information if not out of bounds yet
            IF( S .LE. NMSRC ) THEN
                IFIP  ( S ) = FIP
                IRCLAS( S ) = RCL
                ILINK ( S ) = LNK
                INVYR ( S ) = INY
                TPFLAG( S ) = TPFLGA( J )
                XLOC1 ( S ) = XLOC1A( J )
                XLOC2 ( S ) = XLOC2A( J )
                YLOC1 ( S ) = YLOC1A( J )
                YLOC2 ( S ) = YLOC2A( J )
                VMT   ( S ) = VMTA  ( J )
            ENDIF

122     CONTINUE   ! To head of MOBILE source processing loop

        NSRC = S

        IF( NSRC .NE. NMSRC ) THEN

            EFLAG = .TRUE.
            WRITE( MESG, 94010 )
     &             'Actual Record count                  ', NSRC,
     &              CRLF() // BLANK5 //
     &             'Compiled record-count setting (NMSRC)', NMSRC

            CALL M3MSG2( MESG )
            CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &        'Compiled record count bad for MOBILE SOURCES.', 2 )

        ENDIF

        IF( EFLAG ) THEN
            CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &                   'Error processing RAW MOBILE SOURCE data', 2 )
        ENDIF

C.........  Process the time zone file
 
        CALL M3MSG2( 'Processing TIME ZONE file...' )
        IREC = 0
 
  144   CONTINUE        !  head of time-zone loop
 
            READ( ZDEV,*,END=188,IOSTAT=IOS ) FIP, TZONE
            IREC = IREC + 1
 
            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 )
     &              'Unit number', ZDEV,
     &              'I/O Status ', IOS,
     &              'Line number', IREC
                CALL M3MESG( MESG )
                CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &              'Error reading TIME ZONE file.', 2 )
            END IF
 
            IF ( FIP .EQ. 0 ) THEN              !  fallback -- all sources
 
                DO  155  S = 1, NMSRC
                    TZONES( S ) = TZONE
  155           CONTINUE
 
            ELSE IF ( MOD( FIP, 100 ) .EQ. 0 ) THEN     !  state-specific zone
 
                FIP = FIP / 1000
                DO  166  S = 1, NMSRC
                    IF ( FIP .EQ. IFIP( S ) / 1000 ) THEN
                        TZONES( S ) = TZONE
                    END IF
  166           CONTINUE
 
            ELSE                                        !  county-specific zone
 
                DO  177  S = 1, NMSRC
                    IF ( FIP .EQ. IFIP( S ) ) THEN
                        TZONES( S ) = TZONE
                    END IF
  177           CONTINUE
 
            END IF      !  if fip zero, or nn000, or not.
 
            GO TO  144
 
  188   CONTINUE        !  end of time-zone loop

C.........  Get file name; open output mobile sources file

        FTYPE3D = GRDDED3
        P_ALP3D = DBLE( UZONE )
        P_BET3D = DBLE( BADVAL3 )
        P_GAM3D = DBLE( BADVAL3 )
        XCENT3D = 0.0D0
        YCENT3D = 0.0D0
        XORIG3D = DBLE( BADVAL3 )
        YORIG3D = DBLE( BADVAL3 )
        SDATE3D = 0 !  n/a
        STIME3D = 0 !  n/a
        TSTEP3D = 0             !  time independent
        NVARS3D = 11
        NCOLS3D = 1
        NROWS3D = NMSRC     !  number of rows = # of mobile sources.
        NLAYS3D = 1
        NTHIK3D = 1
        GDTYP3D = UTMGRD3
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        WRITE( GDNAM3D, 94040 ) 'UTM_', UZONE

        FDESC3D( 1 ) = 'Mobile Sources VMT inventory file'

        DO 199 J = 2, MXDESC3
            FDESC3D( J ) = ' '
  199   CONTINUE
           
        J = 1
        VNAME3D( J ) = 'IFIP'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'State and county code'
        VTYPE3D( J ) = M3INT
        J = J + 1
           
        VNAME3D( J ) = 'IRCLAS'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Mobile source road class'
        VTYPE3D( J ) = M3INT
        J = J + 1
           
c       VNAME3D( J ) = 'IVTYPE'
c       UNITS3D( J ) = 'n/a'
c       VDESC3D( J ) = 'Mobile source vehicle type value'
c       VTYPE3D( J ) = M3INT
c       J = J + 1

        VNAME3D( J ) = 'TZONES'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Time zone of inventory data'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'TPFLAG'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Temporal profile category flag'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'INVYR'
        UNITS3D( J ) = 'YYYY'
        VDESC3D( J ) = 'Year of inventory data'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'ILINK'
        UNITS3D( J ) = 'n/a'
        VDESC3D( J ) = 'Mobile source link number for link sources'
        VTYPE3D( J ) = M3INT
        J = J + 1

        VNAME3D( J ) = 'XLOC1'
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'UTM X coordinate for beginning of link'
        VTYPE3D( J ) = M3REAL
        J = J + 1

        VNAME3D( J ) = 'YLOC1'
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'UTM Y coordinate for beginning of link'
        VTYPE3D( J ) = M3REAL
        J = J + 1

        VNAME3D( J ) = 'XLOC2'
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'UTM X coordinate for end of link'
        VTYPE3D( J ) = M3REAL
        J = J + 1

        VNAME3D( J ) = 'YLOC2'
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'UTM Y coordinate for end of link'
        VTYPE3D( J ) = M3REAL
        J = J + 1

        VNAME3D( J ) = 'VMT'
        UNITS3D( J ) = 'mi/day'
        VDESC3D( J ) = 'Vehicle Miles Traveled per typical day'
        VTYPE3D( J ) = M3REAL

        ENAME = PROMPTMFILE( 
     &          'Enter logical name for MOBILE SOURCE output file',
     &          FSUNKN3, 'MOBL', 'RAWMOBIL' )

C.......   Write out the mobile source emissions values:
        
        CALL M3MSG2( 'Writing out MOBILE SOURCES file...' )

        IF ( .NOT. WRITE3( ENAME, 'ALL', 0, 0,  IFIP ) ) THEN
            CALL M3EXIT( 'RAWMOBIL', 0, 0,  
     &                   'Error writing MOBILE OUTPUT file' , 2 )
        END IF

  999   CONTINUE          !  exit program

        CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &               'Normal completion of Program RAWMOBIL', 0 )


1006    MESG = 'Error opening file ' // FNAME( 1:TRIMLEN( FNAME ) )
        CALL M3EXIT( 'RAWMOBIL', 0, 0, MESG, 2 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I10 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93100   FORMAT( 100( A5, 1X ) )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 5( A, :, I8, :, 1X ) )

94020   FORMAT( A, 1X, I5, 1X, A, I7.7, I3.3, 1X, A, I5 )

94040   FORMAT( A, I2.2 )

        END

