
C Version "@(#)$Id$ $Source$ $Date$ 

         PROGRAM RAWMOBIL

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C       
C       
C
C  PRECONDITIONS REQUIRED:
C      
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C	PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL
C
C  REVISION  HISTORY:
C       Prototype  12/95 by MRH
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
       
        INCLUDE 'TMDIMS3.EXT'
        INCLUDE 'MBDIMS3.EXT'
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         ENVYN
        INTEGER         FIND1, FIND2, FIND3
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         LEN2
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        
        EXTERNAL  ENVYN, FIND1, FIND2, FIND3, GETNUM, GETYN, INDEX1,
     &            LEN2, PROMPTFFILE, PROMPTMFILE, TRIMLEN

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

C.........  Unsorted non-link source variables
        INTEGER     INDEXA ( MXNSRC )  ! Non-link sorting index
        INTEGER     IFIPA  ( MXNSRC )
        INTEGER     IRCLSA ( MXNSRC )
        INTEGER     TPFLGA ( MXNSRC )
        INTEGER     INVYRA ( MXNSRC )
        INTEGER     ILINKA ( MXNSRC )
        REAL        XLOC1A ( MXNSRC )
        REAL        XLOC2A ( MXNSRC )
        REAL        YLOC1A ( MXNSRC )
        REAL        YLOC2A ( MXNSRC )
        REAL        VMTA   ( MXNSRC )

C.........  Unsorted link sources variables

        INTEGER     LINDXA ( MXLSRC )         ! Sorting index for links list
        INTEGER     LFIPA  ( MXLSRC )         ! Unsorted LFIP
        INTEGER     LRCLSA ( MXLSRC )         ! Unsorted LRCL
        INTEGER     LLINKA ( MXLSRC )         ! Unsorted LLINK
        INTEGER     LTPFGA ( MXLSRC )         ! Unsorted temporal profile type
        INTEGER     LINVYA ( MXLSRC )         ! Unsorted YEAR
        REAL        LX1A   ( MXLSRC )         ! Unsorted LX1
        REAL        LX2A   ( MXLSRC )         ! Unsorted LX2
        REAL        LY1A   ( MXLSRC )         ! Unsorted LY1
        REAL        LY2A   ( MXLSRC )         ! Unsorted LY2
        REAL        LVMTA  ( MXLSRC )         ! Unsorted LVMT

C.........  File units and logical names

        INTEGER       ADEV              ! Unit for non-link data
        INTEGER       KDEV              ! unit for link data
        INTEGER       LDEV              ! IO/API initialization unit
        INTEGER       ZDEV              ! Zones file unit
        CHARACTER*16  ENAME             ! Mobile source output file name

C.........  Other Local variables

        INTEGER     IREC, J, L, LC, N, NC, S  ! Counters and pointers

        INTEGER     FIP               ! Temporary read in FIP
        INTEGER     IFIL              ! File name counter
        INTEGER     IOS               ! Error status on read
        INTEGER     INYEAR            ! Inventory year for current file
        INTEGER     LFIP              ! Previous FIP in loop
        INTEGER     LRCL              ! Previous road class in loop
        INTEGER     LLNK              ! Previous link ID in loop
        INTEGER     LYR               ! Previous year in loop
        INTEGER     LNK               ! Temporary link number
        INTEGER     LNK_FIP           ! Fip for link source
        INTEGER     LNK_RCL           ! Road class for link source
        INTEGER     NLK_FIP           ! Fip for non-link source
        INTEGER     NLK_RCL           ! Road class for non-link source
        INTEGER     TPF               ! Temporary temporal flag
        INTEGER     IYR               ! Temporary inventory year
        INTEGER     MLNUM             ! Actual number link sources
        INTEGER     MNNUM             ! Actual number non-link sources
        INTEGER     NSRC              ! Actual number of total sources
        INTEGER     RCL               ! Temporary road class
        INTEGER     ROAD              ! Temporary read-in road code
        INTEGER     TZONE             ! Input TZONE for output data file
        INTEGER     UZONE             ! UTM zone to put into file header
        INTEGER     ZONE              ! Temporary UTM zone for link coords

        REAL        TVMT              ! Temporary read-in VMT
        REAL        X1, X2, Y1, Y2    ! Temporary link start and end coord [m]
        REAL        XX, YY            ! Temporary LAT/LON coords

        LOGICAL     EFLAG             ! Error flag
        LOGICAL     LFLAG             ! True if link sources exist
        LOGICAL     NFLAG             ! True if non-link sources exist

        DATA    EFLAG  / .FALSE. /

        CHARACTER*2   RR              ! Temporary area/facility type
        CHARACTER*8   LTYPE           ! The previous type of source in process
        CHARACTER*8   TYPE            ! The current type of source in process
        CHARACTER*16  RNAME           ! Temporary logical raw VMT file name
        CHARACTER*16  LNAME           ! Temporart logical link VMT file name
        CHARACTER*256 MESG
       
C***********************************************************************
C   begin body of program RAWMOBIL
        
        LDEV = INIT3()
        
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
            WRITE( *,92000 ) ' ', 'Ending program RAWMOBIL', ' '
            CALL EXIT( 0 )
        END IF

C.........  Get file name; open input time zone file
        ZDEV = PROMPTFFILE( 
     &           'Enter logical name for TIME ZONE file',
     &           .TRUE., .TRUE., 'ZONES', 'RAWMOBIL' )

C.........  Get the UTM zone to put into file header:

        UZONE = GETNUM( 1, 36, 17, 'Enter UTM zone for this file' )

C.........  Loop over link input files

        WRITE( *,92000 )
     &  ' ',
     &  'Now enter the list of logical file names for',
     &  'RAW LINK-VMT (LNKRAW) MOBILE SOURCE files, and',
     &  'year for each set',
     &  ' '

        IFIL   = 0
        INYEAR = 1989
        LFLAG  = .FALSE.
        S      = 0

   22   CONTINUE

            IFIL = IFIL + 1
            WRITE( LNAME, 94040 ) 'LNKRAW', IFIL

C.............  Get file name; open input link sources file
            KDEV = PROMPTFFILE( 
     &             'Enter the NEXT LNKRAW file or "NONE" to end list',
     &             .TRUE., .TRUE., LNAME, 'RAWMOBIL' )

C.............  End list if NONE entered
            IF( KDEV .LE. 0 ) GO TO 51

C.............  If KDEV ever positive, then we have link data
            LFLAG = .TRUE.
             
C.............  Get year for LNKRAW
            INYEAR = GETNUM( IMISS3, -IMISS3, INYEAR + 1,
     &                        'Enter INVENTORY YEAR for file' )

C.............  Read file for sources with link specific VMT (i.e. LINK sources)

            IREC = 0

   44       CONTINUE  ! Head of loop for reading link specific VMT data

                READ( KDEV, *, END=50, IOSTAT=IOS ) 
     &                FIP, ROAD, LNK, X1, Y1, X2, Y2, ZONE, TVMT

                IREC = IREC + 1

                IF( IOS .NE. 0 ) THEN
                    WRITE( MESG, 94010 ) 'Error ', IOS, 
     &                   'reading ' // LNAME // ' at line', IREC
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

C.................  Ensure no overflow
                IF( S .GE. MXLSRC ) THEN  ! Error

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                     'Maximum number of LINK SOURCES', MXLSRC,
     &                     'exceeded at line', IREC, ', file ' // LNAME
                    CALL M3MESG( MESG )

                ELSE

                    S = S + 1

                    LINDXA( S ) = S     ! Index for sorting

                    LFIPA ( S ) = FIP
                    LRCLSA( S ) = RCL
                    LLINKA( S ) = LNK
                    LTPFGA( S ) = MTPRFAC*WTPRFAC
                    LINVYA( S ) = INYEAR
                    LX1A  ( S ) = X1
                    LX2A  ( S ) = X2
                    LY1A  ( S ) = Y1
                    LY2A  ( S ) = Y2
                    LVMTA ( S ) = TVMT

                ENDIF

            GO TO 44  ! To head of current KDEV read loop

   50   GO TO 22      ! To head of loop on LNKRAW files

   51   CONTINUE

        MLNUM = S

        IF( EFLAG ) THEN   ! Internal dimensions exceeded

            WRITE( *, 92000 ) ' '
            WRITE( *, 92010 ) 
     &             'Records count for LINK sources       :', MLNUM ,
     &             'Max dimensioned record-count (MXLSRC):', MXLSRC,
     &             'Max record-count exceeded.'
            CALL M3MESG( 'Max record-count exceeded in RAW VMT files')

        ENDIF

C.........  Loop over non-link input files

        WRITE( *,92000 )
     &  ' ',
     &  'Now enter the list of logical file names for',
     &  'RAW COUNTY-VMT (MRAW) MOBILE SOURCE files,',
     &  'year for each set',
     &  ' '

        IFIL   = 0
        INYEAR = 1989
        S      = 0

   55   CONTINUE

            IFIL = IFIL + 1
            WRITE( RNAME, 94040 ) 'MRAW', IFIL

C.............  Get file name; open input non-link mobile source file
   65       ADEV = PROMPTFFILE( 
     &               'Enter the NEXT MRAW FILE, or "NONE" to end list',
     &               .TRUE., .TRUE., RNAME, 'RAWMOBIL' )

C.............  End list if NONE entered
            IF( ADEV .LE. 0 ) GO TO 111

C.............  If ADEV ever positive, then we have non-link data
            NFLAG = .TRUE.

C.............  Get year for LNKRAW
            INYEAR = GETNUM( IMISS3, -IMISS3, INYEAR + 1,
     &                        'Enter INVENTORY YEAR for file' )

C.............  Read county-VMT file(s) 

	    IREC   = 0     !  Counter for MRAW records

   66       CONTINUE       !  Head of the ADEV-read loop

		READ( ADEV, *, END=110, IOSTAT=IOS ) FIP, ROAD, TVMT

		IREC = IREC + 1

C.................  Abort if error reading file
		IF( IOS .NE. 0 ) THEN
		    WRITE( MESG, 94010 ) 'Error ', IOS, 
     &                   'reading ' // RNAME // ' at line', IREC
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

                IF( S .GE. MXNSRC ) THEN  ! Error

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                     'Maximum number of NON-LINK SOURCES', MXNSRC,
     &                     'exceeded at line', IREC, ', file ' // RNAME
                    CALL M3MESG( MESG )

                ELSE

                    S = S + 1   

                    INDEXA( S ) = S

                    IFIPA ( S ) = FIP
                    IRCLSA( S ) = RCL
                    ILINKA( S ) = 0
                    TPFLGA( S ) = MTPRFAC*WTPRFAC
                    INVYRA( S ) = INYEAR
                    XLOC1A( S ) = BADVAL3
                    XLOC2A( S ) = BADVAL3
                    YLOC1A( S ) = BADVAL3
                    YLOC2A( S ) = BADVAL3
                    VMTA  ( S ) = TVMT

                ENDIF !  overflow or not

            GO TO 66  ! To head of current ADEV-read loop

  110   GO TO 55      ! To head of loop on MRAW files

  111   CONTINUE

        MNNUM = S 

C.........  Check for overflow and drop sources if necessary
        IF ( EFLAG ) THEN

            WRITE( *, 92000 ) ' '
            WRITE( *,92010 )
     &             'Records count for NON-LINK sources   ', MNNUM ,
     &             'Max dimensioned record-count (MXNSRC)', MXNSRC,
     &             'Max record-count exceeded.'
            CALL M3MESG( 'Max record-count exceeded in RAW VMT files')

        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &                   'Error reading RAW VMT file(s)', 2 )
        END IF

C.........  Sort and merge raw VMT data (both link and non-link)
        WRITE( *,92000 ) ' ', 'Processing the RAW MOBILE SOURCE data...', ' '

C.........  Sort link data
        IF( LFLAG ) THEN
            CALL SORTI3( MLNUM, LINDXA, LFIPA, LRCLSA, LLINKA )
        ENDIF

C.........  Sort non-link data
        IF( NFLAG ) THEN
            CALL SORTI2( MNNUM, INDEXA, IFIPA, IRCLSA )
        ENDIF

C.........  Reorganize non-link data
        LC    = 0       ! Link-source     counter
        NC    = 0       ! Non-link-source counter
        S     = 0       ! Total source    counter
c       OFIP  = 0
c       ORCL  = 0
        LFIP  = -9      ! Previous FIP
        LRCL  = -9      ! Previous road class
        LLNK  = -9      ! Previous link ID
        LYR   = -9      ! Previous inventory year
        LTYPE = 'LINK'  ! Initialize type for last source to link

  122   CONTINUE   ! Continue processing link and non-link sources

            IF( LC + 1 .GT. MLNUM .AND.
     &          NC + 1 .GT. MNNUM       ) THEN

                GO TO 133 ! End loop

            ENDIF

            IF( LC + 1 .GT. MLNUM ) LFLAG = .FALSE.
            IF( NC + 1 .GT. MNNUM ) NFLAG = .FALSE.

C.............  Process for link and non-link both exist
            IF( LFLAG .AND. NFLAG ) THEN

C.................  Set temporary source characteristics for link and non-link
                L = LINDXA( LC + 1 )

                LNK_FIP = LFIPA ( L )
                LNK_RCL = LRCLSA( L )

                N = INDEXA( NC + 1 )

                NLK_FIP = IFIPA ( N )
                NLK_RCL = IRCLSA( N )

C.................  Link sources superceed non-link 
                IF( LNK_FIP .LE. NLK_FIP .AND. 
     &              LNK_RCL .LE. NLK_RCL       ) THEN

                    TYPE = 'LINK'
                    LC   = LC + 1       ! Increase link counter
                    J    = L
                    FIP  = LNK_FIP
                    RCL  = LNK_RCL
                    LNK  = LLINKA( J )
                    TPF  = LTPFGA( J )
                    IYR  = LINVYA( J )
                    X1   = LX1A  ( J )
                    X2   = LX2A  ( J )
                    Y1   = LY1A  ( J )
                    Y2   = LY2A  ( J )
                    TVMT = LVMTA ( J )

C.....................  Skip non-link source if IDs are equal
                    IF( LNK_FIP .EQ. NLK_FIP .AND. 
     &                  LNK_RCL .EQ. NLK_RCL       ) NC = NC + 1

C.................  When no link, go to next non-link source
                ELSE

                    TYPE = 'NON-LINK'
                    NC   = NC + 1
                    J    = N
                    FIP  = NLK_FIP
                    RCL  = NLK_RCL
                    LNK  = ILINKA( J )
                    TPF  = TPFLGA( J )
                    IYR  = INVYRA( J )
                    X1   = XLOC1A( J )
                    X2   = XLOC2A( J )
                    Y1   = YLOC1A( J )
                    Y2   = YLOC2A( J )
                    TVMT = VMTA  ( J )

                ENDIF

C.............  Process for link only
            ELSEIF( LFLAG ) THEN

                TYPE = 'LINK'
                LC   = LC + 1
                J    = LINDXA( LC )
     
                FIP  = LFIPA ( J )
                RCL  = LRCLSA( J )
                LNK  = LLINKA( J )
                TPF  = LTPFGA( J )
                IYR  = LINVYA( J )
                X1   = LX1A  ( J )
                X2   = LX2A  ( J )
                Y1   = LY1A  ( J )
                Y2   = LY2A  ( J )
                TVMT = LVMTA ( J )

C.............  Process for non-link only
            ELSEIF( NFLAG ) THEN

                TYPE = 'NON-LINK'
                NC   = NC + 1
                J    = INDEXA( NC )

                FIP  = IFIPA ( J )
                RCL  = IRCLSA( J )
                LNK  = ILINKA( J )
                TPF  = TPFLGA( J )
                IYR  = INVYRA( J )
                X1   = XLOC1A( J )
                X2   = XLOC2A( J )
                Y1   = YLOC1A( J )
                Y2   = YLOC2A( J )
                TVMT = VMTA  ( J )

C.............  Error
            ELSE

                CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &                       'No raw input data files provided', 2 )
            ENDIF

C.............  Check for duplicates 
            IF( TYPE .EQ. LTYPE .AND.
     &          FIP  .EQ. LFIP  .AND. 
     &          RCL  .EQ. LRCL  .AND.
     &          LNK  .EQ. LLNK        ) THEN  ! Error

                EFLAG = .TRUE.
                CALL M3MESG( 'Duplicate source found:' )
                WRITE( MESG,94010 )
     &                 '    Current  source FIP:', FIP, 'RCL:', RCL,
     &                                     'LNK:', LNK, 'IYR:', IYR
                CALL M3MESG( MESG )

                WRITE( MESG,94010 )
     &                 '    Previous source FIP:', FIP, 'RCL:', RCL,
     &                                     'LNK:', LNK, 'IYR:', IYR
                CALL M3MESG( MESG )

            ELSEIF( S .GE. NMSRC ) THEN       ! Error

                S = S + 1

                EFLAG = .TRUE.
                CALL M3MESG( 
     &               'Maximum number of MOBILE SOURCES exceeded:' )
                WRITE( MESG,94010 )
     &                 '    FIP:', FIP, 'RCL:', RCL,
     &                     'LNK:', LNK, 'IYR:', IYR
                CALL M3MESG( MESG )

C.............  Update sources if no duplicate and no overflow
            ELSE

                LTYPE = TYPE
                LFIP  = FIP
                LRCL  = RCL
                LLNK  = LNK
                LYR   = IYR

                S = S + 1

                IFIP  ( S ) = FIP
                IRCLAS( S ) = RCL
                ILINK ( S ) = LNK
                TPFLAG( S ) = TPF
                XLOC1 ( S ) = X1
                XLOC2 ( S ) = X2
                YLOC1 ( S ) = Y1
                YLOC2 ( S ) = Y2
                INVYR ( S ) = IYR
                TPFLAG( S ) = TPF
                VMT   ( S ) = TVMT

            ENDIF

        GO TO 122  ! To head of MOBILE source processing loop

  133   CONTINUE

        NSRC = S

        IF( NSRC .NE. NMSRC ) THEN

            EFLAG = .TRUE.
            WRITE( *, 92000 ) ' '
            WRITE( *, 92010 )
     &             'Record count                        ', NSRC,
     &             'Max dimensioned record-count (NMSRC)', NMSRC
            CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &        'Max record-count exceeded in POINT SOURCE file.', 2 )

        ENDIF

        IF( EFLAG ) THEN
            CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &                   'Error processing RAW MOBILE SOURCE data', 2 )
        ENDIF

C.........  Process the time zone file
 
        WRITE( *,92000 ) ' ', 'Processing TIME ZONE file...', ' '
        IREC = 0
 
  144   CONTINUE        !  head of time-zone loop
 
            READ( ZDEV,*,END=188,IOSTAT=IOS ) FIP, TZONE
            IREC = IREC + 1
 
            IF ( IOS .NE. 0 ) THEN
                WRITE( *,92010 )
     &              'Unit number', ZDEV,
     &              'I/O Status ', IOS,
     &              'Line number', IREC
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

        FDESC3D( 1 ) = 'NC Mobile Sources VMT values.'
           
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
        
        WRITE( *,92000 ) ' ', 'Writing out MOBILE SOURCES file...', ' '

        IF ( .NOT. WRITE3( ENAME, 'ALL', 0, 0,  IFIP ) ) THEN
            CALL M3EXIT( 'RAWMOBIL', 0, 0,  
     &                   'Error writing MOBILE OUTPUT file' , 2 )
        END IF

  999   CONTINUE          !  exit program

        CALL M3EXIT( 'RAWMOBIL', 0, 0,
     &               'Normal completion of Program RAWMOBIL', 0 )


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

