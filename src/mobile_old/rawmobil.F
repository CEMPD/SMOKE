
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
       
        INCLUDE 'CHDIMS3.EXT' 
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

C.........  Source variables

        INTEGER     IFIP   ( NMSRC )  ! Source FIP
        INTEGER     IRCLAS ( NMSRC )  ! Source Road Class
c       INTEGER     IVTYPE ( NMSRC )  ! Source Vehicle type
        INTEGER     TZONES ( NMSRC )  ! Source time zone number
        INTEGER     INVYR  ( NMSRC )  ! Source inventory year
        INTEGER     ILINK  ( NMSRC )  ! Source Link number
        REAL        XLOC1  ( NMSRC )  ! UTM X coord of link start [m]   
        REAL        XLOC2  ( NMSRC )  ! UTM X coord of link end [m]
        REAL        YLOC1  ( NMSRC )  ! UTM Y coord of link start [m]
        REAL        YLOC2  ( NMSRC )  ! UTM Y coord of link end [m]
        REAL        VMT    ( NMSRC )  ! Source VMT [mi/typical day]

c       COMMON  / EMISREC / IFIP, IRCLAS, IVTYPE, TZONES, INVYR, ILINK, 
c    &                      XLOC1, XLOC2, YLOC1, YLOC2, VMT
        COMMON  / EMISREC / IFIP, IRCLAS, TZONES, INVYR, ILINK, 
     &                      XLOC1, XLOC2, YLOC1, YLOC2, VMT

C.........  Sorted link sources variables

        INTEGER     LFIP   ( MXLSRC )         ! Sorted link sources list FIP
        INTEGER     LIDX1  ( MXLSRC )         ! Sorting index for links list
        INTEGER     LRCL   ( MXLSRC )         ! Sorted link sources list RCL
        INTEGER     LLNK   ( MXLSRC )         ! Sorted user-spec link number
        REAL        LX1    ( MXLSRC )         ! Sorted UTM link sources list X1
        REAL        LX2    ( MXLSRC )         ! Sorted UTM link sources list X2
        REAL        LY1    ( MXLSRC )         ! Sorted UTM link sources list Y1
        REAL        LY2    ( MXLSRC )         ! Sorted UTM link sources list Y2
        REAL        LVMT   ( MXLSRC )         ! Sorted Link specific VMT values

C.........  Unsorted link sources variables

        INTEGER     LFIPA  ( MXLSRC )         ! Unsorted LFIP
        INTEGER     LRCLA  ( MXLSRC )         ! Unsorted LRCL
        INTEGER     LLNKA  ( MXLSRC )         ! Unsirted LLNK
        REAL        LX1A   ( MXLSRC )         ! Unsorted LX1
        REAL        LX2A   ( MXLSRC )         ! Unsorted LX2
        REAL        LY1A   ( MXLSRC )         ! Unsorted LY1
        REAL        LY2A   ( MXLSRC )         ! Unsorted LY2
        REAL        LVMTA  ( MXLSRC )         ! Unsorted LVMT

C......... Other link sources variables

        INTEGER     MLNUM ! Actual number of FIP/RCL combos in link sources file

C.........  Other Non-dimensioned variables

        INTEGER     FIP               ! Temporary read in FIP
        INTEGER     INYEAR            ! Inventory year for current file
        INTEGER     LNK               ! Temporary link number
        INTEGER     RCL               ! Temporary road class
        INTEGER     ROAD              ! Temporary read-in road code
        INTEGER     SRC               ! Source counter
        INTEGER     TZONE             ! Input TZONE for output data file
        INTEGER     UZONE             ! UTM zone to put into file header
        INTEGER     ZONE              ! Temporary UTM zone for link coords

        REAL        TVMT              ! Temporary read-in VMT
        REAL        X1, X2, Y1, Y2    ! Temporary link start and end coord [m]
        REAL        XX, YY            ! Temporary LAT/LON coords

        LOGICAL QFLAG                 ! Flag for quitting the program
        DATA QFLAG  / .FALSE. /

C...........   Other Local variables:

        INTEGER ADEV, KDEV, LDEV, ZDEV
        INTEGER ICNT, IL, IOS, IREC, IS,
     &          J, K, LL, LCNTM1, LNEXT

        LOGICAL SFLAG

        CHARACTER*16  ENAME
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

        SFLAG = ENVYN( 'RAW_SRC_CHECK', 
     &                 'RAWMOBIL check for missing species-records',
     &                 .FALSE.,
     &                 IOS )
        
C.........  Get file name; open input time zone file
        ZDEV = PROMPTFFILE( 
     &           'Enter logical name for TIME ZONE file',
     &           .TRUE., .TRUE., 'ZONES', 'RAWMOBIL' )

C.........  Get the UTM zone to put into file header:

        UZONE = GETNUM( 1, 36, 17, 'Enter UTM zone for this file' )

        INYEAR = 1989
C.........  Loop over input files

C.........  Get file name; open input raw mobile source file
        ADEV = PROMPTFFILE( 
     &           'Enter logical name for RAW MOBILE SOURCES file',
     &           .TRUE., .TRUE., 'MRAW', 'RAWMOBIL' )

C.........  Get file name; open input link sources file
        KDEV = PROMPTFFILE( 
     &           'Enter logical name for RAW LINK SOURCES file',
     &           .TRUE., .TRUE., 'LNKRAW', 'RAWMOBIL' )

C.........  Get year for MRAW and LNKRAW
        INYEAR = GETNUM( IMISS3, -IMISS3, INYEAR + 1,
     &                    'Enter INVENTORY YEAR for these files' )

C.........  Read file for sources with link specific VMT (i.e. LINK sources)

        IL = -1
        ICNT = 0
   44   CONTINUE  ! Head of loop for reading link specific VMT data

            READ( KDEV, *, END=50 ) 
     &            FIP, RCL, LNK, X1, Y1, X2, Y2, ZONE, TVMT

C.............  If link zone is not file zone, convert coords to zone
            IF( ZONE .NE. UZONE ) THEN

                CALL UTM2LL( X1, Y1,  ZONE, XX, YY )  ! X1, Y1 to LAT/LON
                CALL LL2UTM( XX, YY, UZONE, X1, Y1 )  ! X1, Y1 to UZONE
                CALL UTM2LL( X2, Y2,  ZONE, XX, YY )  ! X2, Y2 to LAT/LON
                CALL LL2UTM( XX, YY, UZONE, X2, Y2 )  ! X2, Y2 to UZONE

            ENDIF

C.............  Check if this FIP/RCL/LNK in existing sorted list
            IF( ICNT .GT. 0 ) 
     &          IL = FIND3( FIP, RCL, LNK, ICNT, LFIP, LRCL, LLNK )

            IF( IL .GT. 0 ) THEN  ! Error (shouldn't have been in list)

            ELSE

                ICNT = ICNT + 1

C.................  Ensure no overflow
                IF( ICNT .GT. MXLSRC ) THEN  ! Error
                    QFLAG = .TRUE.

                ELSE

                    LFIPA( ICNT ) = FIP
                    LRCLA( ICNT ) = RCL
                    LLNKA( ICNT ) = LNK
                    LX1A ( ICNT ) = X1
                    LX2A ( ICNT ) = X2
                    LY1A ( ICNT ) = Y1
                    LY2A ( ICNT ) = Y2
                    LVMTA( ICNT ) = TVMT
                    LIDX1( ICNT ) = ICNT  ! Index for sorting

C..................... Sort temporary LFIPA and LRCLA
                    CALL SORTI3( ICNT, LIDX1, LFIPA, LRCLA, LLNKA )

C..................... Keep link records sorted so FIND will work               
                    DO 46 LL = 1, ICNT
                        LFIP ( LL ) = LFIPA( LIDX1( LL ) )
                        LRCL ( LL ) = LRCLA( LIDX1( LL ) )
                        LLNK ( LL ) = LLNKA( LIDX1( LL ) )
                        LX1  ( LL ) = LX1A ( LIDX1( LL ) )
                        LX2  ( LL ) = LX2A ( LIDX1( LL ) )
                        LY1  ( LL ) = LY1A ( LIDX1( LL ) )
                        LY2  ( LL ) = LY2A ( LIDX1( LL ) )
                        LVMT ( LL ) = LVMTA( LIDX1( LL ) )
                        LIDX1( LL ) = LL     ! Re-initialize pointer
   46               CONTINUE 

                ENDIF  

            ENDIF

        GO TO 44           ! End of loop for reading link specific VMT data
   50   MLNUM = ICNT

        IF( QFLAG ) THEN   ! Boundaries exceeded

            WRITE( MESG, 94010 )
     &             'Records count for LINK sources', ICNT,
     &             'Max dimensioned record-count  ', MXLSRC,
     &             'Max record-count exceeded.'
            CALL M3EXIT( 'RAWMOBIL', 0, 0, MESG, 2 )

        ENDIF

C.........  Read raw file and fill in data
        WRITE( *,92000 ) ' ', 'Processing the RAW MOBILE SOURCE file...', ' '

        IREC   = 0     !  Counter for MRAW records
        SRC    = 0     !  Index of sources

   55   CONTINUE       !  Head of the ADEV-read loop

            IREC = IREC + 1
            READ( ADEV, *, END=111, IOSTAT=IOS ) FIP, ROAD, TVMT

            IF( SRC+1 .GT. NMSRC ) THEN
              QFLAG = .TRUE.

            ELSE

C.................  Determine IRCLAS based on Urban/Rural and Facility type
C               IC = FIND1( ROAD, NRCLAS, MROADS3 )
C               IF( IC .GT. 0 ) THEN
C                   RCL = IC
                    RCL = ROAD
C               ELSE      ! Error, road class not found; ABORT 
C                   CALL M3EXIT( )
C               ENDIF

C.................  Ensure source hasn't appearred before
                IL = FIND2( FIP, RCL, SRC, IFIP, IRCLAS ) ! Assume Sorted!

                IF( IL .GT. 0 ) THEN    ! Error

                ELSE
                    SRC = SRC + 1   
    
                    IFIP  ( SRC ) = FIP
                    IRCLAS( SRC ) = RCL
                    ILINK ( SRC ) = 0
                    XLOC1 ( SRC ) = AMISS3
                    XLOC2 ( SRC ) = AMISS3
                    YLOC1 ( SRC ) = AMISS3
                    YLOC2 ( SRC ) = AMISS3
                    INVYR ( SRC ) = INYEAR
                    VMT   ( SRC ) = TVMT
                ENDIF

            ENDIF                ! SRC+1 > NMSRC

        GO TO 55             ! To head of ADEV-read loop
  111   CONTINUE

C.........  Insert link data where necesary
        LL = 0
        DOWHILE( LL+1 .LE. MLNUM )

            LL = LL + 1

C.............  See if main raw list has link FIP/RCL
            IS = FIND2( LFIP( LL ), LRCL( LL ), SRC, IFIP, IRCLAS )

            IF( IS .GT. 0 ) THEN

C.................  Count how many links for this LFIP/LRCL
                LNEXT = LL + 1
                DOWHILE( LLNK( LNEXT ) .GT. LLNK( LL ) )
                    LNEXT = LNEXT + 1
                ENDDO
                LCNTM1 = LNEXT - LL - 1

C.................  Ensure no overflow in source arrays
                IF( SRC+LCNTM1 .GT. NMSRC ) THEN
                    QFLAG = .TRUE.

                ENDIF

C.................  Move non-link data to make room for link data
                DO 125 K = SRC, IS+1, -1
                    IFIP  ( K+LCNTM1 ) = IFIP  ( K ) 
                    IRCLAS( K+LCNTM1 ) = IRCLAS( K ) 
                    ILINK ( K+LCNTM1 ) = ILINK ( K )
                    XLOC1 ( K+LCNTM1 ) = XLOC1 ( K ) 
                    XLOC2 ( K+LCNTM1 ) = XLOC2 ( K ) 
                    YLOC1 ( K+LCNTM1 ) = YLOC1 ( K ) 
                    YLOC2 ( K+LCNTM1 ) = YLOC2 ( K ) 
                    INVYR ( K+LCNTM1 ) = INVYR ( K )
                    VMT   ( K+LCNTM1 ) = VMT   ( K )
  125           CONTINUE

C.................  Insert link data
                DO 133 K = 0,LCNTM1
                    IFIP  ( IS+K ) = IFIP  ( IS   )
                    IRCLAS( IS+K ) = IRCLAS( IS   ) 
                    ILINK ( IS+K ) = LLNK  ( LL+K )
                    XLOC1 ( IS+K ) = LX1   ( LL+K ) 
                    XLOC2 ( IS+K ) = LX2   ( LL+K ) 
                    YLOC1 ( IS+K ) = LY1   ( LL+K )
                    YLOC2 ( IS+K ) = LY2   ( LL+K )
                    INVYR ( IS+K ) = INVYR ( IS   )
                    VMT   ( IS+K ) = LVMT  ( LL+K )
  133           CONTINUE

C.................  Update Link source counter and total number of sources
                LL  = LL  + LCNTM1
                SRC = SRC + LCNTM1

            ELSE  ! Warning

                 WRITE( MESG, 94010 )
     &                  'Warning! Link source FIP:', LFIP( LL ), 
     &                  ' and RCL:', LRCL( LL ), 
     &                  ' do not exist in main sources list.'
                 CALL M3MESG( MESG )

            ENDIF ! Link-specific source found in non-link list

        ENDDO

C.........  Abort for dimensioning problems
        IF( QFLAG .OR. SRC .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &             'Records count for MOBILE sources', SRC,
     &             '. Max dimensioned record-count  ', NMSRC,
     &             '. Dimensions must be equal.'
            CALL M3EXIT( 'RAWMOBIL', 0, 0, MESG, 2 )

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
 
                DO  155  SRC = 1, NMSRC
                    TZONES( SRC ) = TZONE
  155           CONTINUE
 
            ELSE IF ( MOD( FIP, 100 ) .EQ. 0 ) THEN     !  state-specific zone
 
                FIP = FIP / 1000
                DO  166  SRC = 1, NMSRC
                    IF ( FIP .EQ. IFIP( SRC ) / 1000 ) THEN
                        TZONES( SRC ) = TZONE
                    END IF
  166           CONTINUE
 
            ELSE                                        !  county-specific zone
 
                DO  177  SRC = 1, NMSRC
                    IF ( FIP .EQ. IFIP( SRC ) ) THEN
                        TZONES( SRC ) = TZONE
                    END IF
  177           CONTINUE
 
            END IF      !  if fip zero, or nn000, or not.
 
            GO TO  144
 
  188   CONTINUE        !  end of time-zone loop

C.........  Get file name; open output mobile sources file

        FTYPE3D = GRDDED3
        P_ALP3D = DBLE( UZONE )
        P_BET3D = DBLE( AMISS3 )
        P_GAM3D = DBLE( AMISS3 )
        XCENT3D = 0.0D0
        YCENT3D = 0.0D0
        XORIG3D = DBLE( AMISS3 )
        YORIG3D = DBLE( AMISS3 )
        SDATE3D = 0 !  n/a
        STIME3D = 0 !  n/a
        TSTEP3D = 0             !  time independent
c       NVARS3D = 11
        NVARS3D = 10
        NCOLS3D = 1
        NROWS3D = NMSRC     !  number of rows = # of mobile sources.
        NLAYS3D = 1
        NTHIK3D = 1
        GDTYP3D = UTMGRD3
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3
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

94010   FORMAT( 5( A, :, I7, :, 1X ) )

94020   FORMAT( A, 1X, I5, 1X, A, I7.7, I3.3, 1X, A, I5 )

94040   FORMAT( A, I2.2 )

        END

