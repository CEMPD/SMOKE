 
        PROGRAM GETRECS

C***********************************************************************
C  program body starts at line 244
C
C  DESCRIPTION:
C       Searches for a specific source, for all sources in a 
C       specific cell, or for combinations of source keys.  It
C       creates a ascii file which lists all details about the
C       source including source number, grid cell, if found in
C       gridding matrix, temporalization factors, control factors,
C       inventory pollutant emissions, and model species emissions.
C
C  PRECONDITIONS REQUIRED:
C       Completed files for all point source processing stages
C       including PNTS, PTMP, PGMAT, PSMAT, PCMAT 
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, GETNUM, INDEX1, PROMPTFFILE, PROMPTMFILE, GETYN, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  10/95 by M Houyoux
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
C***********************************************************************
 
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'CONST3.EXT'    !  I/O API parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'VGRD3.EXT'     !  vertical grid description

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER       FIND1
        INTEGER       GETMENU
        INTEGER       GETNUM
        LOGICAL       GETYN
        INTEGER       INDEX1
        INTEGER       LBLANK
        CHARACTER*14  MMDDYY
        INTEGER       PROMPTFFILE
        CHARACTER*16  PROMPTMFILE
        INTEGER       TRIMLEN

        EXTERNAL      FIND1, GETMENU, GETNUM, GETYN, INDEX1, LBLANK,
     &                MMDDYY, PROMPTFFILE, PROMPTMFILE, TRIMLEN

C.........  Define PNTS variables
        INTEGER     INDEXA( NPSRC )  !  index for sorting, if needed
        INTEGER     IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER     ISCC  ( NPSRC )  !  source SCC
        INTEGER     ISIC  ( NPSRC )  !  source SIC
        INTEGER     IPLANT( NPSRC )  !  Plant ID
        INTEGER     ISTACK( NPSRC )  !  Stack ID
        INTEGER     TZONES( NPSRC )  !  time zones
        INTEGER     TPFLAG( NPSRC )  !  applicability of temporal profile types
        INTEGER     INVYR ( NPSRC )  !  Inventory year
        REAL        XLOCA ( NPSRC )  !  UTM X-location (m)
        REAL        YLOCA ( NPSRC )  !  UTM Y-location (m)
        REAL        STKHT ( NPSRC )  !  stack height   (m)
        REAL        STKDM ( NPSRC )  !  stack diameter (m)
        REAL        STKTK ( NPSRC )  !  exhaust temperature (deg K)
        REAL        STKVE ( NPSRC )  !  exhaust velocity    (m/s)
        REAL        RULPEN( NPSRC, NIPOL )  !  rule penetration   fraction
        REAL        RULEFF( NPSRC, NIPOL )  !  rule effectiveness fraction
        REAL        CTLEFF( NPSRC, NIPOL )  !  control efficiency fraction
        REAL        EMISV ( NPSRC, NIPOL )  !  emissions values (tons/yr) 

        COMMON / EMISREC / IFIP, ISIC, ISCC, IPLANT, ISTACK, TZONES,
     &                    TPFLAG, INVYR, XLOCA, YLOCA, STKHT, STKDM, 
     &                    STKTK, STKVE, CTLEFF, RULEFF, RULPEN, EMISV

        CHARACTER*15  CPLANT( NPSRC )
        CHARACTER*15  CSTACK( NPSRC )
        CHARACTER*40  CPDESC( NPSRC )

C.........  Define temporary source variable

        REAL        CRITVAL( NPSRC ) !  Criterion for sorting output
        REAL        STKFL  ( NPSRC ) !  Stack flow

C.........  Define PTMP variables
        REAL        EMIST ( NPSRC )  !  ptmp emissions values (tons/yr) 
        REAL        EV

C.........  Define SMAT variables

        REAL            SS( NPSRC, NMPOL )
        INTEGER         NSPCOUT( NIPOL )        !  # of active rows in SPCREF
        INTEGER         SPCREF ( NIPOL, NMPOL )
        CHARACTER*16    INVNAME( NMPOL )        !  inv-spc part of name


C.........  Define GMAT variables

        INTEGER         NS( NGRID )
        INTEGER         IS( NPSRC )
 
        COMMON  / GRIDMAT / NS, IS

C.........  Define CMAT variables

        REAL            CC( NPSRC, NIPOL )

C.........  Define Layer Fractions

        REAL            LFRAC( NPSRC, EMLAYS )

C.........  Local parameters
        INTEGER         MXEOUT            ! Maximum number of output records
        INTEGER         MXTSTEP           ! Maximum number of time steps to proc
        INTEGER         NMETHOD           ! number of menu items 
        INTEGER         NRANGE            ! number of ranges in range menu
        INTEGER         NSTKMTHD          ! no. of ways to sort by stack parms

        PARAMETER     ( MXEOUT  = 2000, 
     &                  MXTSTEP = 144, ! 6 days
     &                  NMETHOD = 5   ,
     &                  NRANGE  = 3   ,
     &                  NSTKMTHD= 5      )

C.........  Local arrays
        INTEGER         EINX ( MXEOUT )   ! index for output sources
        INTEGER         ISLOC( NMATX  )   ! Cell specific IS values
        INTEGER         NOUT             ! actual number of output emis
        INTEGER         SPNT ( NMPOL )   ! pointer to speciation factors
        INTEGER         XNUM ( MXEOUT )   ! X-cell for output sources
        INTEGER         YNUM ( MXEOUT )   ! Y-cell for output sources

        REAL            EISUM( NIPOL )   ! Sum of day's inv pollutants
        REAL            EMSUM( NMPOL )   ! Sum of day's model species
        REAL            EOUT ( MXEOUT, NMPOL, MXTSTEP) ! output emissions
        REAL            ETMP ( MXEOUT, NIPOL, MXTSTEP) ! output temporal inv emission
        REAL            OFRAC( MXEOUT, EMLAYS, MXTSTEP)! output layer fractions

        CHARACTER*7     GSTATE( MXEOUT ) ! state of gridding matrix for output

C.........  Logical names and unit numbers
        INTEGER         LDEV             ! IO/API initialization unit #
        INTEGER         ODEV             ! number for output file
        INTEGER         SDEV             ! PSRC file

        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    CNAME   !  logical name for control    matrix input file
        CHARACTER*16    LNAME   !  logical name for layer fractions   input file
        CHARACTER*16    PNAME   !  logical name for point-source      input file
        CHARACTER*16    SNAME   !  logical name for speciation matrix input file
        CHARACTER*16    TNAME   !  logical name for temporal point source file

C.........  Other variables
        REAL            CFAC, DDX, DDY, DX, DY, SFAC, VAL
        REAL            X0, XMIN, XMAX, Y0, YMIN, YMAX, XX, YY
        REAL            MEANVAL      ! Mean of criterion output values
        REAL            MINVAL       ! Min non-zero of criterion output values
        REAL            MINDIFF      ! Minimum difference between src and mean
        REAL            THISDIFF     ! Current difference between src and mean

        INTEGER         I, ICNT, F, J, K, L, S, T, V   ! pointers and counters

        INTEGER         CELLNUM, RN, SCNT, TSTEP, XN, YN
        INTEGER         EMMTHD       ! method of sorting based on emissions
        INTEGER         ENDSRC       ! Ending source count for output
        INTEGER         FIP, PLANT, SCC, SIC, SID, STACK  ! Temporary vars
        INTEGER         INCSRC       ! Increment for processing sources
        INTEGER         IOS          ! I/O status
        INTEGER         JDATE, JTIME ! Current date and time
        INTEGER         L1, L2
        INTEGER         LDATE        ! Previous date
        INTEGER         MEANCNT      ! count of sources contibuting to mean val
        INTEGER         MEANSRC      ! source ID that has value closest to mean
        INTEGER         METHOD       ! method of source selection
        INTEGER         MINSRC       ! source ID that has min non-zero value
        INTEGER         NINCELL, NSTEPS
        INTEGER         NRPRT        ! number of sources to report
        INTEGER         RANGE        ! method of selecting range of sorted srcs
        INTEGER         SDATE, STIME ! Starting date and time
        INTEGER         SRTSRC       ! Starting source count for output
        INTEGER         STKMTHD      ! method of sorting based on stack parms

        LOGICAL         CTRLFLAG     ! control file input or not
        LOGICAL         EMSFLAG      ! PNTS is EMS-95 origin or not
        LOGICAL         LAYRFLAG     ! layer fractions file input or not
        LOGICAL         SPCFLAG      ! speciation matrix input or not

        CHARACTER*8     TONSUNIT
        CHARACTER*8     MOLEUNIT
        CHARACTER*15    CPLT, CSTK
        CHARACTER*16    SCRBUF  !  scratch point
        CHARACTER*40    CDESC
        CHARACTER*80    MENULST( NMETHOD )
        CHARACTER*80    RNGMENU( NRANGE )
        CHARACTER*80    STKMENU( NSTKMTHD )
        CHARACTER*256   MESG    !  scratch message buffer 

        DATA MENULST / 
     &       'Select all sources in specific grid cell',
     &       'Select sources by FIPs code, SCC, Plant, & Stack',
     &       'Select single source with source ID',
     &       'Select sources by stack parameters',
     &       'Select sources with X emissions of pollutant Y'    /

        DATA RNGMENU / 'Highest X number of sources', 
     &                 'Middle non-zero X number of sources', 
     &                 'Lowest non-zero X number of sources'  /

        DATA STKMENU / 'Stack Height',
     &                 'Stack Diameter',
     &                 'Stack Exit Temperature',
     &                 'Stack Exit Velocity',
     &                 'Stack Exit Flow Rate'   /

C***********************************************************************
C   begin body of program GETRECS

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 )
     &  ' ',
     &  'Program GETRECS searches for a specific source, for all',
     &  'sources in a specific cell, or for combinations of source',
     &  'keys.  It creates a ascii file which lists all details about',
     &  'the source including source number, grid cell, if found in',
     &  'gridding matrix, temporalization factors, control factors,',
     &  'inventory pollutant emissions, and model species emissions.',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are given in',
     &  'brackets [LIKE THIS] and can be accepted by hitting the',
     &  '<RETURN> key.',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            WRITE( *,92000 ) ' ', 'Ending program GETRECS', ' '
            CALL EXIT( 0 )
        END IF


C.........  Prompt for PNTS

        PNAME = PROMPTMFILE(
     &          'Enter logical name for POINT SOURCE input file >> ',
     &          FSREAD3, 'PNTS', 'GETRECS' )

        IF ( .NOT. DESC3( PNAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Could not get description of file "'
     &                  // PNAME( 1:TRIMLEN( PNAME ) ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  PNTS file:', NROWS3D,
     &      'program NPSRC=', NPSRC
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ENDIF

C.........  If file is an EMS-95 file, then prompt for PSRC 
        IF( UPNAM3D .EQ. 'EMSPOINT' ) THEN
            SDEV = PROMPTFFILE(
     &               'Enter name for PSRC file',
     &               .TRUE., .TRUE., 'PSRC', 'GETRECS' )
            EMSFLAG = .TRUE.
        ELSE
            EMSFLAG = .FALSE. 
        ENDIF

C.........  Prompt for PTMP

        TNAME = PROMPTMFILE(
     &          'Enter logical name for POINT SOURCE file',
     &          FSREAD3, 'PTMP', 'GETRECS' )

C.........  Get header information from PTMP

        IF ( .NOT. DESC3( TNAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Could not get description of file "'
     &                  // TNAME( 1:TRIMLEN( TNAME ) ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  PTMP file:', NROWS3D,
     &      'program NPSRC=', NPSRC
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ELSE
            JDATE = SDATE3D
            JTIME = STIME3D
            TSTEP = TSTEP3D
        ENDIF

C.........  Prompt for starting date

        JDATE  = GETNUM( 0, 9999999, JDATE,
     &                  'Enter starting date (YYYYDDD)' )
        SDATE  = JDATE

C.........  Prompt for starting time

        JTIME  = GETNUM( 0, 235959, JTIME,
     &                  'Enter starting time (HHMMSS)' )
        STIME  = JTIME

C.........  Prompt for time number of time steps

        NSTEPS = MXREC3D
        MESG   = 'Enter number of time steps'
        NSTEPS = GETNUM( 1, MXTSTEP, NSTEPS, MESG )
        
C.........  Prompt for PLAY

        LNAME = PROMPTMFILE( 
     &          'Enter logical name for LAYER FRACTIONS file ' //
     &          'or "NONE"', FSREAD3, 'PLAY', 'GETRECS' )

        LAYRFLAG = ( LNAME .NE. 'NONE' ) 

C.........  Prompt for SMAT

        SNAME = PROMPTMFILE(
     &          'Enter logical name for SPECIATION MATRIX file ' //
     &          'or "NONE"', FSREAD3, 'PSMAT', 'GETRECS' )

        SPCFLAG = ( SNAME .NE. 'NONE' )

C.........  Prompt for GMAT

        GNAME = PROMPTMFILE(
     &          'Enter logical name for GRID MATRIX file',
     &          FSREAD3, 'PGMAT', 'GETRECS' )

C.........  Prompt for CMAT

        CNAME = PROMPTMFILE(
     &          'Enter logical name for CONTROL MATRIX file ' //
     &          'or "NONE"', FSREAD3, 'PCMAT', 'GETRECS' )

        CTRLFLAG = ( CNAME .NE. 'NONE' )

C.........  Prompt for Cell or Source or Record number
        METHOD = 1
        METHOD = GETMENU( NMETHOD , METHOD, 
     &                  'Choose source sub-selection method', MENULST )

C.........  If Cell, Prompt for x-dir and y-dir cell numbers

        STKMTHD = 0

        IF( METHOD .EQ. 1 ) THEN

            XN  = 1
            XN  = GETNUM( 1, NCOLS, XN,
     &                    'Enter X-dir cell number' )
            YN  = 1
            YN  = GETNUM( 1, NROWS, YN,
     &                      'Enter Y-dir cell number' )

C.........  Prompt for FIP, SIC, SCC, PLANT, STACK
        ELSEIF( METHOD .EQ. 2 ) THEN

            FIP   = 0
            FIP   = GETNUM( 0, 99999, FIP, 'Enter FIP' )

            SIC   = 0
            SIC   = GETNUM( 0, 9999, SIC, 'Enter SIC' )
           
            SCC   = 0
            SCC   = GETNUM( 0, 99999999, SCC, 'Enter SCC' )

            PLANT = 0
            PLANT = GETNUM( 0, 9999999, PLANT, 'Enter PLANT' )
           
            STACK = 0
            STACK = GETNUM( 0, 9999999, STACK, 'Enter STACK' )

C.........  Prompt for source ID

        ELSEIF( METHOD .EQ. 3 ) THEN

            RN  = 1
            RN  = GETNUM( 1, NPSRC, RN,
     &                    'Enter record number' )

C.........  Prompt for stack parameters

        ELSEIF( METHOD .EQ. 4 ) THEN

            STKMTHD = 1
            STKMTHD = GETMENU( NSTKMTHD, STKMTHD,
     &             'Choose stack feature for ranking sources', STKMENU ) 

C.........  Prompt for emissions based selection

        ELSEIF( METHOD .EQ. 5 ) THEN

            EMMTHD = 1
            EMMTHD = GETMENU( NIPOL, EMMTHD,
     &         'Choose inventory pollutant for ranking sources', EINAM )

        ENDIF

C.........  Check layer fractions file dimensions, if needed
        IF( LAYRFLAG ) THEN
            IF ( .NOT. DESC3( LNAME ) ) THEN
                CALL M3EXIT( 'GETRECS', 0, 0,
     &                      'Could not get description of file "'
     &                      // LNAME( 1:TRIMLEN( LNAME ) ) // '"', 2 )

            ELSEIF( NROWS3D .NE. NPSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  PLAY file:', NROWS3D,
     &          'program NPSRC=', NPSRC
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            ENDIF
        ENDIF

C.........  Get header information from PGMAT

        IF ( .NOT. DESC3( GNAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Could not get description of file "'
     &                  // GNAME( 1:TRIMLEN( GNAME ) ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NGRID ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  PGMAT file:', NROWS3D,
     &      'program NGRID=', NGRID
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ELSEIF( NCOLS3D .NE. NPSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  PGMAT file:', NCOLS3D,
     &      'program NPSRC=', NPSRC
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        END IF

C.........  Read ALL variables in PNTS

        IF ( .NOT. READ3( PNAME, 'ALL', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                   'Error reading from file' // PNAME, 2)
        END IF

C.........  Store parameters
        DX  = SNGL( XCELL3D )
        DY  = SNGL( YCELL3D )
        DDX = 1.0 / DX
        DDY = 1.0 / DY
        X0  = XORIG3D
        Y0  = YORIG3D

C.........  Create unsorted index and preprocess stack parms
        IF( METHOD .EQ. 4 .OR. METHOD .EQ. 5 ) THEN

            MEANVAL = 0.
            MEANCNT = 0
            DO 9 I = 1, NPSRC
                INDEXA( I ) = I

                STKFL( I ) = 0.25 * PI * STKVE( I ) *
     &                                   STKDM( I ) * STKDM( I )

                IF( METHOD .EQ. 5 ) THEN
                    CRITVAL( I ) = EMISV( I, EMMTHD )

                ELSEIF( STKMTHD .EQ. 1 ) THEN
                    CRITVAL( I ) = STKHT( I )

                ELSEIF( STKMTHD .EQ. 2 ) THEN
                    CRITVAL( I ) = STKDM( I )

                ELSEIF( STKMTHD .EQ. 3 ) THEN
                    CRITVAL( I ) = STKTK( I )

                ELSEIF( STKMTHD .EQ. 4 ) THEN
                    CRITVAL( I ) = STKVE( I )

                ELSEIF( STKMTHD .EQ. 5 ) THEN
                    CRITVAL( I ) = STKFL( I )

                ENDIF

                IF( CRITVAL( I ) .GT. 0. ) THEN
                    MEANVAL = MEANVAL + CRITVAL( I )
                    MEANCNT = MEANCNT + 1
                ENDIF

9           CONTINUE  

            MEANVAL = MEANVAL / REAL( MEANCNT )

        ELSE

            DO 10 I = 1, NPSRC
                INDEXA( I ) = I

                STKFL( I ) = 0.25 * PI * STKVE( I ) *
     &                                   STKDM( I ) * STKDM( I )

 10         CONTINUE

        ENDIF

C.........  Calculate cell location if it cell method selected
C.........  and set starting and ending sources
        IF( METHOD .EQ. 1 ) THEN

            XMIN = X0 + ( XN-1 ) * DX
            XMAX = XMIN + DX
            YMIN = Y0 + ( YN-1 ) * DY
            YMAX = YMIN + DY

            SRTSRC = 1
            ENDSRC = NPSRC
            INCSRC = 1

C.........  Set starting and ending sources for source-based methods

        ELSEIF( METHOD .EQ. 2 ) THEN

            SRTSRC = 1
            ENDSRC = NPSRC
            INCSRC = 1

        ELSEIF( METHOD .EQ. 3 ) THEN   ! for single source

            SRTSRC = 1
            ENDSRC = 1
            INCSRC = 1

C.........  If ranked by source characteristic...

        ELSEIF( METHOD .EQ. 4 .OR. METHOD .EQ. 5 ) THEN

C.............  Retrieve source sub-selection characteristics
            RANGE = 1
            RANGE = GETMENU( NRANGE, RANGE,
     &            'Choose range for reporting ranked sources', RNGMENU )

            NRPRT = 50
            NRPRT = GETNUM( 1, NPSRC, NRPRT, 
     &                      'Enter number of sources to report' )

C.............  Sort based on selected category
            CALL SORTR1( NPSRC, INDEXA, CRITVAL )

C.............  Set starting and ending sources and for medium output list, 
C               must first calculate which sources is closest to mean value

            MINDIFF = 1.0E36
            MINVAL  = 1.0E36
            IF( RANGE .EQ. 1 ) THEN      ! Highest X sources

                SRTSRC = NPSRC
                ENDSRC = NPSRC - NRPRT + 1
                INCSRC = -1

            ELSEIF( RANGE .EQ. 2 ) THEN  ! Middle X sources

                DO 11 I = 1, NPSRC

                    J = INDEXA( I )
                    THISDIFF = ABS( MEANVAL - CRITVAL( J ) )
                    IF( THISDIFF .LT. MINDIFF ) THEN
                        MEANSRC = I
                        MINDIFF = THISDIFF
                    ENDIF

11              CONTINUE

                SRTSRC = MEANSRC - NRPRT / 2 + 1
                ENDSRC = MEANSRC + NRPRT / 2 + MOD( NRPRT,2 )
                INCSRC = 1
            
            ELSEIF( RANGE .EQ. 3 ) THEN ! Lowest X sources

                DO 12 I = 1, NPSRC

                    J = INDEXA( I )
                    IF( CRITVAL( J ) .GT. 0      .AND. 
     &                  CRITVAL( J ) .LT. MINVAL       ) THEN
                        MINSRC = I
                        MINVAL = CRITVAL( J )
                    ENDIF

12              CONTINUE

                SRTSRC = MINSRC
                ENDSRC = MINSRC + NRPRT - 1
                INCSRC = 1

            ENDIF

        ENDIF

C.........  Create pointer based on cell or on source parameters

        ICNT = 0
        DO 13 I = SRTSRC, ENDSRC, INCSRC

            J = INDEXA( I )

            IF( ICNT .GE. MXEOUT ) THEN   ! Prevent overflow, but continue count

                ICNT = ICNT + 1

            ELSEIF( METHOD .EQ. 1 ) THEN

                IF( XLOCA( J ) .GE. XMIN .AND.
     &              XLOCA( J ) .LE. XMAX .AND.
     &              YLOCA( J ) .GE. YMIN .AND.
     &              YLOCA( J ) .LE. YMAX       ) THEN

                    ICNT = ICNT + 1
                    EINX( ICNT ) = J

                ENDIF
            
            ELSEIF( METHOD .EQ. 2 ) THEN
                IF( MOD( FIP,1000 ) .EQ. 0 ) THEN

                    SID = ( IFIP( J ) / 1000 ) * 1000 

                    IF( (FIP   .EQ. SID         .OR. FIP   .EQ. 0) .AND.
     &                  (SIC   .EQ. ISIC  ( J ) .OR. SIC   .EQ. 0) .AND.
     &                  (SCC   .EQ. ISCC  ( J ) .OR. SCC   .EQ. 0) .AND.
     &                  (PLANT .EQ. IPLANT( J ) .OR. PLANT .EQ. 0) .AND.
     &                  (STACK .EQ. ISTACK( J ) .OR. STACK .EQ. 0)) THEN
                        ICNT = ICNT + 1
                        EINX( ICNT ) = J
                    ENDIF

                ELSEIF
     &            ( ( FIP   .EQ. IFIP  ( J ) .OR. FIP   .EQ. 0 ) .AND.
     &              ( SIC   .EQ. ISIC  ( J ) .OR. SIC   .EQ. 0 ) .AND.
     &              ( SCC   .EQ. ISCC  ( J ) .OR. SCC   .EQ. 0 ) .AND.
     &              ( PLANT .EQ. IPLANT( J ) .OR. PLANT .EQ. 0 ) .AND.
     &              ( STACK .EQ. ISTACK( J ) .OR. STACK .EQ. 0 ) ) THEN
                    ICNT = ICNT + 1
                    EINX( ICNT ) = J
 
                ENDIF
                
            ELSEIF( METHOD .EQ. 3 ) THEN

                ICNT = 1
                EINX( ICNT ) = RN    

            ELSEIF( METHOD .EQ. 4 .OR. METHOD .EQ. 5 ) THEN

                ICNT = ICNT + 1
                EINX( ICNT ) = J

            ENDIF

   13   CONTINUE

        NOUT = ICNT
        IF( NOUT .GT. MXEOUT ) THEN
            WRITE( MESG, 94010 ) 'Output sources found =', NOUT,
     &                           'but maximum (MXEOUT) =', MXEOUT,
     &                           '. Reset MXEOUT and try again.'
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ENDIF

C.........  Read SMAT

        IF( SPCFLAG ) THEN

            IF( .NOT. READ3( SNAME, 'ALL', 1, 0, 0, SS ) ) THEN
 
            CALL M3ERR( 'GETRECS', 0, 0,
     &           'Could not read SPECIATION matrix from file "'//
     &           SNAME( 1 : TRIMLEN( SNAME ) ) // '".',
     &           .TRUE. )
 
            END IF      !  if read3() failed for speciation matrix

C...........   Build description of how to use speciation matrix:
 
            IF ( .NOT. DESC3( SNAME ) ) THEN
                CALL M3EXIT( 'GETRECS', 0, 0,
     &                       'Could not get description of file ' //
     &                        SNAME, .TRUE. )
            ELSEIF( NROWS3D .NE. NPSRC ) THEN
                WRITE( MESG, 94010 )
     &          'Dimension mismatch.  PSMAT file:', NROWS3D,
     &          'program NPSRC=', NPSRC
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

            END IF
 
            DO  14  V = 1, NIPOL
                NSPCOUT( V ) = 0
14          CONTINUE
 
            DO  15  V = 1, NMPOL
 
                K = INDEX( VNAME3D( V ), '_to_' )
                IF( K .LT. 1 ) THEN
                    CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Bad conversion factor name "' //
     &                  VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) ) //
     &                  '" in ' // SNAME, 2 )
                END IF
                INVNAME( V ) = VNAME3D( V )( 1:K-1 )
 
                I = INDEX1( INVNAME( V ), NIPOL, EINAM )
                IF ( I .LT. 0 ) THEN
                    CALL M3ERR( 'GETRECS', 0, 0,
     &                  'Bad input species "' //
     &                  INVNAME( V )( 1:TRIMLEN( INVNAME( V ) ) ) //
     &                  '" in ' // SNAME,
     &                  .TRUE. )
                END IF
 
                SCRBUF = VNAME3D( V )( K+4:16 )
                J = INDEX1( SCRBUF, NMPOL, EMNAM )
                IF ( J .LT. 0 ) THEN
                    CALL M3ERR( 'GETRECS', 0, 0,
     &                          'Bad output species "' //
     &                          SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &                          '" in ' // SNAME, 2 )
                END IF
 
                F = NSPCOUT( I ) + 1
                SPCREF ( I,F ) = J
                NSPCOUT( I )   = F
 
15          CONTINUE

        ENDIF   ! whether speciation matrix read in or not

C.........  Read GMAT

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NS ) ) THEN
 
            CALL M3ERR( 'GETRECS', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".',
     &           .TRUE. )
 
        END IF      !  if read3() failed for gridding matrix

C......... Check to make sure that the number records in cell matches GMAT info

        IF( METHOD .EQ. 1 ) THEN

            CELLNUM = ( YN-1 )*NCOLS + XN
            NINCELL= NS( CELLNUM )
            IF( NOUT .NE. NINCELL ) THEN

                WRITE( * , 92100 ) 
     &          'Warning: the number of counted records for',
     &          'cell (', XN, ',', YN,') is ', NOUT,
     &          ' but the number in the gridding matrix is ', NINCELL
            ENDIF

C.............  Count up sources up to cell of interest
            SCNT = 1
            DO 17 S = 1, CELLNUM-1
                SCNT = SCNT + NS( S )
   17       CONTINUE
C.............  Create sorted source list for cell
            DO 18 S = 1, NS( CELLNUM )
                ISLOC( S ) = S
   18       CONTINUE

            CALL SORTI1( 
     &           NS( CELLNUM ), ISLOC, IS( SCNT ) ) 

            DO 19 S = 1, NS( CELLNUM )
                ISLOC( S )= IS( SCNT + ISLOC( S )-1 )
   19       CONTINUE

        ENDIF

C.........  Read CMAT

        IF( CTRLFLAG ) THEN

            DO  22 S = 1, NIPOL
 
                SCRBUF = 'CONTROL_' // EINAM( S )( 1:8 )
                IF ( .NOT. READ3( CNAME, SCRBUF, 1, 0, 0, CC( 1,S ) )
     &             ) THEN
 
                    CALL M3ERR( 'GETRECS', 0, 0,
     &                   'Could not read "' // SCRBUF( 1:12 ) //
     &                   '" from CONTROL FILE "' //
     &                   CNAME( 1 : TRIMLEN( CNAME ) ) // '".',
     &                   .TRUE. )
 
                END IF      !  if read3() failed for control matrix

22          CONTINUE

        END IF          !  if controlflag

C.........  Populate plant and stack character arrays (different if have PSRC)
        DO 24 S = 1, NPSRC

            IF( EMSFLAG ) THEN
                READ ( SDEV, 93005, IOSTAT=IOS ) CPLT, CSTK, CDESC

                IF( IOS .NE. 0 ) THEN
                    WRITE( MESG,94010 )
     &                  'Unit number', SDEV,
     &                  'I/O Status ', IOS,
     &                  'Line number', S
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Error reading PSRC file.', 2 )

                ELSE
                    CPLANT( S ) = CPLT
                    CSTACK( S ) = CSTK
                    CPDESC( S ) = CDESC
                ENDIF
                
            ELSE

                WRITE( CPLT, '(I15)' ) IPLANT( S )
                WRITE( CSTK, '(I15)' ) ISTACK( S )

                CPLANT( S ) = CPLT
                CSTACK( S ) = CSTK
                CPDESC( S ) = ' '

            ENDIF

24      CONTINUE

C.........  Loop over time

        DO 66 T = 1, NSTEPS

C.............  Read layer fractions for given hour
            IF( LAYRFLAG ) THEN
                IF( READ3( LNAME, 'LFRAC', ALLAYS3,
     &                     JDATE, JTIME, LFRAC   ) ) THEN

                    DO 27 L = 1, EMLAYS
                        DO 25 I = 1, NOUT
    
                            OFRAC( I, L, T ) = LFRAC( EINX( I ), L )
    
25                      CONTINUE
27                  CONTINUE

                ELSE ! Read failed

                    CALL M3EXIT( 'GETRECS', JDATE, JTIME,
     &                       'Could not read "LFRAC" from '// LNAME, 2 )
 
                ENDIF
            ENDIF

C.............  Loop over inventory species
            DO 55 V = 1, NIPOL    

C.................  Read PTMP for given hour and species
                IF( READ3( TNAME, EINAM( V ), ALLAYS3,
     &                     JDATE, JTIME, EMIST ) ) THEN

C.....................  Loop over list of sources

                    DO 44 I = 1, NOUT

                        IF( METHOD .NE. 1 ) THEN  ! not cell selection

                            XX = DDX * ( XLOCA( EINX( I ) ) -X0 )
                            YY = DDY * ( YLOCA( EINX( I ) ) -Y0 )

                            IF( XX .GT. 0 .AND. YY .GT. 0 ) THEN
                                XNUM( I ) = INT( XX ) + 1
                                YNUM( I ) = INT( YY ) + 1
                            ELSE
                                XNUM( I ) = INT( XX )
                                YNUM( I ) = INT( YY )
                            ENDIF

                        ELSE

                            XNUM( I ) = XN
                            YNUM( I ) = YN

                        ENDIF

                        CELLNUM = ( YNUM( I )-1 )*NCOLS + XNUM( I )

                        IF( XNUM( I ) .LE. 0     .OR. 
     &                      XNUM( I ) .GT. NCOLS .OR.
     &                      YNUM( I ) .LE. 0     .OR.
     &                      YNUM( I ) .GT. NROWS      ) THEN

                            GSTATE( I ) = 'OUTSIDE'
                            GO TO 44
                        ENDIF

                        IF( METHOD .NE. 1 ) THEN  ! not cell selection

C.............................  Count up sources up to cell of interest
                            SCNT = 1
                            DO 29 S = 1, CELLNUM - 1
                                SCNT = SCNT + NS( S )
   29                       CONTINUE

C.............................  Create sorted source list for cell
                            DO 30 S = 1, NS( CELLNUM )
                                ISLOC( S ) = S
   30                       CONTINUE

                            CALL SORTI1( 
     &                           NS( CELLNUM ), ISLOC, IS( SCNT ) ) 

                            DO 31 S = 1, NS( CELLNUM )
                                ISLOC( S )= IS( SCNT + ISLOC( S )-1 )
   31                       CONTINUE

                        ENDIF

                        F = FIND1( EINX( I ),NS( CELLNUM ), ISLOC )
                                 
                        IF( F .GT. 0 ) THEN
                            GSTATE( I ) = 'FOUND'     
                        ELSE
                            GSTATE( I ) = 'UNFOUND'
                        ENDIF

C.........................  Loop over model species per inventory species

                        IF( SPCFLAG ) THEN
                            DO 33 J = 1,  NSPCOUT( V )

                                K = SPCREF( V,J )

C.............................  Multiply PTMP, SMAT, and CMAT for source

                                SFAC = SS( EINX( I ), K ) 
                                IF( CTRLFLAG ) THEN
                                    CFAC = CC( EINX( I ), V )
                                ELSE
                                    CFAC = 1.0
                                ENDIF
                             
                                EOUT( I, K, T ) = 
     &                                EMIST( EINX( I ) ) * SFAC * CFAC

                                EV =  EMISV( EINX( I ), V )

   33                       CONTINUE

                        ENDIF

                        ETMP( I, V, T ) = EMIST( EINX( I ) )

   44               CONTINUE

                ELSE ! Read failed

                    CALL M3EXIT( 'GETRECS', JDATE, JTIME,
     &                           'Could not read "' //
     &                           EINAM( V )( 1:TRIMLEN( EINAM( K ))) //
     &                           '" from ' // TNAME, 2 )
 
                ENDIF
   55       CONTINUE

            CALL NEXTIME( JDATE, JTIME, TSTEP )

   66   CONTINUE

C.........  Prompt for output temporal x-ref file name

        ODEV = PROMPTFFILE(
     &  'Enter name for ASCII OUTPUT file',
     &           .FALSE., .TRUE., 'LISTOUT', 'GETRECS' )
 
        IF ( ODEV .LE. 0 ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                   'Could not open ASCII OUTPUT file', 2 )
        END IF

        TONSUNIT = '[tons]  '
        MOLEUNIT = '[gm mol]'

C.........  Loop over list of sources and output details

        DO 88 I = 1, NOUT
            S = EINX( I )

            L1 = LBLANK( CPLANT( S ) ) + 1
            L2 = LBLANK( CSTACK( S ) ) + 1
            CPLT = CPLANT( S )( L1:TRIMLEN( CPLANT( S ) ) )
            CSTK = CSTACK( S )( L2:TRIMLEN( CSTACK( S ) ) )
            L1 = TRIMLEN( CPLT )
            L2 = TRIMLEN( CSTK )

C.............  Write out all time INdependant data
            WRITE( ODEV, 93000 ) 
     &      I, S, CPDESC( S ),
     &      IFIP( S ), ISIC( S ), ISCC( S ), CPLT, CSTK,
     &      INVYR( S ), TZONES( S ), TPFLAG( S ), 
     &      STKHT( S ), STKDM( S ), STKTK( S ), STKVE( S ), STKFL( S ),
     &      XNUM( I ), YNUM( I ), GSTATE( I ), XLOCA( S ), YLOCA( S )

C.............  Write out control coefficients

            IF( CTRLFLAG )
     &      WRITE( ODEV, 93010 ) 
     &           ( EINAM( V )( 1:TRIMLEN( EINAM( V ) ) ), 
     &             CC( S, V ), V = 1, NIPOL )

C.............  Write out time independant inventory emissions
            WRITE( ODEV, 93020 ) 
     &           ( EINAM( V )( 1:TRIMLEN( EINAM( V ) ) ),
     &             EMISV( EINX( I ), V ), V = 1, NIPOL )

C.............  Create pointer for model species
            IF( SPCFLAG ) THEN
                ICNT = 0
                DO 70 V = 1, NIPOL
                     DO 69 J = 1,  NSPCOUT( V )

                         K = SPCREF( V,J )
                         ICNT = ICNT + 1
                         SPNT( ICNT ) = K

   69                CONTINUE
   70           CONTINUE
            ENDIF

C.............  Write out speciation coefficients
            IF( SPCFLAG ) THEN
                WRITE( ODEV, 93030 ) (  EMNAM( SPNT( J ) ),
     &                                  SS( S, SPNT( J ) ), J=1, NMPOL )
            ENDIF

C.............  Write out header for temporal emissions
            IF( SPCFLAG ) THEN
                WRITE( ODEV, 93040 ) MMDDYY( SDATE ),
     &                           ( EINAM( V )        , V = 1, NIPOL ), 
     &                           ( EMNAM( SPNT( J ) ), J = 1, NMPOL ),
     &                           ( TONSUNIT          , K = 1, NIPOL ),
     &                           ( MOLEUNIT          , L = 1, NMPOL )
            ELSE
                WRITE( ODEV, 93042 ) MMDDYY( SDATE ),
     &                           ( EINAM( V )        , V = 1, NIPOL ), 
     &                           ( TONSUNIT          , K = 1, NIPOL )
            ENDIF

C.............  Write out time dependent emissions
            JDATE = SDATE
            JTIME = STIME
            LDATE = -9
            DO 77 T = 1, NSTEPS

                IF( SPCFLAG ) THEN
                    WRITE( ODEV, 93050 ) 
     &                     JDATE-SDATE+1, JTIME/10000, 
     &                     ( ETMP( I, V, T ), V=1,NIPOL ),
     &                     ( EOUT( I, SPNT( J ), T ), J=1, NMPOL )
                ELSE
                    WRITE( ODEV, 93050 )
     &                     JDATE-SDATE+1, JTIME/10000, 
     &                     ( ETMP( I, V, T ), V=1,NIPOL )
                ENDIF

                IF( JDATE .NE. LDATE ) THEN

                    DO 71 V = 1, NIPOL
                        EISUM( V ) = ETMP( I, V, T )
71                  CONTINUE

                    IF ( SPCFLAG ) THEN
                        DO 72 J = 1, NMPOL
                            EMSUM( J ) = EOUT( I, SPNT( J ), T )
72                      CONTINUE
                    ENDIF

                    LDATE = JDATE

                ELSE

                    DO 75 V = 1, NIPOL
                        EISUM( V ) = EISUM( V ) + ETMP( I, V, T )
75                  CONTINUE

                    IF ( SPCFLAG ) THEN
                        DO 76 J = 1, NMPOL
                          EMSUM( J )= EMSUM( J )+EOUT( I, SPNT( J ), T )
76                      CONTINUE
                    ENDIF

                ENDIF

                CALL NEXTIME( JDATE, JTIME, TSTEP )

                IF( JDATE .NE. LDATE .OR. T .EQ. NSTEPS ) THEN
                    IF( SPCFLAG ) THEN
                        WRITE( ODEV,93055 ) ( EISUM( V ), V=1, NIPOL ),
     &                                      ( EMSUM( J ), J=1, NMPOL )
                    ELSE
                        WRITE( ODEV,93055 ) ( EISUM( V ), V=1, NIPOL )
                    ENDIF
                ENDIF

77          CONTINUE

            IF( LAYRFLAG .AND. GSTATE( I ) .NE. 'OUTSIDE' ) THEN

C.................  Write out header for layer fractions
                WRITE( ODEV, 93060 ) MMDDYY( SDATE ), 
     &                               ( L, L = 1, EMLAYS )

C.................  Write out time dependent layer fractions
C.................  Format will put blanks where there are zeroes
                JDATE = SDATE
                JTIME = STIME
                DO 79 T = 1, NSTEPS

                    WRITE( ODEV, 93070 ) 
     &                     JDATE-SDATE+1, JTIME/10000

                    DO 78 L = 1, EMLAYS - 1

                        VAL = OFRAC( I, L, T ) 
                        IF( VAL .GT. 0. ) THEN
                            WRITE( ODEV, 93080 ) VAL
                        ELSE
                            WRITE( ODEV, 93090 ) 
                        ENDIF

78                  CONTINUE

                    VAL = OFRAC( I, EMLAYS, T )
                    IF( VAL .GT. 0. ) THEN
                        WRITE( ODEV, 93100 ) VAL
                    ELSE 
                        WRITE( ODEV, 93110 )
                    ENDIF

                    CALL NEXTIME( JDATE, JTIME, TSTEP )

79              CONTINUE

            ENDIF

            WRITE( ODEV, '(A)' ) ' '

   88   CONTINUE 

        CALL M3EXIT( 'GETRECS', 0, 0,
     &               'Normal completion  of PROGRAM GETRECS', 0 )

C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Informational (LOG) message formats... 92xxx
 
92000   FORMAT( 5X, A )

92100   FORMAT( 10X, A, / 10X, A, I3, A, I3, A, I5 / 10X, A, I5)
 
C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( I5, ' Src ID:', I7, 1X, A, / ,
     &          5X, ' FIP:', I5.5, ' SIC:', I4.4, ' SCC:', I8.8,
     &              ' PLT:', A<L1>, ' STK:', A<L2>, ' YR:' , I4,
     &              ' ZON:', I2.2, ' TPF:', I1, /,
     &          5X, ' Stack parms...   H[m]:', F7.2, '   D[m]:', F7.2,
     &              '   T[K]:', F7.1, '   V[m/s]:', F7.1, 
     &              '   F[m^3/s]:', F10.1, /,
     &          5X, ' Gridding...      Cell:(', I3,',', I3, 
     &              ')   Status:', A7, 
     &              '   XLOC:', E14.7, '   YLOC:', E14.7 )

93005   FORMAT( 29X, A15, 1X, A12, 14X, A40 )

93010   FORMAT( 5X, ' Control facs...  ',
     &              30( 1X, A, ': ', E11.5 )   )

93020   FORMAT( 5X, ' Inv emissions...',
     &              <NIPOL>( 1X, A, ': ', E11.5 ), '[tons/yr]' )

93030   FORMAT( 5X, ' Speciation coefficients [gm mol/tons]...', /,
     &          5( 9X, 7( 1X, A4, ': ', E11.5 ), / ) )

93040   FORMAT( 5X, ' Hourly emissions starting ', A14, '...' /,
     &          5X, ' Dy Hr', <NIPOL+NMPOL>( 1X, A8 ) , /,
     &          11X, <NIPOL+NMPOL>( 1X, A8 ) )

93042   FORMAT( 5X, ' Hourly emissions starting ', A14, '...' /,
     &          5X, ' Dy Hr', <NIPOL>( 1X, A8 ) , /,
     &          11X, <NIPOL>( 1X, A8 ) )

93050   FORMAT( 6X, I2.2, 1X, I2.2, <NIPOL+NMPOL>( 1X, E8.3 ) )

93055   FORMAT( 5X, ' Total', <NIPOL+NMPOL>( 1X, E8.3 ) )

93060   FORMAT( 5X, ' Layer fractions starting ', A14, '...' /,
     &          5X, ' Dy Hr', <EMLAYS>( 1X, ' Lyr(', I2.2, ')' ) )

93070   FORMAT( 6X, I2.2, 1X, I2.2, $ )

93080   FORMAT( 1X, E8.3, $ )

93090   FORMAT( 9X, $ )

93100   FORMAT( 1X, E8.3 )

93110   FORMAT( 9X )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 1X ) )

        END
