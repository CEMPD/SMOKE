 
        PROGRAM GETRECS

C***********************************************************************
C  program body starts at line 272
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
C       Major update to include area and mobile sources 10/98
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
C***********************************************************************
 
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMDIMS3.EXT'   !  emissions dimensioning params
        INCLUDE 'CHDGEN3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal parameters
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
        INTEGER       SECSDIFF
        INTEGER       TRIMLEN
        REAL          YR2DAY
 
        EXTERNAL      FIND1, GETMENU, GETNUM, GETYN, INDEX1, LBLANK,
     &                MMDDYY, PROMPTFFILE, PROMPTMFILE, SECSDIF,
     &                TRIMLEN, YR2DAY

C.........  Define PNTS variables
        INTEGER     INDEXA( MXESRC )   !  index for sorting, if needed
        INTEGER     IFIP  ( MXESRC )   !  source FIPS (county) ID
        INTEGER     VAR2  ( MXESRC )   !  PT: SCC, AR: ASC7
        INTEGER     VAR3  ( MXESRC )   !  PT: SIC, AR: ASC3 
        INTEGER     VAR4  ( MXESRC )   !  PT: Pltid (dim w/ MXESRC on purpose)
        INTEGER     VAR5  ( MXESRC )   !  PT: Stkid (dim w/ MXESRC on purpose)
        INTEGER     TZONES( MXESRC )   !  time zones
        INTEGER     TPFLAG( MXESRC )   !  applicability of temprl profile types
        INTEGER     INVYR ( MXESRC )   !  Inventory year
        REAL        XLOCA ( MXESRC )   !  stack UTM X-location (m)
        REAL        YLOCA ( MXESRC )   !  stack UTM Y-location (m)
        REAL        XLOC1 ( MXESRC )   !  link start UTM X-location (m)
        REAL        YLOC1 ( MXESRC )   !  link start UTM Y-location (m)
        REAL        XLOC2 ( MXESRC )   !  link end UTM X-location (m)
        REAL        YLOC2 ( MXESRC )   !  link end UTM Y-location (m)
        REAL        STKHT ( MXPSRC )   !  stack height   (m)
        REAL        STKDM ( MXPSRC )   !  stack diameter (m)
        REAL        STKTK ( MXPSRC )   !  exhaust temperature (deg K)
        REAL        STKVE ( MXPSRC )   !  exhaust velocity    (m/s)
        REAL        EMISV ( MXESRC, MXIPOL )  !  emissions values (tons/yr)
        REAL        VMT   ( MXESRC )   !  mobile VMT from inventory 
        REAL        VMTDIS( NMSRC, NVTYPE ) ! VMT per vehicle type

        EQUIVALENCE ( XLOCA, XLOC1 )
        EQUIVALENCE ( YLOCA, YLOC1 )
        EQUIVALENCE ( STKHT, XLOC2 )
        EQUIVALENCE ( STKDM, YLOC2 )
        EQUIVALENCE ( STKTK, VMT   )

        CHARACTER*15   CPLANT( MXPSRC )
        CHARACTER*15   CSTACK( MXPSRC )
        CHARACTER*40   CPDESC( MXPSRC )

C.........  Define temporary source variable

        REAL           CRITVAL( MXESRC ) !  Criterion for sorting output
        REAL           STKFL  ( MXPSRC ) !  Stack flow

C.........  Define hourly emissions variables
        REAL           EMIST ( MXESRC )  !  hourly emissions values (tons/yr) 
        REAL           TMAT  ( NMSRC, NVTYPE, MXTMV )  !  hourly mobile values
        CHARACTER*16   TPRLNAM( MXTMV ) !  temporal emis variable names

        EQUIVALENCE( TMAT, EMIST )

C.........  Define speciation matrix variables

        REAL            SS( MXESRC, MXMPOL )
        REAL            SM( NMSRC*NVTYPE, MXSMV )
        INTEGER         NSPCOUT( MXIPOL*MXEPRC )! no. smat vars per var  
        INTEGER         SPCREF ( MXTMV, MXSMV ) ! which species per tmprl var
        INTEGER         SVARREF( MXTMV, MXSMV ) ! which smat var per tmprl var
        INTEGER         SMREFM ( MXSMV )        ! species indx for each SMAT var
        CHARACTER*16    SPCNAM ( MXSMV )        !  smat variable names

        EQUIVALENCE( SM, SS )

C.........  Gridding matrix variables

        INTEGER         NA( NGRID )   ! area or mobile source gridding matrix
        INTEGER         IA( NMATX )
        REAL            MA( NMATX )

        COMMON  / GAMAT / NA, IA, MA

        INTEGER         NP( NGRID )   ! point source gridding matrix
        INTEGER         IP( MXPSRC )
        REAL            MP( MXPSRC )  ! dummy for rdgridm call
 
        EQUIVALENCE( NA, NP )  ! Will never need to use both
        EQUIVALENCE( IA, IP )  ! Will never need to use both
        EQUIVALENCE( MA, MP )

        INTEGER         PG( NGRID )   ! Points to start of IA/IP for each cell

C.........   Ungridding matrix variables (for mobile only)
 
        INTEGER         NU( NMSRC )
        INTEGER         IU( NMATX )
        REAL            CU( NMATX )
 
        COMMON  / UGRIDMAT / NU, IU, CU

        INTEGER         PU( NMSRC )   ! Points to start of IU for each source

C.........  Define control matrix variables

        REAL            CC( MXESRC, MXIPOL )

C.........  Define Layer Fractions

        REAL            LFRAC1L( MXPSRC )

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
        INTEGER         CCNT ( MXEOUT )   ! Count of cells per output src
        INTEGER         EINX ( MXEOUT )   !   sorted index for output sources
        INTEGER         EINXA( MXEOUT )   ! unsorted index for output sources
        INTEGER         EVIDX( MXTMV  )   ! inventory pol index for TMAT vars  
        INTEGER         IINX ( MXEOUT )   ! index of index for output sources
        INTEGER         ISLOC( NMATX  )   ! Cell specific IS values
        INTEGER         NOUT              ! actual number of output emis
        INTEGER         NSCNT( MXMPOL )   ! no. of emis procs for each species 
        INTEGER         XNUM ( MXEOUT, NGRID )   ! X-cell for output sources
        INTEGER         YNUM ( MXEOUT, NGRID )   ! Y-cell for output sources

        REAL            EISUM( MXIPOL )   ! Sum of day's inv pollutants
        REAL            EMSUM( MXMPOL )   ! Sum of day's model species
        REAL            EOUT ( MXEOUT, MXMPOL, MXTSTEP) ! output emissions
        REAL            ETMP ( MXEOUT, MXIPOL, MXTSTEP) ! out tmprl inv emission
        REAL            GCOEF( MXEOUT, NGRID ) ! gridding coef for output srcs
        REAL            OFRAC( MXEOUT, EMLAYS, MXTSTEP)! output layer fractions
        REAL            SFACS( MXMPOL, NVTYPE )! tmp output mb spec factors

        CHARACTER*7     GSTATE ( MXEOUT  ) ! state of gridding matrix for output
        CHARACTER*16    VNAMINV( MXVARS3 ) ! inventory file var names

C.........  Logical names and unit numbers
        INTEGER         ADEV             ! Area-source ASCT descriptions file
        INTEGER         LDEV             ! IO/API initialization unit #
        INTEGER         ODEV             ! number for output file
        INTEGER         SDEV             ! PSRC file
        INTEGER         VDEV             ! VMT mix file

        CHARACTER*16    CNAME   !  logical name for control    matrix input file
        CHARACTER*16    ENAME   !  logical name for inventory         input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    LNAME   !  logical name for layer fractions   input file
        CHARACTER*16    SNAME   !  logical name for speciation matrix input file
        CHARACTER*16    TNAME   !  logical name for temporal point source file
        CHARACTER*16    UNAME   !  logical name for ungridding matrix

C.........  Other variables
        REAL            CFAC, DDX, DDY, DX, DY, SFAC, VAL
        REAL            X0, Y0, XX, YY
        REAL            MEANVAL      ! Mean of criterion output values
        REAL            MINVAL       ! Min non-zero of criterion output values
        REAL            MINDIFF      ! Minimum difference between src and mean
        REAL            RDUM         ! Dummy real var
        REAL            THISDIFF     ! Current difference between src and mean

        INTEGER         I, C, E, F, J, K, L, N, S, T, V   ! pntrs and counters

        INTEGER         RN, SCNT, TSTEP, XN, YN
        INTEGER         EMMTHD       ! method of sorting based on emissions
        INTEGER         ENDSRC       ! Ending source count for output
        INTEGER         FIP, SID     ! Temporary vars
        INTEGER         IDCELL       ! Cell number from X and Y cells
        INTEGER         INCSRC       ! Increment for processing sources
        INTEGER         ICNT         ! counter for output sources
        INTEGER         IOS          ! I/O status
        INTEGER         JDATE, JTIME ! Current date and time
        INTEGER         KEY2, KEY3, KEY4, KEY5 ! tmp src chars for matching
        INTEGER         L1, L2
        INTEGER         LDATE        ! Previous date
        INTEGER         LPARS        ! String length of NAMPARS
        INTEGER         MEANCNT      ! count of sources contibuting to mean val
        INTEGER         MEANSRC      ! source ID that has value closest to mean
        INTEGER         METHOD       ! method of source selection
        INTEGER         MINSRC       ! source ID that has min non-zero value
        INTEGER         NCINVP       ! number of inventory pols in ctl matrix
        INTEGER         NEVAR        ! number of vars in temporal file
        INTEGER         NIPOL        ! number of inventory pollutants
        INTEGER         NIVAR        ! number of inventory file variables
        INTEGER         NMPOL        ! number of model species
        INTEGER         NSTEPS       ! number of time steps
        INTEGER         NRPRT        ! number of sources to report
        INTEGER         NSRCS        ! actual number of sources in inventory
        INTEGER         NVSMAT       ! actual number of variables in SMAT
        INTEGER         OUTCEL       ! tmp number of output cells per source
        INTEGER         PTR          ! tmp index for gridding matrix
        INTEGER         RANGE        ! method of selecting range of sorted srcs
        INTEGER         SDATE, STIME ! Starting date and time
        INTEGER         SRTSRC       ! Starting source count for output
        INTEGER         STKMTHD      ! method of sorting based on stack parms
        INTEGER         TMPSRC       ! tmp count of sources
        INTEGER         UZONE        ! UTM zone

        LOGICAL         CTRLFLAG     ! control file input or not
        LOGICAL         EFLAG        ! error flag
        LOGICAL         EMSFLAG      ! PNTS is EMS-95 origin or not
        LOGICAL         LAYRFLAG     ! layer fractions file input or not
        LOGICAL         SPCFLAG      ! speciation matrix input or not

        CHARACTER*1     CRL             ! category root letter (A/M/P)
        CHARACTER*3     CTLINV( MXIPOL ) ! inventory pollutants in control matrix
        CHARACTER*3     EINAM( MXIPOL ) ! actual list of inventory pollutants
        CHARACTER*4     NAMPARS         ! string for parsing SMAT variable names
        CHARACTER*6     CATEGORY        ! emissions category (AREA/MOBILE/POINT)
        CHARACTER*8     EMNAM( MXMPOL ) ! actual list of model species
        CHARACTER*8     TONSUNIT
        CHARACTER*8     MOLEUNIT
        CHARACTER*15    CPLT, CSTK
        CHARACTER*16    PRVNAM  !  variable name from previous iteration
        CHARACTER*16    SCRBUF  !  scratch buffer
        CHARACTER*40    CDESC
        CHARACTER*80    MENULST( NMETHOD )
        CHARACTER*80    RNGMENU( NRANGE )
        CHARACTER*80    STKMENU( NSTKMTHD )
        CHARACTER*256   MESG    !  scratch message buffer 

        DATA MENULST / 
     &       'Select all sources in specific grid cell',
     &       'Select sources by St/Co FIPS code',
     &       'Select single source with source ID',
     &       'Select sources with X emissions of pollutant Y',
     &       'Select sources by stack parameters' /

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
     &  'keys.  It creates an ASCII file which lists all details about',
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

C.........  Prompt for inventory file

        ENAME = PROMPTMFILE(
     &          'Enter logical name for INVENTORY input file >> ',
     &          FSREAD3, 'PNTS', 'GETRECS' )

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Error reading header from file "'
     &                  // ENAME( 1:TRIMLEN( ENAME ) ) // '"', 2 )

C.........  Check dimensions of input inventory file
        ELSEIF( NROWS3D .GT. MXESRC ) THEN
            WRITE( MESG, 94010 )
     &          'Dimension overflow.  INVENTORY file:', NROWS3D,
     &          'program (MXESRC):', MXESRC
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

C.........  Determine source category
        ELSE
            I = 0
            I = INDEX1( 'ASC7', NVARS3D, VNAME3D )
            IF( I .GT. 0 ) THEN
                CATEGORY = 'AREA'
            ELSE
                I = INDEX1( 'IRCLAS', NVARS3D, VNAME3D )
                IF( I .GT. 0 ) THEN
                    CATEGORY = 'MOBILE'
                ELSE
                    CATEGORY = 'POINT'
                ENDIF
            ENDIF
        ENDIF

C.........  Store some of the contents of the inventory file header
        NSRCS = NROWS3D
        NIVAR = NVARS3D
        DO 66 I = 1, NIVAR
            VNAMINV( I ) = VNAME3D( I )
66      CONTINUE

C.........  For mobile, source number must be exact and must store UTM zone
        IF( CATEGORY .EQ. 'MOBILE' ) THEN
            IF( NSRCS .NE. NMSRC ) THEN
                WRITE( MESG,94010 )
     &              'Dimension overflow.  INVENTORY file:', NSRCS,
     &              'program (NMSRC):', NMSRC
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            ENDIF
            UZONE = REAL( P_ALP3D )

        ENDIF

C.........  Set flag for EMS-95 point sources inventory input
        EMSFLAG = .FALSE.
        IF( UPNAM3D .EQ. 'EMSPOINT' ) EMSFLAG = .TRUE.
                 
C.........  Set first letter for file names based on category
        CRL = CATEGORY( 1:1 )

C.........  Prompt for temporalized emissions

        TNAME = PROMPTMFILE(
     &          'Enter logical name for TEMPORAL EMISSIONS file',
     &          FSREAD3, CRL//'TMP', 'GETRECS' )
        L = TRIMLEN ( TNAME )

C.........  Get header information from temporalized emissions
C.........  Use this for inventory emissions variable names and no.; default 
C.........  start date, time, and duration; and time step
        IF ( .NOT. DESC3( TNAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Could not get description of file "'
     &                  // TNAME( 1:L ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NSRCS ) THEN
            WRITE( MESG, 94010 )
     &         'Source dimension mismatch. '// CRL// 'TMP file:', 
     &         NROWS3D, 'inventory file= ', NSRCS
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ELSEIF( NVARS3D .GT. MXTMV ) THEN
            WRITE( MESG, 94010 )
     &         'Variables dimension mismatch. '// CRL// 'TMP file:', 
     &         NVARS3D, 'program = ', MXTMV
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ELSE
            JDATE = SDATE3D
            JTIME = STIME3D
            TSTEP = TSTEP3D
            NSTEPS= MXREC3D
            NEVAR = NVARS3D  ! Store number of variables in temporal file

C.............  Set pollutant names. Need to handle mobile inputs that have
C.............  names of the form EXH_CO, and create a pointer for mobile
C.............  input variables to pollutant names
            J = 0
            PRVNAM = '-9'
            DO 88 V = 1, NEVAR

                SCRBUF = VNAME3D( V )
                TPRLNAM( V ) = SCRBUF
                L2 = TRIMLEN( SCRBUF )

                IF( CATEGORY .NE. 'MOBILE' ) THEN ! no pollutants in mobile file
                    I = INDEX1( SCRBUF, NIVAR, VNAMINV )
                    L1 = 0
                ELSE
                    I = 1 
                    L1 = INDEX( SCRBUF, '_' )
                END IF

                IF ( I .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'Pollutant "' //
     &                     SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &                     '" in "' // TNAME( 1:L ) // 
     &                     '" is not in inventory file!'
                    CALL M3MSG2( MESG )
                ELSE
                    SCRBUF =  SCRBUF( L1+1:L2 )
                    IF( SCRBUF .NE. PRVNAM ) THEN
                       J = J + 1                 
                       EINAM( J ) = SCRBUF
                       PRVNAM = SCRBUF
                    ENDIF

                    EVIDX( V ) = J

                END IF

88          CONTINUE

        ENDIF

        IF( EFLAG ) THEN
            MESG = 'Inconsistent hourly and inventory files'
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ELSE
            NIPOL = J
        ENDIF

C.........  Prompt for remaining files
        LNAME = 'NONE'
        IF( CATEGORY .EQ. 'POINT' ) THEN
            LNAME = PROMPTMFILE( 
     &              'Enter logical name for LAYER FRACTIONS file ' //
     &              'or "NONE"', FSREAD3, 'PLAY', 'GETRECS' )
        ENDIF
        LAYRFLAG = ( LNAME .NE. 'NONE' ) 

        SNAME = PROMPTMFILE(
     &          'Enter logical name for SPECIATION MATRIX file ' //
     &          'or "NONE"', FSREAD3, CRL//'SMAT', 'GETRECS' )

        SPCFLAG = ( SNAME .NE. 'NONE' )

        GNAME = PROMPTMFILE(
     &          'Enter logical name for GRIDDING MATRIX file',
     &          FSREAD3, CRL//'GMAT', 'GETRECS' )

C.........  Prompt for control matrix for area and point sources only
        IF( CATEGORY .NE. 'MOBILE' ) THEN
            CNAME = PROMPTMFILE(
     &              'Enter logical name for CONTROL MATRIX file ' //
     &              'or "NONE"', FSREAD3, CRL//'CMAT', 'GETRECS' )
        ELSE
           CNAME = 'NONE'
        ENDIF
        CTRLFLAG = ( CNAME .NE. 'NONE' )

C.........  Get special files for certain source categories
        IF( EMSFLAG .EQ. .TRUE. ) THEN
            SDEV = PROMPTFFILE(
     &             'Enter name for PSRC file',
     &             .TRUE., .TRUE., 'PSRC', 'GETRECS' )
 
        ELSEIF( CATEGORY .EQ. 'AREA' ) THEN
            ADEV = PROMPTFFILE(
     &             'Enter name for ASCT DESCRIPTION file',
     &             .TRUE., .TRUE., 'ASCTNAM', 'GETRECS' )

        ELSEIF( CATEGORY .EQ. 'MOBILE' ) THEN

            UNAME = PROMPTMFILE( 
     &              'Enter logical name for UNGRIDDING MATRIX file',
     &              FSREAD3, 'MUMAT', 'GETRECS' )

            VDEV = PROMPTFFILE(
     &             'Enter logical name for VMT MIX file',
     &             .TRUE., .TRUE., 'MVMTM', 'GETRECS' )

        ENDIF

C.........  Prompt for starting date
        SDATE  = GETNUM( JDATE, 9999999, JDATE,
     &                  'Enter starting date (YYYYDDD)' )

        IF( SDATE .NE. JDATE ) JTIME = 0  ! Reset because now could be 0-24

C.........  Prompt for starting time
        STIME  = GETNUM( JTIME, 235959, JTIME,
     &                  'Enter starting time (HHMMSS)' )

C.........  Prompt for time number of time steps
        I = SECSDIFF( JDATE, JTIME, SDATE, STIME ) / 3600
        NSTEPS = NSTEPS - I
        MESG   = 'Enter number of time steps'
        NSTEPS = GETNUM( 1, NSTEPS, NSTEPS, MESG )
        
C.........  Prompt for Cell or Source or Record number
        N = NMETHOD
        L = TRIMLEN( MENULST( 2 ) )

        IF( CATEGORY .EQ. 'AREA' ) THEN
            N = N - 1  ! Can't select by stack parms
            MENULST( 2 ) = MENULST( 2 )(1:L) // ' and ASCT code'

        ELSEIF( CATEGORY .EQ. 'POINT' ) THEN
            MENULST( 2 ) = MENULST( 2 )(1:L) // 
     &                     ', SCC, Plant ID, and Stack ID'

        ELSEIF( CATEGORY .EQ. 'MOBILE' ) THEN
            N = N - 1  ! Can't select by stack parms
            MENULST( 2 ) = MENULST( 2 )(1:L) // 
     &                     ', road type, and link ID'

        ENDIF
            
        METHOD = 1
        METHOD = GETMENU( N , METHOD, 
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

            FIP  = 0
            FIP  = GETNUM( 0, 99999, FIP, 'Enter FIP' )
            KEY2 = 0
            KEY3 = 0
            KEY4 = 0
            KEY5 = 0

            IF( CATEGORY .EQ. 'AREA' ) THEN
                KEY2 = GETNUM( 0, 9999999, 0, 'Enter ASCT7' )
                KEY3 = GETNUM( 0, 999, 0, 'Enter ASCT3' )

            ELSEIF( CATEGORY .EQ. 'MOBILE' ) THEN
                KEY2 = GETNUM( 0, 9999999, 0, 'Enter ROAD CLASS' )
                KEY3 = GETNUM( 0, 9999999, 0, 'Enter LINK ID (if any)' )

            ELSEIF( CATEGORY .EQ. 'POINT' ) THEN
                KEY2 = GETNUM( 0, 99999999, 0, 'Enter SCC' )
                KEY3 = GETNUM( 0, 9999, 0, 'Enter SIC' )
                KEY4 = GETNUM( 0, 99999999, 0, 'Enter PLANT' )
                KEY5 = GETNUM( 0, 99999999, 0, 'Enter STACK' )

            ENDIF

C.........  Prompt for source ID

        ELSEIF( METHOD .EQ. 3 ) THEN

            RN  = GETNUM( 1, NSRCS, 1,
     &                    'Enter SMOKE source ID number' )

C.........  Prompt for emissions based selection

        ELSEIF( METHOD .EQ. 4 ) THEN

            EMMTHD = 1
            EMMTHD = GETMENU( NIPOL, EMMTHD,
     &         'Choose inventory pollutant for ranking sources', EINAM )

C.........  Prompt for stack parameters

        ELSEIF( METHOD .EQ. 5 ) THEN

            STKMTHD = 1
            STKMTHD = GETMENU( NSTKMTHD, STKMTHD,
     &             'Choose stack feature for ranking sources', STKMENU ) 

        ENDIF

C.........  Check layer fractions file dimensions, if needed
        IF( LAYRFLAG ) THEN
            IF ( .NOT. DESC3( LNAME ) ) THEN
                CALL M3EXIT( 'GETRECS', 0, 0,
     &                      'Could not get description of file "'
     &                      // LNAME( 1:TRIMLEN( LNAME ) ) // '"', 2 )

            ELSEIF( NROWS3D .NE. NSRCS ) THEN
                WRITE( MESG, 94010 )
     &                 'Dimension mismatch.  PLAY file:', NROWS3D,
     &                 'inventory NPSRC=', NSRCS
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            ENDIF
        ENDIF

C.........  Read required variables depending on source category
C.........  For mobile sources, XLOC1 equiv to XLOCA, YLOC1 equiv to YLOCA,
C.........  XLOC2 equiv to STKHT, YLOC2 equiv to STKDM, and VMT equiv to STKTK
        IF( CATEGORY .NE. 'MOBILE' ) THEN
            CALL RDINVEN2( CATEGORY, ENAME, MXESRC, NSRCS, NIPOL, 
     &                     EINAM, IFIP, VAR2, VAR3, VAR4, VAR5,
     &                     TZONES, TPFLAG, INVYR, XLOCA, YLOCA, STKHT, 
     &                     STKDM, STKTK, STKVE, EMISV )
        ELSE
            CALL RDINVEN2( CATEGORY, ENAME, MXESRC, NSRCS, 1, 
     &                     'VMT', IFIP, VAR2, VAR3, VAR4, VAR5,
     &                     TZONES, TPFLAG, INVYR, XLOC1, YLOC1, XLOC2, 
     &                     YLOC2, RDUM, RDUM, VMT )

C.............  Initialize EMISV for mobile because will use it later
            DO 105 V = 1, NIPOL
                DO 101 S = 1, NSRCS
                    EMISV( S,V ) = 0.
101             CONTINUE
105         CONTINUE

        ENDIF

C.........  Get header information from gridding matrix
        SCRBUF = GRDNM
        IF ( .NOT. DESC3( GNAME ) ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                  'Could not get description of file "'
     &                  // GNAME( 1:TRIMLEN( GNAME ) ) // '"', 2 )

        ELSEIF( NROWS3D .NE. NGRID ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  ' // CRL // 'GMAT file:', NROWS3D,
     &      'program NGRID=', NGRID
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ELSEIF( NTHIK3D .NE. NSRCS ) THEN  ! NTHIK3D is used to store sources
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  ' // CRL // 'GMAT file:', NTHIK3D,
     &      'inventory file=', NSRCS
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ELSEIF( GDNAM3D .NE. GRDNM ) THEN 
            WRITE( MESG, 94010 )
     &      'Grid name mismatch.  ' // CRL // 'GMAT file: "' //
     &      SCRBUF ( 1:TRIMLEN( SCRBUF  ) ) // '" inventory file: "' //
     &      GDNAM3D( 1:TRIMLEN( GDNAM3D ) ) // '"'
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        END IF

C.........  Store parameters
        DX  = SNGL( XCELL3D )
        DY  = SNGL( YCELL3D )
        DDX = 1.0 / DX
        DDY = 1.0 / DY
        X0  = XORIG3D
        Y0  = YORIG3D

C.........  Read gridding matrix
        IF( CATEGORY .EQ. 'POINT' ) THEN

            CALL RDGRIDM( CATEGORY, GNAME,
     &                    NGRID, NSRCS, NSRCS, NP, IP, MP )

        ELSE  ! for AREA and MOBILE

            CALL RDGRIDM( CATEGORY, GNAME,
     &                    NGRID, NMATX, NSRCS, NA, IA, MA )

        ENDIF

C.........  Pre-process the gridding matrix. For all sources, create pointer
C.........    to exact position of IA/IP from NA/NP.
        K = 0
        DO 111 I = 1, NGRID
            DO 109 J = 1, NA( I ) ! Use ar source array, but pt source is equiv
                K = K + 1
                IF( J .EQ. 1 ) PG( I ) = K
109         CONTINUE
111     CONTINUE

C.........  Read mobile-specific files
        IF( CATEGORY .EQ. 'MOBILE' ) THEN

C.........  Read vehicle mix fractions tables, and multiply with
C.........  source VMT totals
            CALL RDVMIX( VDEV, INVYR, IFIP, VAR2, VAR3, VMT, VMTDIS )

C.............  Get header information from ungridding matrix
            IF( .NOT. DESC3( UNAME ) ) THEN
                MESG = 'Could not get description for file ' // UNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
 
            ELSEIF ( NROWS3D .NE. NMSRC ) THEN
                WRITE( MESG, 94010 )
     &           'Dimension mismatch.  Ungridding matrix:', NROWS3D,
     &           'program (NMSRC)=', NMSRC
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            ENDIF
            IF( NCOLS3D .NE. NMATX ) THEN
                WRITE( MESG, 94010 )
     &           'Dimension mismatch.  Ungridding matrix:', NCOLS3D,
     &           'program (NMATX)=', NMATX
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            ENDIF
            IF( NTHIK3D .NE. NGRID ) THEN
                WRITE( MESG, 94010 )
     &           'Dimension mismatch.  Ungridding matrix:', NTHIK3D,
     &           'program (NROWS * NCOLS)=', NGRID
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            ENDIF
 
            IF( EFLAG ) THEN
              MESG = 'Bad file dimensions.'
              CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            END IF

C.............  Read ungridding matrix
            IF ( .NOT. READ3( UNAME, ALLVAR3, ALLAYS3, 0, 0, NU ) ) THEN
                MESG = 'Could not read ungridding matrix from "' // 
     &                 UNAME( 1:TRIMLEN( UNAME ) ) // '"'
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            END IF

C.............  Pre-process the ungridding matrix.  Create pointer
C.............    to exact position of IU from NU
            K = 0
            DO 115 I = 1, NSRCS
                DO 113 J = 1, NU( I )
                    K = K + 1
                    IF( J .EQ. 1 ) PU( I ) = K
113             CONTINUE
115         CONTINUE
       

C.............  In order to sort mobile by emissions, and report total 
C.............  inventory emissions, must sum up emissions from MTMP and
C.............  aggregate across pollutants over all of the hours in the time
C.............  range.  Do this and store in EMISV array.

            JDATE = SDATE
            JTIME = STIME
            DO 121 T = 1, NSTEPS

C.................  Read all temporal emissions for given
                IF( .NOT. READ3( TNAME, ALLVAR3, ALLAYS3,
     &                           JDATE, JTIME, TMAT )    ) THEN
                    MESG = 'Could not read hourly emissions from "' //
     &                     TNAME( 1:TRIMLEN( TNAME ) ) // '"'
                    CALL M3EXIT( 'GETRECS', JDATE, JTIME, MESG, 2 )
                ENDIF

C.................  For all sources,sum over vehicle classes and emissions processes
                DO 119 E = 1, NEVAR
                    J = EVIDX( E )
                    DO 118 V = 1, NVTYPE
                        DO 117 S = 1, NSRCS
                            EMISV( S,J ) = EMISV( S,J ) + TMAT ( S,V,E )
           
117                     CONTINUE
118                 CONTINUE
119             CONTINUE

                CALL NEXTIME( JDATE, JTIME, TSTEP )

121         CONTINUE

        ENDIF  ! for mobile only

C.........  Create unsorted index and preprocess stack parms
        IF( METHOD .EQ. 4 .OR. METHOD .EQ. 5 ) THEN

            MEANVAL = 0.
            MEANCNT = 0
            DO 125 I = 1, NSRCS  ! If area sources, remember METHOD won't = 5
                INDEXA( I ) = I

                IF( CATEGORY .EQ. 'POINT' ) THEN
                    STKFL( I ) = 0.25 * PI * STKVE( I ) *
     &                                       STKDM( I ) * STKDM( I ) 
                ENDIF

                IF( METHOD .EQ. 4 ) THEN
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

125         CONTINUE  

            MEANVAL = MEANVAL / REAL( MEANCNT )

        ELSE ! still need index for area processing

            DO 127 I = 1, NSRCS
                INDEXA( I ) = I
127         CONTINUE

            IF( CATEGORY .EQ. 'POINT' ) THEN
                DO 129 I = 1, NSRCS
                    STKFL( I ) = 0.25 * PI * STKVE( I ) *
     &                                       STKDM( I ) * STKDM( I )
129             CONTINUE
            ENDIF

        ENDIF

        SRTSRC = 1  ! Defaults set to have no loop over sources (for METHOD=1)
        ENDSRC = 0
        INCSRC = 1

C.........  Set starting and ending sources for source-based methods
        IF( METHOD .EQ. 2 ) THEN

            SRTSRC = 1
            ENDSRC = NSRCS
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
            NRPRT = GETNUM( 1, NSRCS, NRPRT, 
     &                      'Enter number of sources to report' )

C.............  Sort based on selected category
            CALL SORTR1( NSRCS, INDEXA, CRITVAL )

C.............  Set starting and ending sources and for medium output list, 
C               must first calculate which sources is closest to mean value

            MINDIFF = 1.0E36
            MINVAL  = 1.0E36
            IF( RANGE .EQ. 1 ) THEN      ! Highest X sources

                SRTSRC = NSRCS
                ENDSRC = NSRCS - NRPRT + 1
                INCSRC = -1

            ELSEIF( RANGE .EQ. 2 ) THEN  ! Middle X sources

                DO 133 I = 1, NSRCS

                    J = INDEXA( I )
                    THISDIFF = ABS( MEANVAL - CRITVAL( J ) )
                    IF( THISDIFF .LT. MINDIFF ) THEN
                        MEANSRC = I
                        MINDIFF = THISDIFF
                    ENDIF

133             CONTINUE

                SRTSRC = MEANSRC - NRPRT / 2 + 1
                ENDSRC = MEANSRC + NRPRT / 2 + MOD( NRPRT,2 )
                INCSRC = 1
            
            ELSEIF( RANGE .EQ. 3 ) THEN ! Lowest X sources

                DO 136 I = 1, NSRCS

                    J = INDEXA( I )
                    IF( CRITVAL( J ) .GT. 0      .AND. 
     &                  CRITVAL( J ) .LT. MINVAL       ) THEN
                        MINSRC = I
                        MINVAL = CRITVAL( J )
                    ENDIF

136             CONTINUE

                SRTSRC = MINSRC
                ENDSRC = MINSRC + NRPRT - 1
                INCSRC = 1

            ENDIF

        ENDIF

        ICNT = 0  ! Set for do two loops

C.........  Create list of source IDs based on source parameters
C.........  Process using source-based loop
        DO 144 I = SRTSRC, ENDSRC, INCSRC

            J = INDEXA( I )

            IF( ICNT .GE. MXEOUT ) THEN   ! Prevent overflow, but continue count

                ICNT = ICNT + 1

            ELSEIF( METHOD .EQ. 2 ) THEN  !  By source characeristics
                IF( MOD( FIP,1000 ) .EQ. 0 ) THEN

                    SID = ( IFIP( J ) / 1000 ) * 1000 

                    IF( (FIP  .EQ. SID       .OR. FIP  .EQ. 0) .AND.
     &                  (KEY2 .EQ. VAR2( J ) .OR. KEY2 .EQ. 0) .AND.
     &                  (KEY3 .EQ. VAR3( J ) .OR. KEY3 .EQ. 0) .AND.
     &                  (KEY4 .EQ. VAR4( J ) .OR. KEY4 .EQ. 0) .AND.
     &                  (KEY5 .EQ. VAR5( J ) .OR. KEY5 .EQ. 0)) THEN
                        ICNT = ICNT + 1
                        EINX( ICNT ) = J
                    ENDIF

                ELSEIF( (FIP  .EQ. IFIP( J ) .OR. FIP  .EQ. 0) .AND.
     &                  (KEY2 .EQ. VAR2( J ) .OR. KEY2 .EQ. 0) .AND.
     &                  (KEY3 .EQ. VAR3( J ) .OR. KEY3 .EQ. 0) .AND.
     &                  (KEY4 .EQ. VAR4( J ) .OR. KEY4 .EQ. 0) .AND.
     &                  (KEY5 .EQ. VAR5( J ) .OR. KEY5 .EQ. 0)) THEN
                    ICNT = ICNT + 1
                    EINX( ICNT ) = J
 
                ENDIF
                
            ELSEIF( METHOD .EQ. 3 ) THEN  ! Single source

                ICNT = 1
                EINX( ICNT ) = RN    

            ELSEIF( METHOD .EQ. 4 .OR. METHOD .EQ. 5 ) THEN  ! By emis/stk param

                ICNT = ICNT + 1
                EINX( ICNT ) = J

            ENDIF

144     CONTINUE

C.........  Create pointer based on cell
C.........  Process using cell-based loop
C.........  NOTE: Using area & mobile source arrays NA and IA, but these will
C.........  map to point source arrays NP and IP because of equivalence.

        IF( METHOD .EQ. 1 ) THEN

            IDCELL = ( YN-1 ) * NCOLS + XN  ! Cell selected by user

            DO 151 N = 1, NA( IDCELL )

                PTR = PG( IDCELL ) + N - 1  ! Pointer created to position in IA
                ICNT = ICNT + 1
                IINX( ICNT ) = ICNT         ! pointer for sorting
                EINXA( ICNT ) = IA( PTR )   ! unsorted list of source IDs

151         CONTINUE

            CALL SORTI1( ICNT, IINX, EINXA )
            DO 155 I = 1, ICNT 
                J = IINX( I )
                EINX ( I ) = EINXA ( J )
155         CONTINUE

        ENDIF ! if method 1

        NOUT = ICNT
        IF( NOUT .GT. MXEOUT ) THEN
            WRITE( MESG, 94010 ) 'Output sources found =', NOUT,
     &                           'but maximum (MXEOUT) =', MXEOUT,
     &                           '. Reset MXEOUT and try again.'
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ENDIF

C.........  Read and process speciation matrix file

        IF( SPCFLAG ) THEN

C.............  Category-specific pre-read setup
            IF( CATEGORY .NE. 'MOBILE' ) THEN  ! maximum dims different 
                V = MXMPOL
                N = MXESRC
            ELSE 
                V = MXSMV
                N = NMSRC * NVTYPE
            ENDIF

C.............   Read file
            CALL RDSPECM( CATEGORY, SNAME, N, TMPSRC, V,
     &                    NVSMAT, SS, SPCNAM )

C.............   Check if actual number of sources is consistent
            IF( TMPSRC .NE. NSRCS ) THEN
                WRITE( MESG, 94010 )
     &            'Source number mismatch. Speciation matrix:', 
     &            TMPSRC, '. Inventory file=', NSRCS
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            END IF

C.............  Category-specific post-read setup
            IF( CATEGORY .NE. 'MOBILE' ) THEN 
                NAMPARS = '_to_'
                LPARS = 4
            ELSE 
                NAMPARS = '_'
                LPARS = 1
            END IF


            DO 161 V = 1, NIPOL
                NSPCOUT( V ) = 0
161         CONTINUE

C.............  Check speciation variable names to make sure they are consisent
            L2 = TRIMLEN( SNAME )
            NMPOL = 0
            DO 166 V = 1, NVSMAT

                K = INDEX( SPCNAM( V ), NAMPARS(1:LPARS) ) ! postn of 1st parser
                L = TRIMLEN( SPCNAM( V ) )
                IF( K .LT. 1 ) THEN
                   EFLAG = .TRUE.
                    MESG = 'Bad conversion factor name "' //
     &                     SPCNAM( V )( 1:L ) // '" in file "' //
     &                     SNAME( 1:L2 ) // '"'
                    CALL M3MSG2( MESG )
                END IF

                IF ( CATEGORY .EQ. 'MOBILE' ) THEN
                    K = K + 
     &                  INDEX( SPCNAM( V )( K+1:L ), NAMPARS(1:LPARS) )
                END IF

                SCRBUF = SPCNAM( V )( 1:K-1 )   ! get emis mode (e.g., EXH_CO)
                I = INDEX1( SCRBUF, NEVAR, TPRLNAM ) ! create index

                SCRBUF = SPCNAM( V )( K+LPARS:L )        ! get species name
                J = INDEX1( SCRBUF, NMPOL, EMNAM )
                IF ( J .LE. 0 ) THEN
                    NMPOL = NMPOL + 1      ! add 1 to list of species
                    EMNAM( NMPOL ) = SCRBUF ! store new species
                    J = NMPOL
                ENDIF

                F = NSPCOUT( I ) + 1 ! count no of smat vars per tmprl variable
                SVARREF( I,F ) = V   ! store which smat vars per tmprl variable
                SPCREF ( I,F ) = J   ! store which species for tmprl variable
                SMREFM ( V )   = J   ! store which species for smat variable
                NSPCOUT( I )   = F
 
166         CONTINUE

            IF( EFLAG ) THEN
                MESG = 'Bad speciation matrix'
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            ENDIF

        ENDIF   ! whether speciation matrix read in or not

C.........  Read control matrix

        IF( CTRLFLAG ) THEN

            CALL RDCNTLM( CATEGORY, CNAME, MXESRC, TMPSRC, MXIPOL,
     &                    NCINVP, CC, CTLINV )

C.............   Check if actual number of sources is consistent
            EFLAG = .FALSE.
            IF( TMPSRC .NE. NSRCS ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &            'Source number mismatch. ' // CRL // 'CMAT file:', 
     &            TMPSRC, '. Inventory file=', NSRCS
                CALL M3MSG2( MESG )
            END IF

            IF( NCINVP .NE. NIPOL ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &            'Pollutant number mismatch. ' // CRL // 'CMAT file:',
     &            NCINVP, '. Inventory file=', NIPOL
                CALL M3MSG2( MESG )
            END IF

C.............  Check pollutants
            DO 177 V = 1, NCINVP

                I = INDEX1( CTLINV( V ), MXIPOL, INVSPEC )
                IF ( I .LT. 0 ) THEN
                   EFLAG = .TRUE.
                   MESG = 'Bad control pollutant "' //
     &                     CTLINV( V )( 1:TRIMLEN( CTLINV(V) ) ) //
     &                     '" in ' // SNAME
                    CALL M3MSG2( MESG )
                END IF

177         CONTINUE

            IF( EFLAG ) THEN
                MESG = 'Bad control matrix'
                CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )
            ENDIF

        END IF          !  if control flag

C.........  Populate plant and stack character arrays (different if have PSRC)
        IF( CATEGORY .EQ. 'POINT' ) THEN
            DO 181 S = 1, NSRCS

                IF( EMSFLAG ) THEN
                    READ ( SDEV, 93005, IOSTAT=IOS ) CPLT, CSTK, CDESC

                    IF( IOS .NE. 0 ) THEN
                        WRITE( MESG,94010 )
     &                      'Unit number', SDEV,
     &                      'I/O Status ', IOS,
     &                      'Line number', S
                        CALL M3MSG2( MESG )
                        CALL M3EXIT( 'GETRECS', 0, 0,
     &                      'Error reading PSRC file.', 2 )
                    ELSE
                        CPLANT( S ) = CPLT
                        CSTACK( S ) = CSTK
                        CPDESC( S ) = CDESC
                    ENDIF
                
                ELSE
                    WRITE( CPLT, '(I15)' ) VAR4( S )
                    WRITE( CSTK, '(I15)' ) VAR5( S )

                    CPLANT( S ) = CPLT
                    CSTACK( S ) = CSTK
                    CPDESC( S ) = ' '

                ENDIF

181         CONTINUE
        ENDIF

C.........  Preprocess gridding matrix for area and mobile source processing

        IF( CATEGORY .NE. 'POINT' ) THEN

C.............  Set range of cells to check for sources based on method
            IF( METHOD .EQ. 1 ) THEN
                L1 = IDCELL
                L2 = IDCELL
            ELSE
                L1 = 1
                L2 = NGRID
            ENDIF

C.............  Loop through gridding matrix looking for sources in EINX.  For
C.............  method one, this will just be looping through srcs set by cell.
C.............  Count up cells per source in CCNT and store gridding coeffs.
C.............  Cannot assume that EINX is sorted

C.............  Loop through cells of interest
            DO 188 C = L1, L2
                K = PG( C ) - 1

C.................  Loop through sources in each cell of interest
                DO 186 N = 1, NA( C )
                    K = K + 1

C.....................  Try to find source in the EINX list
                    ICNT = 0
                    DO 183 I = 1, NOUT
                        IF( IA( K ) .EQ. EINX( I ) ) ICNT = I
183                 CONTINUE

C.....................  If this cell has a source in the EINX list store cell
                    IF( ICNT .GT. 0 ) THEN
                        I = ICNT
                        CCNT( I ) = CCNT( I ) + 1
                        J = CCNT( I )

                        YNUM ( I,J ) = INT( C / NCOLS ) + 1
                        XNUM ( I,J ) = C - ( YNUM(I,J) - 1 ) * NCOLS
                        GCOEF( I,J ) = MA( K )

                    ENDIF
186             CONTINUE
188         CONTINUE

        ENDIF

C.........  Preprocess output by output source
        DO 199 I = 1, NOUT

            GSTATE( I ) = 'FOUND'

            IF( CATEGORY .EQ. 'POINT' ) THEN

C.............  Determine cell IDs for point sources

                IF( METHOD .NE. 1 ) THEN  ! not cell selection

                    XX = DDX * ( XLOCA( EINX(I) ) - X0 )
                    YY = DDY * ( YLOCA( EINX(I) ) - Y0 )

                    CCNT( I ) = 1
                    IF( XX .GT. 0 .AND. YY .GT. 0 ) THEN
                        XNUM( I,1 ) = INT( XX ) + 1
                        YNUM( I,1 ) = INT( YY ) + 1
                    ELSE
                        XNUM( I,1 ) = INT( XX )
                        YNUM( I,1 ) = INT( YY )
                    ENDIF

                ELSE ! cell selection
                    CCNT( I ) = 1
                    XNUM( I,1 ) = XN
                    YNUM( I,1 ) = YN

                ENDIF

                IF( XNUM( I,1 ) .LE. 0     .OR. 
     &              XNUM( I,1 ) .GT. NCOLS .OR.
     &              YNUM( I,1 ) .LE. 0     .OR.
     &              YNUM( I,1 ) .GT. NROWS      ) THEN
                    GSTATE( I ) = 'OUTSIDE'
                    GO TO 199 ! to head of loop
                ENDIF

                IDCELL = ( YNUM( I,1 )-1 )*NCOLS + XNUM( I,1 )

C.................  Count up sources up to cell of interest
                SCNT = PG( IDCELL )

C.................  Create sorted source list for cell
                DO 191 S = 1, NP( IDCELL )
                    ISLOC( S ) = S
191             CONTINUE

                CALL SORTI1( NP( IDCELL ), ISLOC, IP( SCNT ) )

                DO 195 S = 1, NP( IDCELL )
                    ISLOC( S )= IP( SCNT + ISLOC( S )-1 )
195             CONTINUE

                F = FIND1( EINX( I ),NP( IDCELL ), ISLOC )
                                
                IF( F .LE. 0 ) THEN
                    GSTATE( I ) = 'UNFOUND'
                ENDIF

            ELSEIF( CATEGORY .NE. 'POINT' .AND. CCNT(I) .EQ. 0 ) THEN
               GSTATE( I ) = 'OUTSIDE'

            ENDIF ! if point sources or not

199     CONTINUE

C.........  Loop over time
        JDATE = SDATE
        JTIME = STIME
        DO 244 T = 1, NSTEPS

C.............  Read layer fractions for given hour
            IF( LAYRFLAG ) THEN

                DO 205 L = 1, EMLAYS

                    IF( READ3( LNAME, 'LFRAC', L,
     &                         JDATE, JTIME, LFRAC1L ) ) THEN

                        DO 201 I = 1, NOUT
                            OFRAC( I, L, T ) = LFRAC1L( EINX( I ) )
201                     CONTINUE

                    ELSE ! Read failed
                        MESG = 'Could not read "LFRAC" from '// LNAME
                        CALL M3EXIT( 'GETRECS', JDATE, JTIME, MESG, 2 )
                    ENDIF
    
205             CONTINUE

            ENDIF

C............  Initialize output emissions
            DO 218 I = 1, NOUT

                DO 216 S = 1, NMPOL
                    EOUT( I,S,T ) = 0.
216             CONTINUE

                DO 217 J = 1, NIPOL
                    ETMP( I,J,T ) = 0.
217             CONTINUE

218         CONTINUE

C.............  Loop over inventory species present in inventory file
            DO 241 V = 1, NEVAR 

C.................  Read temporal emissions for given hour and species
C.................  NOTE: EMIST equivalenced to TMAT
                IF( READ3( TNAME, TPRLNAM( V ), ALLAYS3,
     &                     JDATE, JTIME, EMIST ) ) THEN

C.....................  Loop over list of sources
                    DO 233 I = 1, NOUT

C.........................  Loop over model species for this variable
                        IF( SPCFLAG ) THEN
                            DO 225 J = 1,  NSPCOUT( V )

                                K = SVARREF( V,J ) ! spc variable index
                                S = SPCREF ( V,J ) ! species index

C.................................  Multiply PTMP, SMAT, and CMAT for source

                                IF( CATEGORY .NE. 'MOBILE' ) THEN
                                    SFAC = SS( EINX( I ), K ) 
                                ELSE
                                    SFAC = SM( EINX( I ), K ) 
                                ENDIf

                                IF( CTRLFLAG ) THEN  ! always false for mobile
                                    CFAC = CC( EINX( I ), V )
                                ELSE
                                    CFAC = 1.0
                                ENDIF
                             
                                IF( CATEGORY .NE. 'MOBILE' ) THEN
                                    EOUT( I,S,T ) = EMIST( EINX(I) ) * 
     &                                              SFAC * CFAC
                                ELSE
                                    DO 222 L = 1, NVTYPE
                                      EOUT( I,S,T )=EOUT( I,S,T ) +
     &                                              TMAT( EINX(I),L,1 )*
     &                                              SFAC * CFAC
222                                 CONTINUE
                                ENDIF

225                         CONTINUE

                        ENDIF

                        J = EVIDX( V )
                        IF( CATEGORY .NE. 'MOBILE' ) THEN
                            ETMP( I,J,T ) = EMIST( EINX( I ) )
                        ELSE
                            DO 231 L = 1, NVTYPE
                                ETMP( I,J,T ) = ETMP( I,J,T ) +
     &                                          TMAT( EINX( I ),L,1 )
231                         CONTINUE
                        ENDIF

233                 CONTINUE

                ELSE ! Read failed

                    MESG = 'Could not read "' //
     &                     TPRLNAM( V )( 1:TRIMLEN( TPRLNAM(V) ) ) //
     &                     '" from file "' // TNAME // '"'
                    CALL M3EXIT( 'GETRECS', JDATE, JTIME, MESG, 2 ) 

                ENDIF 

241         CONTINUE

            CALL NEXTIME( JDATE, JTIME, TSTEP )

244     CONTINUE

C.........  Prompt for output temporal x-ref file name

        ODEV = PROMPTFFILE(
     &  'Enter name for ASCII OUTPUT file',
     &           .FALSE., .TRUE., 'LISTOUT', 'GETRECS' )
 
        IF ( ODEV .LE. 0 ) THEN
            CALL M3EXIT( 'GETRECS', 0, 0,
     &                   'Could not open ASCII OUTPUT file', 2 )
        END IF

C.........  Write header to give file a title

        L  = TRIMLEN( MENULST( METHOD ) )
        SCRBUF = GRDNM
        WRITE( ODEV, 93200 ) MENULST( METHOD )(1:L), 
     &                       SCRBUF( 1:TRIMLEN( SCRBUF ) )

        IF( METHOD .EQ. 1 ) THEN
            WRITE( ODEV, 93210 ) XN, YN

        ELSEIF( METHOD .EQ. 2 .AND. CATEGORY .EQ. 'AREA' ) THEN
            WRITE( ODEV, 93220 ) FIP, KEY2, KEY3

        ELSEIF( METHOD .EQ. 2 .AND. CATEGORY .EQ. 'MOBILE' ) THEN
            WRITE( ODEV, 93223 ) FIP, KEY2, KEY3

        ELSEIF( METHOD .EQ. 2 .AND. CATEGORY .EQ. 'POINT' ) THEN
            WRITE( ODEV, 93225 ) FIP, KEY2, KEY3, KEY4, KEY5

        ELSEIF( METHOD .EQ. 4 ) THEN
            L1 = TRIMLEN( RNGMENU( RANGE ) )
            WRITE( ODEV, 93230 ) EINAM( EMMTHD ), 
     &                           RNGMENU( RANGE )(1:L1), NRPRT

        ELSEIF( METHOD .EQ. 5 ) THEN
            L1 = TRIMLEN( RNGMENU( RANGE ) )
            L2 = TRIMLEN( STKMENU( STKMTHD ) )
            WRITE( ODEV, 93240 ) STKMENU( STKMTHD )(1:L2), 
     &                           RNGMENU( RANGE )(1:L1), NRPRT

        ENDIF
 
        WRITE( ODEV, * ) ' '

        TONSUNIT = '[tons]  '
        MOLEUNIT = '[gm mol]'

C.........  Loop over list of sources and output details

        DO 333 I = 1, NOUT
            S = EINX( I )

C.............  Write out all time INdependant data

            IF( CATEGORY .NE. 'POINT' ) THEN
                IF( CATEGORY .EQ. 'AREA' ) THEN
                    WRITE( ODEV, 93002 ) 
     &                     I, S, IFIP( S ), VAR2( S ), VAR3( S ),
     &                     INVYR( S ), TZONES( S ), TPFLAG( S )

                ELSEIF( VAR3( S ) .LE. 0 ) THEN ! mobile non-link
                    WRITE( ODEV, 93003 ) 
     &                     I, S, IFIP( S ), VAR2( S ),
     &                     INVYR( S ), TZONES( S ), TPFLAG( S )

                ELSE ! mobile link
                    WRITE( ODEV, 93004 ) 
     &                     I, S, IFIP( S ), VAR2( S ), VAR3( S ),
     &                     INVYR( S ), TZONES( S ), TPFLAG( S ),
     &                     XLOC1( S ), YLOC1( S ), 
     &                     XLOC2( S ), YLOC2( S ), UZONE
                ENDIF

                OUTCEL = CCNT( I )

                IF( OUTCEL .EQ. 0 ) THEN
                    WRITE( ODEV, 93007 ) GSTATE( I )

                ELSE
                    WRITE( ODEV, 93008 ) GSTATE( I ),
     &                                 ( XNUM ( I,C ), YNUM( I,C ), 
     &                                   GCOEF( I,C ), C = 1,OUTCEL )

C                    IF( CATEGORY .EQ. 'MOBILE' ) THEN
C For ungridding matrix
c                    ENDIF

                ENDIF

            ELSE

                L1 = LBLANK( CPLANT( S ) ) + 1
                L2 = LBLANK( CSTACK( S ) ) + 1
                CPLT = CPLANT( S )( L1:TRIMLEN( CPLANT( S ) ) )
                CSTK = CSTACK( S )( L2:TRIMLEN( CSTACK( S ) ) )
                L1 = TRIMLEN( CPLT )
                L2 = TRIMLEN( CSTK )

                WRITE( ODEV, 93000 )
     &             I, S, CPDESC( S ),
     &             IFIP( S ), VAR3( S ), VAR2( S ), CPLT, CSTK,
     &             INVYR( S ), TZONES( S ), TPFLAG( S ), 
     &             STKHT( S ), STKDM( S ), STKTK( S ), STKVE ( S ),
     &             STKFL( S ), XNUM( I,1 ), YNUM( I,1 ), GSTATE( I ),
     &             XLOCA( S ), YLOCA( S )
            ENDIF

C.............  Write out control coefficients

            IF( CTRLFLAG )
     &      WRITE( ODEV, 93010 ) 
     &           ( CTLINV( V )( 1:TRIMLEN( CTLINV( V ) ) ), 
     &             CC( S, V ), V = 1, NCINVP )

C.............  Convert inventory emissions based on type of input temporal
C.............  Write out time independent inventory emissions
C.............  Week-normal data - input was on day-basis
            IF ( CATEGORY .NE. 'MOBILE' .AND.
     &           MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 .OR.
     &           MOD( TPFLAG( S ), WDTPFAC ) .EQ. 0      ) THEN

                DO 255 V = 1, NIPOL
                    EMISV( S,V ) = EMISV( S,V ) * YR2DAY( INVYR( S ) )
255             CONTINUE

                WRITE( ODEV, 93020 ) 
     &               ( EINAM( V )( 1:TRIMLEN( EINAM( V ) ) ),
     &                 EMISV( S,V ), V = 1, NIPOL ), 'tons/day'

            ELSEIF( CATEGORY .NE. 'MOBILE' ) THEN  

                WRITE( ODEV, 93020 ) 
     &               ( EINAM( V )( 1:TRIMLEN( EINAM( V ) ) ),
     &                 EMISV( S,V ), V = 1, NIPOL ), 'tons/yr'

            ELSE

                WRITE( ODEV,93023 ) 
     &               ( EINAM( V )( 1:TRIMLEN( EINAM( V ) ) ),
     &                 EMISV( S,V ), V = 1, NIPOL ), 'tons/period'

                WRITE( ODEV, 93025 ) 'mi/day', VMT( S ), 
     &               ( VTYPE3( L   ), L=1,NVTYPE ),
     &               ( VMTDIS( S,L ), L=1,NVTYPE )

            ENDIF

C.............  For speciation...
            IF( SPCFLAG ) THEN

C.................  Average mobile speciation matrix over all processes for all
C.................     nehicle types
C.................  NOTE: It is a bit weird to average the speciation factors
C.................        because they do not all get equal weighting when 
C.................        actually applied (weight depends on emissions from 
C.................        each process).  Leave for future improvement.
                IF( CATEGORY .EQ. 'MOBILE' ) THEN

                    DO 277 V = 1, NVTYPE

                        DO 272 F = 1, NMPOL   ! Initialize
                            SFACS( F,V ) = 0.
                            NSCNT( F   ) = 0
272                     CONTINUE

                        DO 273 J = 1, NVSMAT  ! Loop over smat vars

                            F = SMREFM( J )   ! F is model species ID
                            K = ( V-1 )*NSRCS + S
                            SFACS( F,V ) = SFACS( F,V ) + SM( K,J )
                            NSCNT( F ) = NSCNT( F ) + 1

273                     CONTINUE

                        DO 275 F = 1, NMPOL
                            SFACS( F,V ) = SFACS( F,V ) / NSCNT( F )
275                     CONTINUE

277                 CONTINUE 

C.....................  Write out speciation coefficients for mobile
                    WRITE( ODEV,93030 ) ( VTYPE3( V ), V=1, NVTYPE )

                    DO 279 F = 1, NMPOL
                        WRITE( ODEV,93032 )
     &                       EMNAM(F),(  SFACS( F,V ), V=1,NVTYPE )
279                 CONTINUE 
                ELSE

C.....................  Write out speciation coefficients for area/point
                    WRITE( ODEV,93035 ) 
     &                   ( EMNAM( F ), SS( S,F ), F=1,NMPOL )
                ENDIF

            ENDIF

C.............  If the source is outside the domain, don't output time-stepped
            IF( GSTATE( I ) .NE. 'OUTSIDE' ) THEN

C.................  Write out header for temporal emissions
                IF ( SPCFLAG ) THEN
                    WRITE( ODEV, 93040 ) MMDDYY( SDATE ),
     &                           ( EINAM( V ), V = 1, NIPOL ), 
     &                           ( EMNAM( J ), J = 1, NMPOL ),
     &                           ( TONSUNIT  , K = 1, NIPOL ),
     &                           ( MOLEUNIT  , L = 1, NMPOL )

                ELSE
                    WRITE( ODEV, 93042 ) MMDDYY( SDATE ),
     &                           ( EINAM( V )        , V = 1, NIPOL ), 
     &                           ( TONSUNIT          , K = 1, NIPOL )
                ENDIF

C.................  Write out time dependent emissions
                JDATE = SDATE
                JTIME = STIME
                LDATE = -9
                DO 301 T = 1, NSTEPS
   
                    IF( SPCFLAG ) THEN
                        WRITE( ODEV, 93050 ) 
     &                         JDATE-SDATE+1, JTIME/10000, 
     &                         ( ETMP( I, V, T ), V=1,NIPOL ),
     &                         ( EOUT( I, J, T ), J=1,NMPOL )
                    ELSE
                        WRITE( ODEV, 93050 )
     &                         JDATE-SDATE+1, JTIME/10000, 
     &                         ( ETMP( I, V, T ), V=1,NIPOL )
                    ENDIF

                    IF( JDATE .NE. LDATE ) THEN

                        DO 288 V = 1, NIPOL     ! Initialize
                            EISUM( V ) = ETMP( I, V, T )
288                     CONTINUE

                        IF ( SPCFLAG ) THEN
                            DO 291 J = 1, NMPOL ! Initialize
                                EMSUM( J ) = EOUT( I, J, T )
291                         CONTINUE
                        ENDIF

                        LDATE = JDATE

                    ELSE

                        DO 295 V = 1, NIPOL
                            EISUM( V ) = EISUM( V ) + ETMP( I,V,T )
295                     CONTINUE

                        IF ( SPCFLAG ) THEN
                            DO 299 J = 1, NMPOL
                                EMSUM( J ) = EMSUM( J ) + EOUT( I,J,T )
299                         CONTINUE
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

301             CONTINUE

C.............  Write out alternative header for temporal
            ELSE
                WRITE( ODEV, 93044 )

            ENDIF     ! If source is in grid or not

            IF( LAYRFLAG .AND. GSTATE( I ) .NE. 'OUTSIDE' ) THEN

C.................  Write out header for layer fractions
                WRITE( ODEV, 93060 ) MMDDYY( SDATE ), 
     &                               ( L, L = 1, EMLAYS )

C.................  Write out time dependent layer fractions
C.................  Format will put blanks where there are zeroes
                JDATE = SDATE
                JTIME = STIME
                DO 322 T = 1, NSTEPS

                    WRITE( ODEV, 93070 ) 
     &                     JDATE-SDATE+1, JTIME/10000

                    DO 311 L = 1, EMLAYS - 1

                        VAL = OFRAC( I, L, T ) 
                        IF( VAL .GT. 0. ) THEN
                            WRITE( ODEV, 93080 ) VAL
                        ELSE
                            WRITE( ODEV, 93090 ) 
                        ENDIF

311                 CONTINUE

                    VAL = OFRAC( I, EMLAYS, T )
                    IF( VAL .GT. 0. ) THEN
                        WRITE( ODEV, 93100 ) VAL
                    ELSE 
                        WRITE( ODEV, 93110 )
                    ENDIF

                    CALL NEXTIME( JDATE, JTIME, TSTEP )

322             CONTINUE

            ENDIF

            WRITE( ODEV, '(A)' ) ' '

333     CONTINUE 

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

93002   FORMAT( I5, ' Src ID:', I7, / ,
     &          5X, ' FIP:', I5.5, ' ASCT:', I7.7, I3.3,
     &              ' YR:' , I4, ' ZON:', I2.2, ' TPF:', I1 )

93003   FORMAT( I5, ' Src ID:', I7, / ,
     &          5X, ' FIP:', I5.5, ' RDCLAS:', I8,
     &              ' YR:' , I4, ' ZON:', I2.2, ' TPF:', I1 )

93004   FORMAT( I5, ' Src ID:', I7, / ,
     &          5X, ' FIP:', I5.5, ' RDCLAS:', I8, ' ILINK:', I8,
     &              ' YR:' , I4, ' ZON:', I2.2, ' TPF:', I1, /,
     &          5X, ' Link info...   (X1,Y1): (', F10.1, ',', F10.1, 
     &              ')  (X2,Y2): (', F10.1, ',', F10.1, ') UTM Zone= ',
     &              I2 )

93005   FORMAT( 29X, A15, 1X, A12, 14X, A40 )

93007   FORMAT( 5X, ' Gridding... Status:', A7 )

93008   FORMAT( 5X, ' Gridding... Status:', A7, 
     &              ' Cell (X,Y) / gridding factor', / , <OUTCEL>
     &              ( 18X, 3( '(', I3,',', I3, ') / ',E11.5, 3X, : )/ ))

93010   FORMAT( 5X, ' Control facs...  ',
     &              30( 1X, A, ': ', E11.5 ) )

93020   FORMAT( 5X, ' Inv emissions...',
     &              <NIPOL>( 1X, A, ': ', E11.5 ), ' [', A, ']' )

93023   FORMAT( 5X, ' Total emissions...',
     &              <NIPOL>( 1X, A, ': ', E11.5 ) ' [', A, ']' )

93025   FORMAT( 5X, ' Inventory VMT [', A, ']...' / 
     &          7X, ' Total:', 1X, E11.5, /, 
     &          15X, <NVTYPE>( 4X, A5, :, 3X ), /,
     &          7X, ' Disag:' <NVTYPE>( 1X, E11.5, : ) )

93030   FORMAT( 5X, ' Speciation coefficients, VOC values are means',
     &              ' over exhaust, evap, diurnal processes',
     &              ' [gm mol/tons]...', /,
     &              5X, '   VTYPE: ', <NVTYPE>( 4X, A5, 3X ) )

93032   FORMAT( 8X, A5, ':', <NVTYPE>( 1X, E11.5 ) )

93035   FORMAT( 5X, ' Speciation coefficients [gm mol/tons]...', /,
     &          5( 9X, 7( 1X, A4, ': ', E11.5 ), / ) )

93040   FORMAT( 5X, ' Hourly emissions starting ', A14, '...' /,
     &          5X, ' Dy Hr', <NIPOL>(1X, A3, 5X), <NMPOL>( 1X, A8 ), /,
     &          11X, <NIPOL+NMPOL>( 1X, A8 ) )

93042   FORMAT( 5X, ' Hourly emissions starting ', A14, '...' /,
     &          5X, ' Dy Hr', <NIPOL>(1X, A3, 5X) , /,
     &          11X, <NIPOL>( 1X, A8 ) )

93044   FORMAT( 5X, ' Hourly emissions not written because',   
     &              ' source is outside domain.' )

93050   FORMAT( 6X, I2.2, 1X, I2.2, <NIPOL+NMPOL>( 1X, E8.3 ) )

93055   FORMAT( 5X, ' Total', <NIPOL+NMPOL>( 1X, E8.3 ) )

93060   FORMAT( 5X, ' Layer fractions starting ', A14, '...' /,
     &          5X, ' Dy Hr', <EMLAYS>( 1X, ' Lyr(', I2.2, ')' ) )

93070   FORMAT( 6X, I2.2, 1X, I2.2, $ )

93080   FORMAT( 1X, E8.3, $ )

93090   FORMAT( 9X, $ )

93100   FORMAT( 1X, E8.3 )

93110   FORMAT( 9X )

93200   FORMAT( A, ' with grid "', A, '"' )

93210   FORMAT( 5X, 'Selected cell: (', I3, ',', I3, ')' )

93220   FORMAT( 5X, 'St/Co FIPS code= ', I5.5, ' and ASCT= ', 
     &          I7.7, I3.3 )

93223   FORMAT( 5X, 'St/Co FIPS code= ', I5.5, ', Road class= ', 
     &          I8, ', and Link= ', I8 )

93225   FORMAT( 5X, 'St/Co FIPS code= ', I5.5, ', SCC= ', I8.8,
     &          ', SIC= ', I4.4, ', Plant ID= ', I8, 
     &          ', and Stack ID= ', I8 )

93230   FORMAT( 5X, 'Pollutant "', A, '" used to rank ', A, /,
     &          5X, 'with X= ', I8 )

93240   FORMAT( 5X, A, ' used to rank ', A, /,
     &          5X, 'with X= ', I8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 1X ) )

        END
