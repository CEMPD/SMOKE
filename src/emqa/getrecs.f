
        PROGRAM GETRECS

C***********************************************************************
C  program body starts at line
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
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, GETNUM, INDEX1, PROMPTFFILE, PROMPTMFILE, GETYN, LEN_TRIM
C
C  REVISION  HISTORY:
C       Copied from getrecs.F version 1.12 on 3/99 by M Houyoux
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
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER       FIND1
        INTEGER       GETIFDSC
        INTEGER       GETMENU
        INTEGER       GETNUM
        LOGICAL       GETYN
        INTEGER       INDEX1
        CHARACTER*14  MMDDYY
        INTEGER       PROMPTFFILE
        CHARACTER*16  PROMPTMFILE
        INTEGER       SECSDIFF
        REAL          YR2DAY
 
        EXTERNAL      FIND1, GETIFDSC GETMENU, GETNUM, GETYN, INDEX1,
     &                MMDDYY, PROMPTFFILE, PROMPTMFILE, SECSDIF,
     &                YR2DAY

C.........  Local parameters
        INTEGER, PARAMETER :: MXEOUT   = 2000 ! Max no. of output records
        INTEGER, PARAMETER :: NMETHOD  = 5    ! no. of menu items 
        INTEGER, PARAMETER :: NMRANGE  = 3    ! no. of MRANGEs in MRANGE menu
        INTEGER, PARAMETER :: NSTKMTHD = 5    ! no. stack parm sort types

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C.........  Non-module Source arrays
        REAL   , ALLOCATABLE :: STKFL  ( : ) ! stack flow
        REAL   , ALLOCATABLE :: CRITVAL( : ) ! criteria for sorting output
        REAL   , ALLOCATABLE :: LFRAC1L( : ) ! 1-layer fractions info
        REAL   , ALLOCATABLE :: EMIST  ( : ) ! hourly emis values (tons/hr) 
        REAL   , ALLOCATABLE :: EMISV  ( :,: ) ! inven emis values (tons/yr) 

! NOTE: For now, all inventory pollutant allocated at once

c        INTEGER     VAR2  ( MXESRC )   !  PT: SCC, AR: ASC7
c        INTEGER     VAR3  ( MXESRC )   !  PT: SIC, AR: ASC3 
c        INTEGER     VAR4  ( MXESRC )   !  PT: Pltid (dim w/ MXESRC on purpose)
c        INTEGER     VAR5  ( MXESRC )   !  PT: Stkid (dim w/ MXESRC on purpose)
c        REAL        XLOC1 ( MXESRC )   !  link start UTM X-location (m)
c        REAL        YLOC1 ( MXESRC )   !  link start UTM Y-location (m)
c        REAL        XLOC2 ( MXESRC )   !  link end UTM X-location (m)
c        REAL        YLOC2 ( MXESRC )   !  link end UTM Y-location (m)
c        REAL        VMT   ( MXESRC )   !  mobile VMT from inventory 
c        REAL        VMTDIS( NMSRC, NVTYPE ) ! VMT per vehicle type

C.........  Define hourly emissions variables
c        REAL           TMAT  ( NMSRC, NVTYPE, MXTMV )  !  hourly mobile values

C.........  Define speciation matrix variables

C.........   Speciation matrix, either mole-based or mass-based

        REAL, ALLOCATABLE :: SMATX( :,: )    !  speciation coefficients

C...........   Gridding Matrix and helper array

        INTEGER, ALLOCATABLE :: GMATX( : ) ! Contiguous gridding matrix
        INTEGER, ALLOCATABLE ::    PG( : ) ! Points to position in gridding 
                                           !    matrix for each cell

C.........   Ungridding matrix variables (for mobile only)
 
c        INTEGER         NU( NMSRC )
c        INTEGER         IU( NMATX )
c        REAL            CU( NMATX )
 
c        COMMON  / UGRIDMAT / NU, IU, CU

c        INTEGER         PU( NMSRC )   ! Points to start of IU for each source

C.........  Define control matrix variables

        REAL, ALLOCATABLE :: CC( :,: )

C.........  Define Layer Fractions

c        REAL            LFRAC1L( MXPSRC )

C..........  LOCAL ALLOCATABLE ARRAYS...

C..........  Pollutant names and species names
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: EMNAM ( : ) ! species names
        CHARACTER(LEN=PLSLEN3), ALLOCATABLE :: SVDESC( : ) ! spec var descs
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: CTLINV( : ) ! mult cntl names

C..........  Arrays for processing inventory file
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: VNAMINV( : ) ! i/o api var names

C..........  Arrays for processing gridding matrix
        INTEGER, ALLOCATABLE :: ISLOC( : )    ! cell specific grid mat sources
        INTEGER, ALLOCATABLE :: CCNT ( : )    ! count of cells per output src
        INTEGER, ALLOCATABLE :: XNUM ( :,: )  ! x-cell for output sources
        INTEGER, ALLOCATABLE :: YNUM ( :,: )  ! y-cell for output sources

        REAL   , ALLOCATABLE :: GCOEF( :,: )  ! gridding coef for output srcs

        CHARACTER*7, ALLOCATABLE :: GSTATE ( : ) ! state of GMATX for output

C..........  Arrays for processing speciation matrix
        INTEGER, ALLOCATABLE :: NSPCOUT( : )    ! no. smat vars per var  
        INTEGER, ALLOCATABLE :: SMREFM ( : )    ! spcs indx for each SMAT var
        INTEGER, ALLOCATABLE :: SPCREF ( :, : ) ! which species per tmprl var
        INTEGER, ALLOCATABLE :: SVARREF( :, : ) ! which smat var per tmprl var
        INTEGER, ALLOCATABLE :: NSCNT  ( : )    ! no. emis procs for each spcs

        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: SDESCV( : ) ! smat variable names

C..........  Arrays for generating reports
        REAL   , ALLOCATABLE :: EISUM( : )     ! sum of day's inv pollutants
        REAL   , ALLOCATABLE :: EMSUM( : )     ! sum of day's model species
        REAL   , ALLOCATABLE :: EOUT ( :,:,: ) ! output emissions
        REAL   , ALLOCATABLE :: ETMP ( :,:,: ) ! out tmprl inv emission
        REAL   , ALLOCATABLE :: OFRAC( :,:,: ) ! output layer fractions

C..........  LOCAL FIXED-SIZE ARRAYS

        INTEGER         EINX ( MXEOUT )   !   sorted index for output sources
        INTEGER         EINXA( MXEOUT )   ! unsorted index for output sources
        INTEGER         IINX ( MXEOUT )   ! index of index for output sources

c        REAL            SFACS( MXMPOL, NVTYPE )! tmp output mb spec factors
        REAL            OUTL ( MXLAYS3 )  ! tmp output layer fractions

        LOGICAL         LCHR ( 9 )     ! indicator for valid source chars
        CHARACTER*300   CHARS( 9 )     ! source chars in array form

        CHARACTER(LEN=IOVLEN3) IVARNAMS( 50 ) ! all inven file var names

C.........  Logical names and unit numbers
        INTEGER         ADEV             ! Area-source ASCT descriptions file
        INTEGER         LDEV             ! IO/API initialization unit #
        INTEGER         ODEV             ! number for output file
        INTEGER         SDEV             ! PSRC file
        INTEGER         VDEV             ! VMT mix file

        CHARACTER*16    ANAME   !  logical name for additive control matrix
        CHARACTER*16    ENAME   !  logical name for inventory         input file
        CHARACTER*16    GNAME   !  logical name for grid       matrix input file
        CHARACTER*16    INAME   !  logical name for SDEV file (default)
        CHARACTER*16    LNAME   !  logical name for layer fractions   input file
        CHARACTER*16    NNAME   !  logical name for ungridding matrix
        CHARACTER*16    RNAME   !  logical name for reactivity control matrix
        CHARACTER*16    SNAME   !  logical name for speciation matrix input file
        CHARACTER*16    TNAME   !  logical name for temporal point source file
        CHARACTER*16    UNAME   !  logical name for multiplicative control mat

C.........  Other variables
        REAL            CFAC, DDX, DDY, DX, DY, SFAC, VAL
        REAL            X0, Y0, XX, YY
        REAL            MEANVAL      ! Mean of criterion output values
        REAL            MINVALU       ! Min non-zero of criterion output values
        REAL            MINDIFF      ! Minimum difference between src and mean
        REAL            RDUM         ! Dummy real var
        REAL            THISDIFF     ! Current difference between src and mean

        INTEGER         I, C, E, F, J, K, L, N, S, T, V   ! pntrs and counters

        INTEGER         NOUT              ! actual number of output emis

        INTEGER         RN, SCNT, TSTEP, XN, YN
        INTEGER      :: EMLAYS = 0   ! numb;er of layers
        INTEGER         EMMTHD       ! method of sorting based on emissions
        INTEGER         ENDSRC       ! Ending source count for output
        INTEGER         FIP, SID     ! Temporary vars
        INTEGER         GMATDIM      ! Dimension for gridding matrix
        INTEGER         IDCELL       ! Cell number from X and Y cells
        INTEGER         INCSRC       ! Increment for processing sources
        INTEGER         ICNT         ! counter for output sources
        INTEGER         IOS          ! I/O status
        INTEGER         JDATE, JTIME ! Current date and time
        INTEGER         KEY2, KEY3, KEY4, KEY5 ! tmp src chars for matching
        INTEGER         L1, L2
        INTEGER         LDATE        ! Previous date
        INTEGER         MEANCNT      ! count of sources contibuting to mean val
        INTEGER         MEANSRC      ! source ID that has value closest to mean
        INTEGER         METHOD       ! method of source selection
        INTEGER         MINSRC       ! source ID that has min non-zero value
        INTEGER         MXCPSRC      ! max cells per source
        INTEGER         NCHAR        ! number of point-source chars for source
        INTEGER      :: NCINVP = 0   ! no of mult control matrix variables
        INTEGER         NCOLS        ! number of columns
        INTEGER         NGRID        ! number of grid cells
        INTEGER         NGMAT        ! size of gridding matrix
        INTEGER      :: NMSPC  = 0   ! number of model species
        INTEGER      :: NSTEPS = 0   ! number of time steps
        INTEGER         NROWS        ! number of rows
        INTEGER         NRPRT        ! number of sources to report
        INTEGER         NS           ! tmp number of sources in a cell
        INTEGER         NSMATV       ! number of speciation variables
        INTEGER         OUTCEL       ! tmp number of output cells per source
        INTEGER         PTR          ! tmp index for gridding matrix
        INTEGER         MRANGE       ! method of selecting MRANGE of sorted srcs
        INTEGER         SDATE, STIME ! Starting date and time
        INTEGER         SRTSRC       ! Starting source count for output
        INTEGER         STKMTHD      ! method of sorting based on stack parms
        INTEGER         TMPSRC       ! tmp count of sources
        INTEGER         UZONE        ! UTM zone

        LOGICAL      :: AFLAG    = .FALSE. ! true: area sources
        LOGICAL      :: DFLAG    = .FALSE. ! true: use additive controls
        LOGICAL      :: EFLAG    = .FALSE. ! treu: error occured
        LOGICAL      :: LFLAG    = .FALSE. ! true: use layer fractions file
        LOGICAL      :: MFLAG    = .FALSE. ! true: mobile sources
        LOGICAL      :: PFLAG    = .FALSE. ! true: point sources
        LOGICAL      :: RFLAG    = .FALSE. ! true: use reactivity controls
        LOGICAL      :: SFLAG    = .FALSE. ! true: use speciation matrix 
        LOGICAL      :: TFLAG    = .FALSE. ! true: use temporal emissions
        LOGICAL      :: UFLAG    = .FALSE. ! true: use mulitplicative controls

        CHARACTER*8     TONSUNIT
        CHARACTER*8     MOLEUNIT
        CHARACTER*9     CNTLBUF
        CHARACTER*15    CPLT, CSTK
        CHARACTER*16    SCRBUF  !  scratch buffer
        CHARACTER*80    MENULST( NMETHOD )
        CHARACTER*80    RNGMENU( NMRANGE )
        CHARACTER*80    STKMENU( NSTKMTHD )
        CHARACTER*300   FMTBUF    !  scratch format buffer 
        CHARACTER*300   BUFFER    !  scratch output buffer 
        CHARACTER*300   MESG      !  scratch message buffer 

        CHARACTER(LEN=IOVLEN3) GRDNM  !  name for grid from gridding matrix
        CHARACTER(LEN=IOVLEN3) PBUF   !  pollutant name buffer
        CHARACTER(LEN=PLSLEN3) SVBUF  !  species description buffer


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

        CHARACTER*16  :: PROGNAME = 'GETRECS'   !  program name

C***********************************************************************
C   begin body of program GETRECS

        LDEV = INIT3()

! NOTE: For now, I have not segmented this code into subroutines, but should
!       go back later and do this

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Retrieve environment variable for controlling which source category
C           to process and which parts to process.  This will control which
C           input files will be prompted for
        MESG = 'Control source categories and processing parts'
        CALL ENVSTR( 'GETRECS_CONTROL', MESG, ' ', CNTLBUF, IOS  )
        
        AFLAG = ( INDEX( CNTLBUF, 'A' ) .GT. 0 ) ! area sources
        MFLAG = ( INDEX( CNTLBUF, 'M' ) .GT. 0 ) ! mobile sources
        PFLAG = ( INDEX( CNTLBUF, 'P' ) .GT. 0 ) ! point sources
        TFLAG = ( INDEX( CNTLBUF, 'T' ) .GT. 0 ) ! temporal
        SFLAG = ( INDEX( CNTLBUF, 'S' ) .GT. 0 ) ! speciation
        UFLAG = ( INDEX( CNTLBUF, 'U' ) .GT. 0 ) ! multiplicative controls
        DFLAG = ( INDEX( CNTLBUF, 'D' ) .GT. 0 ) ! additive controls
        RFLAG = ( INDEX( CNTLBUF, 'R' ) .GT. 0 ) ! reactivity controls
        LFLAG = ( INDEX( CNTLBUF, 'L' ) .GT. 0 ) ! layer fractions

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, INAME )

        TNAME = CRL // 'TMP'
        GNAME = CRL // 'GMAT'
        SNAME = CRL // 'SMAT_L'
        UNAME = CRL // 'CMAT'
        ANAME = CRL // 'AMAT'
        RNAME = CRL // 'RMAT'
        LNAME = 'PLAY'
        NNAME = ' '     ! mobile ungridding file

C.........  Prompt for inventory file
 
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII INVENTORY file',
     &           .TRUE., .TRUE., INAME, PROGNAME )

C.........  Retrieve header information
        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0,
     &                   'Error reading header from file "'
     &                   // ENAME( 1:LEN_TRIM( ENAME ) ) // '"', 2 )
        ELSE

            CALL GETSINFO

        END IF

C.........  For mobile, source number must be exact and must store UTM zone
c        IF( CATEGORY .EQ. 'MOBILE' ) THEN
c            IF( NSRC .NE. NMSRC ) THEN
c                WRITE( MESG,94010 )
c     &              'Dimension overflow.  INVENTORY file:', NSRC,
c     &              'program (NMSRC):', NMSRC
c                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
c            ENDIF
c            UZONE = REAL( P_ALP3D )

c        ENDIF

        IF( TFLAG ) THEN

C.............  Prompt for temporalized emissions
            MESG = 'Enter logical name for the POINT HOURLY ' //
     &             'EMISSIONS file'
            TNAME = PROMPTMFILE( MESG,FSREAD3,TNAME,PROGNAME )
                                
C.............  Check header information from temporalized emissions
            CALL RETRIEVE_IOAPI_HEADER( TNAME )
            CALL CHKSRCNO( CATDESC, TNAME, NROWS3D, NSRC, EFLAG )

C.............  Store header information
            JDATE = SDATE3D
            JTIME = STIME3D
            TSTEP = TSTEP3D
            NSTEPS= MXREC3D

! NOTE: Removed mobile-specific stuff for variables w/ names like EXH_CO
C............. Reset pollutant information if temporal file input
            J = 0
            DO V = 1, NVARS3D

                SCRBUF = VNAME3D( V )
     
                I = INDEX1( SCRBUF, NIPPA, EANAM )

                IF( I .GT. 0 ) THEN
                    J = J + 1                 
                    EINAM( J ) = SCRBUF

                ELSE
                    MESG = 'Skipping temporal emissions pollutant "' //
     &                     SCRBUF( 1:LEN_TRIM( SCRBUF ) ) //
     &                     '" because it is not in the inventory file.'
                    CALL M3MSG2( MESG ) 

                END IF

            END DO

            NIPOL = J

        END IF  ! End of temporal inputs processing or not

C.........  Open gridding matrix
        GNAME = PROMPTMFILE(
     &          'Enter logical name for GRIDDING MATRIX',
     &          FSREAD3, GNAME, PROGNAME )

C.........  Check dimensions of gridding matrix
        CALL RETRIEVE_IOAPI_HEADER( GNAME )
        CALL CHKSRCNO( CATDESC, GNAME, NTHIK3D, NSRC, EFLAG )

C.........  Store grid setting from the gridding matrix
        NCOLS = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
        NROWS = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
        NGRID = NROWS3D
        GRDNM = GDNAM3D
        DX    = SNGL( XCELL3D )
        DY    = SNGL( YCELL3D )
        DDX   = 1.0 / DX
        DDY   = 1.0 / DY
        X0    = XORIG3D
        Y0    = YORIG3D

C.........  Set size of gridding matrix
        IF( AFLAG .OR. MFLAG ) THEN
            NGMAT = NCOLS3D
        ELSEIF( PFLAG ) THEN
            NGMAT = NSRC
        END IF

        GMATDIM = NGRID + 2 * NGMAT

        IF( LFLAG ) THEN

C.............  Prompt for layer fractions
            LNAME = PROMPTMFILE( 
     &              'Enter logical name for LAYER FRACTIONS MATRIX',
     &              FSREAD3, LNAME, PROGNAME )

C.............  Check header information
            CALL RETRIEVE_IOAPI_HEADER( LNAME )
            CALL CHKSRCNO( CATDESC, LNAME, NROWS3D, NSRC, EFLAG )

C.............  Store header info
            EMLAYS = NLAYS3D

C.............  Allocate memory for reading later
            ALLOCATE( LFRAC1L( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LFRAC1L', PROGNAME )

        END IF

        IF( SFLAG ) THEN

C.............  Prompt for speciation matrix 
            SNAME = PROMPTMFILE(
     &              'Enter logical name for SPECIATION MATRIX',
     &              FSREAD3, SNAME, PROGNAME )

C.............  Check header of speciation matrix 
            CALL RETRIEVE_IOAPI_HEADER( SNAME )
            CALL CHKSRCNO( CATDESC, SNAME, NROWS3D, NSRC, EFLAG )

C.............  Store variable descriptions from header of speciation matrix
            NSMATV = NVARS3D
            NMSPC  = NVARS3D

            ALLOCATE( SVDESC( NSMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SVDESC', PROGNAME )
! NOTE: For now, assume one model species per speciation factor
            ALLOCATE( EMNAM( NSMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMNAM', PROGNAME )

            DO V = 1, NSMATV
                SVDESC( V ) = VDESC3D( V )

                L = INDEX( SVDESC( V ), SPJOIN )
                EMNAM ( V ) = ADJUSTL( VDESC3D( V )( L+1:PLSLEN3 ) )
            END DO

        END IF

C.........  Prompt for multiplicative control matrix 
        IF( UFLAG ) THEN
            UNAME = PROMPTMFILE(
     &              'Enter logical name for MULTIPLICATIVE CONTROL ' //
     &              'MATRIX', FSREAD3, UNAME, PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( UNAME )
            CALL CHKSRCNO( CATDESC, UNAME, NROWS3D, NSRC, EFLAG )

C.............  Store header information from file as needed
            NCINVP = NVARS3D

            ALLOCATE( CTLINV( NCINVP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CTLINV', PROGNAME )

            CTLINV( 1:NCINVP ) = VNAME3D( 1:NCINVP )

        END IF

C.........  Prompt for additive control matrix 
        IF( AFLAG ) THEN
            ANAME = PROMPTMFILE(
     &              'Enter logical name for ADDITIVE CONTROL ' //
     &              'MATRIX', FSREAD3, ANAME, PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( ANAME )
            CALL CHKSRCNO( CATDESC, ANAME, NROWS3D, NSRC, EFLAG )

        END IF

C.........  Prompt for reactivity control matrix 
        IF( RFLAG ) THEN
            RNAME = PROMPTMFILE(
     &              'Enter logical name for REACTIVITY CONTROL ' //
     &              'MATRIX', FSREAD3, RNAME, PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( RNAME )
            CALL CHKSRCNO( CATDESC, RNAME, NROWS3D, NSRC, EFLAG )

        END IF
 
        IF( AFLAG ) THEN
            ADEV = PROMPTFFILE(
     &             'Enter name for ASCT DESCRIPTION file',
     &             .TRUE., .TRUE., 'ASCTNAM', PROGNAME )

! Use for point and mobile sources also??

        ELSEIF( MFLAG ) THEN

            NNAME = PROMPTMFILE( 
     &              'Enter logical name for UNGRIDDING MATRIX file',
     &              FSREAD3, NNAME, PROGNAME )

            VDEV = PROMPTFFILE(
     &             'Enter logical name for VMT MIX file',
     &             .TRUE., .TRUE., 'MVMTM', PROGNAME )

        END IF

        IF( EFLAG ) THEN
! Error because of bad dimensions of inputs
        END IF

! Update this later to be consistent with SMKMERGE.  Note that for point 
! sources, we can be smart and adjust based on PTMP and PLAY
        IF( TFLAG ) THEN

C.............  Prompt for starting date
            SDATE  = GETNUM( JDATE, 9999999, JDATE,
     &                      'Enter starting date (YYYYDDD)' )

            IF( SDATE .NE. JDATE ) JTIME = 0  ! Reset because now could be 0-24

C.............  Prompt for starting time
            STIME  = GETNUM( JTIME, 235959, JTIME,
     &                      'Enter starting time (HHMMSS)' )

C.............  Prompt for time number of time steps
            I = SECSDIFF( JDATE, JTIME, SDATE, STIME ) / 3600
            NSTEPS = NSTEPS - I
            MESG   = 'Enter number of time steps'
            NSTEPS = GETNUM( 1, NSTEPS, NSTEPS, MESG )

        END IF  ! End for temporal inputs

C.........  Allocate arrays

        ALLOCATE( INDEXA( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
        ALLOCATE( STKFL( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKFL', PROGNAME )
        ALLOCATE( CRITVAL( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CRITVAL', PROGNAME )
        ALLOCATE( EMISV( NSRC,NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISV', PROGNAME )
        ALLOCATE( GMATX( GMATDIM ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GMATX', PROGNAME )
        ALLOCATE( PG( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PG', PROGNAME )

        IF( TFLAG ) THEN
            ALLOCATE( EMIST( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMIST', PROGNAME )
            ALLOCATE( EISUM( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EISUM', PROGNAME )
            ALLOCATE( EMSUM( NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMSUM', PROGNAME )
        END IF
 
C.........  Read gridding matrix

        GMATX = 1  ! array

        IF( AFLAG .OR. MFLAG ) THEN

            CALL RDGMAT( GNAME, NGRID, NGMAT, NGMAT, GMATX( 1 ), 
     &                   GMATX( NGRID+1 ), GMATX( NGRID+NGMAT+1 ) )

        ELSEIF( PFLAG ) THEN

            CALL RDGMAT( GNAME, NGRID, NSRC, 1, GMATX( 1 ), 
     &                   GMATX( NGRID+1 ), RDUM            )

        END IF

C.........  Prompt for Cell or Source or Record number
        N = NMETHOD
        L = LEN_TRIM( MENULST( 2 ) )

        IF( AFLAG ) THEN
            N = N - 1  ! Can't select by stack parms
            MENULST( 2 ) = MENULST( 2 )(1:L) // ' and ASCT code'

        ELSEIF( MFLAG ) THEN
            N = N - 1  ! Can't select by stack parms
            MENULST( 2 ) = MENULST( 2 )(1:L) // 
     &                     ', road type, and link ID'

        ELSEIF( PFLAG ) THEN
            MENULST( 2 ) = MENULST( 2 )(1:L) // 
     &                     ', SCC, Plant ID, and Stack ID'

        ENDIF

! NOTE: Make this a subroutine?
! NOTE: Overide N for now b/c only have grid-cell part working for pt sources
        METHOD = 1
        N      = 3
        MESG = 'Choose source sub-selection method'
        METHOD = GETMENU( N , METHOD, MESG, MENULST )

C.........  If Cell, Prompt for x-dir and y-dir cell numbers

        STKMTHD = 0

        IF( METHOD .EQ. 1 ) THEN

            XN  = 1
            XN  = GETNUM( 1, NCOLS, XN,
     &                    'Enter X-dir cell number' )
            YN  = 1
            YN  = GETNUM( 1, NROWS, YN,
     &                    'Enter Y-dir cell number' )

C.........  Prompt for FIP, SIC, SCC, PLANT, STACK
        ELSEIF( METHOD .EQ. 2 ) THEN

            MESG = 'METHOD NOT SUPPORTED!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            FIP  = 0
            FIP  = GETNUM( 0, 99999, FIP, 'Enter FIP' )
            KEY2 = 0
            KEY3 = 0
            KEY4 = 0
            KEY5 = 0

            IF( AFLAG ) THEN
                KEY2 = GETNUM( 0, 9999999, 0, 'Enter ASCT7' )
                KEY3 = GETNUM( 0, 999, 0, 'Enter ASCT3' )

            ELSEIF( MFLAG ) THEN
                KEY2 = GETNUM( 0, 9999999, 0, 'Enter ROAD CLASS' )
                KEY3 = GETNUM( 0, 9999999, 0, 'Enter LINK ID (if any)' )

            ELSEIF( PFLAG ) THEN
                KEY2 = GETNUM( 0, 99999999, 0, 'Enter SCC' )
                KEY3 = GETNUM( 0, 9999, 0, 'Enter SIC' )
                KEY4 = GETNUM( 0, 99999999, 0, 'Enter PLANT' )
                KEY5 = GETNUM( 0, 99999999, 0, 'Enter STACK' )

            ENDIF

C.........  Prompt for source ID

        ELSEIF( METHOD .EQ. 3 ) THEN

            RN  = GETNUM( 1, NSRC, 1,
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

C.........  Read required variables from inventory file depending on source 
C           category
        IF( AFLAG ) THEN

c            CALL RAREACHR( ENAME, SDEV, NSRC, NVAR, IVARNAMS )

        ELSEIF( MFLAG ) THEN

c            CALL RMOBLCHR( ENAME, SDEV, NSRC, NVAR, IVARNAMS )

        ELSEIF( PFLAG ) THEN

! Make this a subroutine that depends on source category, later:
! NOTE: These should be set for only those vars that exist. For example, boiler
!       ID might not be there, so don't want to include if it is not there
C.............  Set names of input inventory variables:
            IVARNAMS( 1 ) = 'IFIP'
            IVARNAMS( 2 ) = 'ISIC'
            IVARNAMS( 3 ) = 'CORIS'
            IVARNAMS( 4 ) = 'TZONES'
            IVARNAMS( 5 ) = 'TPFLAG'
            IVARNAMS( 6 ) = 'INVYR'
            IVARNAMS( 7 ) = 'XLOCA'
            IVARNAMS( 8 ) = 'YLOCA'
            IVARNAMS( 9 ) = 'STKHT'
            IVARNAMS( 10 ) = 'STKDM'
            IVARNAMS( 11 ) = 'STKTK'
            IVARNAMS( 12 ) = 'STKVE'
            IVARNAMS( 13 ) = 'CSCC'
            IVARNAMS( 14 ) = 'CPDESC'
            IVARNAMS( 15 ) = 'CSOURC'

C.............  Read source properties
            CALL RDINVCHR( 'POINT', ENAME, SDEV, NSRC, 15, IVARNAMS )

C.............  Read inventory emissions
            DO V = 1, NIPOL 

                PBUF = EINAM( V )

                IF( .NOT. READ3( ENAME, PBUF, ALLAYS3, 0, 0, 
     &                           EMISV( 1,V )                ) ) THEN
                    EFLAG = .TRUE.
                    L1 = LEN_TRIM( PBUF )
                    L2 = LEN_TRIM( ENAME )
                    MESG = 'Error reading "' // PBUF( 1:L1 ) //
     &                     '" from file "' // ENAME( 1:L2 ) // '."'
                    CALL M3MSG2( MESG )

                ENDIF

            END DO

C.............  Set indicator for valid point source characteristics
            LCHR = .FALSE.  ! array
            DO I = 1, NCHARS
                LCHR( I ) = .TRUE.
            END DO
 
        ENDIF

        IF( EFLAG ) THEN
! ERror reading inventory file

        END IF

C.............  Initialize EMISV for mobile because will use it later
c            DO 105 V = 1, NIPOL
c                DO 101 S = 1, NSRC
c                    EMISV( S,V ) = 0.
c101             CONTINUE
c105         CONTINUE
        
C.........  Pre-process the gridding matrix. For all sources, create index
C           to exact position in source-arrays for 
        K = 0
        DO I = 1, NGRID
            DO J = 1, GMATX( I ) 
                K = K + 1
                IF( J .EQ. 1 ) PG( I ) = K
            END DO
        END DO

! NOTE: Removed big section for mobile-only stuff.  See old code to
C       see what was being done here.

C.........  Create unsorted index and preprocess stack parms
        IF( METHOD .EQ. 4 .OR. METHOD .EQ. 5 ) THEN

            MEANVAL = 0.
            MEANCNT = 0
            DO S = 1, NSRC  ! If area sources, remember METHOD won't = 5
                INDEXA( S ) = I

                IF( PFLAG ) THEN
                    STKFL( S ) = 0.25 * PI * STKVE( S ) *
     &                                       STKDM( S ) * STKDM( S ) 
                ENDIF

                IF( METHOD .EQ. 4 ) THEN
                    CRITVAL( S ) = EMISV( I, EMMTHD )

                ELSEIF( STKMTHD .EQ. 1 ) THEN
                    CRITVAL( S ) = STKHT( S )

                ELSEIF( STKMTHD .EQ. 2 ) THEN
                    CRITVAL( S ) = STKDM( S )

                ELSEIF( STKMTHD .EQ. 3 ) THEN
                    CRITVAL( S ) = STKTK( S )

                ELSEIF( STKMTHD .EQ. 4 ) THEN
                    CRITVAL( S ) = STKVE( S )

                ELSEIF( STKMTHD .EQ. 5 ) THEN
                    CRITVAL( S ) = STKFL( S )

                ENDIF

                IF( CRITVAL( S ) .GT. 0. ) THEN
                    MEANVAL = MEANVAL + CRITVAL( S )
                    MEANCNT = MEANCNT + 1
                ENDIF

            END DO  

            MEANVAL = MEANVAL / REAL( MEANCNT )

        ELSE ! still need index for area processing

            DO S = 1, NSRC
                INDEXA( S ) = S
            END DO

            IF( PFLAG ) THEN
                DO S = 1, NSRC
                    STKFL( S ) = 0.25 * PI * STKVE( S ) *
     &                                       STKDM( S ) * STKDM( S )
                END DO
            ENDIF

        ENDIF

        SRTSRC = 1  ! Defaults set to have no loop over sources (for METHOD=1)
        ENDSRC = 0
        INCSRC = 1

C.........  Set starting and ending sources for source-based methods
        IF( METHOD .EQ. 2 ) THEN

            SRTSRC = 1
            ENDSRC = NSRC
            INCSRC = 1

        ELSEIF( METHOD .EQ. 3 ) THEN   ! for single source

            SRTSRC = 1
            ENDSRC = 1
            INCSRC = 1

C.........  If ranked by source characteristic...

        ELSEIF( METHOD .EQ. 4 .OR. METHOD .EQ. 5 ) THEN

C.............  Retrieve source sub-selection characteristics
            MRANGE = 1
            MRANGE = GETMENU( NMRANGE, MRANGE,
     &           'Choose MRANGE for reporting ranked sources', RNGMENU )

            NRPRT = 50
            NRPRT = GETNUM( 1, NSRC, NRPRT, 
     &                      'Enter number of sources to report' )

C.............  Sort based on selected category
            CALL SORTR1( NSRC, INDEXA, CRITVAL )

C.............  Set starting and ending sources and for medium output list, 
C               must first calculate which sources is closest to mean value

            MINDIFF = 1.0E36
            MINVALU = 1.0E36
            IF( MRANGE .EQ. 1 ) THEN      ! Highest X sources

                SRTSRC = NSRC
                ENDSRC = NSRC - NRPRT + 1
                INCSRC = -1

            ELSEIF( MRANGE .EQ. 2 ) THEN  ! Middle X sources

                DO S = 1, NSRC

                    J = INDEXA( S )
                    THISDIFF = ABS( MEANVAL - CRITVAL( J ) )
                    IF( THISDIFF .LT. MINDIFF ) THEN
                        MEANSRC = S
                        MINDIFF = THISDIFF
                    ENDIF

                END DO

                SRTSRC = MEANSRC - NRPRT / 2 + 1
                ENDSRC = MEANSRC + NRPRT / 2 + MOD( NRPRT,2 )
                INCSRC = 1
            
            ELSEIF( MRANGE .EQ. 3 ) THEN ! Lowest X sources

                DO S = 1, NSRC

                    J = INDEXA( S )
                    IF( CRITVAL( J ) .GT. 0      .AND. 
     &                  CRITVAL( J ) .LT. MINVALU       ) THEN
                        MINSRC = S
                        MINVALU = CRITVAL( J )
                    ENDIF

                END DO

                SRTSRC = MINSRC
                ENDSRC = MINSRC + NRPRT - 1
                INCSRC = 1

            ENDIF

        ENDIF

        ICNT = 0  ! Set for do two loops

C.........  Create list of source IDs based on source parameters
C.........  Process using source-based loop

! NOTE: This part has not been updated to work yet!
!       It will need to be different, & perhaps use character-based 
!       matching some way
        DO S = SRTSRC, ENDSRC, INCSRC

            J = INDEXA( S )

            IF( ICNT .GE. MXEOUT ) THEN   ! Prevent overflow, but continue count

                ICNT = ICNT + 1

            ELSEIF( METHOD .EQ. 2 ) THEN  !  By source characeristics
c                IF( MOD( FIP,1000 ) .EQ. 0 ) THEN

c                    SID = ( IFIP( J ) / 1000 ) * 1000 

c                    IF( (FIP  .EQ. SID       .OR. FIP  .EQ. 0) .AND.
c     &                  (KEY2 .EQ. VAR2( J ) .OR. KEY2 .EQ. 0) .AND.
c     &                  (KEY3 .EQ. VAR3( J ) .OR. KEY3 .EQ. 0) .AND.
c     &                  (KEY4 .EQ. VAR4( J ) .OR. KEY4 .EQ. 0) .AND.
c     &                  (KEY5 .EQ. VAR5( J ) .OR. KEY5 .EQ. 0)) THEN
c                        ICNT = ICNT + 1
c                        EINX( ICNT ) = J
c                    ENDIF

c                ELSEIF( (FIP  .EQ. IFIP( J ) .OR. FIP  .EQ. 0) .AND.
c     &                  (KEY2 .EQ. VAR2( J ) .OR. KEY2 .EQ. 0) .AND.
c     &                  (KEY3 .EQ. VAR3( J ) .OR. KEY3 .EQ. 0) .AND.
c     &                  (KEY4 .EQ. VAR4( J ) .OR. KEY4 .EQ. 0) .AND.
c     &                  (KEY5 .EQ. VAR5( J ) .OR. KEY5 .EQ. 0)) THEN
c                    ICNT = ICNT + 1
c                    EINX( ICNT ) = J
c 
c                ENDIF
                
            ELSEIF( METHOD .EQ. 3 ) THEN  ! Single source

                ICNT = 1
                EINX( ICNT ) = RN    

            ELSEIF( METHOD .EQ. 4 .OR. METHOD .EQ. 5 ) THEN  ! By emis/stk param

                ICNT = ICNT + 1
                EINX( ICNT ) = J

            ENDIF

        END DO

C.........  Create pointer based on cell.  This loop adds to index ICNT when
C           there are more sources as a result of the cell selected.

! NOTE: For now, IINX, EINXA, and EINX are compile-time dimensioned.  However,
!       there are not overflow checks, so could get core dump

        IF( METHOD .EQ. 1 ) THEN

            IDCELL = ( YN-1 ) * NCOLS + XN  ! Cell selected by user

            DO N = 1, GMATX( IDCELL )

                PTR = PG( IDCELL ) + N - 1 + NGRID  ! Pointer created to position in IA
                ICNT = ICNT + 1
                IINX( ICNT ) = ICNT          ! pointer for sorting
                EINXA( ICNT ) = GMATX( PTR ) ! unsorted list of source IDs

            END DO

            CALL SORTI1( ICNT, IINX, EINXA )
            DO I = 1, ICNT 
                J = IINX( I )
                EINX ( I ) = EINXA ( J )
            END DO

        END IF ! if method 1

        NOUT = ICNT
        IF( NOUT .GT. MXEOUT ) THEN
            WRITE( MESG, 94010 ) 'Output sources found =', NOUT,
     &                           'but maximum (MXEOUT) =', MXEOUT,
     &                           '. Reset MXEOUT and try again.'
            CALL M3EXIT( 'GETRECS', 0, 0, MESG, 2 )

        ENDIF
c        WRITE( MESG, 94010 ) 'The number of output sources is', NOUT
c        CALL M3MSG2( MESG )

C.........  Read and process speciation matrix file
! NOTE: Still can't handle mass-based (shouldn't be hard)
        IF( SFLAG ) THEN

C.............   Allocate memory for entire speciation matrix
            ALLOCATE( SMATX( NSRC,NSMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SMATX', PROGNAME )
           
            ALLOCATE( NSPCOUT( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NSPCOUT', PROGNAME )
            ALLOCATE( SVARREF( NIPOL,NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SVARREF', PROGNAME )
            ALLOCATE( SPCREF( NIPOL,NMSPC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPCREF', PROGNAME )
            ALLOCATE( SMREFM( NSMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SMREFM', PROGNAME )

            NSPCOUT = 0   ! array
            SVARREF = 0   ! array
            SPCREF  = 0   ! array
            SMREFM  = 0   ! array

C.............   Read file
            DO V = 1, NSMATV

                SVBUF = SVDESC( V )
                CALL RDSMAT( SNAME, SVBUF, SMATX( 1,V ) )
  
            END DO

C.............  Create index from speciation variables to inventory pollutants
C               that have been read in already
            DO  V = 1, NSMATV
                
                SCRBUF = SVDESC( V )( 1:IOVLEN3 )
                I = INDEX1( SCRBUF, NIPOL, EINAM )   ! create index

! NOTE: For now, this is a placeholder, because EMNAM are SVDESC match up
                L = INDEX( SVDESC( V ), SPJOIN )
                SCRBUF = SVDESC( V )( L+1:PLSLEN3 )  ! get species name
                J = INDEX1( SCRBUF, NMSPC, EMNAM )

                F = NSPCOUT( I ) + 1 ! count no of smat vars per pollutant
                SVARREF( I,F ) = V   ! store which smat vars per pollutant
                SPCREF ( I,F ) = J   ! store which species for pollutant
                SMREFM ( V )   = J   ! store which species for smat variable
                NSPCOUT( I )   = F
 
            END DO

        END IF   ! whether speciation matrix read in or not
 
! Insert section for reading and processing control matrix
        IF( UFLAG ) THEN
            ALLOCATE( CC( NSRC,NCINVP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CC', PROGNAME )

            DO V = 1, NCINVP
                IF( .NOT. READ3( UNAME, CTLINV( V ), 1, 0, 0,
     &                           CC( 1,V ) ) ) THEN
                    L = LEN_TRIM( UNAME )
                    MESG= 'Could not read ' // CATEGORY // 
     &                    ' CONTROL MATRIX'//
     &                    ' from file "' // UNAME( 1:L ) // '"'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
 
                END IF
            END DO
        END IF

! NOTE: Later, set the maximum cells per source based on area/mobile
        MXCPSRC = 1
        ALLOCATE( XNUM( NOUT,MXCPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XNUM', PROGNAME )
        ALLOCATE( YNUM( NOUT,MXCPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'YNUM', PROGNAME )
        ALLOCATE( CCNT( NOUT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CCNT', PROGNAME )
        ALLOCATE( GSTATE( NOUT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GSTATE', PROGNAME )

C.........  Preprocess gridding matrix for area and mobile source processing

        IF( .NOT. PFLAG ) THEN

            ALLOCATE( GCOEF( NOUT,MXCPSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GCOEF', PROGNAME )

C.............  Set MRANGE of cells to check for sources based on method
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
            DO C = L1, L2
                K = PG( C ) - 1 + NGRID

C.................  Loop through sources in each cell of interest
                DO N = 1, GMATX( C )
                    K = K + 1

C.....................  Try to find source in the EINX list
                    ICNT = 0
                    DO I = 1, NOUT
                        IF( GMATX( K ) .EQ. EINX( I ) ) ICNT = I
                    END DO

C.....................  If this cell has a source in the EINX list store cell
                    IF( ICNT .GT. 0 ) THEN
                        I = ICNT
                        CCNT( I ) = CCNT( I ) + 1
                        J = CCNT( I )

                        YNUM ( I,J ) = INT( C / NCOLS ) + 1
                        XNUM ( I,J ) = C - ( YNUM(I,J) - 1 ) * NCOLS
                        L = NGMAT + K
!                        GCOEF( I,J ) = CONVERT_GVAL( GMATX( K ) )
! NOTE: Later, add CONVERT_GVAL to get real value from memory and return it
                    ENDIF
                END DO
            END DO

        ENDIF

C.........  Preprocess output by output source
        DO I = 1, NOUT

            GSTATE( I ) = 'FOUND'

            IF( PFLAG ) THEN

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
                    CYCLE    ! to head of loop
                ENDIF

                IDCELL = ( YNUM( I,1 )-1 )*NCOLS + XNUM( I,1 )

C.................  Count up sources up to cell of interest
                SCNT = PG( IDCELL ) + NGRID

C.................  Create sorted source list for cell
                NS = GMATX( IDCELL )

C.................  NOTE: The following is total overkill, but it ensures
C                   each source is in the correct place in the gridding matrix
                ALLOCATE( ISLOC( NS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ISLOC', PROGNAME )
    
                DO J = 1, NS
                    ISLOC( J ) = J   ! create sorting index
                END DO

                CALL SORTI1( NS, ISLOC, GMATX( SCNT ) )

                DO J = 1, NS
                    ISLOC( J )= GMATX( SCNT + ISLOC( J )-1 )
                END DO

                F = FIND1( EINX( I ), NS, ISLOC )
                                
                IF( F .LE. 0 ) THEN
                    GSTATE( I ) = 'UNFOUND'
                ENDIF

                DEALLOCATE( ISLOC )

            ELSEIF( .NOT. PFLAG .AND. CCNT(I) .EQ. 0 ) THEN
               GSTATE( I ) = 'OUTSIDE'

            ENDIF ! if point sources or not

        END DO

C.........  Allocate memory for and initialize output emissions
        ALLOCATE( EOUT( MXEOUT, NMSPC+NIPOL, NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOUT', PROGNAME )
        ALLOCATE( ETMP( MXEOUT, NMSPC+NIPOL, NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOUT', PROGNAME )
        ALLOCATE( OFRAC( MXEOUT, EMLAYS, NSTEPS ), STAT=IOS ) ! EMLAYS might = 0
        CALL CHECKMEM( IOS, 'OFRAC', PROGNAME )

        EOUT = 0.  ! array
        ETMP = 0.  ! array

C.........  Loop over time steps, if any 

        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS

C.............  Read layer fractions for given hour
            IF( LFLAG ) THEN

                DO L = 1, EMLAYS

                    IF( READ3( LNAME, 'LFRAC', L,
     &                         JDATE, JTIME, LFRAC1L ) ) THEN

                        DO I = 1, NOUT
                            OFRAC( I, L, T ) = LFRAC1L( EINX( I ) )
                        END DO

                    ELSE ! Read failed
                        MESG = 'Could not read "LFRAC" from '// LNAME
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF
    
                END DO

            END IF

C.............  Loop over inventory species present in inventory file
            DO V = 1, NIPOL 

C.................  Read temporal emissions for given hour and species
C.................  NOTE: EMIST equivalenced to TMAT
                IF( READ3( TNAME, EINAM( V ), ALLAYS3,
     &                     JDATE, JTIME, EMIST ) ) THEN

C.....................  Loop over list of sources
                    DO I = 1, NOUT

C.........................  Loop over model species for this variable
                        IF( SFLAG ) THEN
                            DO J = 1, NSPCOUT( V )

                                K = SVARREF( V,J ) ! spc variable index
                                S = SPCREF ( V,J ) ! species index

C.................................  Multiply PTMP, SMAT, and CMAT for source

                                SFAC = SMATX( EINX( I ), K ) 

c                                IF( CTRLFLAG ) THEN  ! always false for mobile
c                                    CFAC = CC( EINX( I ), V )
c                                ELSE
                                    CFAC = 1.0
c                                ENDIF
                             
c                                IF( CATEGORY .NE. 'MOBILE' ) THEN
                                    EOUT( I,S,T ) = EMIST( EINX(I) ) * 
     &                                              SFAC * CFAC
c                                ELSE
c                                    DO L = 1, NVTYPE
c                                      EOUT( I,S,T )=EOUT( I,S,T ) +
c     &                                              TMAT( EINX(I),L,1 )*
c     &                                              SFAC * CFAC
c                                    END DO
c                                ENDIF

                            END DO

                        END IF

c                        IF( CATEGORY .NE. 'MOBILE' ) THEN

                            ETMP( I,V,T ) = EMIST( EINX( I ) )

c                        ELSE
c                            DO L = 1, NVTYPE
c                                ETMP( I,J,T ) = ETMP( I,J,T ) +
c     &                                          TMAT( EINX( I ),L,1 )
c                            END DO
c                        ENDIF

                    END DO

                ELSE ! Read failed

                    MESG = 'Could not read "' //
     &                     EINAM( V )( 1:LEN_TRIM( EINAM(V) ) ) //
     &                     '" from file ' // TNAME
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 ) 

                ENDIF 

            END DO

            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO

C.........  Prompt for output report file name

        ODEV = PROMPTFFILE( 'Enter name for ASCII OUTPUT REPORT',
     &                      .FALSE., .TRUE., 'LISTOUT', PROGNAME )
 
C.........  Write header to give file a title

        L  = LEN_TRIM( MENULST( METHOD ) )
        SCRBUF = GRDNM
        WRITE( ODEV, 93200 ) MENULST( METHOD )(1:L), 
     &                       SCRBUF( 1:LEN_TRIM( SCRBUF ) )

        IF( METHOD .EQ. 1 ) THEN
            WRITE( ODEV, 93210 ) XN, YN

        ELSEIF( METHOD .EQ. 2 .AND. AFLAG ) THEN
            WRITE( ODEV, 93220 ) FIP, KEY2, KEY3

        ELSEIF( METHOD .EQ. 2 .AND. MFLAG ) THEN
            WRITE( ODEV, 93223 ) FIP, KEY2, KEY3

        ELSEIF( METHOD .EQ. 2 .AND. PFLAG ) THEN
            WRITE( ODEV, 93225 ) FIP, KEY2, KEY3, KEY4, KEY5

        ELSEIF( METHOD .EQ. 4 ) THEN
            L1 = LEN_TRIM( RNGMENU( MRANGE ) )
            WRITE( ODEV, 93230 ) EINAM( EMMTHD ), 
     &                           RNGMENU( MRANGE )(1:L1), NRPRT

        ELSEIF( METHOD .EQ. 5 ) THEN
            L1 = LEN_TRIM( RNGMENU( MRANGE ) )
            L2 = LEN_TRIM( STKMENU( STKMTHD ) )
            WRITE( ODEV, 93240 ) STKMENU( STKMTHD )(1:L2), 
     &                           RNGMENU( MRANGE )(1:L1), NRPRT

        ENDIF

        WRITE( ODEV, * ) ' '

        TONSUNIT = '[tons]  '
        MOLEUNIT = '[gm mol]'

C.........  Loop over list of sources and output details

        DO I = 1, NOUT

C.............  Set source number from sorted index of output sources
            S = EINX( I )  

C.............  Write out all time INdependant data
C.............  For non-point sources...
! NOTE: This has not been updated for new SMOKE
            IF( .NOT. PFLAG ) THEN
c                IF( AFLAG ) THEN
c                    WRITE( ODEV, 93002 ) 
c     &                     I, S, IFIP( S ), VAR2( S ), VAR3( S ),
c     &                     INVYR( S ), TZONES( S ), TPFLAG( S )

c                ELSEIF( VAR3( S ) .LE. 0 ) THEN ! mobile non-link
c                    WRITE( ODEV, 93003 ) 
c     &                     I, S, IFIP( S ), VAR2( S ),
c     &                     INVYR( S ), TZONES( S ), TPFLAG( S )

c                ELSE ! mobile link
c                    WRITE( ODEV, 93004 ) 
c     &                     I, S, IFIP( S ), VAR2( S ), VAR3( S ),
c     &                     INVYR( S ), TZONES( S ), TPFLAG( S ),
c     &                     XLOC1( S ), YLOC1( S ), 
c     &                     XLOC2( S ), YLOC2( S ), UZONE
c                ENDIF

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

C.............  For point sources...
            ELSE

C.................  Build plant + characteristic list
                CALL PARSCSRC( CSOURC( S ), NCHARS, SC_BEGP, SC_ENDP, 
     &                         LCHR, NCHAR, CHARS )

                BUFFER = ' PLT: ' // CHARS( 2 )
                DO J = 3, NCHAR
                    L1 = LEN_TRIM( BUFFER )
                    L2 = LEN_TRIM( CHARS( J ) )
                    WRITE( BUFFER, '(A,I1,A)' ) BUFFER( 1:L1 ) // 
     &                     ' CHAR', J-2, ': ' // CHARS( J )( 1:L2 )
                END DO

                L1 = LEN_TRIM( CPDESC( S ) )
                L2 = LEN_TRIM( BUFFER )
                WRITE( ODEV, 93000 )
     &             I, S, CPDESC( S )( 1:L1 ),
     &             IFIP( S ), ISIC( S ), CSCC( S ), CORIS( S ),
     &             INVYR( S ), TZONES( S ), TPFLAG( S ), BUFFER( 1:L2 ),
     &             STKHT( S ), STKDM( S ), STKTK( S ), STKVE ( S ),
     &             STKFL( S ), XNUM( I,1 ), YNUM( I,1 ), GSTATE( I ),
     &             XLOCA( S ), YLOCA( S )

            ENDIF

C.............  Write out control coefficients
C NOTE: this is set up for multiplicative controls only at this time
            IF( UFLAG ) THEN

                WRITE( ODEV, 93010 ) 
     &               ( CTLINV( V )( 1:LEN_TRIM( CTLINV( V ) ) ), 
     &                 CC( S, V ), V = 1, NCINVP )

            END IF

! Insert here:

C.............  Convert inventory emissions based on type of input temporal
C.............  Write out time independent inventory emissions
C.............  Week-normal data - input was on day-basis
            IF ( .NOT. MFLAG .AND.
     &           MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 .OR.
     &           MOD( TPFLAG( S ), WDTPFAC ) .EQ. 0      ) THEN

                DO V = 1, NIPOL
                    EMISV( S,V ) = EMISV( S,V ) * YR2DAY( INVYR( S ) )
                END DO

                WRITE( FMTBUF, 93020 ) NIPOL
                WRITE( ODEV, FMTBUF ) 
     &               ( EINAM( V )( 1:LEN_TRIM( EINAM( V ) ) ),
     &                 EMISV( S,V ), V = 1, NIPOL ), 'tons/day'

            ELSEIF( .NOT. MFLAG ) THEN  

                WRITE( FMTBUF, 93020 ) NIPOL
                WRITE( ODEV, FMTBUF ) 
     &               ( EINAM( V )( 1:LEN_TRIM( EINAM( V ) ) ),
     &                 EMISV( S,V ), V = 1, NIPOL ), 'tons/yr'

            ELSE

                WRITE( FMTBUF, 93023 ) NIPOL
                WRITE( ODEV,FMTBUF ) 
     &               ( EINAM( V )( 1:LEN_TRIM( EINAM( V ) ) ),
     &                 EMISV( S,V ), V = 1, NIPOL ), 'tons/period'

c                WRITE( ODEV, 93025 ) 'mi/day', VMT( S ), 
c     &               ( VTYPE3( L   ), L=1,NVTYPE ),
c     &               ( VMTDIS( S,L ), L=1,NVTYPE )

            ENDIF

C.............  For speciation...
            IF( SFLAG ) THEN

C.................  Average mobile speciation matrix over all processes for all
C.................     nehicle types
C.................  NOTE: It is a bit weird to average the speciation factors
C.................        because they do not all get equal weighting when 
C.................        actually applied (weight depends on emissions from 
C.................        each process).  Leave for future improvement.
                IF( MFLAG ) THEN

! NOTE: Removed mobile-specific section

C.................  Write out speciation coefficients for area/point
                ELSE

!  NOTE: Redo format for 80-columns wide no matter what
                    WRITE( ODEV,93035 ) 
     &                   ( EMNAM( F ), SMATX( S,F ), F = 1, NSMATV )

                END IF

            END IF

C.............  If the source is outside the domain, don't output time-stepped
            IF( GSTATE( I ) .NE. 'OUTSIDE' ) THEN

C.................  Write out header for temporal emissions
                IF ( TFLAG .AND. SFLAG ) THEN

                    WRITE( FMTBUF, 93040 ) NIPOL, NMSPC, NIPOL+NMSPC
                    WRITE( ODEV, FMTBUF ) MMDDYY( SDATE ),
     &                           ( EINAM( V ), V = 1, NIPOL ), 
     &                           ( EMNAM( J ), J = 1, NMSPC ),
     &                           ( TONSUNIT  , K = 1, NIPOL ),
     &                           ( MOLEUNIT  , L = 1, NMSPC )

! NOTE: 93040 and 93042 (maybe others) use A3 to print out EINAM

                ELSEIF( TFLAG ) THEN
                    WRITE( FMTBUF, 93042 ) NIPOL, NIPOL
                    WRITE( ODEV, FMTBUF ) MMDDYY( SDATE ),
     &                           ( EINAM( V )        , V = 1, NIPOL ), 
     &                           ( TONSUNIT          , K = 1, NIPOL )
                ENDIF

C.................  Write out time dependent emissions
                JDATE = SDATE
                JTIME = STIME
                LDATE = -9
                DO T = 1, NSTEPS  ! NSTEPS could be zero
   
                    IF( SFLAG ) THEN
                        WRITE( FMTBUF, 93050 ) NIPOL+NMSPC
                        WRITE( ODEV, FMTBUF ) 
     &                         JDATE-SDATE+1, JTIME/10000, 
     &                         ( ETMP( I, V, T ), V=1,NIPOL ),
     &                         ( EOUT( I, J, T ), J=1,NMSPC )
                    ELSE
                        WRITE( FMTBUF, 93050 ) NIPOL
                        WRITE( ODEV, FMTBUF )
     &                         JDATE-SDATE+1, JTIME/10000, 
     &                         ( ETMP( I, V, T ), V=1,NIPOL )
                    END IF

                    IF( JDATE .NE. LDATE ) THEN

                        DO V = 1, NIPOL     ! Initialize
                            EISUM( V ) = ETMP( I, V, T )
                        END DO

                        IF ( SFLAG ) THEN
                            DO J = 1, NMSPC ! Initialize
                                EMSUM( J ) = EOUT( I, J, T )
                            END DO
                        ENDIF

                        LDATE = JDATE

                    ELSE

                        DO V = 1, NIPOL
                            EISUM( V ) = EISUM( V ) + ETMP( I,V,T )
                        END DO

                        IF ( SFLAG ) THEN
                            DO J = 1, NMSPC
                                EMSUM( J ) = EMSUM( J ) + EOUT( I,J,T )
                            END DO
                        ENDIF

                    ENDIF

                    CALL NEXTIME( JDATE, JTIME, TSTEP )

                    IF( JDATE .NE. LDATE .OR. T .EQ. NSTEPS ) THEN
                      IF( SFLAG ) THEN
                        WRITE( FMTBUF, 93055 ) NIPOL + NMSPC
                        WRITE( ODEV,FMTBUF ) ( EISUM( V ), V=1, NIPOL ),
     &                                      ( EMSUM( J ), J=1, NMSPC )
                      ELSE
                        WRITE( FMTBUF, 93055 ) NIPOL
                        WRITE( ODEV,FMTBUF ) ( EISUM( V ), V=1, NIPOL )
                      ENDIF
                    ENDIF

                END DO  ! End loop on time steps for emissions

C.............  Write out alternative header for temporal
            ELSEIF( TFLAG ) THEN

                WRITE( ODEV, 93044 )

            ENDIF     ! If source is in grid or not

            IF( LFLAG .AND. GSTATE( I ) .NE. 'OUTSIDE' ) THEN

C.................  Write out header for layer fractions
                WRITE( FMTBUF, 93060 ) EMLAYS
                WRITE( ODEV, FMTBUF ) MMDDYY( SDATE ), 
     &                               ( L, L = 1, EMLAYS )

C.................  Write out time dependent layer fractions
C.................  Format will put blanks where there are zeroes
                JDATE = SDATE
                JTIME = STIME
                DO T = 1, NSTEPS

                    WRITE( ODEV, 93070 ) 
     &                     JDATE-SDATE+1, JTIME/10000

                    J = 0
                    FMTBUF = '('
                    DO L = 1, EMLAYS

                        VAL = OFRAC( I, L, T ) 
                        L2  = LEN_TRIM( FMTBUF )
                        IF( VAL .GT. 0. ) THEN
                            WRITE( FMTBUF, 93080 ) FMTBUF( 1:L2 )
                            J = J + 1
                            OUTL( J ) = VAL
                        ELSE
                            WRITE( FMTBUF, 93090 ) FMTBUF( 1:L2 )
                        ENDIF

                    END DO

                    L2   = LEN_TRIM( FMTBUF )
                    WRITE( FMTBUF, '(A)' ) FMTBUF( 1:L2 ) // '1X)'

                    WRITE( ODEV,FMTBUF ) ( OUTL( K ), K=1, J )

                    CALL NEXTIME( JDATE, JTIME, TSTEP )

                END DO   ! End loop on time steps for layer fractions

            ENDIF

            WRITE( ODEV, '(A)' ) ' '

        END DO   ! End loop on output sources

C.........  Normal completion
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Informational (LOG) message formats... 92xxx
 
92000   FORMAT( 5X, A )

92100   FORMAT( 10X, A, / 10X, A, I3, A, I3, A, I5 / 10X, A, I5)
 
C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( I5, ' Src ID:', I7, 1X, A, / ,
     &          5X, ' FIP:', I6.6, ' SIC:', I4.4, ' SCC:', A10,
     &              ' ORIS:', A, ' YR:' , I4, ' ZON:', I2.2, 
     &              ' TPF:', I1, /,
     &          5X, A, /, 
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
     &              ' Cell (X,Y) / gridding factor', / , 500
     &              ( 18X, 3( '(', I3,',', I3, ') / ',E11.5, 3X, : )/ ))

93010   FORMAT( 5X, ' Control facs...  ',
     &              30( 1X, A, ': ', E11.5 ) )

93020   FORMAT( '( 5X, " Inv emissions...", ', I4,
     &              '( 1X, A, ": ", E11.5 ), " [", A, "]" )' )

93023   FORMAT( '( 5X, " Total emissions...", ', I4,
     &              '( 1X, A, ": ", E11.5 ), " [", A, "]" )' )

c93025   FORMAT( 5X, ' Inventory VMT [', A, ']...' / 
c     &          7X, ' Total:', 1X, E11.5, /, 
c     &          15X, <NVTYPE>( 4X, A5, :, 3X ), /,
c     &          7X, ' Disag:' <NVTYPE>( 1X, E11.5, : ) )

c93030   FORMAT( 5X, ' Speciation coefficients, VOC values are means',
c     &              ' over exhaust, evap, diurnal processes',
c     &              ' [gm mol/tons]...', /,
c     &              5X, '   VTYPE: ', <NVTYPE>( 4X, A5, 3X ) )

93035   FORMAT( 5X, ' Speciation coefficients [gm mol/tons]...', /,
     &          5( 9X, 7( 1X, A4, ': ', E11.5 ), / ) )

93040   FORMAT( '( 5X, " Hourly emissions starting ", A14, "..." /,',
     &          '5X, " Dy Hr",', I4, '(1X, A3, 5X), ', I4, 
     &          '( 1X, A8 ), /,11X, ', I4, '( 1X, A8 ) ) ' )

93042   FORMAT( '( 5X, " Hourly emissions starting ", A14, "..." /,',
     &          '5X, " Dy Hr",', I4, '(1X, A3, 5X), /, 11X,', I4,
     &          '( 1X, A8 ) )' )

93044   FORMAT( 5X, ' Hourly emissions not written because',   
     &              ' source is outside domain.' )

93050   FORMAT( '( 6X, I2.2, 1X, I2.2, ', I4, '( 1X, E8.3 ) )' )

93055   FORMAT( '( 5X, " Total", ', I4, '( 1X, E8.3 ) )' )

93060   FORMAT( '( 5X, " Layer fractions starting ", A14, "..." /,',
     &          '5X, " Dy Hr", ', I4, '( 1X, " Lyr(", I2.2, ")" ) )' )

93070   FORMAT( 6X, I2.2, 1X, I2.2, $ )

93080   FORMAT( A, '1X, E8.3,' )

93090   FORMAT( A, '9X,' )

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

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal subprogram tries to retrieve the I/O API header
C               and aborts if it was not successful
            SUBROUTINE RETRIEVE_IOAPI_HEADER( FILNAM )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM

C----------------------------------------------------------------------

            IF ( .NOT. DESC3( FILNAM ) ) THEN

                MESG = 'Could not get description of file "' //
     &                 FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF
 
            END SUBROUTINE RETRIEVE_IOAPI_HEADER

        END PROGRAM GETRECS
