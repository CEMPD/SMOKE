
        SUBROUTINE OPENREPIN( ENAME, ANAME, GNAME, LNAME, SLNAME, 
     &                        SSNAME, TNAME, SDEV, GDEV, PDEV, TDEV, 
     &                        EDEV, YDEV, NDEV )

C***********************************************************************
C  subroutine OPENREPIN body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to open all of the necessary
C      files for the Smkreport routine and set the episode information 
C      for the calling program.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 7/2000 by M. Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains the temporal profile tables
        USE MODTMPRL

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

C...........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        CHARACTER*50    GETCFDSC  
        INTEGER         GETIFDSC  
        INTEGER         PROMPTFFILE  
        CHARACTER*16    PROMPTMFILE  
        INTEGER         SECSDIFF  

        EXTERNAL  CRLF, GETCFDSC, GETIFDSC, PROMPTFFILE, 
     &            PROMPTMFILE, SECSDIFF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(OUT) :: ENAME  ! name for I/O API inven input
        CHARACTER(*), INTENT(OUT) :: ANAME  ! name for ASCII inven input 
        CHARACTER(*), INTENT(OUT) :: GNAME  ! gridding matrix name
        CHARACTER(*), INTENT(OUT) :: LNAME  ! layer fractions file name
        CHARACTER(*), INTENT(OUT) :: SLNAME ! speciation matrix name
        CHARACTER(*), INTENT(OUT) :: SSNAME ! speciation matrix name
        CHARACTER(*), INTENT(OUT) :: TNAME  ! hourly emissions file
        INTEGER     , INTENT(OUT) :: SDEV   ! unit no.: ASCII inven file
        INTEGER     , INTENT(OUT) :: GDEV   ! gridding supplemental file
        INTEGER     , INTENT(OUT) :: PDEV   ! speciation supplemental file
        INTEGER     , INTENT(OUT) :: TDEV   ! temporal supplemental file
        INTEGER     , INTENT(OUT) :: EDEV   ! unit no.: elevated ID file (PELV)
        INTEGER     , INTENT(OUT) :: YDEV   ! unit no.: cy/st/co file
        INTEGER     , INTENT(OUT) :: NDEV   ! unit no.: SCC descriptions

C.........  Temporary array for speciation variable names
        CHARACTER(LEN=IODLEN3) SLVNAMS( MXVARS3 )

C.........  Local units and logical file names
        INTEGER      :: MDEV = 0     ! unit no. emission processes file

C.........  Other local variables

        INTEGER         I, J, L, L1, L2, N, V       ! counters and indices

        INTEGER         IOS           ! tmp I/O status

        LOGICAL      :: EFLAG = .FALSE.  ! true: error found

        CHARACTER*16    NAMBUF       ! tmp file name buffer
        CHARACTER*300   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'OPENREPIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENREPIN

C.........  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open input I/O API and ASCII files
C.........  Use NAMBUF for using on the HP
        NAMBUF = PROMPTMFILE( 
     &          'Enter logical name for the I/O API INVENTORY file',
     &          FSREAD3, ENAME, PROGNAME )
        ENAME = NAMBUF

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for the ASCII INVENTORY file',
     &           .TRUE., .TRUE., ANAME, PROGNAME )

C.........  Get source category information from the inventory files
C.........  Get header description of inventory file
C.........  Exit if getting the description fails
        IF( .NOT. DESC3( ENAME ) ) THEN

            L = LEN_TRIM( ENAME )
            MESG = 'Could not get description of file "' //
     &             ENAME( 1:L ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Otherwise, store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API head info is passed by include file and the
C           results are stored in module MODINFO.
        ELSE

            CALL GETSINFO

C.............  Store non-category-specific header information
            NSRC   = NROWS3D
            TSTEP  = 000000
            NSTEPS = 1

        END IF        

C.........  Reset the maximum input data if any reports did not select
C           specific data values.  MXINDAT might get larger than needed.
        IF( DATAMISS ) THEN
            MXINDAT = MAX( NIPPA, MXINDAT )
        END IF

C.........  Determine the year and projection status of the inventory
c       CALL CHECK_INVYEAR( ENAME, APRJFLAG, FDESC3D )

C.........  For temporal inputs, prompt for hourly file
        IF( TFLAG ) THEN

            MESG = 'Enter logical name for the HOURLY ' //
     &             'EMISSIONS file'
            TNAME = PROMPTMFILE( MESG, FSREAD3, CRL//'TMP', PROGNAME )

C.............  Set parameters and pollutants from hourly file
            CALL RETRIEVE_IOAPI_HEADER( TNAME )
            CALL CHKSRCNO( CATDESC, TNAME, NROWS3D, NSRC, EFLAG )
            CALL UPDATE_TIME_INFO( TNAME )

C.............  Determine ozone-season emissions status from hourly file
            INVPIDX = GETIFDSC( FDESC3D, '/OZONE SEASON/', .FALSE. )
            IF( INVPIDX .EQ. 1 ) THEN
                MESG = 'NOTE: Ozone-season emissions in hourly ' //
     &                 'emissions file'
                CALL M3MSG2( MESG )
            END IF

C.............  Store variable number, names, and units from the hourly 
C               emissions file
            NTPDAT = NVARS3D
            ALLOCATE( TPNAME( NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TPNAME', PROGNAME )
            ALLOCATE( TPUNIT( NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TPUNIT', PROGNAME )
            ALLOCATE( TPDESC( NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TPDESC', PROGNAME )

            TPNAME = VNAME3D( 1:NTPDAT )  ! array
            TPUNIT = UNITS3D( 1:NTPDAT )  ! array
            TPDESC = VDESC3D( 1:NTPDAT )  ! array

C.............  Determine the year and projection status of the hourly
c           CALL CHECK_INVYEAR( TNAME, PRJFLAG, FDESC3D )

        END IF

        IF( TSFLAG ) THEN
    
            MESG = 'Enter logical name for the TEMPORAL '//
     &             'SUPPLEMENTAL file'
            TDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                          CRL//'TSUP', PROGNAME )

        END IF

C.........  Open gridding matrix and compare number of sources
        IF( GFLAG ) THEN

            GNAME = PROMPTMFILE( 
     &         'Enter logical name for the GRIDDING MATRIX',
     &         FSREAD3, CRL//'GMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( GNAME )
            CALL CHKSRCNO( CATDESC, GNAME, NTHIK3D, NSRC, EFLAG )

C.............  Initialize grid description
            CALL CHKGRID( CATDESC, 'GMAT', 0, EFLAG )

C.............  Store gridding matrix size
            NMATX = NCOLS3D

        END IF

        IF( GSFLAG ) THEN
    
            MESG = 'Enter logical name for the GRIDDING SUPPLEMENTAL '//
     &             'file'
            GDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                          CRL//'GSUP', PROGNAME )

        END IF

C.........  Open mole speciation matrix, compare number of sources, store
C           speciation variable descriptions, and store mass or moles.
        IF( SLFLAG ) THEN

            SLNAME = PROMPTMFILE( 
     &           'Enter logical name for the MOLE SPECIATION MATRIX',
     &           FSREAD3, CRL//'SMAT_L', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( SLNAME )
            CALL CHKSRCNO( CATDESC, SLNAME, NROWS3D, NSRC, EFLAG )

            NSVARS  = NVARS3D
            SLVNAMS = VDESC3D  ! array

        END IF  ! end of mole speciation open

        IF( PSFLAG ) THEN
    
            MESG = 'Enter logical name for the SPECIATION '//
     &             'SUPPLEMENTAL file'
            PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                          CRL//'SSUP', PROGNAME )

        END IF

C.........  Open mass speciation matrix, compare number of sources, store
C           speciation variable descriptions, and store mass or moles.
        IF( SSFLAG ) THEN

            SSNAME = PROMPTMFILE( 
     &           'Enter logical name for the MASS SPECIATION MATRIX',
     &           FSREAD3, CRL//'SMAT_S', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( SSNAME )
            CALL CHKSRCNO( CATDESC, SSNAME, NROWS3D, NSRC, EFLAG )

C.............  Compare matrix header with mole-based, if available
            IF( SLFLAG ) THEN

C.................  Check the number of variables
                IF( NSVARS .NE. NVARS3D ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Inconsistent number '//
     &                'of speciation variables.'// CRLF()// BLANK10// 
     &                'Mole file:', NSVARS,'; Mass file:', NVARS3D
                    CALL M3MSG2( MESG )

                END IF

C.................  Check the dscriptions of variables
                DO V = 1, NSVARS

                    IF( SLVNAMS( V ) .NE. VDESC3D( V ) ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Inconsistent '//
     &                    'variable descriptions in speciation '//
     &                    'matrices for variable', V, CRLF() // 
     &                    BLANK10 //'Mole file: ', SLVNAMS( V ) //
     &                    CRLF() // BLANK10 //'Mass file: ', VDESC3D(V)
                        CALL M3MSG2( MESG )

                    END IF

                END DO

C.............  Otherwise, set number of speciation variables
            ELSE

        	NSVARS  = NVARS3D

            END IF 

        END IF  ! end of mass speciation open

C.............  Open multiplicative control matrix, compare number of sources, 
C               and store control variable names.
c        IF( CUFLAG ) THEN

c            MESG = 'Enter logical name for the ' //
c     &             'MULTIPLICATIVE CONTROL MATRIX'
c            UNAME = PROMPTMFILE( MESG, FSREAD3, CRL//'CMAT', PROGNAME )

c            CALL RETRIEVE_IOAPI_HEADER( UNAME )
c            CALL CHKSRCNO( 'area', UNAME, NROWS3D, NASRC, EFLAG )
c            NUMATV = NVARS3D
c            ALLOCATE( UVNAMS( NUMATV ), STAT=IOS )
c            CALL CHECKMEM( IOS, 'UVNAMS', PROGNAME )
c            CALL STORE_VNAMES( 1, 1, NUMATV, UVNAMS )

c        END IF  ! end of multiplicative control open

C.............  Open additive control matrix, compare number of sources, 
C               and store control variable names.
c        IF( CAFLAG ) THEN
c            MESG= 'INTERNAL ERROR: Area additive controls not ' //
c     &            'yet implemented in ' // PROGNAME
c            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

c        END IF  ! end of additive control open

C.............  Open reactivity control matrix, compare number of sources, and
C               store control variable descriptions, and store mass or moles.
c        IF( CRFLAG ) THEN
c            RNAME = PROMPTMFILE( 
c     &           'Enter logical name for the REACTIVITY MATRIX',
c     &           FSREAD3, CRL//'RMAT', PROGNAME )

c            CALL RETRIEVE_IOAPI_HEADER( RNAME )
c            CALL CHKSRCNO( CATDESC, RNAME, NTHIK3D, NSRC, EFLAG )
c            NRMATV = NVARS3D
c            NSREAC = NROWS3D
c            ALLOCATE( RVDESC( NRMATV ), STAT=IOS )
c            CALL CHECKMEM( IOS, 'RVDESC', PROGNAME )
c            CALL STORE_VDESCS( 1, 1, NRMATV, RVDESC )

C.................  Retrieve the number of speciation factors 
c            RNMSPC = GETIFDSC( FDESC3D, '/SPECIES VARS/', .TRUE. )

C.................  Check the year and projection year of the matrix
c           CALL CHECK_INVYEAR( ARNAME, APRJFLAG, FDESC3D )

c        END IF  ! end of reactivity control open

C.........  Open layer fractions file, compare number of sources, check 
C           met information, and store the vertical coordinates info
        IF( LFLAG ) THEN

            MESG= 'Enter logical name for the POINT LAYER ' //
     &            'FRACTIONS MATRIX'
            LNAME = PROMPTMFILE( MESG, FSREAD3, 'PLAY', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( LNAME )
            CALL CHKSRCNO( CATDESC, LNAME, NROWS3D, NSRC, EFLAG )
            CALL UPDATE_TIME_INFO( LNAME )
            EMLAYS = NLAYS3D

        END IF  ! End of layer fractions open

C.........  Open elevated/low-level 
        IF( VFLAG ) THEN

C.............  Open elevated/plume-in-grid file 
            MESG = 'Enter logical name for the ELEVATED/PING file'
            EDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                          'PELV', PROGNAME      )

        END IF

C.........  Get country, state, and county names, if needed
        IF( YFLAG ) THEN

            MESG = 'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file'
            YDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'COSTCY',PROGNAME )

        END IF

C.........  Get SCC descriptions, if needed
        IF( NFLAG ) THEN

            MESG = 'Enter logical name for SCC DESCRIPTIONS'
            NDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'SCCDESC',PROGNAME )

        END IF

C.........  If there were any errors inputing files or while comparing
C           with one another, then abort
        IF( EFLAG ) THEN

           MESG = 'Problems opening input files. See ERROR(S) above.'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  If we are using temporalized emissions, then update date/time and
C           duration using environment variable settings, then prompt.
        IF( TFLAG ) THEN

C.............  Write explanation 
            MESG = 'For time-based reports, enter the starting date, '//
     &             'starting time, and output' // CRLF() // BLANK10 //
     &             'duration.  Defaults are set based on the input ' //
     &             'files.'
            CALL M3MSG2( MESG )

C.............  Subselect dates and times
            CALL GETM3EPI( TZONE, SDATE, STIME, NSTEPS )
            EDATE = SDATE
            ETIME = STIME
            CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

        END IF   !  if have temporalized inputs and outputs

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

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

C----------------------------------------------------------------------
C----------------------------------------------------------------------

C.............  This subprogram updates the time (episode) information
C               and compares to the existing information, if it has been
C               previously set.
            SUBROUTINE UPDATE_TIME_INFO( FILNAM )

C.............  Subprogram arguments
            CHARACTER(*) FILNAM

C.............  Local variables
            INTEGER ISECS   ! number of seconds different between dates/times
            INTEGER ED      ! tmp ending date
            INTEGER ET      ! tmp ending time            
            INTEGER LOCZONE ! tmp time zone

            INTEGER, SAVE :: EDATE = 0    ! Ending date
            INTEGER, SAVE :: ETIME = 0    ! Ending time

            LOGICAL, SAVE :: IFLAG = .FALSE.  ! true: episode settings init
            LOGICAL, SAVE :: ZFLAG = .FALSE.  ! true: time zone init

            CHARACTER*300    MESG             ! message buffer

C----------------------------------------------------------------------

C.............  If time information has already been initialized...
            IF( IFLAG ) THEN
                ISECS = SECSDIFF( SDATE, STIME, SDATE3D, STIME3D )

                IF( ISECS .GT. 0 ) THEN  ! SDATE3D/STIME3D are later
                    SDATE = SDATE3D
                    STIME = STIME3D
                END IF

                ED = SDATE3D
                ET = STIME3D
                CALL NEXTIME( ED, ET, ( MXREC3D-1 ) * TSTEP3D )
        
                ISECS = SECSDIFF( EDATE, ETIME, ED, ET )

                IF( ISECS .LT. 0 ) THEN  ! ED/ET are earlier
                    EDATE = ED
                    ETIME = ET
                END IF

                NSTEPS = 1+ SECSDIFF( SDATE, STIME, EDATE, ETIME )/ 3600

                IF( TFLAG .AND. NSTEPS .LE. 0 ) THEN
                    MESG = 'Because of file ' // FILNAM // 
     &                     ', dates and times do not overlap at all!'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C.................  Check time step 
                IF( TSTEP3D .NE. TSTEP ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Time step is not one hour in ' // 
     &                     FILNAM // ' file!'
                    CALL M3MSG2( MESG )
                END IF

C.............  If time information needs to be initialized...
            ELSE
                SDATE  = SDATE3D
                STIME  = STIME3D
                NSTEPS = MXREC3D
                TSTEP  = TSTEP3D

                EDATE  = SDATE
                ETIME  = STIME
                CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

                IFLAG = .TRUE.

            END IF

C.............  Make sure that time step is one hour
            IF( TSTEP .NE. 10000 ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: Time step is not one hour in ' // 
     &                 FILNAM // ' file!'
                CALL M3MSG2( MESG )

            END IF

C.............  Use layers to screen for non-layer-fractions files
C.............  If not layer-fractions file, retrieve and compare time zone
            IF( NLAYS3D .LE. 1 ) THEN

                LOCZONE = GETIFDSC( FDESC3D, '/TZONE/', .TRUE. )

                IF( ZFLAG .AND. LOCZONE .NE. TZONE ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                 'Time zone ', LOCZONE, 'in ' // FILNAM // 
     &                 ' hourly emissions file is not consistent ' //
     &                 'with initialized value of', TZONE
                    CALL M3MSG2( MESG )

                ELSE IF( .NOT. ZFLAG ) THEN
                    ZFLAG = .TRUE.
                    TZONE = LOCZONE

                    MESG = 'NOTE: Time zone initialized using ' // 
     &                     FILNAM // ' hourly emissions file.'

                    CALL M3MSG2( MESG )
                END IF
            END IF

C------------------  FORMAT  STATEMENTS   -----------------------------

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE UPDATE_TIME_INFO

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram initializes and checks the inventory year
C               of the emissions and the projection status
            SUBROUTINE CHECK_INVYEAR( FNAME, PRJFLAG, IODESC )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN)     :: FNAME
            LOGICAL     , INTENT (IN OUT) :: PRJFLAG
            CHARACTER(*), INTENT (IN)     :: IODESC( * )

C.............  Local variables
            INTEGER           YY      ! tmp year

            INTEGER , SAVE :: BYEAR   ! base year
            INTEGER , SAVE :: FLEN    ! name length of savnam

            LOGICAL           STRICT  ! flag for strict checks or not

            CHARACTER*20      BUFFER  ! program name buffer
            CHARACTER(LEN=IOVLEN3), SAVE :: SAVNAM  ! name of file used to init

C----------------------------------------------------------------------

            STRICT = .TRUE.

C.............  First determine whether to abort when projected year does not
C               match.  This is used for reactivity matrices, which will
C               always have a projection year, even if the inventory isn't
C               projected.
            IF( .NOT. PRJFLAG ) THEN
                BUFFER = GETCFDSC( FDESC3D, '/FROM/', .FALSE. )
                IF( BUFFER .EQ. 'OPENRMAT' ) STRICT = .FALSE.
            END IF

C.............  If time information has already been initialized...
            IF( YFLAG ) THEN

                YY = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )
                IF( YY .LE. 0 ) THEN

                    YY = GETIFDSC( IODESC, '/BASE YEAR/', .FALSE. ) 
                    IF( YY .NE. BYEAR ) THEN
                        WRITE( MESG,94010 ) 
     &                        'Base year of ' // FNAME // ' file:', YY,
     &                        CRLF() // BLANK10 //
     &                        ', does not equal emissions year of ' //
     &                        SAVNAM( 1:FLEN ) // ' file:', BYEAR
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                ELSE IF ( STRICT            .AND. 
     &                    YY     .GT. 0     .AND. 
     &                    YY     .NE. PYEAR      ) THEN

                    WRITE( MESG,94010 ) 
     &                    'Projected year of ' // FNAME // ' file:', YY,
     &                    CRLF() // BLANK10 //
     &                    ', does not equal emissions year of ' //
     &                    SAVNAM( 1:FLEN ) // ' file:', PYEAR
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C.............  If year information needs to be initialized...
            ELSE
                
                BYEAR = GETIFDSC( IODESC, '/BASE YEAR/', .FALSE. ) 
                PYEAR = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )

                IF( YY .GT. 0 ) THEN
                    PRJFLAG = .TRUE.
                ELSE
                    PYEAR = BYEAR
                END IF

                SAVNAM = FNAME
                FLEN   = LEN_TRIM( SAVNAM )
                YFLAG  = .TRUE.

            END IF

C------------------  FORMAT  STATEMENTS   -----------------------------

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_INVYEAR

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores I/O API NetCDF variable names into
C               a local array based on indices in subprogram call.
            SUBROUTINE STORE_VNAMES( ISTART, INCRMT, NNAM, NAMES )

C.............  Subprogram arguments
            INTEGER      ISTART        ! starting position in VNAMES of names
            INTEGER      INCRMT        ! increment of VNAMES for names
            INTEGER      NNAM          ! number of names
            CHARACTER(*) NAMES( NNAM ) ! stored variable names

C.............  Local variables
            INTEGER  I, J

C----------------------------------------------------------------------

            J = ISTART
            DO I = 1, NNAM

                NAMES( I ) = VNAME3D( J )
                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VNAMES

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores I/O API NetCDF variable descriptions into
C               a local array based on indices in subprogram call.
            SUBROUTINE STORE_VDESCS( ISTART, INCRMT, NDESC, DESCS )

C.............  Subprogram arguments
            INTEGER      ISTART        ! starting position in VDESCS of names
            INTEGER      INCRMT        ! increment of VDESCS for names
            INTEGER      NDESC         ! number of descriptions
            CHARACTER(*) DESCS( NDESC )! stored variable descriptions

C.............  Local variables
            INTEGER  I, J, L

C----------------------------------------------------------------------

            DESCS = ' '

            J = ISTART
            DO I = 1, NDESC

                L = LEN_TRIM( VDESC3D( J ) )
                DESCS( I ) = VDESC3D( J )( 1:L )
                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VDESCS

        END SUBROUTINE OPENREPIN
