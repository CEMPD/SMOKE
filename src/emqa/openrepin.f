
        SUBROUTINE OPENREPIN( ENAME, ANAME, CUNAME, GNAME, LNAME, 
     &                        PRNAME, SLNAME, SSNAME, TNAME, RDEV, 
     &                        SDEV, GDEV, PDEV, TDEV, EDEV, YDEV, NDEV,
     &                        NIDEV, NPDEV, ADEV, NMDEV, NNDEV, NODEV )

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
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: AFLAG, TSTEP, NSTEPS, DATAMISS, MXINDAT,
     &                      TSFLAG, TFLAG, GFLAG, GSFLAG, SLFLAG,
     &                      PSFLAG, SSFLAG, PRFLAG, PRRPTFLG, NMATX,
     &                      PRBYR, PRPYR, PYEAR, CHKPFX, CUFLAG,
     &                      LFLAG, EMLAYS, VFLAG, YFLAG, NFLAG,
     &                      ASCREC, ASCDATA, STIME, SDATE, ETIME,
     &                      EDATE, TZONE, NIFLAG, NMFLAG, NNFLAG,
     &                      NOFLAG, SDFLAG

C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: NTPDAT, TPNAME, TPUNIT, TPDESC

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NSVARS

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: NVPROJ, PNAMPROJ, NVCMULT, PNAMMULT

C...........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CRL, CATDESC, NIPPA, NCHARS, JSCC,
     &                     JSTACK, PLTIDX, MXCHRS, NSRC, INVPIDX, BYEAR,
     &                     EANAM, EAUNIT, SC_BEGP, SC_ENDP

C.........  This module is required for the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        CHARACTER(50)   GETCFDSC  
        INTEGER         GETIFDSC  
        INTEGER         PROMPTFFILE  
        CHARACTER(16)   PROMPTMFILE  
        INTEGER         SECSDIFF  
        LOGICAL         USEEXPGEO

        EXTERNAL  CRLF, GETCFDSC, GETIFDSC, PROMPTFFILE, 
     &            PROMPTMFILE, SECSDIFF, USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(OUT) :: ENAME  ! name for I/O API inven input
        CHARACTER(*), INTENT(OUT) :: ANAME  ! name for ASCII inven input 
        CHARACTER(*), INTENT(OUT) :: CUNAME ! multiplicative control matrix name
        CHARACTER(*), INTENT(OUT) :: GNAME  ! gridding matrix name
        CHARACTER(*), INTENT(OUT) :: LNAME  ! layer fractions file name
        CHARACTER(*), INTENT(OUT) :: PRNAME ! projection matrix name
        CHARACTER(*), INTENT(OUT) :: SLNAME ! speciation matrix name
        CHARACTER(*), INTENT(OUT) :: SSNAME ! speciation matrix name
        CHARACTER(*), INTENT(OUT) :: TNAME  ! hourly emissions file
        INTEGER     , INTENT(OUT) :: RDEV(3)! control report files
        INTEGER     , INTENT(OUT) :: SDEV   ! unit no.: ASCII inven file
        INTEGER     , INTENT(OUT) :: GDEV   ! gridding supplemental file
        INTEGER     , INTENT(OUT) :: PDEV   ! speciation supplemental file
        INTEGER     , INTENT(OUT) :: TDEV   ! temporal supplemental file
        INTEGER     , INTENT(OUT) :: EDEV   ! unit no.: elevated ID file (PELV)
        INTEGER     , INTENT(OUT) :: YDEV   ! unit no.: cy/st/co file
        INTEGER     , INTENT(OUT) :: NDEV   ! unit no.: SCC descriptions
        INTEGER     , INTENT(OUT) :: NIDEV  ! unit no.: SIC descriptions
        INTEGER     , INTENT(OUT) :: NPDEV  ! unit no.: GSPRO descriptions
        INTEGER     , INTENT(OUT) :: NMDEV  ! unit no.: MACT descriptions
        INTEGER     , INTENT(OUT) :: NNDEV  ! unit no.: NAICS descriptions
        INTEGER     , INTENT(OUT) :: NODEV  ! unit no.: ORIS descriptions
        INTEGER     , INTENT(OUT) :: ADEV   ! unit no.: ASCII elevated file

C.........  Temporary array for speciation variable names
        CHARACTER(IODLEN3), ALLOCATABLE :: SLVNAMS( : )

C.........  Local units and logical file names
        INTEGER         IDEV      ! tmp unit number if ENAME is map file
        INTEGER      :: MDEV = 0  ! unit no. emission processes file

        CHARACTER(16)   INAME     ! tmp name for inven file of unknown fmt

C.........  Other local variables

        INTEGER         I, J, L, L1, L2, N, V       ! counters and indices

        INTEGER         IOS           ! tmp I/O status
        INTEGER         ISD           ! start time of ASCII elevated file
        INTEGER         IED           ! end time of ASCII elevated file
        INTEGER         TSTEP_T       ! unused time step from environment

        LOGICAL      :: EFLAG = .FALSE.  ! true: error found
        LOGICAL      :: TIMEFLAG = .FALSE.  ! true: time info already init

        CHARACTER(16)   NAMBUF       ! tmp file name buffer
        CHARACTER(16)   UNITS        ! units of ASCII elevated file
        CHARACTER(256)  MESG         ! message buffer
        CHARACTER(300)  LINE         ! tmp line buffer

        CHARACTER(IOULEN3) GRDENV      ! gridded output units from envrmt

        CHARACTER(16) :: PROGNAME = 'OPENREPIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENREPIN

        IF( .NOT. AFLAG ) THEN
C.........  Get inventory file names given source category
            CALL GETINAME( CATEGORY, ENAME, ANAME )

C.........  Prompt for and open inventory file 
            INAME = ENAME
            MESG = 'Enter logical name for the MAP INVENTORY file'
            IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Open and read map file
            CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

C.........  Store source-category-specific header information, 
C           including the inventory pollutants in the file (if any).  Note that 
C           the I/O API header info is passed by include file and the
C           results are stored in module MODINFO.
            CALL GETSINFO( ENAME )

C.........  Store non-category-specific header information
            TSTEP  = 000000
            NSTEPS = 1

C.........  Reset the maximum input data if any reports did not select
C           specific data values.  MXINDAT might get larger than needed.
            IF( DATAMISS ) THEN
                MXINDAT = MAX( NIPPA, MXINDAT )
            END IF

C.........  Determine the year and projection status of the inventory
c           CALL CHECK_INVYEAR( ENAME, APRJFLAG, FDESC3D )

        ELSE
            NCHARS = 3
            JSCC = 0
            JSTACK = 3

            ALLOCATE( SC_BEGP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SC_BEGP', PROGNAME )
            ALLOCATE( SC_ENDP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SC_ENDP', PROGNAME )

            PLTIDX = 2
            MXCHRS = MXPTCHR3

            DO I = 1, NCHARS
                SC_BEGP( I ) = PTBEGL3( I )
                SC_ENDP( I ) = PTENDL3( I )
            END DO

        END IF

C.........  For temporal inputs, prompt for hourly file
        IF( TFLAG ) THEN

            MESG = 'Enter logical name for the HOURLY ' //
     &             'EMISSIONS file'
            TNAME = PROMPTSET( MESG, FSREAD3, CRL//'TMP', PROGNAME )

C.............  Set parameters and pollutants from hourly file
            CALL RETRIEVE_SET_HEADER( TNAME )
            CALL CHKSRCNO( CATDESC, TNAME, NROWS3D, NSRC, EFLAG )
            CALL UPDATE_TIME_INFO( TNAME )

C.............  Determine average day emissions status from hourly file
            INVPIDX = GETIFDSC( FDESC3D, '/AVERAGE DAY/', .FALSE. )
            IF( INVPIDX .EQ. 1 ) THEN
                MESG = 'NOTE: Average day emissions in hourly ' //
     &                 'emissions file'
                CALL M3MSG2( MESG )
            END IF

C.............  Store variable number, names, and units from the hourly 
C               emissions file
            NTPDAT = NVARSET
            ALLOCATE( TPNAME( NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TPNAME', PROGNAME )
            ALLOCATE( TPUNIT( NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TPUNIT', PROGNAME )
            ALLOCATE( TPDESC( NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TPDESC', PROGNAME )

            TPNAME = VNAMESET( 1:NTPDAT )  ! array
            TPUNIT = VUNITSET( 1:NTPDAT )  ! array
            TPDESC = VDESCSET( 1:NTPDAT )  ! array

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

            SLNAME = PROMPTSET( 
     &           'Enter logical name for the MOLE SPECIATION MATRIX',
     &           FSREAD3, CRL//'SMAT_L', PROGNAME )

            CALL RETRIEVE_SET_HEADER( SLNAME )
            CALL CHKSRCNO( CATDESC, SLNAME, NROWS3D, NSRC, EFLAG )

            NSVARS = NVARSET
            
            ALLOCATE( SLVNAMS( NSVARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SLVNAMS', PROGNAME )
            
            SLVNAMS = VDESCSET  ! array

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

            SSNAME = PROMPTSET( 
     &           'Enter logical name for the MASS SPECIATION MATRIX',
     &           FSREAD3, CRL//'SMAT_S', PROGNAME )

            CALL RETRIEVE_SET_HEADER( SSNAME )
            CALL CHKSRCNO( CATDESC, SSNAME, NROWS3D, NSRC, EFLAG )

C.............  Compare matrix header with mole-based, if available
            IF( SLFLAG ) THEN

C.................  Check the number of variables
                IF( NSVARS .NE. NVARSET ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Inconsistent number '//
     &                'of speciation variables.'// CRLF()// BLANK10// 
     &                'Mole file:', NSVARS,'; Mass file:', NVARSET
                    CALL M3MSG2( MESG )

                END IF

C.................  Check the dscriptions of variables
                DO V = 1, NSVARS

                    IF( SLVNAMS( V ) .NE. VDESCSET( V ) ) THEN

                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Inconsistent '//
     &                    'variable descriptions in speciation '//
     &                    'matrices for variable', V, CRLF() // 
     &                    BLANK10 //'Mole file: ', SLVNAMS( V ) //
     &                    CRLF() // BLANK10 //'Mass file: ', VDESCSET(V)
                        CALL M3MSG2( MESG )

                    END IF

                END DO

C.............  Otherwise, set number of speciation variables
            ELSE

                NSVARS  = NVARSET

            END IF 

        END IF  ! end of mass speciation open

C.............  Open projection matrix, compare number of sources, 
C               and store projection variable names.
        IF( PRFLAG ) THEN

            MESG = 'Enter logical name for the ' //
     &             'PROJECTION MATRIX'
            PRNAME = PROMPTSET( MESG, FSREAD3, CRL//'PMAT', PROGNAME )

            CALL RETRIEVE_SET_HEADER( PRNAME )
            CALL CHKSRCNO( CATDESC, PRNAME, NROWS3D, NSRC, EFLAG )
            NVPROJ = NVARS3D

C...........  Set allocation size depending on whether report is also
C             read in.
            I = NVPROJ
            IF( PRRPTFLG ) I = 2 * I

C...........  Allocate memory for variables (including possible 
C             other variables for report check) and store names
            ALLOCATE( PNAMPROJ( I ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PNAMPROJ', PROGNAME )
            PNAMPROJ = ' '  ! array

            CALL STORE_VNAMES( 1, 1, NVPROJ, PNAMPROJ )

            IF( NVPROJ .NE. 1 ) THEN
                MESG = 'INTERNAL ERROR: Smkreport is not set up to ' //
     &                 'support more than 1 variable in ' //
     &                 CRLF() // BLANK10 // 'the projection matrix.'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            PRBYR = GETIFDSC( FDESC3D, '/BASE YEAR/', .TRUE. )
            PRPYR = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .TRUE. )

            IF( PYEAR .LE. 0 ) THEN
                IF ( PRBYR .NE. BYEAR ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Inventory base'//
     &                     ' year', BYEAR, 'is not consistent '//
     &                     'with projection base year', PRBYR
                    CALL M3MSG2( MESG )

                ELSE
                    WRITE( MESG,94010 ) 'NOTE: Base year', BYEAR, 
     &                     'is consistent between the inventory and '//
     &                     CRLF() // BLANK10 // 'projection matrix.'
                    CALL M3MSG2( MESG )
                END IF

            ELSE 

                IF ( PRBYR .NE. PYEAR ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Inventory projected'//
     &                     ' data year', PYEAR, 'is not consistent '//
     &                     'with projection base year', PRBYR
                    CALL M3MSG2( MESG )

                ELSE
                    WRITE( MESG,94010 ) 'WARNING: Inventory projected'//
     &                     ' year', PYEAR, 'is being projected to',
     &                     PRPYR, 'by projection matrix.'
                    CALL M3MSG2( MESG )

                END IF
            END IF
        END IF  ! end of multiplicative control open

C.............  Open projections report
        IF( PRRPTFLG ) THEN

            MESG = 'Enter logical name for input PROJECTION REPORT '//
     &             'from Cntlmat'
            RDEV(1) = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                             CRL // 'PROJRPT', PROGNAME )

            DO V = NVPROJ+1, 2*NVPROJ
                J = V - NVPROJ
                PNAMPROJ( V ) = CHKPFX // PNAMPROJ( J )
            END DO

        END IF

C.............  Open multiplicative control matrix, compare number of sources, 
C               and store control variable names.
        IF( CUFLAG ) THEN

            MESG = 'Enter logical name for the ' //
     &             'MULTIPLICATIVE CONTROL MATRIX'
            CUNAME = PROMPTSET( MESG, FSREAD3, CRL//'CMAT', PROGNAME )

            CALL RETRIEVE_SET_HEADER( CUNAME )
            CALL CHKSRCNO( CATDESC, CUNAME, NROWS3D, NSRC, EFLAG )
            NVCMULT = NVARS3D
            ALLOCATE( PNAMMULT( NVCMULT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PNAMMULT', PROGNAME )
            CALL STORE_VNAMES( 1, 1, NVCMULT, PNAMMULT )

        END IF  ! end of multiplicative control open

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
        IF( YFLAG .AND. .NOT. USEEXPGEO() ) THEN

            MESG = 'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file'
            YDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'COSTCY',PROGNAME )

        END IF

C.........  Get SCC descriptions, if needed
        IF( NFLAG ) THEN

            MESG = 'Enter logical name for SCC DESCRIPTIONS'
            NDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'SCCDESC',PROGNAME )

        END IF

C.........  Get SIC descriptions, if needed
        IF( NIFLAG ) THEN

            MESG = 'Enter logical name for SIC DESCRIPTIONS'
            NIDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'SICDESC',PROGNAME )

        END IF

C.........  Get GSPRO descriptions, if needed
        IF( SDFLAG ) THEN

            MESG = 'Enter logical name for GSPRO DESCRIPTIONS'
            NPDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'GSPRODESC',PROGNAME )

        END IF

C.........  Get MACT descriptions, if needed
        IF( NMFLAG ) THEN

            MESG = 'Enter logical name for MACT DESCRIPTIONS'
            NMDEV = PROMPTFFILE( 
     &                  MESG,.TRUE.,.TRUE.,'MACTDESC',PROGNAME )

        END IF

C.........  Get NAICS descriptions, if needed
        IF( NNFLAG ) THEN

            MESG = 'Enter logical name for NAICS DESCRIPTIONS'
            NNDEV = PROMPTFFILE( 
     &                  MESG,.TRUE.,.TRUE.,'NAICSDESC',PROGNAME )

        END IF

C.........  Get ORIS descriptions, if needed
        IF( NOFLAG ) THEN

            MESG = 'Enter logical name for ORIS DESCRIPTIONS'
            NODEV = PROMPTFFILE( 
     &                  MESG,.TRUE.,.TRUE.,'ORISDESC',PROGNAME )

        END IF

C.........  Open ASCII elevation file output by SMKMERGE, if needed
        IF( AFLAG ) THEN

            CALL ENVSTR( 'MRG_GRDOUT_UNIT', ' ', ' ', GRDENV, IOS)

            IF( GRDENV( 1:1 ) .EQ. 'm' ) THEN

                ADEV = PROMPTFFILE(
     &              'Enter name for ASCII ELEVATED SOURCES file', 
     &              .TRUE., .TRUE., 'ELEVTS_L', PROGNAME )

            ELSE

                ADEV = PROMPTFFILE(
     &              'Enter name for ASCII ELEVATED SOURCES file',
     &              .TRUE., .TRUE., 'ELEVTS_S', PROGNAME )

            END IF

C.........  Read ASCII elevated file
            MESG = 'Reading ASCII elevated file...'
            CALL M3MSG2( MESG )

            ASCREC = 0 

C............  Read in units
            ASCREC = ASCREC + 1
            READ( ADEV, '(10X,A)') UNITS

C............  Skip header lines
            DO I = 1, 2
                ASCREC = ASCREC + 1
                READ( ADEV, '(A)' ) LINE
            END DO

C.............  Get number of species and point stacks from file
            ASCREC = ASCREC + 1
            READ( ADEV, '(I10,10X,I10)' ) ASCDATA, NSRC

            NIPPA = ASCDATA
            NSVARS = NIPPA
            MXINDAT = MAX( NIPPA, MXINDAT )

            ALLOCATE( EANAM( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
            ALLOCATE( EAUNIT( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EAUNIT', PROGNAME )
            EANAM = ''
            EAUNIT = TRIM( UNITS )

            DO I = 1, 4
                ASCREC = ASCREC + 1
                READ( ADEV, '(A)' ) LINE
            END DO

C..............  Read in list of species
            DO I = 1, ASCDATA
                ASCREC = ASCREC + 1
                READ( ADEV, '(A)' ) LINE
                EANAM( I ) = TRIM( LINE )
            END DO

C..............  Read in start and end dates and times
            ASCREC = ASCREC + 1
            READ( ADEV, '(4I10)' ) ISD, STIME, IED, ETIME
            SDATE = ISD + 1900000
            EDATE = IED + 1900000
            BYEAR = INT( SDATE/1000 )
            TSTEP = 10000
            NSTEPS = 1 + SECSDIFF( SDATE, STIME, EDATE, ETIME )/3600

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
            CALL GETM3EPI( TZONE, SDATE, STIME, TSTEP_T, NSTEPS )
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
            CHARACTER(*), INTENT (IN) :: FILNAM
            
C----------------------------------------------------------------------

            IF ( .NOT. DESC3( FILNAM ) ) THEN

                MESG = 'Could not get description of file "' //
     &                 FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF
 
            END SUBROUTINE RETRIEVE_IOAPI_HEADER

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram tries to retrieve the description for a file
C               set and aborts if it was not successful
            SUBROUTINE RETRIEVE_SET_HEADER( FILNAM )

            INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C.............  Subprogram arguments
            CHARACTER(*) FILNAM

C----------------------------------------------------------------------

            IF ( .NOT. DESCSET( FILNAM, ALLFILES ) ) THEN

                MESG = 'Could not get description of file set "' //
     &                 FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF
 
            END SUBROUTINE RETRIEVE_SET_HEADER

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

            CHARACTER(300)   MESG             ! message buffer

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

            CHARACTER(20)     BUFFER  ! program name buffer
            CHARACTER(IOVLEN3), SAVE :: SAVNAM  ! name of file used to init

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
            IF( TIMEFLAG ) THEN

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
                TIMEFLAG  = .TRUE.

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

                NAMES( I ) = VNAMESET( J )
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

                L = LEN_TRIM( VDESCSET( J ) )
                DESCS( I ) = VDESCSET( J )( 1:L )
                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VDESCS

        END SUBROUTINE OPENREPIN
