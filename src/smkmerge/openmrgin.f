
        SUBROUTINE OPENMRGIN

C***********************************************************************
C  subroutine OPENMRGIN body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to open all of the necessary
C      files for the merge routine and set the episode information 
C      for the calling program.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures
        INCLUDE 'FLTERR.EXT'    !  error filter statement function

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        CHARACTER*50    GETCFDSC  
        INTEGER         GETIFDSC  
        INTEGER         PROMPTFFILE  
        CHARACTER*16    PROMPTMFILE  
        INTEGER         SECSDIFF  

        EXTERNAL  CRLF, GETCFDSC, GETIFDSC, PROMPTFFILE, PROMPTMFILE,
     &            SECSDIFF

C.........  Other local variables

        INTEGER         J, N, V       ! counters and indices

        INTEGER         EDATE         ! episode end date (YYYYDDD)
        INTEGER         ETIME         ! episode end time (HHMMSS)
        INTEGER         IOS           ! tmp I/O status
        INTEGER         ISECS         ! tmp duration in seconds

        LOGICAL      :: CFLAG = .FALSE.  ! true: speciation type has been init
        LOGICAL      :: EFLAG = .FALSE.  ! true: error in routine
        LOGICAL      :: GFLAG = .FALSE.  ! true: grid settings have been init
        LOGICAL      :: IFLAG = .FALSE.  ! true: episode settings have been init
        LOGICAL      :: OFLAG = .FALSE.  ! true: met info has been init
        LOGICAL      :: YFLAG = .FALSE.  ! true: year/projection info been init
        LOGICAL      :: ZFLAG = .FALSE.  ! true: time zone has been init

        CHARACTER*4     SPCTYPE      ! type of speciation matrix (mass|mole)
        CHARACTER*50    METSCENR     ! met scenario name
        CHARACTER*50    METCLOUD     ! met cloud scheme name
        CHARACTER*50    METTMP       ! temporary buffer for met info
        CHARACTER*300   MESG         ! message buffer

        CHARACTER*16 :: PROGNAME = 'OPENMRGIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENMRGIN

C.........  For area sources... 
        IF( AFLAG ) THEN

C.............  Prompt for inventory files
            AENAME = PROMPTMFILE( 
     &       'Enter logical name for the I/O API AREA INVENTORY file',
     &       FSREAD3, 'AREA', PROGNAME )

            ASDEV = PROMPTFFILE( 
     &       'Enter logical name for the ASCII AREA INVENTORY file',
     &       .TRUE., .TRUE., 'ASRC', PROGNAME )

C.............  Get number of sources
            CALL RETRIEVE_IOAPI_HEADER( AENAME )
            NASRC = NROWS3D

C.............  Determine the year and projection status of the inventory
            CALL CHECK_INVYEAR( PENAME, APRJFLAG, FDESC3D )

C.............  For temporal inputs, prompt for hourly file
            IF( TFLAG ) THEN

                MESG = 'Enter logical name for the AREA HOURLY ' //
     &                 'EMISSIONS file'
                ATNAME = PROMPTMFILE( MESG, FSREAD3, 'ATMP', PROGNAME )

C.................  Set parameters and pollutants from hourly file
                CALL RETRIEVE_IOAPI_HEADER( ATNAME )
                CALL CHKSRCNO( 'area', 'ATMP', NROWS3D, NASRC, EFLAG )
                CALL UPDATE_TIME_INFO( 'ATMP' )
                INVUNIT = UNITS3D( 1 )

                ANIPOL = NVARS3D
                ALLOCATE( AEINAM( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AEINAM', PROGNAME )
                CALL STORE_VNAMES( 1, 1, ANIPOL, AEINAM )

C.................  Determine the year and projection status of the hourly
                CALL CHECK_INVYEAR( ATNAME, APRJFLAG, FDESC3D )

C.............  Otherwise, just set parameters and pollutants from inven file
            ELSE
                ATNAME = AENAME
                ANIPOL = ( NVARS3D - NARVAR3 ) / NARPPOL3
                ALLOCATE( AEINAM( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AEINAM', PROGNAME )
                CALL STORE_VNAMES( NARVAR3+1, NARPPOL3, ANIPOL, AEINAM )
                INVUNIT = UNITS3D( NARVAR3 + 1 )

            ENDIF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            AGNAME = PROMPTMFILE( 
     &       'Enter logical name for the AREA GRIDDING MATRIX',
     &       FSREAD3, 'AGMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( AGNAME )
            CALL CHKSRCNO( 'area', 'AGMAT', NTHIK3D, NASRC, EFLAG )
            CALL CHECK_GRID_INFO( 'area', 'GMAT' )
            ANGMAT = NCOLS3D

C.............  Open speciation matrix, compare number of sources, store
C               speciation variable descriptions, and store mass or moles.
            IF( SFLAG ) THEN
                ASNAME = PROMPTMFILE( 
     &           'Enter logical name for the AREA SPECIATION MATRIX',
     &           FSREAD3, 'ASMAT', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( ASNAME )
                CALL CHKSRCNO( 'area', 'ASMAT', NROWS3D, NASRC, EFLAG )
                ANSMATV = NVARS3D
                ALLOCATE( ASVDESC( ANSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ASVDESC', PROGNAME )
                CALL STORE_VDESCS( 1, 1, ANSMATV, ASVDESC )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'area' )
                SPCUNIT = UNITS3D( 1 )

            END IF  ! end of speciation open

C.............  Open multiplicative control matrix, compare number of sources, 
C               and store control variable names.
            IF( AUFLAG ) THEN
                MESG = 'Enter logical name for the AREA ' //
     &                 'MULTIPLICATIVE CONTROL MATRIX'
                AUNAME = PROMPTMFILE( MESG, FSREAD3, 'AXMAT', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( AUNAME )
                CALL CHKSRCNO( 'area', 'AXMAT', NROWS3D, NASRC, EFLAG )
                ANUMATV = NVARS3D
                ALLOCATE( AUVNAMS( ANUMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AUVNAMS', PROGNAME )
                CALL STORE_VNAMES( 1, 1, ANUMATV, AUVNAMS )

            END IF  ! end of multiplicative control open

C.............  Open additive control matrix, compare number of sources, 
C               and store control variable names.
            IF( AAFLAG ) THEN
                MESG= 'INTERNAL ERROR: Area additive controls not ' //
     &                'yet implemented in ' // PROGNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF  ! end of additive control open

C.............  Open reactivity control matrix, compare number of sources, and
C               store control variable descriptions, and store mass or moles.
            IF( ARFLAG ) THEN
                ARNAME = PROMPTMFILE( 
     &           'Enter logical name for the AREA REACTIVITY MATRIX',
     &           FSREAD3, 'ARMAT', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( ARNAME )
                CALL CHKSRCNO( 'area', 'ARMAT', NTHIK3D, NASRC, EFLAG )
                ANRMATV = NVARS3D
                ANSREAC = NROWS3D
                ALLOCATE( ARVDESC( ANRMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ARVDESC', PROGNAME )
                CALL STORE_VDESCS( 1, 1, ANRMATV, ARVDESC )

C.................  Retrieve the number of speciation factors 
                ARNMSPC = GETIFDSC( FDESC3D, '/SPECIES VARS/', .TRUE. )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'area' )
                SPCUNIT = UNITS3D( NVARS3D )

C.................  Check the year and projection year of the matrix
                CALL CHECK_INVYEAR( ARNAME, APRJFLAG, FDESC3D )

            END IF  ! end of reactivity control open

        END IF  ! End of section for area sources

C.........  If we have biogenic sources 
        IF( BFLAG ) THEN

            MESG= 'INTERNAL ERROR: Biogenics not yet implemented in '//
     &            PROGNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

c            BTNAME = PROMPTMFILE(
c     &         'Enter logical name for TIME-STEPPED BIOGENIC EMIS file',
c     &          FSREAD3, 'BGTS', 'SMKMERGE' )

c            IF ( .NOT. DESC3( BTNAME ) ) THEN
c                MESG = 'Could not get description of file ' // BTNAME
c                CALL M3EXIT( 'SMKMERGE', 0, 0, MESG, 2 )
c            END IF

c            BMETHEAD( 1 ) = FDESC3D( 2 ) 
c            BMETHEAD( 2 ) = FDESC3D( 3 ) 
c            DO  101  V = 1, NMPOL	! Just set vname3d(*)
c                BSPCFLAG( V ) = 
c     &          ( INDEX1( EMNAM(V), NVARS3D, VNAME3D ) .GT. 0 )
c101         CONTINUE

        END IF  ! End of section for biogenic sources

C.........  If we have mobile sources 
        IF( MFLAG ) THEN

            MESG= 'INTERNAL ERROR: Mobile not yet implemented in '//
     &            PROGNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF  ! End of section for mobile sources

C.........  If we have point sources 
        IF( PFLAG ) THEN

C.............  Prompt for inventory files
            PENAME = PROMPTMFILE( 
     &       'Enter logical name for the I/O API POINT INVENTORY file',
     &       FSREAD3, 'PNTS', PROGNAME )

            PSDEV = PROMPTFFILE( 
     &       'Enter logical name for the ASCII POINT INVENTORY file',
     &       .TRUE., .TRUE., 'PSRC', PROGNAME )

C.............  Get number of sources
            CALL RETRIEVE_IOAPI_HEADER( PENAME )
            NPSRC = NROWS3D

C.............  Determine the year and projection status of the inventory
            CALL CHECK_INVYEAR( PENAME, PPRJFLAG, FDESC3D )

C.............  For temporal inputs, prompt for hourly file
            IF( TFLAG ) THEN

                MESG = 'Enter logical name for the POINT HOURLY ' //
     &                 'EMISSIONS file'
                PTNAME = PROMPTMFILE( MESG, FSREAD3, 'PTMP', PROGNAME )

C.................  Set parameters and pollutants from hourly file
                CALL RETRIEVE_IOAPI_HEADER( PTNAME )
                CALL CHKSRCNO( 'point', 'PTMP', NROWS3D, NPSRC, EFLAG )
                CALL UPDATE_TIME_INFO( 'PTMP' )

                PNIPOL = NVARS3D
                ALLOCATE( PEINAM( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PEINAM', PROGNAME )
                CALL STORE_VNAMES( 1, 1, PNIPOL, PEINAM )
                INVUNIT = UNITS3D( 1 )

C.................  Determine the year and projection status of the hourly 
                CALL CHECK_INVYEAR( PTNAME, PPRJFLAG, FDESC3D )

C.............  Otherwise, just set parameters and pollutants from inven file
            ELSE
                PTNAME = PENAME
                PNIPOL = ( NVARS3D - NPTVAR3 ) / NPTPPOL3
                ALLOCATE( PEINAM( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PEINAM', PROGNAME )
                CALL STORE_VNAMES( NPTVAR3+1, NPTPPOL3, PNIPOL, PEINAM )
                INVUNIT = UNITS3D( NPTVAR3 + 1 )

            ENDIF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            PGNAME = PROMPTMFILE( 
     &       'Enter logical name for the POINT GRIDDING MATRIX',
     &       FSREAD3, 'PGMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( PGNAME )
            CALL CHKSRCNO( 'point', 'PGMAT', NCOLS3D, NPSRC, EFLAG )
            CALL CHECK_GRID_INFO( 'point', 'GMAT' )

C.............  Open speciation matrix, compare number of sources, store
C               speciation variable names, and store mass or moles.
            IF( SFLAG ) THEN
                PSNAME = PROMPTMFILE( 
     &           'Enter logical name for the POINT SPECIATION MATRIX',
     &           FSREAD3, 'PSMAT', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( PSNAME )
                CALL CHKSRCNO( 'point','PSMAT',NROWS3D,NPSRC,EFLAG )
                PNSMATV = NVARS3D
                ALLOCATE( PSVDESC( PNSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PSVDESC', PROGNAME )
                CALL STORE_VDESCS( 1, 1, PNSMATV, PSVDESC )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'point' )
                SPCUNIT = UNITS3D( 1 )

            END IF  ! end of speciation open

C.............  Open multiplicative control matrix, compare number of sources, 
C               and store control variable names.
            IF( PUFLAG ) THEN
                MESG = 'Enter logical name for the POINT ' //
     &                 'MULTIPLICATIVE CONTROL MATRIX'
                PUNAME = PROMPTMFILE( MESG, FSREAD3, 'PXMAT', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( PUNAME )
                CALL CHKSRCNO( 'point', 'PXMAT', NROWS3D, NPSRC, EFLAG )
                PNUMATV = NVARS3D
                ALLOCATE( PUVNAMS( PNUMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PUVNAMS', PROGNAME )
                CALL STORE_VNAMES( 1, 1, PNUMATV, PUVNAMS )

            END IF  ! end of multiplicative control open

C.............  Open additive control matrix, compare number of sources, 
C               and store control variable names.
            IF( PAFLAG ) THEN
                MESG= 'INTERNAL ERROR: Point additive controls not ' //
     &                'yet implemented in ' // PROGNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF  ! end of additive control open

C.............  Open reactivity control matrix, compare number of sources, and
C               store control variable descriptions, and store mass or moles.
            IF( PRFLAG ) THEN
                PRNAME = PROMPTMFILE( 
     &           'Enter logical name for the POINT REACTIVITY MATRIX',
     &           FSREAD3, 'PRMAT', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( PRNAME )
                CALL CHKSRCNO( 'point', 'PRMAT', NTHIK3D, NPSRC, EFLAG )
                PNRMATV = NVARS3D
                PNSREAC = NROWS3D
                ALLOCATE( PRVDESC( PNRMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PRVDESC', PROGNAME )
                CALL STORE_VDESCS( 1, 1, PNRMATV, PRVDESC )

C.................  Retrieve the number of speciation factors 
                PRNMSPC = GETIFDSC( FDESC3D, '/SPECIES VARS/', .TRUE. )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'point' )
                SPCUNIT = UNITS3D( NVARS3D )

C.................  Check the year and projection year of the matrix
                CALL CHECK_INVYEAR( PRNAME, PPRJFLAG, FDESC3D )

            END IF  ! end of reactivity control open

C.............  Open layer fractions file, compare number of sources, check 
C               met information, and store the vertical coordinates info
            IF( LFLAG ) THEN
                MESG= 'Enter logical name for the POINT LAYER ' //
     &                'FRACTIONS MATRIX'
                PLNAME = PROMPTMFILE( MESG, FSREAD3, 'PLAY', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( PLNAME )
                CALL CHKSRCNO( 'point', 'PLAY', NROWS3D, NPSRC, EFLAG )
                CALL UPDATE_TIME_INFO( 'PLAY' )

                IF( LMETCHK ) CALL CHECK_MET_INFO( 'point' ) 

                EMLAYS = NLAYS3D
                VGTYP  = VGTYP3D
                VGTOP  = VGTOP3D

C.................  Deal with vertical coordinate info, but be adaptive
C                   to the potential for 0-based or 1-based VGLVS3D
                ALLOCATE( VGLVS( 0:EMLAYS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'VGLVS', PROGNAME )
              
                J = LBOUND( VGLVS3D,1 )
                DO V = 0, EMLAYS
                    VGLVS( V ) = VGLVS3D( J )
                    J = J + 1
                END DO

            END IF  ! End of layer fractions open

        END IF      ! End of section for point sources

C.........  Get master pollutants list so we will be able to output in the
C           proper order in case different source categories have different
C           pollutants.

        PDEV = PROMPTFFILE( 
     &         'Enter logical name for POLLUTANT CODES & NAMES file',
     &         .TRUE., .TRUE., 'SIPOLS', PROGNAME )

C.........  If there were any errors inputing files or while comparing
C           with one another, then abort
        IF( EFLAG ) THEN
           MESG = 'Problems opening input files. See ERROR(S) above.'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  If we are using temporalized emissions, then update date/time and
C           duration using environment variable settings, then prompt.
        IF( TFLAG ) THEN

             CALL GETM3EPI( TZONE, SDATE, STIME, NSTEPS )

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

C----------------------------------------------------------------------

C.............  If time information has already been initialized...
            IF( IFLAG ) THEN
                ISECS = SECSDIFF( SDATE, STIME, SDATE3D, STIME3D )

                IF( ISECS .GT. 0 ) THEN  ! SDATE3D/SDATE3D are later
                    SDATE = SDATE3D
                    STIME = STIME3D
                END IF

                ED = SDATE3D
                ET = STIME3D
                CALL NEXTIME( ED, ET, MXREC3D * TSTEP3D )
        
                ISECS = SECSDIFF( EDATE, ETIME, ED, ET )

                IF( ISECS .LT. 0 ) THEN  ! ED/ET are earlier
                    EDATE = ED
                    ETIME = ET
                END IF

                NSTEPS = SECSDIFF( SDATE, STIME, EDATE, ETIME ) / 3600

                IF( TFLAG .AND. NSTEPS .LE. 0 ) THEN
                    MESG = 'Because of file ' // FILNAM // 
     &                     ', dates and times do not overlap at all!'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C.............  If time information needs to be initialized...
            ELSE
                SDATE  = SDATE3D
                STIME  = STIME3D
                NSTEPS = MXREC3D

                EDATE  = SDATE
                ETIME  = STIME
                CALL NEXTIME( EDATE, ETIME, NSTEPS * TSTEP3D )

                IFLAG = .TRUE.

            END IF

C.............  Make sure that time step is one hour
            IF( TSTEP3D .NE. 10000 ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: Time step is not one hour in ' // 
     &                 FILNAM // 'file!'
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

                ELSE
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
            LOGICAL           STRICT  ! flag for strict checks or not
            CHARACTER*20      BUFFER  ! program name buffer
            INTEGER , SAVE :: FLEN    ! name length of savnam
            CHARACTER(LEN=IOVLEN3), SAVE :: SAVNAM  ! name of file used to init

C----------------------------------------------------------------------

            STRICT = .TRUE.

C.............  First determine whether to abort when projected year does not
C               match.  This is used for reactivity matrices, which will
C               always have a projection year, even if the inventory isn't
C               projected.
            IF( .NOT. PRJFLAG ) THEN
                BUFFER = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
                IF( BUFFER .EQ. 'OPENRMAT' ) STRICT = .FALSE.
            END IF

C.............  If time information has already been initialized...
            IF( YFLAG ) THEN

                YY = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )
                IF( YY .LE. 0 ) THEN

                    YY = GETIFDSC( IODESC, '/BASE YEAR/', .TRUE. ) 
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
                
                BYEAR = GETIFDSC( IODESC, '/BASE YEAR/', .TRUE. ) 
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
C.............  This subprogram updates the grid information and compares to 
C               the existing information, if it has been previously set.
            SUBROUTINE CHECK_GRID_INFO( CATDESC, FTYPE )

C.............  Subprogram arguments
            CHARACTER(*) CATDESC  ! category descriptions
            CHARACTER(*) FTYPE    ! file type (GMAT or GRID) to get grid info

C.............  Local variables
            INTEGER       L       ! length of file description
            INTEGER       NC      ! tmp number of columns
            INTEGER       NR      ! tmp number of rwos
            CHARACTER*20  FILDESC ! description of input file

C----------------------------------------------------------------------

C.............  Set tmp rows, columns, and total cells depending on file type
            IF( FTYPE .EQ. 'GMAT' ) THEN
                NC = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
                NR = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
                FILDESC = 'gridding matrix'

            ELSEIF( FTYPE .EQ. 'GRID' ) THEN
                NC = NCOLS3D
                NR = NROWS3D
                FILDESC = 'gridded file'

            ELSE
                MESG = 'INTERNAL ERROR: File type "' // FTYPE // 
     &                 '" not known in call to CHECK_GRID_INFO!'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF

            L = LEN_TRIM( FILDESC )

C.............  If grid information has already been initialized, then compare
C               existing to this file.
            IF( GFLAG ) THEN

                IF ( NCOLS .NE. NC      .OR.
     &               NROWS .NE. NR      .OR.
     &               GDTYP .NE. GDTYP3D .OR.
     &               FLTERR( XCELL, SNGL( XCELL3D ) ) .OR.
     &               FLTERR( YCELL, SNGL( YCELL3D ) ) .OR.
     &               FLTERR( XORIG, SNGL( XORIG3D ) ) .OR.
     &               FLTERR( YORIG, SNGL( YORIG3D ) ) .OR.
     &               FLTERR( XCENT, SNGL( XCENT3D ) ) .OR.
     &               FLTERR( YCENT, SNGL( YCENT3D ) ) .OR.
     &               FLTERR( P_ALP, SNGL( P_ALP3D ) ) .OR.
     &               FLTERR( P_BET, SNGL( P_BET3D ) ) .OR.
     &               FLTERR( P_GAM, SNGL( P_GAM3D ) )      ) THEN

                    EFLAG = .TRUE.
                    MESG = 'Grid parameters in ' // CATDESC // ' ' //
     &                     FILDESC( 1:L ) // ' are not consistent ' //
     &                     'with initialized values.'
                    CALL M3MSG2( MESG )

                END IF

C.............  Initialize grid information
            ELSE

                GFLAG = .TRUE.
                GRDNM = GDNAM3D
                GDTYP = GDTYP3D
                P_ALP = SNGL( P_ALP3D )
                P_BET = SNGL( P_BET3D )
                P_GAM = SNGL( P_GAM3D )
                XCENT = SNGL( XCENT3D )
                YCENT = SNGL( YCENT3D )
                XORIG = SNGL( XORIG3D )
                YORIG = SNGL( YORIG3D )
                XCELL = SNGL( XCELL3D )
                YCELL = SNGL( YCELL3D )
                NCOLS = NC
                NROWS = NR
                NGRID = NCOLS * NROWS

                MESG = 'NOTE: Grid settings initialized using ' // 
     &                 CATDESC // ' ' // FILDESC( 1:L ) // '.'

                CALL M3MSG2( MESG )

            ENDIF

            END SUBROUTINE CHECK_GRID_INFO

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram sets the speciation type and compares to 
C               the existing information, if it has been previously set.
            SUBROUTINE CHECK_SPEC_TYPE( CATDESC )

C.............  Subprogram arguments
            CHARACTER(*) CATDESC  ! category descriptions

C.............  Local variables
            CHARACTER*4  LOCTYPE  ! tmp speciation type 

C----------------------------------------------------------------------

            LOCTYPE = GETCFDSC( FDESC3D, '/SMATTYPE/', .TRUE. )

C.............  If speciation type has already been initialized, then compare
C               existing to this file.
            IF( CFLAG ) THEN

                IF ( LOCTYPE .NE. SPCTYPE ) THEN

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Speciation type "' // LOCTYPE // 
     &                     '" in ' // CATDESC // ' speciation matrix '//
     &                     'is inconsistent with initialized type "' //
     &                     SPCTYPE // '"'
                    CALL M3MSG2( MESG )

                END IF

C.............  Initialize speciation type information
            ELSE

                CFLAG   = .TRUE.
                SPCTYPE = LOCTYPE

                MESG = 'NOTE: Speciation type initialized '//
     &                 'using '// CATDESC // ' speciation matrix.'
                CALL M3MSG2( MESG )

            ENDIF

            END SUBROUTINE CHECK_SPEC_TYPE

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram updates the met information and compares to 
C               the existing information, if it has been previously set.
            SUBROUTINE CHECK_MET_INFO( CATDESC )

C.............  Subprogram arguments
            CHARACTER(*) CATDESC  ! category descriptions

C.............  Local variables
            INTEGER       L, L1, L2  ! length of strings
            CHARACTER*20  FILDESC    ! description of input file

C----------------------------------------------------------------------

C.............  Set tmp rows, columns, and total cells depending on file type
            IF( CATDESC .EQ. 'biogenic' ) THEN
                FILDESC = 'gridded emissions file'

            ELSEIF( CATDESC .EQ. 'mobile' ) THEN
                FILDESC = 'hourly emissions file'

            ELSEIF( CATDESC .EQ. 'point' ) THEN
                FILDESC = 'layer fractions file'

            ELSE
                MESG= 'INTERNAL ERROR: Category description "' // 
     &                CATDESC// '" not known in call to CHECK_MET_INFO!'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ENDIF

            L = LEN_TRIM( FILDESC )

C.............  If met information has already been initialized, then compare
C               existing to this file.
            IF( OFLAG ) THEN

                METTMP = GETCFDSC( FDESC3D, '/MET SCENARIO/', .TRUE. )
                IF ( METTMP .NE. METSCENR ) THEN

                    L1 = LEN_TRIM( METTMP )
                    L2 = LEN_TRIM( METSCENR )

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Meteorology scenario name "' // 
     &                     METTMP( 1:L1 ) // '" in ' // CATDESC //
     &                     FILDESC( 1:L ) // ' is inconsistent with '//
     &                     'initialized value "'// METSCENR(1:L2)// '"'
                    CALL M3MSG2( MESG )

                END IF

                METTMP = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .TRUE. )
                IF ( METTMP .NE. METCLOUD ) THEN

                    L1 = LEN_TRIM( METTMP )
                    L2 = LEN_TRIM( METCLOUD )

                    EFLAG = .TRUE.
                    MESG = 'ERROR: Meteorology cloud scheme "' // 
     &                     METTMP( 1:L1 ) // '" in ' // CATDESC //
     &                     FILDESC( 1:L ) // ' is inconsistent with '//
     &                     'initialized value "'// METCLOUD(1:L2)// '"'
                    CALL M3MSG2( MESG )

                END IF

C.............  Initialize meteorology information
            ELSE

                OFLAG    = .TRUE.
                METSCENR = GETCFDSC( FDESC3D, '/MET SCENARIO/', .TRUE. )
                METCLOUD = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .TRUE. )

                MESG = 'NOTE: Meteorology description initialized '//
     &                 'using '// CATDESC// ' '// FILDESC( 1:L )// '.'
                CALL M3MSG2( MESG )

            ENDIF

            END SUBROUTINE CHECK_MET_INFO

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

        END SUBROUTINE OPENMRGIN
