
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
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2     CRLF
        LOGICAL         DSCM3GRD
        LOGICAL         ENVYN
        CHARACTER*50    GETCFDSC  
        INTEGER         GETIFDSC  
        INTEGER         PROMPTFFILE  
        CHARACTER*16    PROMPTMFILE  
        INTEGER         SECSDIFF  

        EXTERNAL  CRLF, ENVYN, GETCFDSC, GETIFDSC, PROMPTFFILE, 
     &            PROMPTMFILE, SECSDIFF

C.........  Other local variables

        INTEGER         I, J, N, V       ! counters and indices

        INTEGER         IOS           ! tmp I/O status
        INTEGER         ISECS         ! tmp duration in seconds
        INTEGER         NPACT         ! no. variables per activity
        INTEGER         NPPOL         ! no. variables per pollutant
        INTEGER         NDIM          ! tmp dimensioning variable 
        INTEGER         NVAR          ! tmp no. variables 

        LOGICAL      :: CFLAG = .FALSE.  ! true: speciation type has been init
        LOGICAL      :: DFLAG = .FALSE.  ! true: use pollutants list
        LOGICAL      :: EFLAG = .FALSE.  ! true: error in routine
        LOGICAL      :: IFLAG = .FALSE.  ! true: episode settings have been init
        LOGICAL      :: KFLAG = .FALSE.  ! true: use activies list
        LOGICAL      :: OFLAG = .FALSE.  ! true: met info has been init
        LOGICAL      :: YFLAG = .FALSE.  ! true: year/projection info been init
        LOGICAL      :: ZFLAG = .FALSE.  ! true: time zone has been init

        CHARACTER*4     SPCTYPE      ! type of speciation matrix (mass|mole)
        CHARACTER*16    DUMNAME      ! tmp file name
        CHARACTER*50    METSCENR     ! met scenario name
        CHARACTER*50    METCLOUD     ! met cloud scheme name
        CHARACTER*50    METTMP       ! temporary buffer for met info
        CHARACTER*80    GDESC        ! grid description
        CHARACTER*300   MESG         ! message buffer
        CHARACTER(LEN=IOVLEN3) COORD3D    ! coordinate system name 
        CHARACTER(LEN=IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(LEN=IOVLEN3) PROJTYPE   ! projection type

        CHARACTER*16 :: PROGNAME = 'OPENMRGIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENMRGIN

C.........  Set controls for reading the pollutants and activities files
C.........  Default is for mobile to read in activities and not pollutants
C           and for other source categories to read in pollutants and not
C           activities
        IF( MFLAG ) THEN
            KFLAG = .TRUE.
        ELSE
            KFLAG = .FALSE.
        END IF

        IF( AFLAG .OR. PFLAG .OR. BFLAG .OR. ( MFLAG .AND. TFLAG )) THEN
            DFLAG = .TRUE.
        ELSE
            DFLAG = .FALSE.
        END IF

C.........  Get value of these controls from the environment
        MESG = 'Indicator for using pollutants list'
        DFLAG = ENVYN( 'SMK_USE_SIPOLS', MESG, DFLAG, IOS )

        MESG = 'Indicator for using activities list'
        KFLAG = ENVYN( 'SMK_USE_ACTVNAMS', MESG, KFLAG, IOS )

C.........  Initialize gridded information with grid description file
        IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL CHKGRID( 'general', 'GRIDDESC', EFLAG )   ! Initialization call
        
C.........  For area sources... 
        IF( AFLAG ) THEN

C.............  Get inventory file names given source category
            CALL GETINAME( 'AREA', AENAME, DUMNAME )

C.............  Prompt for inventory files
            AENAME = PROMPTMFILE( 
     &       'Enter logical name for the I/O API AREA INVENTORY file',
     &       FSREAD3, AENAME, PROGNAME )

            ASDEV = PROMPTFFILE( 
     &       'Enter logical name for the ASCII AREA INVENTORY file',
     &       .TRUE., .TRUE., DUMNAME, PROGNAME )

C.............  Get number of sources
            CALL RETRIEVE_IOAPI_HEADER( AENAME )
            NASRC = NROWS3D

C.............  Determine the year and projection status of the inventory
            CALL CHECK_INVYEAR( AENAME, APRJFLAG, FDESC3D )

C.............  For temporal inputs, prompt for hourly file
            IF( TFLAG ) THEN

                MESG = 'Enter logical name for the AREA HOURLY ' //
     &                 'EMISSIONS file'
                ATNAME = PROMPTMFILE( MESG, FSREAD3, 'ATMP', PROGNAME )

C.................  Set parameters and pollutants from hourly file
                CALL RETRIEVE_IOAPI_HEADER( ATNAME )
                CALL CHKSRCNO( 'area', ATNAME, NROWS3D, NASRC, EFLAG )
                CALL UPDATE_TIME_INFO( ATNAME )

                ANIPOL = NVARS3D
                ALLOCATE( AEINAM( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AEINAM', PROGNAME )
                ALLOCATE( AONAMES( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AONAMES', PROGNAME )
                ALLOCATE( AOUNITS( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AOUNITS', PROGNAME )

                CALL STORE_VNAMES( 1, 1, ANIPOL, AEINAM )
                CALL STORE_INVINFO( 1, 1, ANIPOL, 1, INVPIDX,
     &                              AONAMES, AOUNITS )

C.................  Determine the year and projection status of the hourly
                CALL CHECK_INVYEAR( ATNAME, APRJFLAG, FDESC3D )

C.............  Otherwise, just set parameters and pollutants from inven file
            ELSE
                ATNAME = AENAME
        	NVAR   = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
        	ANIPOL = GETIFDSC( FDESC3D, '/POLLUTANTS/', .FALSE. )
        	NPPOL  = GETIFDSC( FDESC3D, '/PER POLLUTANT/', .FALSE. )

                ALLOCATE( AEINAM ( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AEINAM', PROGNAME )
                ALLOCATE( AONAMES( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AONAMES', PROGNAME )
                ALLOCATE( AOUNITS( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AOUNITS', PROGNAME )

                CALL STORE_VNAMES( NVAR+1, NPPOL, ANIPOL, AEINAM )
                CALL STORE_INVINFO( NVAR+1, NPPOL, ANIPOL, 1, INVPIDX,
     &                              AONAMES, AOUNITS )

            END IF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            AGNAME = PROMPTMFILE( 
     &         'Enter logical name for the AREA GRIDDING MATRIX',
     &         FSREAD3, 'AGMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( AGNAME )
            CALL CHKSRCNO( 'area', 'AGMAT', NTHIK3D, NASRC, EFLAG )
            CALL CHKGRID( 'area', 'GMAT', EFLAG )
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
                ALLOCATE( ASVUNIT( ANSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ASVUNIT', PROGNAME )
                CALL STORE_VDESCS( 1, 1, ANSMATV, ASVDESC )
                CALL STORE_VUNITS( 1, 1, ANSMATV, ASVUNIT )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'area' )

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

C.................  Check the year and projection year of the matrix
                CALL CHECK_INVYEAR( ARNAME, APRJFLAG, FDESC3D )

            END IF  ! end of reactivity control open

        END IF  ! End of section for area sources

C.........  If we have biogenic sources 
        IF( BFLAG ) THEN

            MESG = 'Enter logical name for TIME-STEPPED BIOGENIC ' //
     &             'EMISSIONS file'
            BTNAME = PROMPTMFILE( MESG, FSREAD3, 'BGTS', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( BTNAME )
            CALL UPDATE_TIME_INFO( BTNAME )
            CALL CHKGRID( 'biogenics', 'GRID', EFLAG )

            IF( LMETCHK ) CALL CHECK_MET_INFO( 'biogenics' ) 

C.............  Store biogenic species names as if they are stored as
C               speciation matrix variable descriptions
            BIOUNIT = UNITS3D( 1 )
            BNSMATV = NVARS3D
            ALLOCATE( BSVDESC( BNSMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BSVDESC', PROGNAME )

            DO I = 1, BNSMATV
                IF( VNAME3D( I ) .NE. 'NO' ) THEN
                    DUMNAME = 'VOC'
                    BSVDESC( I ) = DUMNAME // SPJOIN // VNAME3D( I )
                ELSE
                    DUMNAME = 'NOX'
                    BSVDESC( I ) = DUMNAME // SPJOIN // VNAME3D( I )
                END IF

            END DO

C.............  Store biogenic pollutant names and units
            BNIPOL  = 2
            BEINAM ( 1 ) = 'NOX'
            BEINAM ( 2 ) = 'VOC'
            BONAMES      = BEINAM
            BOUNITS      = UNITS3D( 1 )

        END IF  ! End of section for biogenic sources

C.........  If we have mobile sources 
        IF( MFLAG ) THEN

C.............  Get inventory file names given source category
            CALL GETINAME( 'MOBILE', MENAME, DUMNAME )

C.............  Prompt for inventory files
            MENAME = PROMPTMFILE( 
     &       'Enter logical name for the I/O API MOBILE INVENTORY file',
     &       FSREAD3, MENAME, PROGNAME )

            MSDEV = PROMPTFFILE( 
     &       'Enter logical name for the ASCII MOBILE INVENTORY file',
     &       .TRUE., .TRUE., DUMNAME, PROGNAME )

C.............  Get number of sources
            CALL RETRIEVE_IOAPI_HEADER( MENAME )
            NMSRC = NROWS3D

C.............  Determine the year and projection status of the inventory
            CALL CHECK_INVYEAR( MENAME, MPRJFLAG, FDESC3D )

C.............  For temporal inputs, prompt for hourly file
            IF( TFLAG ) THEN

                MESG = 'Enter logical name for the MOBILE HOURLY ' //
     &                 'EMISSIONS file'
                MTNAME = PROMPTMFILE( MESG, FSREAD3, 'MTMP', PROGNAME )

C.................  Set parameters and pollutants from hourly file
                CALL RETRIEVE_IOAPI_HEADER( MTNAME )
                CALL CHKSRCNO( 'mobile', MTNAME, NROWS3D, NMSRC, EFLAG )
                CALL UPDATE_TIME_INFO( MTNAME )

                MNIPPA = NVARS3D
                ALLOCATE( MEANAM( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MEANAM', PROGNAME )
                ALLOCATE( MONAMES( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MONAMES', PROGNAME )
                ALLOCATE( MOUNITS( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MOUNITS', PROGNAME )

                CALL STORE_VNAMES( 1, 1, MNIPPA, MEANAM )
                CALL STORE_INVINFO( 1, 1, MNIPPA, 1, INVPIDX,
     &                              MONAMES, MOUNITS )

C.................  Determine the year and projection status of the hourly
                CALL CHECK_INVYEAR( MTNAME, MPRJFLAG, FDESC3D )

C.............  Otherwise, just set parameters and pollutants from inven file
            ELSE
                MTNAME = MENAME

        	NVAR   = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
        	MNIPOL = GETIFDSC( FDESC3D, '/POLLUTANTS/', .FALSE. )
        	NPPOL  = GETIFDSC( FDESC3D, '/PER POLLUTANT/', .FALSE. )
        	MNIACT = GETIFDSC( FDESC3D, '/ACTIVITIES/', .FALSE. )
        	NPACT  = GETIFDSC( FDESC3D, '/PER ACTIVITY/', .FALSE. )

        	MNIPOL = MAX( 0, MNIPOL )
        	NPPOL  = MAX( 0, NPPOL )
        	MNIACT = MAX( 0, MNIACT )
        	NPACT  = MAX( 0, NPACT )
        	MNIPPA = MNIPOL + MNIACT

                ALLOCATE( MEANAM( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MEANAM', PROGNAME )
                ALLOCATE( MONAMES( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MONAMES', PROGNAME )
                ALLOCATE( MOUNITS( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MOUNITS', PROGNAME )

C.................  Store pollutant names and other information
                CALL STORE_VNAMES( NVAR+1, NPPOL, MNIPOL, MEANAM )
                CALL STORE_INVINFO( NVAR+1, NPPOL, MNIPOL, 1, INVPIDX,
     &                              MONAMES, MOUNITS )

C.................  Store activity names and other information
                I    = MNIPOL * NPPOL
                NVAR = NVAR + I
                CALL STORE_VNAMES( NVAR+1, NPACT, MNIACT, 
     &                             MEANAM( MNIPOL+1 )     )
                CALL STORE_INVINFO( NVAR+1, NPACT, MNIACT, I + 1, 1,
     &                              MONAMES, MOUNITS )

            END IF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            MGNAME = PROMPTMFILE( 
     &       'Enter logical name for the MOBILE GRIDDING MATRIX',
     &       FSREAD3, 'MGMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( MGNAME )
            CALL CHKSRCNO( 'mobile', 'MGMAT', NTHIK3D, NMSRC, EFLAG )
            CALL CHKGRID( 'mobile', 'GMAT', EFLAG )
            MNGMAT = NCOLS3D

C.............  Open speciation matrix, compare number of sources, store
C               speciation variable descriptions, and store mass or moles.
            IF( SFLAG ) THEN
                MSNAME = PROMPTMFILE( 
     &           'Enter logical name for the MOBILE SPECIATION MATRIX',
     &           FSREAD3, 'MSMAT', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( MSNAME )
                CALL CHKSRCNO( 'mobile', 'MSMAT', NROWS3D, NMSRC, EFLAG)
                MNSMATV = NVARS3D
                ALLOCATE( MSVDESC( MNSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MSVDESC', PROGNAME )
                ALLOCATE( MSVUNIT( MNSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MSVUNIT', PROGNAME )
                CALL STORE_VDESCS( 1, 1, MNSMATV, MSVDESC )
                CALL STORE_VUNITS( 1, 1, MNSMATV, MSVUNIT )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'mobile' )

            END IF  ! end of speciation open

C.............  Open reactivity control matrix, compare number of sources, and
C               store control variable descriptions, and store mass or moles.
            IF( MRFLAG ) THEN
                MRNAME = PROMPTMFILE( 
     &           'Enter logical name for the MOBILE REACTIVITY MATRIX',
     &           FSREAD3, 'MRMAT', PROGNAME )

                CALL RETRIEVE_IOAPI_HEADER( ARNAME )
                CALL CHKSRCNO( 'mobile', 'MRMAT', NTHIK3D, NMSRC, EFLAG)
                MNRMATV = NVARS3D
                MNSREAC = NROWS3D
                ALLOCATE( MRVDESC( MNRMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MRVDESC', PROGNAME )
                CALL STORE_VDESCS( 1, 1, MNRMATV, MRVDESC )

C.................  Retrieve the number of speciation factors 
                MRNMSPC = GETIFDSC( FDESC3D, '/SPECIES VARS/', .TRUE. )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'mobile' )

C.................  Check the year and projection year of the matrix
                CALL CHECK_INVYEAR( MRNAME, MPRJFLAG, FDESC3D )

            END IF  ! end of reactivity control open

        ENDIF  ! End of section for mobile sources

C.........  If we have point sources 
        IF( PFLAG ) THEN

C.............  Get inventory file names given source category
            CALL GETINAME( 'POINT', PENAME, DUMNAME )

C.............  Prompt for inventory files
            PENAME = PROMPTMFILE( 
     &       'Enter logical name for the I/O API POINT INVENTORY file',
     &       FSREAD3, PENAME, PROGNAME )

            PSDEV = PROMPTFFILE( 
     &       'Enter logical name for the ASCII POINT INVENTORY file',
     &       .TRUE., .TRUE., DUMNAME, PROGNAME )

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
                CALL CHKSRCNO( 'point', PTNAME, NROWS3D, NPSRC, EFLAG )
                CALL UPDATE_TIME_INFO( PTNAME )

                PNIPOL = NVARS3D
                ALLOCATE( PEINAM( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PEINAM', PROGNAME )
                ALLOCATE( PONAMES( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PONAMES', PROGNAME )
                ALLOCATE( POUNITS( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'POUNITS', PROGNAME )

                CALL STORE_VNAMES( 1, 1, PNIPOL, PEINAM )
                CALL STORE_INVINFO( 1, 1, PNIPOL, 1, INVPIDX,
     &                              PONAMES, POUNITS )

C.................  Determine the year and projection status of the hourly 
                CALL CHECK_INVYEAR( PTNAME, PPRJFLAG, FDESC3D )

C.............  Otherwise, just set parameters and pollutants from inven file
            ELSE
                PTNAME = PENAME
        	NVAR   = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
        	PNIPOL = GETIFDSC( FDESC3D, '/POLLUTANTS/', .FALSE. )
        	NPPOL  = GETIFDSC( FDESC3D, '/PER POLLUTANT/', .FALSE. )

                ALLOCATE( PEINAM( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PEINAM', PROGNAME )
                ALLOCATE( PONAMES( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PONAMES', PROGNAME )
                ALLOCATE( POUNITS( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'POUNITS', PROGNAME )

                CALL STORE_VNAMES( NPTVAR3+1, NPTPPOL3, PNIPOL, PEINAM )
                CALL STORE_INVINFO( NVAR+1, NPPOL, PNIPOL, 1, 1,
     &                              PONAMES, POUNITS )

            END IF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            PGNAME = PROMPTMFILE( 
     &       'Enter logical name for the POINT GRIDDING MATRIX',
     &       FSREAD3, 'PGMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( PGNAME )
            CALL CHKSRCNO( 'point', 'PGMAT', NCOLS3D, NPSRC, EFLAG )
            CALL CHKGRID( 'point', 'GMAT', EFLAG )

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
                ALLOCATE( PSVUNIT( PNSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PSVUNIT', PROGNAME )
                CALL STORE_VDESCS( 1, 1, PNSMATV, PSVDESC )
                CALL STORE_VUNITS( 1, 1, PNSMATV, PSVUNIT )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'point' )

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
                CALL CHKSRCNO( 'point', PLNAME, NROWS3D, NPSRC, EFLAG )
                CALL UPDATE_TIME_INFO( PLNAME )

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

C.............  Open elevated/plume-in-grid file - for plume-in-grid outputs
            IF( PINGFLAG ) THEN
                MESG = 'Enter logical name for the ELEVATED/PING ' //
     &                 'file'
                EDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                              'PELV', PROGNAME      )
            END IF

        END IF      ! End of section for point sources

C.........  Get master pollutants list so we will be able to output in the
C           proper order in case different source categories have different
C           pollutants.
        IF( DFLAG ) THEN

            PDEV = PROMPTFFILE( 
     &         'Enter logical name for POLLUTANT CODES & NAMES file',
     &         .TRUE., .TRUE., 'SIPOLS', PROGNAME )
        END IF

C.........  Get master activities list 
        IF( KFLAG ) THEN

            VDEV = PROMPTFFILE( 
     &         'Enter logical name for ACTIVITY CODES & NAMES file',
     &         .TRUE., .TRUE., 'ACTVNAMS', PROGNAME )

        END IF

C.........  Get country, state, and county names no matter what, because it is
C           needed to allocate memory for the state and county totals, even
C           when they aren't going to be output
        CDEV = PROMPTFFILE( 
     &             'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )

C.........  If reporting state and/or county emissions, and processing for
C           biogenic sources, get gridding surrogates
        IF( LREPANY ) THEN

            IF( BFLAG ) THEN
                GDEV = PROMPTFFILE( 
     &             'Enter logical name for SURROGATE COEFFICIENTS file',
     &             .TRUE., .TRUE., 'AGPRO', PROGNAME )
            END IF

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

C.............  If time information needs to be initialized...
            ELSE
                SDATE  = SDATE3D
                STIME  = STIME3D
                NSTEPS = MXREC3D

                EDATE  = SDATE
                ETIME  = STIME
                CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP3D )

                IFLAG = .TRUE.

            END IF

C.............  Make sure that time step is one hour
            IF( TSTEP3D .NE. 10000 ) THEN

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
            CHARACTER*30  FILDESC    ! description of input file

C----------------------------------------------------------------------

C.............  Set tmp rows, columns, and total cells depending on file type
            IF( CATDESC .EQ. 'biogenics' ) THEN
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
C.............  This subprogram stores I/O API NetCDF variable names and
C               units from the inventory file into a local array based 
C               on indices in subprogram call.
            SUBROUTINE STORE_INVINFO( ISTART, NPER, NPOA, IDX1, IDX2,
     &                                NAMES, UNITS )

C.............  Subprogram arguments
            INTEGER      ISTART        ! starting position in VNAMES of names
            INTEGER      NPER          ! no. variables per pollutant
            INTEGER      NPOA          ! number of pollutants or activities
            INTEGER      IDX1          ! start index for output variables
            INTEGER      IDX2          ! index to which pol-assoc variable
            CHARACTER(*) NAMES( NPOA ) ! stored variable names
            CHARACTER(*) UNITS( NPOA ) ! stored variable units

C.............  Local variables
            INTEGER  I, J, K

C----------------------------------------------------------------------

            J = ISTART + IDX2 - NPER - 1
            K = IDX1 - 1
            DO I = 1, NPOA

                J = J + NPER
                K = K + 1
                NAMES( K ) = VNAME3D( J )
                UNITS( K ) = UNITS3D( J )

            END DO
 
            END SUBROUTINE STORE_INVINFO

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

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores I/O API NetCDF variable units into
C               a local array based on indices in subprogram call.
            SUBROUTINE STORE_VUNITS( ISTART, INCRMT, NUNIT, UNITS )

C.............  Subprogram arguments
            INTEGER      ISTART        ! starting position in VDESCS of names
            INTEGER      INCRMT        ! increment of VDESCS for names
            INTEGER      NUNIT         ! number of units
            CHARACTER(*) UNITS( NUNIT )! stored variable units

C.............  Local variables
            INTEGER  I, J, L

C----------------------------------------------------------------------

            UNITS = ' '

            J = ISTART
            DO I = 1, NUNIT

                L = LEN_TRIM( UNITS3D( J ) )
                UNITS( I ) = UNITS3D( J )( 1:L )
                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VUNITS

        END SUBROUTINE OPENMRGIN
