
        SUBROUTINE OPENMRGIN( SRGNROWS, SRGNCOLS, SRGGRDNM, SRGFMT )

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
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: LREPANY, GDEV, AFLAG, AENAME, ASDEV, ANMAP, 
     &          AMAPNAM, AMAPFIL, NASRC, APRJFLAG, TFLAG, AFLAG_BD, 
     &          ATNAME, ASDATE, ANIPOL, AEINAM, AONAMES, 
     &          AOUNITS, AGNAME, ANGMAT, SFLAG, ASNAME, ANSMATV, 
     &          ASVDESC, ASVUNIT, AUFLAG, AUNAME, ANUMATV, AUVNAMS, 
     &          ARFLAG, ARNAME, ANRMATV, ANSREAC, ARVDESC, ARNMSPC,
     &          BTNAME, LMETCHK, BIOUNIT, BNSMATV, BNIPOL, BONAMES,
     &          BEINAM, BOUNITS, MFLAG, MENAME, MSDEV, MNMAP, MMAPNAM,
     &          MMAPFIL, NMSRC, MPRJFLAG, MFLAG_BD, MTNAME, MSDATE,
     &          MNIPPA, MEANAM, MONAMES, MOUNITS, MNIPOL, MNIACT,
     &          MGNAME, MNGMAT, MSNAME, MNSMATV, MSVDESC, MSVUNIT,
     &          MUFLAG, MUNAME, MNUMATV, MUVNAMS, MRFLAG, MRNAME,
     &          MNRMATV, MNSREAC, MRVDESC, MRNMSPC, PFLAG, PENAME,
     &          PSDEV, PNMAP, PMAPNAM, PMAPFIL, NPSRC, ELEVFLAG, JSTACK,
     &          PPRJFLAG, PFLAG_BD, PTNAME, PSDATE, PNIPOL, PEINAM, 
     &          PONAMES, POUNITS, PGNAME, PSNAME, PNSMATV, PSVDESC, 
     &          PSVUNIT, PUFLAG, PUNAME, PNUMATV, PUVNAMS, PRFLAG,
     &          PRNAME, PNRMATV, PNSREAC, PRVDESC, PRNMSPC, LFLAG,
     &          PLNAME, EXPLFLAG, PHNAME, EMLAYS, PINGFLAG,
     &          INLINEFLAG, SRCGRPFLAG, SGDEV, EDEV, SUBSECFLAG,
     &          PVNAME, PVSDATE, PVSTIME, PDEV, CDEV, TZONE, SDATE, 
     &          STIME, TSTEP, NSTEPS, EDATE, ETIME, BYEAR, PYEAR,
     &          BSVDESC, BFLAG, VARFLAG, PFACFLAG

C...........  This module contains the information about the source category
        USE MODINFO, ONLY: NMAP, MAPNAM, MAPFIL, INVPIDX

C.........  This module contains arrays for plume-in-grid and major sources
        USE MODELEV, ONLY: NGROUP, NHRSRC, SGFIREFLAG

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM, NCOLS, NROWS, VGTYP, VGTOP, VGLVS

C.........  This module is required for the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C.........  INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C.........  EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        LOGICAL         DSCM3GRD
        LOGICAL         ENVYN
        CHARACTER(50)   GETCFDSC  
        INTEGER         GETIFDSC  
        INTEGER         PROMPTFFILE  
        CHARACTER(16)   PROMPTMFILE  
        INTEGER         SECSDIFF  
        LOGICAL         SETENVVAR
        LOGICAL         USEEXPGEO

        EXTERNAL  CRLF, ENVYN, GETCFDSC, GETIFDSC, PROMPTFFILE, 
     &            PROMPTMFILE, SECSDIFF, SETENVVAR, USEEXPGEO

C...........   Subroutine arguments

        INTEGER      , INTENT (OUT) :: SRGNROWS  ! no. rows in surrogates file
        INTEGER      , INTENT (OUT) :: SRGNCOLS  ! no. columns in surrogates file
        CHARACTER(*) , INTENT (OUT) :: SRGGRDNM  ! name of srgs grid
        CHARACTER(*) , INTENT (OUT) :: SRGFMT    ! gridding surrogates format

C.........  Other local variables

        INTEGER         I, J, K, N, V       ! counters and indices

        INTEGER         IDEV          ! tmp unit number if ENAME is map file
        INTEGER         IOS           ! tmp I/O status
        INTEGER         ISECS         ! tmp duration in seconds
        INTEGER         NPACT         ! no. variables per activity
        INTEGER         NPPOL         ! no. variables per pollutant
        INTEGER         NDIM          ! tmp dimensioning variable 
        INTEGER         NVAR          ! tmp no. variables 

        LOGICAL      :: CFLAG = .FALSE.  ! true: speciation type has been init
        LOGICAL      :: EFLAG = .FALSE.  ! true: error in routine
        LOGICAL      :: IFLAG = .FALSE.  ! true: episode settings have been init
        LOGICAL      :: OFLAG = .FALSE.  ! true: met info has been init
        LOGICAL      :: YFLAG = .FALSE.  ! true: year/projection info been init
        LOGICAL      :: ZFLAG = .FALSE.  ! true: time zone has been init

        CHARACTER(4)    SPCTYPE      ! type of speciation matrix (mass|mole)
        CHARACTER(16)   DUMNAME      ! tmp file name
        CHARACTER(16)   INAME        ! tmp name for inven file of unknown fmt
        CHARACTER(50)   METSCENR     ! met scenario name
        CHARACTER(50)   METCLOUD     ! met cloud scheme name
        CHARACTER(50)   METTMP       ! temporary buffer for met info
        CHARACTER(80)   GDESC        ! grid description
        CHARACTER(256)  MESG         ! message buffer
        CHARACTER(IOVLEN3) COORD3D    ! coordinate system name 
        CHARACTER(IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(IOVLEN3) PROJTYPE   ! projection type

        CHARACTER(16) :: PROGNAME = 'OPENMRGIN' ! program name

C***********************************************************************
C   begin body of subroutine OPENMRGIN

C.........  If reporting state and/or county emissions, and processing for
C           biogenic sources, get gridding surrogates
        IF( ( LREPANY .OR. SRCGRPFLAG ) .AND. BFLAG ) THEN
        
            IF( VARFLAG ) THEN
                MESG = 'Cannot report biogenic emissions when ' //
     &                 'using a variable grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            GDEV = PROMPTFFILE( 
     &             'Enter logical name for SURROGATE COEFFICIENTS file',
     &             .TRUE., .TRUE., 'BGPRO', PROGNAME )

C.............  Read surrogate file header    
            CALL RDSRGHDR( .FALSE., GDEV, SRGFMT )   ! CHKGRID may be initialized here
            SRGGRDNM = GRDNM
            SRGNCOLS = NCOLS
            SRGNROWS = NROWS

        END IF

C.........  Initialize gridded information with grid description file
        IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Check or initialize the grid; do not allow subgrids 
C           when using a variable grid
        IF( VARFLAG ) THEN
            CALL CHKGRID( 'general', 'GRIDDESC', 0, EFLAG )
        ELSE
            CALL CHKGRID( 'general', 'GRIDDESC', 1, EFLAG )
        END IF
        
C.........  For area sources... 
        IF( AFLAG ) THEN

C.............  Get inventory file names given source category
            CALL GETINAME( 'AREA', AENAME, DUMNAME )

C.........  Prompt for and open inventory file 
            MESG= 'Enter logical name for the MAP ' //
     &            'AREA INVENTORY file'
            INAME = AENAME
            IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Read map-formatted inventory file
            CALL RDINVMAP( INAME, IDEV, AENAME, DUMNAME, ASDEV )

C.............  Transfer from MODINFO arrays to MODMERGE arrays
            ALLOCATE( AMAPNAM( NMAP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AMAPNAM', PROGNAME )
            ALLOCATE( AMAPFIL( NMAP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AMAPFIL', PROGNAME )
                
            ANMAP = NMAP
            AMAPNAM = MAPNAM   ! array
            AMAPFIL = MAPFIL   ! array

C.............  Get number of sources
            CALL RETRIEVE_SET_HEADER( AENAME )
            NASRC = NROWS3D

C.............  Determine the year and projection status of the inventory
            CALL CHECK_INVYEAR( AENAME, APRJFLAG, FDESC3D )

C.............  For temporal inputs, prompt for hourly file
            IF( TFLAG ) THEN

C.................  Open all temporal files for either by-day or standard
C                   processing. 
C.................  Compare headers to make sure files are consistent.
                CALL OPEN_TMP_FILES( 'AREA', AFLAG_BD, ATNAME, ASDATE )

C.................  Set pollutants from hourly file
                ANIPOL = NVARSET
                ALLOCATE( AEINAM( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AEINAM', PROGNAME )
                ALLOCATE( AONAMES( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AONAMES', PROGNAME )
                ALLOCATE( AOUNITS( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AOUNITS', PROGNAME )

                CALL STORE_VNAMES( 1, 1, ANIPOL, .TRUE., AEINAM )
                CALL STORE_INVINFO( 1, 1, ANIPOL, 1, INVPIDX, .TRUE.,
     &                              AONAMES, AOUNITS )

C.................  Determine the year and projection status of the hourly
                CALL CHECK_INVYEAR( ATNAME( 1 ), APRJFLAG, FDESC3D )

C.............  Otherwise, just set parameters and pollutants from inven file
            ELSE
                ATNAME = AENAME  ! array
                NVAR   = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
                ANIPOL = GETIFDSC( FDESC3D, '/POLLUTANTS/', .FALSE. )
                NPPOL  = GETIFDSC( FDESC3D, '/PER POLLUTANT/', .FALSE. )

                ALLOCATE( AEINAM ( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AEINAM', PROGNAME )
                ALLOCATE( AONAMES( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AONAMES', PROGNAME )
                ALLOCATE( AOUNITS( ANIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AOUNITS', PROGNAME )

C...............  For map-formatted inventories
                IF( ANMAP .GT. 0 ) THEN
                    CALL STORE_INVEN_VARS( ANMAP,ANIPOL,NPPOL,1+INVPIDX,
     &                     AMAPNAM, AMAPFIL, AEINAM, AONAMES, AOUNITS ) 

C...............  For old format of inventories
                ELSE
                    CALL STORE_VNAMES( NVAR+1, NPPOL, ANIPOL, .TRUE., 
     &                                 AEINAM )
                    CALL STORE_INVINFO( NVAR+1, NPPOL, ANIPOL, 1, 
     &                              INVPIDX, .TRUE. , AONAMES, AOUNITS ) 
                END IF

            END IF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            AGNAME = PROMPTMFILE( 
     &         'Enter logical name for the AREA GRIDDING MATRIX',
     &         FSREAD3, 'AGMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( AGNAME )
            
            IF( VARFLAG ) THEN
                DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', .TRUE. )
            END IF
                
            CALL CHKSRCNO( 'area', 'AGMAT', NTHIK3D, NASRC, EFLAG )

C.............  Check the grid definition; do not allow subgrids if using
C               a variable grid
            IF( VARFLAG ) THEN
                CALL CHKGRID( 'area', 'GMAT', 0, EFLAG )
            ELSE
                CALL CHKGRID( 'area', 'GMAT', 1, EFLAG )
            END IF
            
            ANGMAT = NCOLS3D

C.............  Open speciation matrix, compare number of sources, store
C               speciation variable descriptions, and store mass or moles.
            IF( SFLAG ) THEN
                ASNAME = PROMPTSET( 
     &           'Enter logical name for the AREA SPECIATION MATRIX',
     &           FSREAD3, 'ASMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( ASNAME )
                CALL CHKSRCNO( 'area', 'ASMAT', NROWS3D, NASRC, EFLAG )
                ANSMATV = NVARSET
                ALLOCATE( ASVDESC( ANSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ASVDESC', PROGNAME )
                ALLOCATE( ASVUNIT( ANSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ASVUNIT', PROGNAME )
                CALL STORE_VDESCS( 1, 1, ANSMATV, .TRUE., ASVDESC )
                CALL STORE_VUNITS( 1, 1, ANSMATV, .TRUE., ASVUNIT )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'area' )

            END IF  ! end of speciation open

C.............  Open multiplicative control matrix, compare number of sources, 
C               and store control variable names.
            IF( AUFLAG ) THEN
                MESG = 'Enter logical name for the AREA ' //
     &                 'MULTIPLICATIVE CONTROL MATRIX'
                AUNAME = PROMPTSET( MESG, FSREAD3, 'ACMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( AUNAME )
                CALL CHKSRCNO( 'area', 'ACMAT', NROWS3D, NASRC, EFLAG )
                ANUMATV = NVARS3D

C.................  Supporting "pfac" to apply matrix to all pollutants (group)
                IF( VNAME3D( 1 ) == 'pfac' ) ANUMATV = ANMAP

                ALLOCATE( AUVNAMS( ANUMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'AUVNAMS', PROGNAME )

C.................  Store all pollutant namse to AUVNAMS to support "pfac"
                IF( VNAME3D( 1 ) == 'pfac' ) THEN 
                    AUVNAMS = AMAPNAM
                    PFACFLAG = .TRUE.
                ELSE
                    CALL STORE_VNAMES( 1, 1, ANUMATV, .TRUE., AUVNAMS )
                ENDIF

            END IF  ! end of multiplicative control open

C.............  Open reactivity control matrix, compare number of sources, and
C               store control variable descriptions, and store mass or moles.
            IF( ARFLAG ) THEN
                ARNAME = PROMPTSET( 
     &           'Enter logical name for the AREA REACTIVITY MATRIX',
     &           FSREAD3, 'ARMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( ARNAME )
                CALL CHKSRCNO( 'area', 'ARMAT', NTHIK3D, NASRC, EFLAG )
                ANRMATV = NVARS3D
                ANSREAC = NROWS3D
                ALLOCATE( ARVDESC( ANRMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ARVDESC', PROGNAME )
                CALL STORE_VDESCS( 1, 1, ANRMATV, .FALSE., ARVDESC )

C.................  Retrieve the number of speciation factors 
                ARNMSPC = GETIFDSC( FDESC3D, '/SPECIES VARS/', .TRUE. )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'area' )

C.................  Check the year and projection year of the matrix
                CALL CHECK_INVYEAR( ARNAME, APRJFLAG, FDESC3D )

            END IF  ! end of reactivity control open

        ELSE

            ALLOCATE( AEINAM ( ANIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'AEINAM', PROGNAME )

        END IF  ! End of section for area sources

C.........  If we have biogenic sources 
        IF( BFLAG ) THEN

            MESG = 'Enter logical name for TIME-STEPPED BIOGENIC ' //
     &             'EMISSIONS file'
            BTNAME = PROMPTMFILE( MESG, FSREAD3, 'BGTS', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( BTNAME )

            IF( VARFLAG ) THEN
                DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', .TRUE. )
            END IF
            
            CALL UPDATE_TIME_INFO( BTNAME )
            
C.............  Check the grid definition; do not allow subgrids if using
C               a variable grid
            IF( VARFLAG ) THEN
                CALL CHKGRID( 'biogenics', 'GRID', 0, EFLAG )
            ELSE
                CALL CHKGRID( 'biogenics', 'GRID', 1, EFLAG )
            END IF

            IF( LMETCHK ) CALL CHECK_MET_INFO( 'biogenics' ) 

C.............  Store biogenic species names as if they are stored as
C               speciation matrix variable descriptions
            BIOUNIT = UNITS3D( 1 )
            BNSMATV = NVARS3D
            ALLOCATE( BSVDESC( BNSMATV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'BSVDESC', PROGNAME )

            DO I = 1, BNSMATV
                IF( VNAME3D( I ) .NE. 'NO' ) THEN
                    DUMNAME = 'BIO' // ETJOIN // 'VOC'
                    BSVDESC( I ) = DUMNAME // SPJOIN // VNAME3D( I )
                ELSE
                    DUMNAME = 'BIO' // ETJOIN // 'NOX'
                    BSVDESC( I ) = DUMNAME // SPJOIN // VNAME3D( I )
                END IF

            END DO

C.............  Store biogenic pollutant names and units
            BNIPOL  = 2
            BEINAM ( 1 ) = 'BIO' // ETJOIN // 'NOX'
            BEINAM ( 2 ) = 'BIO' // ETJOIN // 'VOC'
            BONAMES      = BEINAM
            BOUNITS      = UNITS3D( 1 )  ! array

        END IF  ! End of section for biogenic sources

C.........  If we have mobile sources 
        IF( MFLAG ) THEN

C.............  Get inventory file names given source category
            CALL GETINAME( 'MOBILE', MENAME, DUMNAME )

C.........  Prompt for and open inventory file 
            MESG= 'Enter logical name for the MAP ' //
     &            'MOBILE INVENTORY file'
            INAME = MENAME
            IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Read map-formatted inventory file
            CALL RDINVMAP( INAME, IDEV, MENAME, DUMNAME, MSDEV )

C.............  Transfer from MODINFO arrays to MODMERGE arrays
            ALLOCATE( MMAPNAM( NMAP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MMAPNAM', PROGNAME )
            ALLOCATE( MMAPFIL( NMAP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MMAPFIL', PROGNAME )
                
            MNMAP = NMAP
            MMAPNAM = MAPNAM   ! array
            MMAPFIL = MAPFIL   ! array

C.............  Get number of sources
            CALL RETRIEVE_SET_HEADER( MENAME )
            NMSRC = NROWS3D

C.............  Determine the year and projection status of the inventory
            CALL CHECK_INVYEAR( MENAME, MPRJFLAG, FDESC3D )

C.............  For temporal inputs, prompt for hourly file
            IF( TFLAG ) THEN

C.................  Open all temporal files for either by-day or standard
C                   processing. 
C.................  Compare headers to make sure files are consistent.
                CALL OPEN_TMP_FILES( 'MOBILE', MFLAG_BD, MTNAME, MSDATE)

                MNIPPA = NVARSET
                ALLOCATE( MEANAM( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MEANAM', PROGNAME )
                ALLOCATE( MONAMES( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MONAMES', PROGNAME )
                ALLOCATE( MOUNITS( MNIPPA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MOUNITS', PROGNAME )

                CALL STORE_VNAMES( 1, 1, MNIPPA, .TRUE., MEANAM )
                CALL STORE_INVINFO( 1, 1, MNIPPA, 1, INVPIDX, .TRUE., 
     &                              MONAMES, MOUNITS )

C.................  Determine the year and projection status of the hourly
                CALL CHECK_INVYEAR( MTNAME( 1 ), MPRJFLAG, FDESC3D )

C.............  Otherwise, just set parameters and pollutants from inven file
            ELSE
                MTNAME = MENAME  ! array

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

C...............  For map-formatted inventories
                IF( MNMAP .GT. 0 ) THEN
                    CALL STORE_INVEN_VARS( MNMAP-MNIACT, MNIPOL, NPPOL,
     &                             1+INVPIDX, MMAPNAM, MMAPFIL, MEANAM,
     &                                               MONAMES, MOUNITS ) 

                    K = MNIPOL + 1
                    CALL STORE_INVEN_VARS( MNMAP-MNIPOL, MNIACT,NPACT,2,
     &                     MMAPNAM( K ), MMAPFIL( K ), MEANAM( K ), 
     &                     MONAMES( K ), MOUNITS( K ) ) 

C...............  For old format of inventories
                ELSE

C.................  Store pollutant names and other information
                    CALL STORE_VNAMES( NVAR+1, NPPOL, MNIPOL, .TRUE., 
     &                                 MEANAM )
                    CALL STORE_INVINFO( NVAR+1, NPPOL, MNIPOL, 1, 
     &                               INVPIDX, .TRUE., MONAMES, MOUNITS )

C.................  Store activity names and other information
                    I    = MNIPOL * NPPOL
                    NVAR = NVAR + I
                    CALL STORE_VNAMES( NVAR+1, NPACT, MNIACT, .TRUE., 
     &                                 MEANAM( MNIPOL+1 )     )
                    CALL STORE_INVINFO( NVAR+1, NPACT, MNIACT, 1, 1, 
     &                                  .TRUE.,  MONAMES( MNIPOL+1 ), 
     &                                  MOUNITS( MNIPOL+1 )           ) 
                END IF

            END IF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            MGNAME = PROMPTMFILE( 
     &       'Enter logical name for the MOBILE GRIDDING MATRIX',
     &       FSREAD3, 'MGMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( MGNAME )
            
            IF( VARFLAG ) THEN
                DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', .TRUE. )
            END IF
            
            CALL CHKSRCNO( 'mobile', 'MGMAT', NTHIK3D, NMSRC, EFLAG )

C.............  Check the grid definition; do not allow subgrids if using
C               a variable grid
            IF( VARFLAG ) THEN
                CALL CHKGRID( 'mobile', 'GMAT', 0, EFLAG )
            ELSE
                CALL CHKGRID( 'mobile', 'GMAT', 1, EFLAG )
            END IF
            
            MNGMAT = NCOLS3D

C.............  Open speciation matrix, compare number of sources, store
C               speciation variable descriptions, and store mass or moles.
            IF( SFLAG ) THEN
                MSNAME = PROMPTSET( 
     &           'Enter logical name for the MOBILE SPECIATION MATRIX',
     &           FSREAD3, 'MSMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( MSNAME )
                CALL CHKSRCNO( 'mobile', 'MSMAT', NROWS3D, NMSRC, EFLAG)
                MNSMATV = NVARSET
                ALLOCATE( MSVDESC( MNSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MSVDESC', PROGNAME )
                ALLOCATE( MSVUNIT( MNSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MSVUNIT', PROGNAME )
                CALL STORE_VDESCS( 1, 1, MNSMATV, .TRUE., MSVDESC )
                CALL STORE_VUNITS( 1, 1, MNSMATV, .TRUE., MSVUNIT )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'mobile' )

            END IF  ! end of speciation open

C.............  Open multiplicative control matrix, compare number of sources, 
C               and store control variable names.
            IF( MUFLAG ) THEN
                MESG = 'Enter logical name for the MOBILE ' //
     &                 'MULTIPLICATIVE CONTROL MATRIX'
                MUNAME = PROMPTSET( MESG, FSREAD3, 'MCMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( MUNAME )
                CALL CHKSRCNO( 'mobile', 'MCMAT', NROWS3D, NMSRC, EFLAG)
                MNUMATV = NVARS3D

C.................  Supporting "pfac" to apply matrix to all pollutants (group)
                IF( VNAME3D( 1 ) == 'pfac' ) MNUMATV = MNMAP

                ALLOCATE( MUVNAMS( MNUMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MUVNAMS', PROGNAME )

C.................  Store all pollutant namse to MUVNAMS to support "pfac"
                IF( VNAME3D( 1 ) == 'pfac' ) THEN 
                    MUVNAMS = MMAPNAM
                    PFACFLAG = .TRUE.
                ELSE
                    CALL STORE_VNAMES( 1, 1, MNUMATV, .TRUE., MUVNAMS )
                ENDIF

            END IF  ! end of multiplicative control open

C.............  Open reactivity control matrix, compare number of sources, and
C               store control variable descriptions, and store mass or moles.
            IF( MRFLAG ) THEN
                MRNAME = PROMPTSET( 
     &           'Enter logical name for the MOBILE REACTIVITY MATRIX',
     &           FSREAD3, 'MRMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( ARNAME )
                CALL CHKSRCNO( 'mobile', 'MRMAT', NTHIK3D, NMSRC, EFLAG)
                MNRMATV = NVARS3D
                MNSREAC = NROWS3D
                ALLOCATE( MRVDESC( MNRMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'MRVDESC', PROGNAME )
                CALL STORE_VDESCS( 1, 1, MNRMATV, .FALSE., MRVDESC )

C.................  Retrieve the number of speciation factors 
                MRNMSPC = GETIFDSC( FDESC3D, '/SPECIES VARS/', .TRUE. )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'mobile' )

C.................  Check the year and projection year of the matrix
                CALL CHECK_INVYEAR( MRNAME, MPRJFLAG, FDESC3D )

            END IF  ! end of reactivity control open

        ELSE

            ALLOCATE( MEANAM ( MNIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MEANAM', PROGNAME )

        ENDIF  ! End of section for mobile sources

C.........  If we have point sources 
        IF( PFLAG ) THEN

C.............  Get inventory file names given source category
            CALL GETINAME( 'POINT', PENAME, DUMNAME )

C.........  Prompt for and open inventory file 
            MESG= 'Enter logical name for the MAP ' //
     &            'POINT INVENTORY file'
            INAME = PENAME
            IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

C.........  Read map-formatted inventory file
            CALL RDINVMAP( INAME, IDEV, PENAME, DUMNAME, PSDEV )

C.............  Transfer from MODINFO arrays to MODMERGE arrays
            ALLOCATE( PMAPNAM( NMAP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PMAPNAM', PROGNAME )
            ALLOCATE( PMAPFIL( NMAP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PMAPFIL', PROGNAME )
                
            PNMAP = NMAP
            PMAPNAM = MAPNAM   ! array
            PMAPFIL = MAPFIL   ! array

C.............  Get number of sources
            CALL RETRIEVE_SET_HEADER( PENAME )
            NPSRC = NROWS3D

C.............  If outputing ASCII elevated sources, retrieve the position
C               of the stack in the source characteristics
            IF( ELEVFLAG ) 
     &        JSTACK   = GETIFDSC( FDESC3D, '/STACK POSITION/', .TRUE. )

C.............  Determine the year and projection status of the inventory
            CALL CHECK_INVYEAR( PENAME, PPRJFLAG, FDESC3D )

C.............  For temporal inputs, prompt for hourly file
            IF( TFLAG ) THEN

C.................  Open all temporal files for either by-day or standard
C                   processing. 
C.................  Compare headers to make sure files are consistent.
                CALL OPEN_TMP_FILES( 'POINT', PFLAG_BD, PTNAME, PSDATE)

                PNIPOL = NVARSET
                ALLOCATE( PEINAM( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PEINAM', PROGNAME )
                ALLOCATE( PONAMES( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PONAMES', PROGNAME )
                ALLOCATE( POUNITS( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'POUNITS', PROGNAME )

                CALL STORE_VNAMES( 1, 1, PNIPOL, .TRUE., PEINAM )
                CALL STORE_INVINFO( 1, 1, PNIPOL, 1, INVPIDX, .TRUE., 
     &                              PONAMES, POUNITS )

C.................  Determine the year and projection status of the hourly 
                CALL CHECK_INVYEAR( PTNAME( 1 ), PPRJFLAG, FDESC3D )

C.............  Otherwise, just set parameters and pollutants from inven file
            ELSE
                PTNAME = PENAME  ! array
                NVAR   = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
                PNIPOL = GETIFDSC( FDESC3D, '/POLLUTANTS/', .FALSE. )
                NPPOL  = GETIFDSC( FDESC3D, '/PER POLLUTANT/', .FALSE. )

                ALLOCATE( PEINAM( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PEINAM', PROGNAME )
                ALLOCATE( PONAMES( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PONAMES', PROGNAME )
                ALLOCATE( POUNITS( PNIPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'POUNITS', PROGNAME )

C...............  For map-formatted inventories
                IF( PNMAP .GT. 0 ) THEN
                    CALL STORE_INVEN_VARS( PNMAP,PNIPOL,NPPOL,1+INVPIDX,
     &                     PMAPNAM, PMAPFIL, PEINAM, PONAMES, POUNITS ) 

C...............  For old format of inventories
                ELSE

                    CALL STORE_VNAMES( NPTVAR3+1, NPTPPOL3, PNIPOL,  
     &                                 .TRUE., PEINAM )
                    CALL STORE_INVINFO( NVAR+1, NPPOL, PNIPOL, 1,
     &                               INVPIDX, .TRUE., PONAMES, POUNITS )
                END IF

            END IF

C.............  Open gridding matrix, compare number of sources, and 
C               compare or initialize grid information.
            PGNAME = PROMPTMFILE( 
     &       'Enter logical name for the POINT GRIDDING MATRIX',
     &       FSREAD3, 'PGMAT', PROGNAME )

            CALL RETRIEVE_IOAPI_HEADER( PGNAME )
            
            IF( VARFLAG ) THEN
                DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', .TRUE. )
            END IF
            
            CALL CHKSRCNO( 'point', 'PGMAT', NTHIK3D, NPSRC, EFLAG )

C.............  Check the grid definition; do not allow subgrids if using
C               a variable grid
            IF( VARFLAG ) THEN
                CALL CHKGRID( 'point', 'GMAT', 0, EFLAG )
            ELSE
                CALL CHKGRID( 'point', 'GMAT', 1, EFLAG )
            END IF

C.............  Open speciation matrix, compare number of sources, store
C               speciation variable names, and store mass or moles.
            IF( SFLAG ) THEN
                PSNAME = PROMPTSET( 
     &           'Enter logical name for the POINT SPECIATION MATRIX',
     &           FSREAD3, 'PSMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( PSNAME )
                CALL CHKSRCNO( 'point','PSMAT',NROWS3D,NPSRC,EFLAG )
                PNSMATV = NVARSET
                ALLOCATE( PSVDESC( PNSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PSVDESC', PROGNAME )
                ALLOCATE( PSVUNIT( PNSMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PSVUNIT', PROGNAME )
                CALL STORE_VDESCS( 1, 1, PNSMATV, .TRUE., PSVDESC )
                CALL STORE_VUNITS( 1, 1, PNSMATV, .TRUE., PSVUNIT )

C.................  Ensure consistent spec matrix type for all source categories
                CALL CHECK_SPEC_TYPE( 'point' )

            END IF  ! end of speciation open

C.............  Open multiplicative control matrix, compare number of sources, 
C               and store control variable names.
            IF( PUFLAG ) THEN
                MESG = 'Enter logical name for the POINT ' //
     &                 'MULTIPLICATIVE CONTROL MATRIX'
                PUNAME = PROMPTSET( MESG, FSREAD3, 'PCMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( PUNAME )
                CALL CHKSRCNO( 'point', 'PCMAT', NROWS3D, NPSRC, EFLAG )
                PNUMATV = NVARS3D

C.................  Supporting "pfac" to apply matrix to all pollutants (group)
                IF( VNAME3D( 1 ) == 'pfac' ) PNUMATV = PNMAP

                ALLOCATE( PUVNAMS( PNUMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PUVNAMS', PROGNAME )

C.................  Store all pollutant namse to PUVNAMS to support "pfac"
                IF( VNAME3D( 1 ) == 'pfac' ) THEN 
                    PUVNAMS = PMAPNAM
                    PFACFLAG = .TRUE.
                ELSE
                    CALL STORE_VNAMES( 1, 1, PNUMATV, .TRUE., PUVNAMS )
                ENDIF

            END IF  ! end of multiplicative control open

C.............  Open reactivity control matrix, compare number of sources, and
C               store control variable descriptions, and store mass or moles.
            IF( PRFLAG ) THEN
                PRNAME = PROMPTSET( 
     &           'Enter logical name for the POINT REACTIVITY MATRIX',
     &           FSREAD3, 'PRMAT', PROGNAME )

                CALL RETRIEVE_SET_HEADER( PRNAME )
                CALL CHKSRCNO( 'point', 'PRMAT', NTHIK3D, NPSRC, EFLAG )
                PNRMATV = NVARS3D
                PNSREAC = NROWS3D
                ALLOCATE( PRVDESC( PNRMATV ), STAT=IOS )
                CALL CHECKMEM( IOS, 'PRVDESC', PROGNAME )
                CALL STORE_VDESCS( 1, 1, PNRMATV, .FALSE., PRVDESC )

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
                
                IF( VARFLAG ) THEN
                    DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', .TRUE. )
                END IF
                
                CALL CHKSRCNO( 'point', PLNAME, NROWS3D, NPSRC, EFLAG )
                CALL UPDATE_TIME_INFO( PLNAME )

                IF( LMETCHK ) CALL CHECK_MET_INFO( 'point' ) 

C.............  Get file name and open daily input inventory file
            ELSE IF( EXPLFLAG ) THEN

                MESG = 'Enter logical name for EXPLICIT LAYER ' //
     &                 'FRACTIONS MATRIX'
                PHNAME = PROMPTMFILE( MESG,FSREAD3,'PLAY_EX',PROGNAME )

C.................  Check to see if appropriate variable list exists
                CALL RETRIEVE_IOAPI_HEADER( PHNAME )
                
                IF( VARFLAG ) THEN
                    DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', .TRUE. )
                END IF
                
                CALL CHKSRCNO( 'point', PHNAME, NTHIK3D, NPSRC, EFLAG )
                CALL UPDATE_TIME_INFO( PHNAME )

                NHRSRC = NROWS3D

            END IF             ! End of layer fractions open

C.............  If either PLAY or PLAY_EX were just opened, store vertical
C               layer information from the header of these files.
            IF ( LFLAG .OR. EXPLFLAG ) THEN
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

            END IF

C.............  For plume-in-grid outputs or for UAM-style elevated point
C               sources (PTSRCE input file)...
            IF( ELEVFLAG .OR. PINGFLAG .OR. INLINEFLAG ) THEN

C.................  If elevated ASCII and units are grams, print warning
                IF( ELEVFLAG .AND. SPCTYPE .EQ. MASSSTR ) THEN
                    MESG = 'WARNING: Processing with mass-based ' //
     &                     'speciation for elevated ASCII outputs.'
                    CALL M3MSG2( MESG )
                END IF

C.................  Open elevated/plume-in-grid file 
                MESG = 'Enter logical name for the ELEVATED/PING ' //
     &                 'file'
                EDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                              'PELV', PROGNAME      )

C.................  Open stack groups file output from Elevpoint
                MESG = 'Enter logical name for the ELEVATED STACK ' //
     &                 'GROUPS file'
                PVNAME = PROMPTMFILE( MESG, FSREAD3, 'STACK_GROUPS', 
     &                   PROGNAME )
                CALL RETRIEVE_IOAPI_HEADER( PVNAME )

C.................  Check if stack groups file has fire data
                DO I = 1, NVARS3D
                    IF( VNAME3D( I ) .EQ. 'ACRESBURNED' ) THEN
                        SGFIREFLAG = .TRUE.
                        EXIT
                    END IF
                END DO
                  
                IF( VARFLAG ) THEN
                    DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', 
     &                                  .TRUE. )
                END IF
                
                NGROUP = NROWS3D
                PVSDATE = SDATE3D
                PVSTIME = STIME3D
                
C.................  Check grid definition; do not allow subgrids when
C                   using a variable grid
                IF( VARFLAG ) THEN
                    CALL CHKGRID( 'point', 'GROUPS', 0, EFLAG )
                ELSE
                    CALL CHKGRID( 'point', 'GROUPS', 1, EFLAG )
                END IF

            END IF

        ELSE

            ALLOCATE( PEINAM ( PNIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PEINAM', PROGNAME )

        END IF      ! End of section for point sources

C.........  Get file name for inventory pollutants codes/names
        MESG = 'Enter logical name for INVENTORY DATA TABLE file'
        PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'INVTABLE',
     &                      PROGNAME )

C.........  Get country, state, and county names no matter what, because it is
C           needed to allocate memory for the state and county totals, even
C           when they aren't going to be output
        IF( .NOT. USEEXPGEO() ) THEN
            CDEV = PROMPTFFILE( 
     &             'Enter logical name for COUNTRY, STATE, AND ' //
     &             'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )
        END IF

C.........  Open source groups file if needed
        IF( SRCGRPFLAG ) THEN
            MESG = 'Enter logical name for SOURCE GROUPS file'
            SGDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 
     &                           'SOURCE_GROUPS', PROGNAME )
        END IF

C.........  Open sub-sector source groups file if needed
        IF( SUBSECFLAG ) THEN
            MESG = 'Enter logical name for SUB-SECTOR SOURCE GROUPS file'
            SGDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,
     &                           'SUB_SEC_SOURCES', PROGNAME )
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

            CALL GETM3EPI( TZONE, SDATE, STIME, TSTEP, NSTEPS )
            TSTEP = 10000   ! only 1-hour time steps supported
            EDATE = SDATE
            ETIME = STIME
            CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

        END IF   !  if have temporalized inputs and outputs

C.........  Compare base year with episode and warn if not consistent
        IF( BYEAR .NE. 0 .AND. SDATE / 1000 .NE. BYEAR ) THEN

            WRITE( MESG,94010 ) 'WARNING: Inventory base year ', BYEAR, 
     &             'is inconsistent with year ' // CRLF() // BLANK10 //
     &             'of episode start date', SDATE/1000
            CALL M3MSG2( MESG )

        ENDIF

C.........  Give a note if running for a projected year
        IF( PYEAR .NE. BYEAR ) THEN

            WRITE( MESG,94010 ) 'NOTE: Emissions based on projected '//
     &             'year', PYEAR
            CALL M3MSG2( MESG )

        END IF

C.........  Write message stating grid name and description
        N = LEN_TRIM( GRDNM )
        MESG = 'NOTE: Output grid "' // GRDNM( 1:N ) // 
     &         '" set; described as' // CRLF() // BLANK10 // GDESC
        CALL M3MSG2( MESG )


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
C.............  This subprogram tries to retrieve the description for a file
C               set and aborts if it was not successful
            SUBROUTINE RETRIEVE_SET_HEADER( FILNAM )

            INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

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
C.............  This subprogram opens the temporal emissions files. If their 
C               are multiple files, it compares the files to make sure that they
C               are consistent with each other.  The number of sources
C               are compared to the master number of sources.
            SUBROUTINE OPEN_TMP_FILES( LOCCAT, LBDSTAT, FNAME, SDATE )

C.............  Subprogram arguments
            CHARACTER(*), INTENT (IN) :: LOCCAT
            LOGICAL     , INTENT (IN) :: LBDSTAT
            CHARACTER(*), INTENT(OUT) :: FNAME( 7 )
            INTEGER     , INTENT(OUT) :: SDATE( 7 )

C.............  Local parameters
            CHARACTER(3), PARAMETER :: SUFFIX( 7 ) = 
     &                                ( / 'MON', 'TUE', 'WED', 'THU', 
     &                                    'FRI', 'SAT', 'SUN'        / )

C.............  Local allocatable arrays
            CHARACTER(IOVLEN3), ALLOCATABLE :: LOCVNAM ( : )
            CHARACTER(IOULEN3), ALLOCATABLE :: LOCVUNIT( : )

C.............  Local arrays
            INTEGER        IDX( 7 )     ! index for per-file arrays

C.............  Local variables
            INTEGER        D, L, N      ! counters and indices

            INTEGER        INVPIDX   ! tmp index for average day or not
            INTEGER        LOCZONE   ! tmp time zone
            INTEGER        LOCNVAR   ! tmp local number of variables in file 
            INTEGER        NFILE     ! no. hourly emission files

            LOGICAL     :: NFLAG = .FALSE.  ! true: no. vars inconsistent
            LOGICAL     :: VFLAG = .FALSE.  ! true: var names inconsistent
            LOGICAL     :: UFLAG = .FALSE.  ! true: var units inconsistent

            CHARACTER      CRL      ! 1-letter src category indicator
            CHARACTER(16)  TMPNAM   ! temporary logical file name
            CHARACTER(300) MESG     ! message buffer

C----------------------------------------------------------------------

            IF( LOCCAT .EQ. 'AREA'   ) CRL = 'A'
            IF( LOCCAT .EQ. 'MOBILE' ) CRL = 'M'
            IF( LOCCAT .EQ. 'POINT'  ) CRL = 'P'

C.............  Set the number of files and open the files...
C.............  For by-day processing...
            IF( LBDSTAT ) THEN
                NFILE = 7

                DO D = 1, NFILE

                    MESG = 'Enter logical name for the ' // SUFFIX( D )
     &                     // ' ' // LOCCAT // ' HOURLY EMISSIONS file'
                    TMPNAM = CRL // 'TMP_' // SUFFIX( D )

                    FNAME( D ) = PROMPTSET( MESG,FSREAD3,
     &                                        TMPNAM,PROGNAME )
                    IDX( D ) = D
                END DO

C.............  For standard processing...
            ELSE
                NFILE = 1

                MESG = 'Enter logical name for the ' // LOCCAT // 
     &                 ' HOURLY EMISSIONS file'
                TMPNAM = CRL // 'TMP'

                FNAME = PROMPTSET( MESG,FSREAD3,TMPNAM,PROGNAME ) ! array
                IDX( NFILE ) = 1

            END IF

C.............  Loop through each file and ensure they are consistent
            DO D = 1, NFILE

                TMPNAM = FNAME( IDX( D ) )

C.................  Get header and compare source number and time range
                CALL RETRIEVE_SET_HEADER( TMPNAM )

C.................  Store the starting date
                SDATE( IDX( D ) ) = SDATE3D

C.................  Check the number of sources
                SELECT CASE( LOCCAT )
                CASE( 'AREA' )
                    CALL CHKSRCNO( 'area', TMPNAM, NROWS3D, 
     &                             NASRC, EFLAG )
                CASE( 'MOBILE' ) 
                    CALL CHKSRCNO( 'mobile', TMPNAM, NROWS3D, 
     &                             NMSRC, EFLAG )

                CASE( 'POINT' ) 
                    CALL CHKSRCNO( 'point', TMPNAM, NROWS3D, 
     &                             NPSRC, EFLAG )

                END SELECT

C.................  Determine average day emissions status from hourly file
                INVPIDX = GETIFDSC( FDESC3D, '/AVERAGE DAY/', .FALSE. )
                IF( INVPIDX .EQ. 1 ) THEN
                    MESG = 'NOTE: Average day emissions in ' //
     &                     LOCCAT // ' hourly emissions file'
                    CALL M3MSG2( MESG )
                END IF

C.................  For standard processing, compare time info to master
                IF( .NOT. LBDSTAT .AND. D .EQ. 1 ) THEN
                    CALL UPDATE_TIME_INFO( TMPNAM )
                END IF

C.................  For by-day files, make sure that the file starts at hour 0
                IF( LBDSTAT .AND. STIME3D .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    L = LEN_TRIM( TMPNAM )
                    WRITE( MESG,94010 ) 'ERROR: Start time of', STIME3D,
     &                     'in file "'// TMPNAM( 1:L ) // 
     &                     '" is invalid.' // CRLF() // BLANK10 //
     &                     'Only start time of 000000 is valid for' //
     &                     'processing by day.'
                    CALL M3MSG2( MESG )

                END IF

C.................  Make sure that the file has at least 24 hours 
                IF( LBDSTAT .AND. MXREC3D .LT. 24 ) THEN
                    EFLAG = .TRUE.
                    L = LEN_TRIM( TMPNAM )
                    WRITE( MESG,94010 ) 'ERROR: Number of hours', 
     &                     MXREC3D, 'in file "'// TMPNAM( 1:L ) // 
     &                     '" is invalid.' // CRLF() // BLANK10 //
     &                     'Minimum number of 24 hours is needed for' //
     &                     'processing by day.'
                    CALL M3MSG2( MESG )

                END IF

                LOCZONE = GETIFDSC( FDESC3D, '/TZONE/', .TRUE. )

                IF( ZFLAG .AND. LOCZONE .NE. TZONE ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                 'Time zone ', LOCZONE, 'in ' // TMPNAM // 
     &                 ' hourly emissions file is not consistent ' //
     &                 'with initialized value of', TZONE
                    CALL M3MSG2( MESG )

                ELSE IF( .NOT. ZFLAG ) THEN
                    ZFLAG = .TRUE.
                    TZONE = LOCZONE

                    MESG = 'NOTE: Time zone initialized using ' // 
     &                     TMPNAM // ' hourly emissions file.'

                    CALL M3MSG2( MESG )
                END IF

C.................  For first file, store the pollutant names and units for
C                   making comparisons with other files.
                IF( D .EQ. 1 ) THEN

                    LOCNVAR = NVARS3D
                    ALLOCATE( LOCVNAM( LOCNVAR ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'LOCVNAM', PROGNAME )
                    ALLOCATE( LOCVUNIT( LOCNVAR ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'LOCVUNIT', PROGNAME )

                    LOCVNAM ( 1:LOCNVAR ) = VNAMESET( 1:LOCNVAR )
                    LOCVUNIT( 1:LOCNVAR ) = VUNITSET( 1:LOCNVAR )

C.................  Compare the pollutant names and units
                ELSE

C.....................  Check to make sure the number is consistent first
                    IF( NVARSET .NE. LOCNVAR ) NFLAG = .TRUE.

C.....................  Make sure no overflows                    
                    N = MIN( NVARSET, LOCNVAR )

C.....................  compare variable names and units among files
                    DO V = 1, N
                        IF( LOCVNAM( V ) .NE. VNAMESET( V ) ) THEN
                            VFLAG = .TRUE.
                        END IF

                        IF( LOCVUNIT( V ) .NE. VUNITSET( V ) ) THEN
                            UFLAG = .TRUE.
                        END IF
                    END DO

                END IF

            END DO

C.............  Write message and set error if any inconsistencies
            IF( NFLAG ) THEN
c bbh               EFLAG = .TRUE.  ! removed to prevent false errer of odd nubmer
c                                     of species for tmp files
                MESG = 'WARNING: ' // LOCCAT // ' source hourly ' //
     &                 'emission files have inconsistent ' //
     &                 CRLF() // BLANK10 // 'number of variables.'
                CALL M3MSG2( MESG )
            END IF

            IF( VFLAG ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: ' // LOCCAT // ' source hourly ' //
     &                 'emission files have inconsistent ' //
     &                 CRLF() // BLANK10 // 'variable names.'
                CALL M3MSG2( MESG )
            END IF

            IF( UFLAG ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: ' // LOCCAT // ' source hourly ' //
     &                 'emission files have inconsistent ' //
     &                 CRLF() // BLANK10 // 'variable units.'
                CALL M3MSG2( MESG )
            END IF

C.............  Deallocate local memory
            DEALLOCATE( LOCVNAM, LOCVUNIT )

            RETURN

C------------------  FORMAT  STATEMENTS   -----------------------------

C...........   Internal buffering formats.............94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE OPEN_TMP_FILES

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
            INTEGER           L
            INTEGER           YY      ! tmp year
            LOGICAL           STRICT  ! flag for strict checks or not
            CHARACTER(20)     BUFFER  ! program name buffer
            INTEGER,  SAVE :: FLEN    ! name length of savnam
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

C.........................  If there is projection, abort
                        IF ( PRJFLAG ) THEN
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........................  Otherwise, make it a warning
                        ELSE
                            L = LEN_TRIM( MESG )
                            MESG = 'WARNING: ' // MESG( 1:L )
                            CALL M3MSG2( MESG )
                        END IF

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

                IF( PYEAR .GT. 0 ) THEN
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

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_INVYEAR

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram sets the speciation type and compares to 
C               the existing information, if it has been previously set.
            SUBROUTINE CHECK_SPEC_TYPE( CATDESC )

C.............  Subprogram arguments
            CHARACTER(*) CATDESC  ! category descriptions

C.............  Local variables
            CHARACTER(4) LOCTYPE  ! tmp speciation type 

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
            CHARACTER(30) FILDESC    ! description of input file

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
            SUBROUTINE STORE_VNAMES( ISTART, INCRMT, NNAM, LFSET, NAMES)

C.............  Subprogram arguments
            INTEGER,      INTENT (IN) :: ISTART ! starting position in VNAMES of names
            INTEGER,      INTENT (IN) :: INCRMT ! increment of VNAMES for names
            INTEGER,      INTENT (IN) :: NNAM   ! number of names
            LOGICAL,      INTENT (IN) :: LFSET  ! true: input file is FileSetAPI
            CHARACTER(*), INTENT(OUT) :: NAMES( NNAM ) ! stored variable names

C.............  Local variables
            INTEGER  I, J

C----------------------------------------------------------------------

            J = ISTART
            DO I = 1, NNAM

                IF( LFSET ) THEN    ! From FileSetAPI
                    NAMES( I ) = VNAMESET( J )
                ELSE                ! From standard I/O API
                    NAMES( I ) = VNAME3D( J )
                END IF

                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VNAMES

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores I/O API NetCDF variable names and
C               units from the inventory file into a local array based 
C               on indices in subprogram call.
            SUBROUTINE STORE_INVINFO( ISTART, NPER, NPOA, IDX1, IDX2,
     &                                LFSET, NAMES, UNITS )

C.............  Subprogram arguments
            INTEGER,       INTENT (IN) :: ISTART ! starting position in VNAMES of names
            INTEGER,       INTENT (IN) :: NPER   ! no. variables per pollutant
            INTEGER,       INTENT (IN) :: NPOA   ! number of pollutants or activities
            INTEGER,       INTENT (IN) :: IDX1   ! start index for output variables
            INTEGER,       INTENT (IN) :: IDX2   ! index to which pol-assoc variable
            LOGICAL,       INTENT (IN) :: LFSET  ! true: input file is FileSetAPI
            CHARACTER(*), INTENT (OUT) :: NAMES( NPOA ) ! stored variable names
            CHARACTER(*), INTENT (OUT) :: UNITS( NPOA ) ! stored variable units

C.............  Local variables
            INTEGER  I, J, K

C----------------------------------------------------------------------

            J = ISTART + IDX2 - NPER - 1
            K = IDX1 - 1
            DO I = 1, NPOA

                J = J + NPER
                K = K + 1

C................  If input file is from the FileSetAPI
                IF( LFSET ) THEN
                    NAMES( K ) = VNAMESET( J )
                    UNITS( K ) = VUNITSET( J )

C................  If input file is from the standard I/O
                ELSE
                    NAMES( K ) = VNAME3D( J )
                    UNITS( K ) = UNITS3D( J )
                END IF

            END DO
 
            END SUBROUTINE STORE_INVINFO

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores I/O API NetCDF variable descriptions into
C               a local array based on indices in subprogram call.
            SUBROUTINE STORE_VDESCS( ISTART,INCRMT,NDESC,LFSET,DESCS )

            INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: ISTART   ! starting position in VDESCS of names
            INTEGER     , INTENT (IN) :: INCRMT   ! increment of VDESCS for names
            INTEGER     , INTENT (IN) :: NDESC    ! number of descriptions
            LOGICAL     , INTENT (IN) :: LFSET    ! number of descriptions
            CHARACTER(*), INTENT(OUT) :: DESCS( NDESC )! stored variable descriptions

C.............  Local variables
            INTEGER  I, J

C----------------------------------------------------------------------

            DESCS = ' '

            J = ISTART
            DO I = 1, NDESC

                IF( LFSET ) THEN  ! From FileSetAPI
                    DESCS( I ) = TRIM( VDESCSET( J ) )
                ELSE              ! From standard I/O API
                    DESCS( I ) = TRIM( VDESC3D( J ) )
                END IF

                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VDESCS

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores I/O API NetCDF variable units into
C               a local array based on indices in subprogram call.
            SUBROUTINE STORE_VUNITS( ISTART,INCRMT,NUNIT,LFSET,UNITS )

            INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: ISTART        ! starting position in VDESCS of names
            INTEGER     , INTENT (IN) :: INCRMT        ! increment of VDESCS for names
            INTEGER     , INTENT (IN) :: NUNIT         ! number of units
            LOGICAL     , INTENT (IN) :: LFSET         ! number of descriptions
            CHARACTER(*), INTENT(OUT) :: UNITS( NUNIT )! stored variable units

C.............  Local variables
            INTEGER  I, J, L

C----------------------------------------------------------------------

            UNITS = ' '

            J = ISTART
            DO I = 1, NUNIT

                IF( LFSET ) THEN  ! From FileSetAPI
                    UNITS( I ) = TRIM( VUNITSET( J ) )
                ELSE              ! From standard I/O API
                    UNITS( I ) = TRIM( UNITS3D( J ) )
                END IF 

                J = J + INCRMT

            END DO
 
            END SUBROUTINE STORE_VUNITS

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C.............  This subprogram stores pollutant and activity names and
C               units using the map-formatted inventory information 
            SUBROUTINE STORE_INVEN_VARS( NMAPVAR, NHDRVAR, NPPOL, IDX2, 
     &                                   MAPVARS, MAPFILES, NAMES1, 
     &                                   NAMES2, UNITS )

            INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.............  Subprogram arguments
            INTEGER,      INTENT (IN) :: NMAPVAR    ! no. data in map file
            INTEGER,      INTENT (IN) :: NHDRVAR    ! no. data in inven header
            INTEGER,      INTENT (IN) :: NPPOL      ! no. variables per data var
            INTEGER,      INTENT (IN) :: IDX2       ! index to which pol-assoc var
            CHARACTER(*), INTENT (IN) :: MAPVARS ( NMAPVAR )   ! map file var names
            CHARACTER(*), INTENT (IN) :: MAPFILES( NMAPVAR )   ! map file file names
            CHARACTER(*), INTENT(OUT) :: NAMES1  ( NHDRVAR )   ! var names list
            CHARACTER(*), INTENT(OUT) :: NAMES2  ( NHDRVAR )   ! var names read list
            CHARACTER(*), INTENT(OUT) :: UNITS   ( NHDRVAR )   ! var units

C.............  Local variables
            INTEGER  I, J

            CHARACTER(16)   :: TMPNAME = ' '   ! tmp logical file name

C----------------------------------------------------------------------

C............  Ensure that the number of map variables and the number of 
C              variables listed in the input file are consistent
            IF( NMAPVAR .NE. NHDRVAR ) THEN
                WRITE( MESG,94010 ) 'Number of map-formatted ' //
     &             'variables is inconsistent with header of  '//
     &             'inventory file.' // CRLF()// BLANK10 // 
     &             'I/O API header has', NHDRVAR,'variables indicated'//
     &              ', but map-formated inventory file has', NMAPVAR
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C............  Loop through variables requested
            DO I = 1, NMAPVAR

                TMPNAME = 'IOAPI_DAT'
                IF( .NOT. SETENVVAR( TMPNAME, MAPFILES( I ) ) ) THEN
                    MESG = 'ERROR: Could not set logical file name '
     &                          // CRLF() // BLANK10 // 'for file ' //
     &                          TRIM( MAPFILES( I ) )
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

                ELSE IF( .NOT. OPENSET( TMPNAME,FSREAD3,PROGNAME )) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not open "'//TRIM(MAPVARS(I))//
     &                      '" file listed in ' // CRLF() // BLANK10 // 
     &                      'map-formatted inventory for file name: ' //
     &                      CRLF() // BLANK10 // TRIM( MAPFILES( I ) )
                    CALL M3MSG2( MESG )
                    CYCLE

                ELSE IF( .NOT. DESCSET( TMPNAME,ALLFILES ) ) THEN
                    MESG = 'Could not get description of file "' //
     &                     TRIM( TMPNAME ) // '".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C................  Ensure that map variable name and file variable 
C                  name are the same
                IF( MAPVARS( I ) .NE. VNAMESET( 2 ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Variable name "'// TRIM(MAPVARS(I))//
     &                    '" in map-formatted inventory is mapped to '//
     &                     CRLF()// BLANK10// 'a file with variable "'//
     &                     TRIM( VNAMESET( 2 ) )//'". Map-formatted '//
     &                     'inventory has been corrupted.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

                NAMES1( I ) = VNAMESET( 2 )
                NAMES2( I ) = VNAMESET( IDX2 )
                UNITS ( I ) = VUNITSET( IDX2 )

C...............  Close output file for this variable
                IF( .NOT. CLOSESET( TMPNAME ) ) THEN
                    MESG = 'Could not close file:'//CRLF()//BLANK10//
     &                     TRIM( MAPFILES( I ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END DO

C...........   Internal buffering formats.............94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE STORE_INVEN_VARS

        END SUBROUTINE OPENMRGIN
