
        PROGRAM LNKMERGE 

C***********************************************************************
C  program LNKMERGE body starts at line
C
C  DESCRIPTION:
C      This program allocates link-level houlry/min/sec emissions into
C      the grid cells and compute houlry emissions. The grided houlry
C      emissions are merged with the speciation matrix.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     07/17: Created by B.H. Baek
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
        USE MODMERGE, ONLY: TSVDESC,
     &          MFLAG_BD, LREPANY, LGRDOUT, LREPSRC,
     &          LREPSTA, LREPCNY, LREPSCC,         ! report flags, gridded output
     &          CDEV, MONAME, NMSRC,               ! costcy
     &          MCFIP, EMNAM, EANAM, NIPPA,        ! inv counties and poll names
     &          NSMATV, NMSPC, SIINDEX, SPINDEX,   ! EANAM & EMNAM idx
     &          SDATE, STIME, NSTEPS, TSTEP, GRDFAC, ! episode information
     &          EDATE, ETIME, MSDATE, EMGGRD, CUTOFF,! dates for by-day hrly emis
     &          APRT_ELEV, APRT_CODE, NAPRT, APFLAG  ! airport height information

C.........  This module contains data structures and flags specific to lnkmerge
        USE MODMVSMRG, ONLY: METNAME, MGRNAME, GRDENV, TOTENV,
     &                       MSMATX_L, MSMATX_S, MSNAME_L, MSNAME_S,
     &                       MNSMATV_L, MNSMATV_S
        
C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, NINVSCC, INVSCC

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: MICNY, NCOUNTY, NSTATE

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID, NROWS, NCOLS, XORIG, YORIG, XOFF, YOFF,
     &                     GDTYP, XCELL, YCELL, XCENT, YCENT,
     &                     P_ALP, P_BET, P_GAM, GRDNM, COORD,
     &                     VGTYP, VGLVS, VGTOP, NLAYS

C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: CSCC, CIFIP, TZONES, CLINK, CDPTID, CARRID, CSOURC

C.........  MODULES for I/O API INTERFACEs, geo-transform codes:
        USE M3UTILIO
        USE MODGCTP

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        INCLUDE 'MVSCNST3.EXT'   !  MOVES constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         USEEXPGEO, INGRID, GETFLINE, CHKINT, BLKORCMT
        EXTERNAL        USEEXPGEO, INGRID, GETFLINE, CHKINT, BLKORCMT

C.........  LOCAL PARAMETERS and their descriptions:
        INTEGER, PARAMETER :: MXSEG = 50    ! number of segments in line
        CHARACTER(30)      :: SEGMENT( MXSEG ) = ''

        CHARACTER(50), PARAMETER :: CVSW = '$Name $' ! CVS release tag

C...........   Local arrays for per-source information
        INTEGER, ALLOCATABLE :: DAYBEGT( : )   ! daily start time for each source
        INTEGER, ALLOCATABLE :: DAYENDT( : )   ! daily end time for each source
        INTEGER, ALLOCATABLE :: POLIDX ( : )   ! index for processing pollutants
        LOGICAL, ALLOCATABLE :: LDAYSAV( : )   ! true: src uses DST
        REAL,    ALLOCATABLE :: LFRAC ( : )       ! model layer fractions
        REAL,    ALLOCATABLE :: TERRAIN( : )      ! terrain ENDHGT (m)
        REAL,    ALLOCATABLE :: SFCHGT( : )       ! surface ENDHGT (m)
        REAL,    ALLOCATABLE :: VGLVLS( : )       ! gridded mask values to be output
        REAL,    ALLOCATABLE :: VGLVSXG( : )      ! gridded mask values to be output
        REAL,    ALLOCATABLE :: ZZF   ( :,: )     ! layer's full ENDHGT (m)
        REAL,    ALLOCATABLE :: TMP3D ( :,:,:,: ) ! tmp emissions
        REAL,    ALLOCATABLE :: POLVAL( : )           ! array for poll values
        CHARACTER(IOVLEN3),ALLOCATABLE:: POLNAM( : )  ! array for poll names
        CHARACTER (256),ALLOCATABLE :: LINKLIST( : )  ! list of hourly link input files

C.........   Local arrays dimensioned by subroutine arguments
C.........   Note that the NGRID dimension could conceivably be too small if
C            a link winds through the whole domain, but this is a case that
C            is not worth going to extra trouble for since it is not realistic
        INTEGER,       ALLOCATABLE :: ACEL( : )    ! number of cell intersections per src
        REAL,          ALLOCATABLE :: AFAC( : )    ! fraction of link in cell

C.........  File units and logical names
        INTEGER      :: FDEV = 0            ! unit no for a list of of link inventory files
        INTEGER      :: FLDEV= 0            ! unit no for the list link inventory files
        INTEGER      :: SGDEV= 0            ! unit no for an individual link inventory file
        INTEGER      :: RDEV = 0            ! unit no for report file
        INTEGER      :: LDEV = 0            ! unit no for log file

        CHARACTER*16    GNAME               ! grid-point layered met file
        CHARACTER*16    PNAME               ! cross-point met file for surface pressure values
     
C...........   Other local variables
        INTEGER         C, I, II, J, K, L, N, NC, NL, NS, ES, P, S, SS, T, NV, V  ! counters and indices
        INTEGER         IOS                 ! i/o status
        INTEGER         IREC                ! line counter
        INTEGER         NCEL                ! tmp number of cells
        INTEGER         MXLAYS              ! max output layers
        INTEGER         ORG_CELLID          ! origin cell id
        INTEGER         END_CELLID          ! ending cell id
        INTEGER         NFILES              ! no of hourly link input files
        INTEGER         MXWARN              ! maximum number of warnings
        INTEGER      :: NWARN = 0           ! current number of warnings
        INTEGER         ROW, COL            ! grid cell row/col
        INTEGER         NLNK                ! location ID
        INTEGER         NVARS, NPOL         ! no. output variables, input pollutans
        INTEGER         SEGHOUR             ! segment processing hour
        INTEGER         SEGTIME             ! segment duration
        INTEGER         JDATE, JTIME, TDATE, TTIME, HTIME ! Processing date/time
        INTEGER         STR, DSTR, LST, DLST ! Processing start to end time
        INTEGER         YEAR, MON, DAY, HOUR, MIN, SEC

        INTEGER         LTOP, LBOT, DDP     ! layer# for top/bottom cells
        INTEGER         STRID, ENDID, INCID ! start/end grid cells and increment
        INTEGER         ZONE, DZONE         ! county-specific time zone

        REAL            Po, Ph, Z, Zo, Zh, ZBOT, ZTOP, PDIFF

        REAL            ZFRAC, PFRAC, LTOT

        REAL            TMPVAL
        REAL         :: LTOALT = 0.0        ! LTO operations altitude (ft)
        REAL         :: ALEN   = 0.0        ! link length
        REAL         :: STRLAT = 0.0        ! absolute latitude
        REAL         :: STRLON = 0.0        ! absolute longitude
        REAL         :: ENDLAT = 0.0        ! previous absolute latitude
        REAL         :: ENDLON = 0.0        ! previous absolute longitude
        REAL         :: STRHGT = 0.0        ! previous altitude
        REAL         :: ENDHGT = 0.0        ! altitude
        REAL         :: DELTAZ = 0.0        ! delta altitude
        REAL         :: RATIO  = 1.0        ! model_length/total_length

        LOGICAL      :: FIRSTIME  = .TRUE.  ! true: first time
        LOGICAL      :: EFLAG = .FALSE.     ! true: ERROR

        CHARACTER(3)         TZN                 ! tmp time zone
        CHARACTER(80)        GDESC               ! grid description
        CHARACTER(1056)      LINE                ! input line buffer
        CHARACTER(256)       MESG                ! message buffer
        CHARACTER(SRCLEN3)   CSRC                ! tmp source chars string
        CHARACTER(SCCLEN3)   TSCC                ! tmp scc code
        CHARACTER(IOVLEN3)   VBUF                ! tmp variable name buffer
        CHARACTER(FIPLEN3)   CFIP, LFIP          ! tmp FIPS code
        CHARACTER(IOVLEN3)   TMPOUT              ! temporal resolution
        CHARACTER(LNKLEN3)   DPRTID              ! tmp departure airport code
        CHARACTER(LNKLEN3)   ARRVID              ! tmp arriving airport code
        CHARACTER(LNKLEN3)   LNKID               ! aircraft flight identifier
        CHARACTER(LNKLEN3)   SEGID               ! aircraft flight segment identifie
        CHARACTER(IODLEN3)   TDESC               ! tmp combined pollutant & species

        CHARACTER(16) :: PROGNAME = 'LNKMERGE' ! program name

C***********************************************************************
C   begin body of program LNKMERGE 
        
        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Retrieve control environment variables and set logical control
C           flags. Use a local module to pass the control flags.
        CALL GETMRGEV

C.........  Open input files, read inventory data, and retrieve episode information
        CALL OPENMRGIN

C.........  Create array of sorted unique pol-to-species
C.........  Create array of sorted unique pollutants
C.........  Create array of sorted unique species
        CALL MRGVNAMS

C.........  Determine units conversion factors
        CALL MRGUNITS

C.........  Read the state and county names file and store for the 
C           states and counties in the grid
C.........  Use FIPS list from the inventory for limiting state/county list
        IF( USEEXPGEO() ) THEN
            CALL RDGEOCODES( NINVIFIP, INVCFIP )
        ELSE
            CALL RDSTCY( CDEV, NINVIFIP, INVCFIP )
        END IF

C.........  Allocate memory for fixed-size arrays...        
        ALLOCATE( LDAYSAV( NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LDAYSAV', PROGNAME )
        LDAYSAV = .FALSE.  ! array
        
        ALLOCATE( DAYBEGT( NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYBEGT', PROGNAME )
        DAYBEGT = 0   ! array
        
        ALLOCATE( DAYENDT( NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAYENDT', PROGNAME )
        DAYENDT = 0   ! array

C.........  Determine sources that observe DST
        CALL GETDYSAV( NMSRC, CIFIP, LDAYSAV )

C.........  Store local layer info
        ALLOCATE( VGLVLS( 0:MXLAYS3 ), STAT= IOS)
        CALL CHECKMEM( IOS, 'VGLVLS', PROGNAME )
        ALLOCATE( VGLVSXG( 0:MXLAYS3 ), STAT= IOS)
        CALL CHECKMEM( IOS, 'VGLVSXG', PROGNAME )

        NLAYS  = NLAYS3D
        VGTYP  = VGTYP3D
        VGTOP  = VGTOP3D
        VGLVLS = 1.0 - VGLVS3D   ! array

C.........  Store local layer information
        J = LBOUND( VGLVS3D, 1 )
        VGLVSXG( 0 ) = VGLVS3D( J )
        DO I = 1, NLAYS
            J = J + 1
            VGLVSXG( I ) = VGLVS3D( J )
        END DO

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.........  Build indicies for pollutant/species groups
        CALL BLDMRGIDX

C.........  Intialize state/county summed emissions to zero
        CALL INITSTCY

C.........  Open NetCDF output files, open ASCII report files, and write headers
        CALL OPENMRGOUT

C.........  Buid speciation matrices
        ALLOCATE( POLIDX( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLIDX', PROGNAME )
        ALLOCATE( POLNAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLNAM', PROGNAME )
        ALLOCATE( POLVAL( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )
        POLIDX = 0
        POLNAM = ''
        POLVAL = 0.0

        NVARS = MNSMATV_L
        ALLOCATE( MSMATX_L( NMSRC, MNSMATV_L ), STAT=IOS )    ! mole speciation matrix
        CALL CHECKMEM( IOS, 'MSMATX_L', PROGNAME )

        ALLOCATE( MSMATX_S( NMSRC, MNSMATV_S ), STAT=IOS )    ! mass speciation matrix
        CALL CHECKMEM( IOS, 'MSMATX_S', PROGNAME )

C.........  Read speciation matrices for current variable
         DO V = 1, NSMATV

C.............  Extract name of variable
            TDESC = TSVDESC( V )
            CALL RDSMAT( MSNAME_L, TDESC, MSMATX_L( 1,V ) )
            CALL RDSMAT( MSNAME_S, TDESC, MSMATX_S( 1,V ) )

C.................  Switch SPC matrix (mole/mass) based on MRG_GRDOUT_UNIT, MRG_TOTOUT_UNIT
            IF( INDEX( GRDENV, 'mole' ) < 1 ) THEN
                MSMATX_L( 1,V ) = MSMATX_S( 1,V )
            END IF

        END DO

C.........  Allocate arrays
        ALLOCATE( TERRAIN( NGRID ), STAT= IOS)
        CALL CHECKMEM( IOS, 'TERRAIN', PROGNAME )
        ALLOCATE( ZZF( NGRID,NLAYS ), STAT= IOS)
        CALL CHECKMEM( IOS, 'ZZF', PROGNAME )
        ALLOCATE( TMP3D( NGRID,NLAYS,NVARS,NSTEPS ), STAT= IOS)
        CALL CHECKMEM( IOS, 'TMP3D', PROGNAME )
        ALLOCATE( LFRAC( NLAYS ), STAT= IOS)
        CALL CHECKMEM( IOS, 'LFRAC', PROGNAME )
        ALLOCATE( ACEL( NGRID ), STAT= IOS)
        CALL CHECKMEM( IOS, 'ACEL', PROGNAME )
        ALLOCATE( AFAC( NGRID ), STAT= IOS)
        CALL CHECKMEM( IOS, 'AFAC', PROGNAME )
        ACEL   = 0
        AFAC   = 0.0
        TERRAIN= 0.0
        ZZF    = 0.0
        LFRAC  = 0.0
        TMP3D  = 0.0

C.........  Open actual Segment input file
        MESG = 'Enter logical name for hourly link inventory list file'
        FLDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'LINKHOUR', PROGNAME )

        NFILES = GETFLINE( FLDEV, 'Opening houlry LINK inventory list file' )

C.........  Allocate array for storing the list of link inventory files
        ALLOCATE( LINKLIST( NFILES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LINKLIST', PROGNAME )
        LINKLIST = ' '

C.........  Store lines of LINKHOUR list file
        MESG = 'Opening hourly link inventory file'
        CALL RDLINES( FLDEV, MESG, NFILES, LINKLIST )

        DO II = 1, NFILES

          IF( BLKORCMT( LINKLIST(II) ) ) CYCLE

          OPEN( SGDEV, FILE=LINKLIST( II ), STATUS='OLD', IOSTAT=IOS )

C...........  Check for errors while opening file
          IF( IOS /= 0 ) THEN
              MESG = 'ERROR: Could not open file:' //
     &             CRLF()// BLANK5 // TRIM( LINKLIST( (II) ) )
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          ELSE
              MESG = 'Successful OPEN for inventory file(s):' //
     &             CRLF()// BLANK5 // TRIM( LINKLIST( (II) ) )
              CALL M3MSG2( MESG )
          END IF

C.........  Reading individual hourly link inventory files         
          IREC  = 0
          LFIP = ''
          DO
            
            READ( SGDEV, 93000, END=299 ) LINE

            IREC = IREC + 1

C.............  Check the header line first
            IF( IREC == 1 ) THEN
                I = INDEX( LINE, 'FF10_LINK_HOURLY' )
                IF( I < 1 ) THEN
                    MESG = 'ERROR: Missing header "FF10_LINK_HOURLY": '//
     &                  'Define the header of houlry link inventory file' 
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )
                END IF
            END IF

C.............   Store the list of inventory pollutatns
            IF( LINE( 1:6 ) == '#POLID' ) THEN
                L = LEN_TRIM( LINE )
                CALL PARSLINE( LINE(7:L), MXSEG, SEGMENT )

                NPOL = 0
                DO I = 1, NIPPA
                    IF( LEN_TRIM( SEGMENT( I ) ) < 1 ) CYCLE
                    IF( INDEX1( SEGMENT( I ), NIPPA, EANAM ) < 1 ) CYCLE
                    POLIDX( I ) = I
                    NPOL = NPOL + 1
                    POLNAM( NPOL ) = SEGMENT( I )
                END DO
            END IF

            IF( LINE( 1:1 ) == '#' ) CYCLE

            IF( MOD( IREC,100000 ) == 0 ) THEN
                WRITE( MESG, 94010 ) 'Processing line at', IREC
                CALL M3MSG2( MESG )
            END IF

C..............  Read data from LINE
            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C..............  Skip lines
            IF( .NOT. CHKINT( SEGMENT( 7 ) ) ) CYCLE

C.............  Read FIPS code
            IF( USEEXPGEO() ) THEN
                CFIP(  1: 3 ) = ADJUSTR( SEGMENT( 1 )( 1:3 ) )
                CFIP(  4: 9 ) = ADJUSTR( SEGMENT( 2 )( 1:6 ) )
                CFIP( 10:12 ) = ADJUSTR( SEGMENT( 3 )( 1:3 ) )
            ELSE
                CFIP( FIPEXPLEN3+2:FIPLEN3 ) = ADJUSTR( SEGMENT( 2 )( 1:5 ) )  ! country/state/county code
            END IF

C.............  Replace blanks with zeros
            DO I = 1,FIPLEN3
                IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
            END DO

            LNKID = TRIM( SEGMENT( 4 ) )       ! Flight ID
            TSCC  = TRIM( SEGMENT( 5 ) )       ! SCC code
            CALL PADZERO( TSCC )

C.............  Build source characteristics field for searching inventory
            CALL BLDCSRC( CFIP, LNKID, TSCC, CHRBLNK3,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                    CHRBLNK3, CSRC )

C.............  Look up the source ID for this record
            S = FINDC( CSRC, NMSRC, CSOURC )
            IF( S < 1 ) THEN
                WRITE( MESG,94010 ) 'WARNING: Can not find a matched source: '
     &               // 'Skipping line at ',IREC
                NWARN = NWARN + 1
                IF( NWARN < MXWARN ) CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Compute processing dates from LINK_HOURLY file 
            YEAR  = STR2INT( SEGMENT( 7 ) )    ! integer current year
            MON   = STR2INT( SEGMENT( 8 ) )    ! integer current month
            DAY   = STR2INT( SEGMENT( 9 ) )    ! integer current day
            HOUR  = STR2INT( SEGMENT( 10 ) )   ! integer current hour
            MIN   = STR2INT( SEGMENT( 11 ) )   ! integer current mins

            JDATE = YEAR * 1000 + JULIAN( YEAR, MON, DAY )  ! segment julian date
            JTIME = HOUR * 10000 + MIN * 100

C.............  local to output time zone shift 
            TZN = TRIM( SEGMENT( 12 ) )
            K = INDEX1( TZN, MXTZONE, TZONNAM )
            IF( K < 1 ) THEN
                WRITE( MESG,94010 ) 
     &              'ERROR: Time zone ('//TRIM(SEGMENT(12))//') at line',IREC,
     &              ' is not supported'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            ZONE = TZONNUM( K ) 
            CALL NEXTIME( JDATE, JTIME, ZONE * 10000 )

C.............  Check segment time
            SEGTIME =  STR2INT( SEGMENT( 13 ) )  ! Duration (seconds)
            IF( SEGTIME == 0 ) THEN
                WRITE( MESG,94010 )
     &              'WARNING: Can NOT process zero segment duration: Skipping line at ',IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
            HTIME = SEC2TIME( SEGTIME )

C............. compute ending date and time
            TDATE = JDATE
            TTIME = JTIME
            CALL NEXTIME( TDATE, TTIME, HTIME )

C.............  Skipping hours
            IF( SECSDIFF( SDATE, STIME, TDATE, TTIME ) < 0 ) CYCLE
            IF( SECSDIFF( EDATE, ETIME, JDATE, JTIME ) > 0 ) CYCLE

C.............  Calculate processing start hour (T)
C.............  Determine time step pointer based on reference time
            STR = 1 + SECSDIFF( SDATE, STIME, JDATE, JTIME ) / 3600
            LST = 1 + SECSDIFF( SDATE, STIME, TDATE, TTIME ) / 3600
            SEGHOUR = 1 + SECSDIFF( JDATE, JTIME, TDATE, TTIME ) / 3600

C.............  Adjust starting/ending time modeling period
            DSTR = 0
            DLST = 0
            IF( STR < 1 ) THEN
                DSTR = 1 - STR
                STR  = 1
            END IF
            IF( LST > NSTEPS ) THEN
                DLST = LST - NSTEPS
                LST  = NSTEPS
            END IF

C.............  Store poll name and values and compute emission rate based on segment duration (se
            POLVAL = 0.0

C..................  Compute the ratio of segment based on segment duration and modeling duration
            P = 0
            DO I = 1, NIPPA
                IF( POLIDX( I ) < 1 ) CYCLE
                P = P + 1
                RATIO = 1.0 - ( FLOAT( DSTR + DLST ) / FLOAT( SEGHOUR ) )
                POLVAL( P ) = RATIO * STR2REAL( SEGMENT(21 + I) )
            END DO

C.............  Skip when it is out of range of episode dates or non-matched link IDs
            LNKID = TRIM( SEGMENT( 4 ) )       ! Flight ID
            NLNK = INDEX1( LNKID, NMSRC, CLINK ) 
C            IF( NLNK < 1 ) CYCLE

C.............  Convert source coordinates from lat-lon to output grid
            STRLON = STR2REAL( SEGMENT( 16 ) )
            STRLAT = STR2REAL( SEGMENT( 17 ) )
            STRHGT = FT2M * STR2REAL( SEGMENT( 18 ) )   ! altitude in unit of feet

            ENDLON = STR2REAL( SEGMENT( 19 ) )
            ENDLAT = STR2REAL( SEGMENT( 20 ) )
            ENDHGT = FT2M * STR2REAL( SEGMENT( 21 ) )   ! altitude in unit of feet

            Zo = STRHGT                ! origin height
            Zh = ENDHGT                ! end height

            CALL CONVRTXY( 1, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,
     &                     XCENT, YCENT, STRLON, STRLAT )

            CALL CONVRTXY( 1, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,
     &                     XCENT, YCENT, ENDLON, ENDLAT )

C.............  If link source, determine the number of cells for this source
            NCEL = 0
            ACEL = 0
            AFAC = 0.0
            CALL LNK2GRD( NGRID, STRLON, STRLAT, ENDLON, ENDLAT,
     &                    NCEL, ACEL, AFAC, ALEN, EFLAG)

C.............  Make sure that there was enough storage
            IF( EFLAG ) THEN
                WRITE( MESG,94010 )
     &              'INTERNAL ERROR: Overflow for source at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip if there is no intersected grid cell
            IF( NCEL == 0 ) CYCLE

c bbaek        write(*,'(10I10,2f12.5)') s,irec,ncel,jdate,jtime,tdate,ttime,segtime,str,lst,zo,zh
C.............  If source is in the domain, get cell number and store
            ORG_CELLID = 0
            END_CELLID = 0
            IF( INGRID( ENDLON,ENDLAT,NCOLS,NROWS,COL,ROW ) ) THEN
                ORG_CELLID = COL + ( ROW - 1 ) * NCOLS
            END IF

            IF( INGRID( STRLON,STRLAT,NCOLS,NROWS,COL,ROW ) ) THEN
                END_CELLID = COL + ( ROW - 1 ) * NCOLS
            END IF

C.............  Read layers top ENDHGT (meter)
            JTIME = HOUR * 10000
            IF( .NOT. READ3( METNAME, 'ZF', -1,
     &                       JDATE, JTIME, ZZF ) ) THEN
                MESG = 'ERROR : Could not read ZF from file '
     &                 // TRIM( METNAME )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Adjust altitudes for climbing/landing flight trajectory
            IF( APFLAG ) THEN
                IF( .NOT. READ3( MGRNAME, 'HT', -1,
     &                       JDATE, JTIME, TERRAIN ) ) THEN
                    MESG = 'ERROR : Could not read TERRAIN from file '
     &                     // TRIM( MGRNAME )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                SEGID = TRIM( SEGMENT( 6 ) )
                CALL UPCASE( SEGID )
                IF( SEGID < 'TAXI_OUT' ) THEN
                    N = INDEX1( CARRID( S ), NAPRT, APRT_CODE )
                    Zh = Zh - APRT_ELEV( N )   ! departure airport elev
                ELSE IF( SEGID > 'TAXI_IN' ) THEN
                    N = INDEX1( CDPTID( S ), NAPRT, APRT_CODE )
                    Zh = Zh - APRT_ELEV( N )   ! arrival airport elev
                ELSE
                    IF( END_CELLID>0 ) Zh = Zh - TERRAIN(END_CELLID)
                ENDIF
            END IF

            IF( Zo < 0.0 ) Zo = 0.0
            IF( Zh < 0.0 ) Zh = 0.0

C............  Apply CUTOFF method (i.e., 10K ft)
C              if previous ENDHGT > CUTOFF, skip processing
            IF( Zh >= CUTOFF .AND. Zo >= CUTOFF ) CYCLE

C.............  Sort the order of grid cell processing. Starting from org_cellid....
            ZBOT = Zo
            DELTAZ = Zh - Zo   ! delta z (<0:langind, >0:climbing)

            IF( ACEL( 1 ) /= ORG_CELLID ) THEN
                STRID = NCEL
                ENDID = 1
                INCID = -1
            ELSE
                STRID = 1
                ENDID = NCEL
                INCID = 1
            END IF

            FIRSTIME = .TRUE.

C.............  loop over assigned grid cells
            DO NC = STRID, ENDID, INCID

                LFRAC = 0.0         ! Gridded x-y link vertical fractions 

                C     = ACEL( NC )  ! cell id
                ZFRAC = AFAC( NC )  ! cell-fraction value

                IF( C > NGRID ) THEN
                    MESG='ERROR: Incorrect link to grid conversion'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
                 
                IF( ZFRAC < 0.0 ) THEN
                    MESG = 'ERROR: Can not process negative x-y '//
     &                     'link fraction'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

C.................  retrieve col/row from cellid using C=(ROW-1)*NCOLS+COL
                ROW = C / NCOLS
                IF( MOD( C, NCOLS ) .GT. 0. ) ROW = ROW + 1
                COL = C - ( ROW-1 ) * NCOLS

C.................  Update bottom and top layers
                Z = ZFRAC * DELTAZ
                IF( DELTAZ < 0.0 ) THEN   ! aircraft landing mode
                    ZTOP = ZBOT 
                    ZBOT = ZTOP + Z
                ELSE                    ! aircraft climbing mode
                    ZBOT = ZBOT
                    ZTOP = ZBOT + Z
                END IF

c          write(*,'(a,6f10.3)')lnkid,zo,zbot,zh,ztop
C.................  Looping through layers to determine associated layer for each link
                DO L = 1, NLAYS - 1

                    IF ( ZBOT <= ZZF( C,L ) ) THEN
                        LBOT = L
                        GO TO  111   ! end loop and skip reset of LBOT
                    END IF

                END DO

                LBOT = NLAYS           !  fallback

C.................  hard corded to switch vertical allocation method between sigma and pressure
111             CONTINUE                !  loop exit:  bottom found at LBOT
 
                IF ( ZTOP <= ZZF( C,LBOT ) ) THEN  !  plume in this layer
 
                    PFRAC = 1.0
                    LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                    LTOP = LBOT

c        write(*,'(a,2F10.3,5i8,f15.7,a)')lnkid,zbot,ztop,col,row,LBOT,LTOP,LBOT,pfrac,' onelayer'

                ELSE IF( LBOT == NLAYS ) THEN    ! plume above top layer
 
                    PFRAC = 1.0
                    LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                    LTOP = NLAYS
c        write(*,'(a,2F10.3,5i8,f15.7,a)')lnkid,zbot,ztop,col,row,LBOT,LTOP,LTOP,pfrac,' toplayer'
                
                ELSE                               ! plume crosses layers
 
                    DO L = LBOT + 1, NLAYS
                        IF ( ZTOP <= ZZF( C,L ) ) THEN
                            LTOP = L
                            GO TO 222  ! end loop and skip reset of LTOP
                        END IF
                    END DO
                    LTOP = NLAYS
 
222                 CONTINUE
 
C.....................  Calculate between layer 
                    PDIFF = ZTOP - ZBOT
                    
C.....................  Calculate a fraction for the bottom layer
                    PFRAC = ( ( ZZF( C,LBOT ) - ZBOT )
     &                              / PDIFF )
                    LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
c            write(*,'(a,2F10.3,5i8,f15.7)')lnkid,zbot,ztop,col,row,LBOT,LTOP,LBOT,pfrac

C.....................  Calculate a fraction for the top layer
                    PFRAC = ( (ZTOP-ZZF( C,LTOP-1 )) / PDIFF )
                    LFRAC( LTOP ) = LFRAC( LTOP ) + PFRAC

                    DDP = LTOP - LBOT
                    IF( DDP >= 2 ) THEN
                        DO L = LBOT+1, LTOP-1 !  layers in plume
                            
                            PFRAC=( (ZZF(C,L)-ZZF(C,L-1)) / PDIFF )
                            LFRAC( L ) = LFRAC( L ) + PFRAC
c            write(*,'(a,2F10.3,5i8,f15.7)')lnkid,zbot,ztop,col,row,LBOT,LTOP,L,pfrac
                        END DO
                    ENDIF
c            write(*,'(a,2F10.3,5i8,f15.7)')lnkid,zbot,ztop,col,row,LBOT,LTOP,LTOP,lfrac(ltop)
                END IF

C.................  initialize pressure-based vertical allocation by link
                IF( DELTAZ >= 0.0 ) ZBOT = ZTOP   ! climbing mode: 

C................. Before applying layer fractions make sure that they add to 1.0
                LTOT = 0.0
                DO NL = 1, NLAYS
                    LTOT = LTOT + LFRAC( NL )
                ENDDO

                IF( LTOT < 0.99999 ) THEN
                    MESG = 'ERROR: Total of layer fractions '//
     &                     'are less than 1.0.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE IF( LTOT > 1.00001 ) THEN
                    MESG = 'ERROR: Total of layer fractions '//
     &                     'are greater than 1.0.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

C.................  Loop over allocated layers by link
                DO L = LBOT, LTOP

C.....................  Apply speciation matrix and layer fractions to compute
C                       gridded/speciated/hourly model specie emissions rate in output unit
                    DO T = STR, LST
                    DO P = 1, NPOL
                    DO V = 1, NVARS
                        IF( POLNAM( P ) == EANAM( SIINDEX( V,1 ) ) ) THEN
                            TMPVAL = POLVAL( P ) * MSMATX_L( S, SPINDEX(V,1) ) * GRDFAC( SPINDEX(V,1) ) 
                            TMP3D( C,L,V,T ) = TMP3D( C,L,V,T ) + TMPVAL * ZFRAC * LFRAC( L )
                        END IF
                    END DO
                    END DO
                    END DO

                    MXLAYS = MAX( LTOP, MXLAYS )     ! define max layer #

                END DO      ! end of link loop

            ENDDO       ! end of flight loop

          END DO   ! end of loop

299       CONTINUE ! Exit from read loop

        END DO

C......... Error message
        IF( EFLAG ) THEN 
            WRITE( MESG,94010 ) 'I/O error', IOS,
     &          'reading input file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Define top layer for output file
        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS
            DO V = 1, NVARS
                IF ( .NOT. WRITE3( MONAME, EMNAM(V), JDATE, JTIME, 
     &                            TMP3D( :,1:NLAYS,V,T ) ) ) THEN
                    WRITE( MESG, 93000 ) 'Could not write to "'
     &                    // TRIM( MONAME ) // '".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
            END DO
            CALL NEXTIME( JDATE, JTIME, 10000 )
        END DO

        IF ( .NOT. CLOSE3( MONAME ) ) THEN
            CALL M3ERR( PROGNAME, 0, 0, TRIM( MONAME ) // '".', .TRUE. )
        END IF      !  if close3() failed 

C.........  Successful completion of program
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C.........  Formatted file I/O formats...... 93xxx
93000   FORMAT( A )

C.......  Internal buffering formats...... 94xxx
94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

        END PROGRAM LNKMERGE 
