
        PROGRAM MOVESMRG 

C***********************************************************************
C  program MOVESMRG body starts at line
C
C  DESCRIPTION:
C      This program reads emission factors from MOVES and calculates
C      hourly emissions based on VMT or vehicle population data. The
C      hourly emissions are merged with the gridding matrix and
C      speciation matrix.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/10 by C. Seppanen - based on smkmerge.f
C     04/11: Modified by B.H. Baek
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
        USE MODMERGE, ONLY: 
     &          MFLAG_BD,                          ! by-day hourly emis flags
     &          LREPSTA, LREPCNY, LREPSCC, LREPSRC, LGRDOUT,! report flags, gridded output
     &          CDEV,                              ! costcy
     &          MGNAME, MTNAME, MONAME,    ! input files
     &          NMSRC, MNGMAT,                     ! no. of srcs, no. gridding matrix entries
     &          NMSPC,                             ! no. species
     &          EMNAM,                             ! species names
     &          TSVDESC,                           ! var names
     &          SIINDEX, SPINDEX,                  ! EANAM & EMNAM idx
     &          SDATE, STIME, NSTEPS, TSTEP,       ! episode information
     &          MSDATE,                            ! dates for by-day hrly emis
     &          GRDFAC, TOTFAC,                    ! conversion factors
     &          NSMATV,                            ! speciation matrices
     &          MEBCNY, MEBSTA, MEBSUM,            ! cnty/state/src total spec emissions
     &          MEBSCC, MEBSTC,                    ! scc total spec emissions
     &          EANAM, NIPPA                       ! pol/act names

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: RPDFLAG, RPVFLAG, RPPFLAG, 
     &          TVARNAME, METNAME,
     &          NREFSRCS, REFSRCS, NSRCCELLS, SRCCELLS, SRCCELLFRACS,
     &          EMPROCIDX, EMPOLIDX,
     &          NEMTEMPS, EMTEMPS, EMXTEMPS, EMTEMPIDX, AVGMIN, AVGMAX,
     &          RPDEMFACS, RPVEMFACS, RPPEMFACS,
     &          SPDFLAG, SPDPRO, MISCC, 
     &          MSNAME_L, MSMATX_L, MNSMATV_L, 
     &          MSNAME_S, MSMATX_S, MNSMATV_S,
     &          EANAMREP

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVIFIP, NINVSCC, INVSCC

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: MICNY, NCOUNTY, NSTATE

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID

C.........  This module is used for reference county information
        USE MODMBSET, ONLY: NREFC, MCREFIDX,
     &                      NREFF, FMREFSORT, NFUELC, FMREFLIST

C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: SPEED, CSCC, VPOP, IFIP, TZONES

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        INCLUDE 'MVSCNST3.EXT'   !  MOVES constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(10)   HHMMSS
        INTEGER         INDEX1
        INTEGER         FIND1
        INTEGER         FIND1FIRST
        INTEGER         FINDC
        INTEGER         WKDAY

        EXTERNAL    HHMMSS, INDEX1, FIND1, FIND1FIRST, FINDC, WKDAY

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Local arrays for per-source information
        INTEGER, ALLOCATABLE :: DAYBEGT( : )   ! daily start time for each source
        INTEGER, ALLOCATABLE :: DAYENDT( : )   ! daily end time for each source
        LOGICAL, ALLOCATABLE :: LDAYSAV( : )   ! true: src uses DST

C...........   Local arrays for hourly data
        REAL, ALLOCATABLE :: VMT( : )
        REAL, ALLOCATABLE :: TEMPG( : )
        REAL, ALLOCATABLE :: EMGRD( :,: )     ! emissions for each grid cell and species
        REAL, ALLOCATABLE :: TMPEMGRD( :,: )  ! tmp emissions for each grid cell and species

C...........   Local temporary array for input and output variable names
        CHARACTER(IOVLEN3), ALLOCATABLE :: VARNAMES( : )

C...........   Local array for gridding matrix
        REAL, ALLOCATABLE :: MGMATX( : )

C...........   Logical names and unit numbers (not in MODMERGE)
        INTEGER         LDEV
     
C...........   Other local variables
    
        INTEGER          I, J, K, L1, L2, M, N, NG, V, S, T ! counters and indices

        INTEGER          BIN1, BIN2    ! speed bins for current source
        INTEGER          CELL          ! current grid cell
        INTEGER          DAY           ! day-of-week index (monday=1)
        INTEGER          DAYMONTH      ! day-of-month
        INTEGER          DAYIDX        ! current day value index
        INTEGER          FUELMONTH     ! current fuel month
        INTEGER          HOURIDX       ! current hour of the day
        INTEGER          IDX1, IDX2    ! temperature indexes for current cell
        INTEGER          IOS           ! tmp I/O status
        INTEGER          JDATE         ! Julian date (YYYYDDD)
        INTEGER          JTIME         ! time (HHMMSS)
        INTEGER       :: K1 = 0        ! tmp index
        INTEGER       :: K5 = 0        ! tmp index
        INTEGER          KM            ! tmp index to src-category species
        INTEGER          LDATE         ! Julian date from previous iteration
        INTEGER          MJDATE        ! mobile-source Julian date for by-day
        INTEGER          MONTH         ! current month
        INTEGER          OCNT          ! tmp count output variable names
        INTEGER       :: PDAY = 0      ! previous iteration day no.
        INTEGER          POLIDX        ! current pollutant index
        INTEGER          PROCIDX       ! current emission process index
        INTEGER          SCCIDX        ! current SCC index
        INTEGER          SRC           ! current source number
        INTEGER          UUIDX, UOIDX, OUIDX, OOIDX  ! indexes for matching profiles
        INTEGER          USTART, UEND, OSTART, OEND

        REAL             F1, F2, FG0   ! tmp conversion
        REAL             GFRAC         ! grid cell fraction
        REAL             SPEEDVAL      ! average speed value for current source
        REAL             TEMPVAL       ! temperature value for current grid cell
        REAL             VMTVAL        ! hourly VMT value for current source and hour
        REAL             VPOPVAL       ! annual vehicle population value for current source
        REAL             SPDFAC        ! speed interpolation factor
        REAL             TEMPFAC       ! temperature interpolation factor
        REAL             MINVAL, MAXVAL  ! min and max temperature for current source
        REAL             UMIN, OMIN      ! bounding minimum temperature profile values
        REAL             UMAX, OMAX      ! bounding maximum temperature profile values
        REAL             MINFAC, MAXFAC  ! min and max temp interpolation factors
        REAL             PDIFF, TDIFF  ! temperature differences
        REAL             EFVAL1, EFVAL2, EFVALA, EFVALB, EFVAL   ! emission factor values
        REAL             EMVAL         ! emissions value

        LOGICAL       :: NO_INTRPLT = .FALSE.   ! true: single interploation, false: bi-interpolation
        LOGICAL       :: LAST_CNTY  = .FALSE.   ! true: reach the last county
        LOGICAL       :: SKIPSUM    = .FALSE.   ! true: skip MEBSUM calculation

        CHARACTER(300)     MESG    ! message buffer
        CHARACTER(IOVLEN3) LBUF    ! previous species or pollutant name
        CHARACTER(IOVLEN3) PBUF    ! tmp pollutant or emission type name
        CHARACTER(IOVLEN3) SBUF    ! tmp species or pollutant name
        CHARACTER(PLSLEN3) VBUF    ! pol to species or pol description buffer
        CHARACTER(SCCLEN3) SCC     ! current source SCC

        CHARACTER(16) :: PROGNAME = 'MOVESMRG' ! program name

C***********************************************************************
C   begin body of program MOVESMRG 
        
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
        CALL RDSTCY( CDEV, NINVIFIP, INVIFIP )

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

        ALLOCATE( MGMATX( NGRID + 2 * MNGMAT ), STAT=IOS )    ! contiguous gridding matrix
        CALL CHECKMEM( IOS, 'MGMATX', PROGNAME )

        ALLOCATE( EMGRD( NGRID, NMSPC ), STAT=IOS )     ! gridded emissions
        CALL CHECKMEM( IOS, 'EMGRD', PROGNAME )
        ALLOCATE( TMPEMGRD( NGRID, NMSPC ), STAT=IOS )     ! gridded emissions
        CALL CHECKMEM( IOS, 'TMPEMGRD', PROGNAME )
        
        IF( LREPSCC ) THEN
            ALLOCATE( MEBSCC( NINVSCC, NMSPC+NIPPA ), STAT=IOS )    ! SCC totals
            CALL CHECKMEM( IOS, 'MEBSCC', PROGNAME )
        END IF
        
        IF( LREPSTA ) THEN
            ALLOCATE( MEBSTA( NSTATE, NMSPC+NIPPA ), STAT=IOS )    ! state totals
            CALL CHECKMEM( IOS, 'MEBSTA', PROGNAME )
            
            IF( LREPSCC ) THEN
                ALLOCATE( MEBSTC( NSTATE, NINVSCC, NMSPC+NIPPA ), STAT=IOS )     ! state-scc totals
                CALL CHECKMEM( IOS, 'MEBSTC', PROGNAME )
            END IF
        END IF

        IF( LREPCNY ) THEN
            ALLOCATE( MEBCNY( NCOUNTY, NMSPC+NIPPA ), STAT=IOS )    ! county totals
            CALL CHECKMEM( IOS, 'MEBCNY', PROGNAME )
        END IF
        
        ALLOCATE( MEBSUM( NMSRC, NMSPC+NIPPA ), STAT=IOS )    ! source totals
        CALL CHECKMEM( IOS, 'MEBSUM', PROGNAME )

        ALLOCATE( MSMATX_L( NMSRC, MNSMATV_L ), STAT=IOS )    ! mole speciation matrix
        CALL CHECKMEM( IOS, 'MSMATX_L', PROGNAME )

        ALLOCATE( MSMATX_S( NMSRC, MNSMATV_S ), STAT=IOS )    ! mass speciation matrix
        CALL CHECKMEM( IOS, 'MSMATX_S', PROGNAME )
        
        IF( RPDFLAG ) THEN
            ALLOCATE( VMT( NMSRC ), STAT=IOS )     ! hourly VMT
            CALL CHECKMEM( IOS, 'VMT', PROGNAME )
        END IF
        
        ALLOCATE( TEMPG( NGRID ), STAT=IOS )    ! hourly temperatures
        CALL CHECKMEM( IOS, 'TEMPG', PROGNAME )

C.........  Determine sources that observe DST
        CALL GETDYSAV( NMSRC, IFIP, LDAYSAV )

C.........  Read gridding matrix
        CALL RDGMAT( MGNAME, NGRID, MNGMAT, MNGMAT,
     &               MGMATX( 1 ), MGMATX( NGRID + 1 ),
     &               MGMATX( NGRID + MNGMAT + 1 ) )

C.........  Build list of grid cells for each source
        CALL BLDSRCCELL( NMSRC, NGRID, MNGMAT, MGMATX( 1 ), 
     &                   MGMATX( NGRID + 1 ), MGMATX( NGRID + MNGMAT + 1 ) )

C.........  Set up reference county information
        CALL SETREFCNTY

C.........  Build emission process mapping
        CALL BLDPROCIDX

C.........  Build indicies for pollutant/species groups
        CALL BLDMRGIDX

C.........  Open NetCDF output files, open ASCII report files, and write headers
        CALL OPENMRGOUT

C.........  Intialize state/county summed emissions to zero
        CALL INITSTCY

C.........  Allocate memory for temporary list of species and pollutant names
        ALLOCATE( VARNAMES( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VARNAMES', PROGNAME )

C.........  Loop through pollutant-species combos and read speciation matrix
        OCNT = 0
        LBUF = ' '
        VARNAMES = ' '  ! array
        DO V = 1, NSMATV

C.............  Extract name of variable
            VBUF = TSVDESC( V )

C.............  Update list of output species names for message
            SBUF = EMNAM( SPINDEX( V,1 ) )
            M = INDEX1( SBUF, OCNT, VARNAMES )

            IF( M .LE. 0 .AND. SBUF .NE. LBUF ) THEN
                OCNT = OCNT + 1                            
                VARNAMES( OCNT ) = SBUF
                LBUF = SBUF
            END IF

C.............  Read speciation matrices for current variable
            CALL RDSMAT( MSNAME_L, VBUF, MSMATX_L( 1,V ) )
            CALL RDSMAT( MSNAME_S, VBUF, MSMATX_S( 1,V ) )

        END DO

C.........  Write out message with list of species
        CALL POLMESG( OCNT, VARNAMES )

C.........  Loop over reference counties
        DO I = 1, NREFC

C.................  Determine Last county
            If ( I .EQ. NREFC ) THEN
               LAST_CNTY = .TRUE.
            END IF

C.............  Determine fuel month for current time step and reference county
            K = FIND1FIRST( MCREFIDX( I,1 ), NREFF, FMREFSORT( :,1 ) )
            M = FIND1( MCREFIDX( I,1 ), NFUELC, FMREFLIST( :,1 ) )

C.................  Determine month
            CALL DAYMON( SDATE, MONTH, DAYMONTH )
                
            IF( K .LT. 0 .OR. M .LT. 0 ) THEN
                WRITE( MESG, 94010 ) 'No fuel month data for ' //
     &            'reference county', MCREFIDX( I,1 )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF
                
            FUELMONTH = 0
            DO J = K, K + FMREFLIST( M,2 )
                IF( FMREFSORT( J,3 ) == MONTH ) THEN
                    FUELMONTH = FMREFSORT( J,2 )
                    EXIT
                END IF
            END DO
                
            IF( FUELMONTH == 0 ) THEN
                WRITE( MESG, 94010 ) 'Could not determine ' //
     &            'fuel month for reference county', MCREFIDX( I,1 ),
     &            'and episode month', MONTH
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

C.............  Read emission factors for reference county and month
            WRITE( MESG,94010 ) 'Processing MOVES lookup ' //
     &           'tables for reference county', MCREFIDX( I,1 ),  
     &           ' of fuel month:', FUELMONTH 
            CALL M3MSG2( MESG )
 
            IF( RPDFLAG ) THEN
                CALL RDRPDEMFACS( I, FUELMONTH )
            END IF
             
            IF( RPVFLAG ) THEN
                CALL RDRPVEMFACS( I, FUELMONTH )
            END IF
                
            IF( RPPFLAG ) THEN
                CALL RDRPPEMFACS( I, FUELMONTH )
            END IF

C.............  Initializations before main time loop 
            JDATE  = SDATE
            JTIME  = STIME
            LDATE  = 0
            DAY    = 1

C.............  Loop through output time steps
            DO T = 1, NSTEPS
                EMGRD = 0.  ! array
                TMPEMGRD = 0.  ! array

C.................  Determine weekday index (Monday is 1)
                DAY = WKDAY( JDATE )
                IF( DAY .GT. 5 ) THEN
                    DAYIDX = 1
                ELSE
                    DAYIDX = 2
                END IF

C.................  Determine month
                CALL DAYMON( JDATE, MONTH, DAYMONTH )

C.................  Write out message for new day.
                IF( JDATE .NE. LDATE ) THEN
c                    CALL WRDAYMSG( JDATE, MESG )

C.....................  Set start hour of day for all sources
                    CALL SETSRCDY( NMSRC, JDATE, TZONES, LDAYSAV, .TRUE.,
     &                             DAYBEGT, DAYENDT )
                END IF

C.................  Write out files that are being used for by-day treatment
                IF( RPDFLAG .AND. DAY .NE. PDAY ) THEN
                    IF( MFLAG_BD ) THEN
                        MESG = '   with MTMP file ' // MTNAME( DAY )
                        CALL M3MSG2( MESG )
                    END IF
                    PDAY = DAY
                END IF

C.................  For new hour...
C.................  Write to screen because WRITE3 only writes to LDEV
c                WRITE( *, 93020 ) HHMMSS( JTIME )

C.................  Initialize current date
                MJDATE = JDATE

C.................  Reset the date when by-day processing is being done
                IF( MFLAG_BD ) MJDATE = MSDATE( DAY )

C.................  In RPD mode, read VMT for current hour
                IF( RPDFLAG ) THEN
                    IF( .NOT. READSET( MTNAME( DAY ), 'VMT', 1, ALLFILES,
     &                                 MJDATE, JTIME, VMT ) ) THEN
                        MESG = 'Could not read VMT from MTMP file'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF
                END IF

C.................  In RPD and RPV modes, read temperatures for current hour
                IF( RPDFLAG .OR. RPVFLAG ) THEN
                    IF( .NOT. READ3( METNAME, TVARNAME, 1, 
     &                               JDATE, JTIME, TEMPG ) ) THEN
                        MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                         ' from ' // TRIM( METNAME )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF
                END IF

C.................  Loop over sources in reference county
                DO S = 1, NREFSRCS( I )
                
                    SRC = REFSRCS( I,S )
                
                    IF( RPDFLAG ) THEN
                        VMTVAL = VMT( SRC )
                    END IF
                    
                    IF( RPPFLAG .OR. RPVFLAG ) THEN
                        VPOPVAL = VPOP( SRC )
                    END IF

C.....................  Determine hour index based on source's local time
                    HOURIDX = ( JTIME - DAYBEGT( SRC ) ) / 10000
                    IF( HOURIDX < 0 ) THEN
                        HOURIDX = HOURIDX + 24
                    END IF
                    HOURIDX = HOURIDX + 1  ! array index is 1 to 24

C.....................  Determine SCC index for source
                    SCCIDX = MISCC( SRC )

C.....................  Determine speed bins for source
                    IF( RPDFLAG ) THEN
                        SPEEDVAL = BADVAL3
                        IF( SPDFLAG ) THEN
                            SPEEDVAL = SPDPRO( MICNY( SRC ), SCCIDX, DAYIDX, HOURIDX )
                        END IF

C.........................  Fall back to inventory speed if hourly speed isn't available
                        IF( SPEEDVAL .LT. AMISS3 ) THEN
                            SPEEDVAL = SPEED( SRC )
                        END IF

                        BIN1 = 0
                        BIN2 = 0
                        DO K = MXSPDBINS, 1, -1
                            IF( SPEEDVAL < SPDBINS( K ) ) CYCLE
    
                            IF( SPEEDVAL == SPDBINS( K ) ) THEN
                                BIN1 = K
                                BIN2 = K
                            ELSE
                                BIN1 = K
                                BIN2 = K + 1
                            END IF
                            EXIT
                        END DO
                        
                        IF( BIN2 > MXSPDBINS ) BIN2 = MXSPDBINS
                        IF( BIN1 == 0 ) THEN
                            BIN1 = 1
                            BIN2 = 1
                        END IF

C.........................  Calculate speed interpolation factor
                        IF( BIN1 .NE. BIN2 ) THEN
                            SPDFAC = ( SPEEDVAL - SPDBINS( BIN1 ) ) / 
     &                               ( SPDBINS( BIN2 ) - SPDBINS( BIN1 ) )
                        ELSE
                            SPDFAC = 0.
                        END IF
                    END IF

C.....................  Determine profiles for current inventory county
C                       There will be 4 profiles used in total:
C                         UU - both min and max profile temps are under county temps
C                         UO - min profile temp is under county min, max profile temp is over county max
C                         OU - min profile temp is over county min, max profile temp is under county max
C                         OO - both min and max profile temps are over county temps
                    IF( RPPFLAG ) THEN
                        MINVAL = AVGMIN( MICNY( SRC ), MONTH, DAYMONTH )
                        MAXVAL = AVGMAX( MICNY( SRC ), MONTH, DAYMONTH )

C.........................  Check that min and max county temps were found
                        IF( MINVAL .LT. AMISS3 .OR.
     &                      MAXVAL .LT. AMISS3 ) THEN
                            WRITE( MESG, 94010 ) 'Could not find minimum ' //
     &                        'and maximum temperatures for county', IFIP( SRC ),
     &                        'and episode month', MONTH
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF

C.........................  Find indexes of bounding minimum temperatures
                        USTART = 0
                        OSTART = 0
                        PDIFF = 999
                        NO_INTRPLT = .FALSE.     ! no interpolation flag
                        DO K = 1, NEMTEMPS
                            TDIFF = MINVAL - EMTEMPS( EMTEMPIDX( K ) )

C.............................  Once profile min temp is greater than county temp, this loop is done                            
                            IF( TDIFF .LT. 0 ) THEN
                                OSTART = K
                                OMIN = EMTEMPS( EMTEMPIDX( K ) )
                                EXIT
                            END IF

C.............................  If current profile min temp is closer to county temp, store index                            
                            IF( TDIFF .LT. PDIFF ) THEN
                                USTART = K
                                UMIN = EMTEMPS( EMTEMPIDX( K ) )
                            END IF
                        
                            PDIFF = TDIFF
                        END DO

C.........................  Check that appropriate minimum temperatures were found
                        IF( USTART == 0 ) THEN
                            WRITE( MESG, 94040 ) 'ERROR: Lowest profile ' //
     &                        'minimum temperature', 
     &                        EMTEMPS( EMTEMPIDX( 1 ) ),
     &                        'is higher than county minimum temperature',
     &                        MINVAL
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF
                        
                        IF( OSTART == 0 ) THEN
                            WRITE( MESG, 94040 ) 'ERROR: Highest profile ' //
     &                        'minimum temperature', 
     &                        EMTEMPS( EMTEMPIDX( NEMTEMPS ) ),
     &                        'is lower than county minimum temperature',
     &                        MINVAL
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF

C.........................  Find indexes of bounding maximum temperatures
                        UUIDX = 0
                        UOIDX = 0
                        UMAX = EMXTEMPS( EMTEMPIDX( USTART ) )
                        DO K = USTART, NEMTEMPS

C.............................  Check that profile minimum temperature hasn't changed
                            IF( EMTEMPS( EMTEMPIDX( K ) ) .NE. UMIN ) THEN
                                EXIT
                            END IF
                            
                            OMAX = EMXTEMPS( EMTEMPIDX( K ) )
                        
                            IF( OMAX > MAXVAL ) THEN
                                UOIDX = EMTEMPIDX( K )
                                IF( K > USTART ) THEN
                                    UUIDX = EMTEMPIDX( K - 1 )
                                END IF
                                EXIT
                            END IF
                        END DO
C.........................  Check that appropriate maximum temperatures were found
                        IF( UUIDX .EQ. 0 .AND. UOIDX .NE. 0 ) THEN
                            WRITE( MESG, 94040 ) 'ERROR: Lowest profile of ' //
     &                        'max temperature', UMAX,'- min temperature',UMIN,
     &                        ' is higher than county maximum temperature:',
     &                        MAXVAL
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF
                        
                        IF( UUIDX .EQ. 0 .AND. UOIDX .EQ. 0 ) THEN
                            WRITE( MESG, 94040 ) 'ERROR: Highest profile of ' //
     &                        'max temperature', OMAX,'- min temperature',UMIN, 
     &                        ' is lower than county maximum temperature',
     &                        MAXVAL
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF

                        OUIDX = 0
                        OOIDX = 0
                        UMAX = EMXTEMPS( EMTEMPIDX( OSTART ) )
                        DO K = OSTART, NEMTEMPS

C.............................  Check that profile minimum temperature hasn't changed
                            IF( EMTEMPS( EMTEMPIDX( K ) ) .NE. OMIN ) THEN
                                EXIT
                            END IF
                            
                            OMAX = EMXTEMPS( EMTEMPIDX( K ) )
                        
                            IF( OMAX > MAXVAL ) THEN
                                OOIDX = EMTEMPIDX( K )
                                IF( K > OSTART ) THEN
                                    OUIDX = EMTEMPIDX( K - 1 )
                                END IF
                                EXIT
                            END IF
                        END DO

C.........................  Determine whether min/max temp within one temp bin
                        IF( OMAX .EQ. OMIN ) NO_INTRPLT = .TRUE.

C.........................  Check that appropriate maximum temperatures were found
                        IF( .NOT. NO_INTRPLT ) THEN
                          IF( OUIDX .EQ. 0 .AND. OOIDX .NE. 0 ) THEN
                            WRITE( MESG, 94040 ) 'ERROR: Lowest profile of ' //
     &                        'max temperature', UMAX,'- min temperature',OMIN,
     &                        ' is higher than county maximum temperature:',
     &                        MAXVAL
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                          END IF

                          IF( OUIDX .EQ. 0 .AND. OOIDX .EQ. 0 ) THEN
                            WRITE( MESG, 94040 ) 'ERROR: Highest profile of ' //
     &                        'max temperature', OMAX,'- min temperature',OMIN,
     &                        ' is lower than county maximum temperature',
     &                        MAXVAL
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                          END IF

C.........................  Check that maximum temperatures of profiles match
                          IF( EMXTEMPS( UUIDX ) .NE. EMXTEMPS( OUIDX ) .OR.
     &                      EMXTEMPS( UOIDX ) .NE. EMXTEMPS( OOIDX ) ) THEN
                            MESG = 'ERROR: Inconsistent temperature ' //
     &                        'profiles.'
C                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                            CALL M3MESG( MESG )
                          END IF
                        
                          MINFAC = ( MINVAL - EMTEMPS( UUIDX ) ) /
     &                           ( EMTEMPS( OUIDX ) - EMTEMPS( UUIDX ) )

                          MAXFAC = ( MAXVAL - EMXTEMPS( OUIDX ) ) /
     &                           ( EMXTEMPS( OOIDX ) - EMXTEMPS( OUIDX ) )

                      END IF

                    END IF

C.....................  Loop over grid cells for this source
                    DO NG = 1, NSRCCELLS( SRC )

                        CELL = SRCCELLS( SRC, NG )
                        GFRAC = SRCCELLFRACS( SRC, NG )
                        
                        IF( RPDFLAG .OR. RPVFLAG ) THEN
                            TEMPVAL = ( TEMPG( CELL ) - CTOK ) * CTOF + 32.
    
C.............................  Determine temperature indexes for cell
                            IDX1 = 0
                            IDX2 = 0
                            DO K = NEMTEMPS, 1, -1
                                IF( TEMPVAL < EMTEMPS( K ) ) CYCLE
                                
                                IF( TEMPVAL == EMTEMPS( K ) ) THEN
                                    IDX1 = K
                                    IDX2 = K
                                ELSE
                                    IDX1 = K
                                    IDX2 = K + 1
                                END IF
                                EXIT
                            END DO
    
C.............................  Check that an appropriate temperature was found
                            IF( IDX1 == 0 .OR. IDX2 > NEMTEMPS ) THEN
                                WRITE( MESG, 94040 ) 'ERROR: Grid cell ' //
     &                            'temperature', TEMPVAL, 'out of range ' //
     &                            'of emission factor data.'
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                            END IF
                            
                            IF( IDX1 .NE. IDX2 ) THEN
                                TEMPFAC = ( TEMPVAL - EMTEMPS( IDX1 ) ) /
     &                                    ( EMTEMPS( IDX2 ) - EMTEMPS( IDX1 ) )
                            ELSE
                                TEMPFAC = 1.
                            END IF
                        END IF

C.........................  Loop through pollutant-species combos
                        LBUF = ' '
                        DO V = 1, NSMATV
                        
                            PBUF = EANAM( SIINDEX( V,1 ) )  ! process/pollutant

                            PROCIDX = EMPROCIDX( V )
                            POLIDX = EMPOLIDX( V )

C.............................  Check if emission factors exist for this process/pollutant
                            IF( PROCIDX .EQ. 0 .OR. POLIDX .EQ. 0 ) THEN
                                CYCLE
                            END IF

C.............................  Calculate interpolated emission factor if process/pollutant has changed
                            IF( PBUF .NE. LBUF ) THEN
                                IF( RPDFLAG ) THEN
                                    IF( BIN1 .NE. BIN2 ) THEN
                                        EFVAL1 = RPDEMFACS( SCCIDX, BIN1, IDX1, PROCIDX, POLIDX )
                                        EFVAL2 = RPDEMFACS( SCCIDX, BIN2, IDX1, PROCIDX, POLIDX )
                                        EFVALA = SPDFAC * (EFVAL2 - EFVAL1) + EFVAL1
        
                                        EFVAL1 = RPDEMFACS( SCCIDX, BIN1, IDX2, PROCIDX, POLIDX )
                                        EFVAL2 = RPDEMFACS( SCCIDX, BIN2, IDX2, PROCIDX, POLIDX )
                                        EFVALB = SPDFAC * (EFVAL2 - EFVAL1) + EFVAL1
                                    ELSE
                                        EFVALA = RPDEMFACS( SCCIDX, BIN1, IDX1, PROCIDX, POLIDX )
                                        EFVALB = RPDEMFACS( SCCIDX, BIN1, IDX2, PROCIDX, POLIDX )
                                    END IF
                                    
                                    EFVAL = TEMPFAC * (EFVALB - EFVALA) + EFVALA
                                END IF
                                
                                IF( RPVFLAG ) THEN
                                    EFVALA = RPVEMFACS( DAYIDX, SCCIDX, HOURIDX, IDX1, PROCIDX, POLIDX )
                                    EFVALB = RPVEMFACS( DAYIDX, SCCIDX, HOURIDX, IDX2, PROCIDX, POLIDX )

                                    EFVAL = TEMPFAC * (EFVALB - EFVALA) + EFVALA
                                END IF
                                
                                IF( RPPFLAG ) THEN
                                    EFVAL1 = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, UUIDX, PROCIDX, POLIDX )
                                    EFVAL2 = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, UOIDX, PROCIDX, POLIDX )
                                    EFVALA = MAXFAC * (EFVAL2 - EFVAL1) + EFVAL1
                                    
                                    EFVAL1 = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, OUIDX, PROCIDX, POLIDX )
                                    EFVAL2 = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, OOIDX, PROCIDX, POLIDX )
                                    EFVALB = MAXFAC * (EFVAL2 - EFVAL1) + EFVAL1

                                    EFVAL = MINFAC * (EFVALB - EFVALA) + EFVALA
                                    
                                    IF( NO_INTRPLT ) THEN
                                        EFVAL = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, UOIDX, PROCIDX, POLIDX )
                                    ENDIF
                                END IF

                            END IF
                            
                            LBUF = PBUF

C.............................  Set units conversion factor
                            F1 = GRDFAC( SPINDEX( V,1 ) )
                            F2 = TOTFAC( SPINDEX( V,1 ) )

C.............................  Calculate gridded, hourly emissions
                            IF( RPDFLAG ) THEN
                                EMVAL = VMTVAL * EFVAL * GFRAC
                            END IF
                            
                            IF( RPVFLAG .OR. RPPFLAG ) THEN
                                EMVAL = VPOPVAL * EFVAL * GFRAC
                            END IF

                            EMGRD( CELL,SPINDEX( V,1 ) ) = 
     &                          EMGRD( CELL,SPINDEX( V,1 ) ) + 
     &                          EMVAL * MSMATX_L( SRC,V ) * F1

C.............................  Add this cell's emissions to source totals
                            IF( LREPSTA .OR. LREPCNY .OR. LREPSCC .OR. LREPSRC ) THEN
                                IF( .NOT. SKIPSUM ) THEN
                                    MEBSUM( SRC,SPINDEX( V,1 )) =
     &                                  MEBSUM( SRC,SPINDEX( V,1 )) + 
     &                                  EMVAL * MSMATX_S( SRC,V ) * F2
                            
                                    IF( EANAMREP( V ) ) THEN
                                       F2 = TOTFAC( NMSPC+SIINDEX( V,1 ) )
                                        MEBSUM( SRC,NMSPC+SIINDEX( V,1 ) ) =
     &                                      MEBSUM( SRC,NMSPC+SIINDEX( V,1 ) ) +
     &                                      EMVAL * F2
                                    END IF
                                END IF
                            END IF

                        END DO    ! end loop over pollutant-species combos

                    END DO    ! end loop over grid cells for source
                
                END DO    ! end loop over sources in inv. county
                DO V = 1, NMSPC 
                    SBUF = EMNAM( V )
C.................  Read out old data if not first county
                    IF ( I > 1 ) THEN
  		        IF(.NOT. READSET( MONAME, SBUF, 1,  ALLFILES,
     &                            JDATE, JTIME, TMPEMGRD( 1,V ) ) )THEN
                             MESG = 'Could not read "' // SBUF // '" ' //
     &                         'from file "' // MONAME // '"'
                             CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF
C.................  sum old county data with new county
                        EMGRD( :,V ) = EMGRD( :,V ) + TMPEMGRD( :,V )
                    END IF

                    IF( LGRDOUT ) THEN
                        IF( .NOT. WRITESET( MONAME, SBUF, ALLFILES,
     &                              JDATE, JTIME, EMGRD( 1,V ) ) ) THEN
                            MESG = 'Could not write "' // SBUF // '" ' //
     &                       'to file "' // MONAME // '"'
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF
                    END IF
                END DO
           
C.............  Write state, county, and SCC emissions (all that apply) 
C.............  The subroutine will only write for certain hours and 
C               will reinitialize the totals after output
                If ( LAST_CNTY) THEN
                    CALL WRMRGREP( JDATE, JTIME )
                END IF
C...............Do not count the hour after 230000
                SKIPSUM = .FALSE. 
                IF( JTIME .EQ. 230000 ) SKIPSUM = .TRUE.

                LDATE = JDATE

                CALL NEXTIME( JDATE, JTIME, TSTEP )     !  update model clock

            END DO   ! End loop on time steps

        END DO   ! end loop over inventory counties

C.........  Close output file
        IF( .NOT. CLOSESET( MONAME ) ) THEN
            MESG = 'Could not close file:"'// MONAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        DEALLOCATE( DAYBEGT, DAYENDT, LDAYSAV )
        DEALLOCATE( TEMPG, EMGRD, TMPEMGRD, VARNAMES, MGMATX )

C.........  Successful completion of program
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( 8X, 'at time ', A8 )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( A )

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

94020   FORMAT( A, I4, 2X, 10 ( A, :, 1PG14.6, :, 2X ) )

94030   FORMAT( 8X, 'at time ', A8 )

94040   FORMAT( 10 ( A, :, F10.2, :, 2X ) )

        END PROGRAM MOVESMRG 

