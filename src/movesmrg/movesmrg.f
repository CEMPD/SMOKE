
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
C       Created 3/10 by C. Seppanen - based on smkmerge.f
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
     &          LREPSTA, LREPANY, LGRDOUT,         ! state report, any reports, gridded output
     &          CDEV,                              ! costcy
     &          MGNAME, MTNAME, MSNAME, MONAME,    ! input files
     &          NMSRC, MNGMAT,                     ! no. of srcs, no. gridding matrix entries
     &          NMSPC,                             ! no. species
     &          EMNAM,                             ! species names
     &          TSVDESC,                           ! var names
     &          SIINDEX, SPINDEX,                  ! EANAM & EMNAM idx
     &          SDATE, STIME, NSTEPS, TSTEP,       ! episode information
     &          MSDATE,                            ! dates for by-day hrly emis
     &          GRDFAC, TOTFAC,                    ! conversion factors
     &          MSMATX, MNSMATV, NSMATV,           ! speciation matrices
     &          MEBCNY, MEBSTA,                    ! cnty/state total spec emissions
     &          EANAM                              ! pol/act names

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: RPDFLAG, RPVFLAG, RPPFLAG, 
     &          TVARNAME, METNAME,
     &          NREFSRCS, REFSRCS, NSRCCELLS, SRCCELLS, SRCCELLFRACS,
     &          EMPROCIDX, EMPOLIDX,
     &          NEMTEMPS, EMTEMPS, EMXTEMPS, AVGMIN, AVGMAX,
     &          RPDEMFACS, RPVEMFACS, RPPEMFACS

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
        REAL, ALLOCATABLE :: EMGRD( :,: )   ! emissions for each grid cell and species

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
        REAL             MINFAC, MAXFAC  ! min and max temp interpolation factors
        REAL             PDIFF, TDIFF  ! temperature differences
        REAL             EFVAL1, EFVAL2, EFVALA, EFVALB, EFVAL   ! emission factor values
        REAL             EMVAL         ! emissions value

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
        EMGRD = 0.  ! array
        
        IF( LREPSTA ) THEN
            ALLOCATE( MEBSTA( NSTATE, NMSPC ), STAT=IOS )    ! state totals
            CALL CHECKMEM( IOS, 'MEBSTA', PROGNAME )
        END IF

        ALLOCATE( MEBCNY( NCOUNTY, NMSPC ), STAT=IOS )    ! county totals
        CALL CHECKMEM( IOS, 'MEBCNY', PROGNAME )

        ALLOCATE( MSMATX( NMSRC, MNSMATV ), STAT=IOS )    ! speciation matrix
        CALL CHECKMEM( IOS, 'MSMATX', PROGNAME )
        
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

C.............  Read speciation matrix for current variable
            CALL RDSMAT( MSNAME, VBUF, MSMATX( 1,V ) )

        END DO

C.........  Write out message with list of species
        CALL POLMESG( OCNT, VARNAMES )

C.........  Initializations before main time loop 
        JDATE  = SDATE
        JTIME  = STIME
        LDATE  = 0
        DAY    = 1

C.........  Loop through output time steps
        DO T = 1, NSTEPS

C.............  Determine weekday index (Monday is 1)
            DAY = WKDAY( JDATE )
            IF( DAY .GT. 5 ) THEN
                DAYIDX = 1
            ELSE
                DAYIDX = 2
            END IF

C.............  Determine month
            CALL DAYMON( JDATE, MONTH, DAYMONTH )

C.............  Write out message for new day.
            IF( JDATE .NE. LDATE ) THEN
                CALL WRDAYMSG( JDATE, MESG )

C.................  Set start hour of day for all sources
                CALL SETSRCDY( NMSRC, JDATE, TZONES, LDAYSAV, .TRUE.,
     &                         DAYBEGT, DAYENDT )
            END IF

C.............  Write out files that are being used for by-day treatment
            IF( RPDFLAG .AND. DAY .NE. PDAY ) THEN
                IF( MFLAG_BD ) THEN
                    MESG = '   with MTMP file ' // MTNAME( DAY )
                    CALL M3MSG2( MESG )
                END IF

                PDAY = DAY
            END IF

C.............  For new hour...
C.............  Write to screen because WRITE3 only writes to LDEV
            WRITE( *, 93020 ) HHMMSS( JTIME )

C.............  Initialize current date
            MJDATE = JDATE

C.............  Reset the date when by-day processing is being done
            IF( MFLAG_BD ) MJDATE = MSDATE( DAY )

C.............  In RPD mode, read VMT for current hour
            IF( RPDFLAG ) THEN
                IF( .NOT. READSET( MTNAME( DAY ), 'VMT', 1, ALLFILES,
     &                             MJDATE, JTIME, VMT ) ) THEN
                    MESG = 'Could not read VMT from MTMP file'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF
            END IF

C.............  In RPD and RPV modes, read temperatures for current hour
            IF( RPDFLAG .OR. RPVFLAG ) THEN
                IF( .NOT. READ3( METNAME, TVARNAME, 1, 
     &                           JDATE, JTIME, TEMPG ) ) THEN
                    MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                     ' from ' // TRIM( METNAME )
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF
            END IF

C.............  Loop over reference counties
            DO I = 1, NREFC

C.................  Determine fuel month for current time step and reference county
                K = FIND1FIRST( MCREFIDX( I,1 ), NREFF, FMREFSORT( :,1 ) )
                M = FIND1( MCREFIDX( I,1 ), NFUELC, FMREFLIST( :,1 ) )
                
                IF( K .LT. 0 .OR. M .LT. 0 ) THEN
                    WRITE( MESG, 94010 ) 'No fuel month data for ' //
     &                'reference county', MCREFIDX( I,1 )
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
     &                'fuel month for reference county', MCREFIDX( I,1 ),
     &                'and episode month', MONTH
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF

C.................  Read emission factors for reference county and month
                IF( RPDFLAG ) THEN
                    CALL RDRPDEMFACS( I, FUELMONTH )
                END IF
                
                IF( RPVFLAG ) THEN
                    CALL RDRPVEMFACS( I, FUELMONTH )
                END IF
                
                IF( RPPFLAG ) THEN
                    CALL RDRPPEMFACS( I, FUELMONTH )
                END IF

C.................  Loop over sources in reference county
                DO S = 1, NREFSRCS( I )
                
                    SRC = REFSRCS( I,S )
                
                    IF( RPDFLAG ) THEN
                        VMTVAL = VMT( SRC )
                        SPEEDVAL = SPEED( SRC )
                    END IF
                    
                    IF( RPPFLAG .OR. RPVFLAG ) THEN
                        VPOPVAL = VPOP( SRC )

C.........................  Determine hour index based on source's local time
                        HOURIDX = ( JTIME - DAYBEGT( SRC ) ) / 10000
                        IF( HOURIDX < 0 ) THEN
                            HOURIDX = HOURIDX + 24
                        END IF
                        HOURIDX = HOURIDX + 1  ! array index is 1 to 24
                    END IF
                    
                    SCC = CSCC( SRC )

C.....................  Determine SCC index for source
                    SCCIDX = FINDC( SCC, NINVSCC, INVSCC )
                    
C.....................  Determine speed bins for source
                    IF( RPDFLAG ) THEN
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
                    IF( RPPFLAG ) THEN
                        MINVAL = AVGMIN( MICNY( SRC ) )
                        MAXVAL = AVGMAX( MICNY( SRC ) )
                        
                        UUIDX = 0
                        UOIDX = 0
                        OUIDX = 0
                        OOIDX = 0
                        
                        USTART = 0
                        UEND = 0
                        OSTART = 0
                        OEND = 0
                        
                        PDIFF = 999
                        DO K = NEMTEMPS, 1, -1
                            TDIFF = EMTEMPS( K ) - MINVAL

C.............................  Each time the temperature difference is smaller,
C                               this could potentially be the ending index
                            IF( TDIFF .GT. 0 .AND. TDIFF .LT. PDIFF ) THEN
                                OEND = K
                            END IF
                            
                            IF( PDIFF .GT. 0 .AND. TDIFF .LT. 0 ) THEN
                                OSTART = K + 1
                                UEND = K
                            END IF
                            
                            IF( PDIFF .LT. 0 .AND. TDIFF .LT. PDIFF ) THEN
                                USTART = K + 1
                                EXIT
                            END IF
                            
                            PDIFF = TDIFF
                        END DO
                        
                        IF( OEND == 0 ) THEN
                            MESG = 'ERROR: Highest profile minimum ' //
     &                        'temperature is lower than county ' //
     &                        'minimum temperature.'
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF
                        
                        IF( UEND == 0 ) THEN
                            MESG = 'ERROR: Lowes profile minimum ' //
     &                        'temperature is higher than county ' //
     &                        'minimum temperature.'
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF
                        
                        IF( USTART == 0 ) THEN
                            USTART = 1
                        END IF
                        
                        DO K = UEND, USTART, -1
                            IF( EMXTEMPS( K ) > MAXVAL ) CYCLE
                            
                            UOIDX = K + 1
                            UUIDX = K
                        END DO
                        
                        DO K = OEND, OSTART, -1
                            IF( EMXTEMPS( K ) > MAXVAL ) CYCLE
                            
                            OOIDX = K + 1
                            OUIDX = K
                        END DO
                        
                        MINFAC = ( MINVAL - EMTEMPS( UUIDX ) ) /
     &                           ( EMTEMPS( OUIDX ) - EMTEMPS( UUIDX ) )
     
                        MAXFAC = ( MAXVAL - EMXTEMPS( OUIDX ) ) /
     &                           ( EMXTEMPS( OOIDX ) - EMXTEMPS( OUIDX ) )
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
                                END IF
                            END IF
                            
                            LBUF = PBUF

C.............................  Set units conversion factor
                            F1 = GRDFAC( SPINDEX( V,1 ) )
                            F2 = TOTFAC( SPINDEX( V,1 ) )

C.............................  Calculate hourly emissions
                            IF( RPDFLAG ) THEN
                                EMVAL = VMTVAL * EFVAL * GFRAC * MSMATX( SRC,V )
                            END IF
                            
                            IF( RPVFLAG .OR. RPPFLAG ) THEN
                                EMVAL = VPOPVAL * EFVAL * GFRAC * MSMATX( SRC,V )
                            END IF

                            EMGRD( CELL,SPINDEX( V,1 ) ) = 
     &                          EMGRD( CELL,SPINDEX( V,1 ) ) + EMVAL * F1

C.............................  Add this cells emissions to county totals
                            MEBCNY( MICNY( SRC ),SPINDEX( V,1 ) ) =
     &                          MEBCNY( MICNY( SRC ),SPINDEX( V,1 ) ) + EMVAL * F2

                        END DO

                    END DO    ! end loop over grid cells for source
                
                END DO    ! end loop over sources in inv. county

            END DO    ! end loop over inventory counties

C.............  Output gridded emissions for all species
            DO V = 1, NMSPC
            
                SBUF = EMNAM( V )

C.................  Write out gridded data
                IF( LGRDOUT ) THEN
                    IF( .NOT. WRITESET( MONAME, SBUF, ALLFILES,
     &                                  JDATE, JTIME, EMGRD( 1,V ) ) ) THEN
                        MESG = 'Count not write "' // SBUF // '" ' //
     &                    'to file "' // MONAME // '"'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF
                END IF
            END DO

C.............  Initialize gridded emissions
            EMGRD = 0.   ! array

C.............  Write country, state, and county emissions (all that apply) 
C.............  The subroutine will only write for certain hours and 
C               will reinitialize the totals after output
            IF( LREPANY ) THEN
                CALL WRMRGREP( JDATE, JTIME )
            END IF

            LDATE = JDATE

            CALL NEXTIME( JDATE, JTIME, TSTEP )     !  update model clock

        END DO          ! End loop on time steps

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

