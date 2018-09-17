
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
     &          MFLAG_BD, LREPANY,                          ! by-day hourly emis flags
     &          LREPSTA, LREPCNY, LREPSCC, LREPSRC, LGRDOUT,! report flags, gridded output
     &          CDEV,                              ! costcy
     &          MGNAME, MTNAME, MONAME,            ! input files
     &          NMSRC, MNGMAT, MGMATX,             ! no. of srcs, no. gridding matrix entries
     &          NMSPC, MCFIP,                      ! no. species
     &          EMNAM,                             ! species names
     &          TSVDESC, SMATCHK,                  ! var names, use SMAT for model spc calc
     &          SIINDEX, SPINDEX,                  ! EANAM & EMNAM idx
     &          SDATE, STIME, NSTEPS, TSTEP,       ! episode information
     &          MSDATE,                            ! dates for by-day hrly emis
     &          GRDFAC, TOTFAC,                    ! conversion factors
     &          NSMATV,                            ! speciation matrices
     &          MEBCNY, MEBSTA, MEBSUM,            ! cnty/state/src total spec emissions
     &          MEBSCC, MEBSTC,                    ! scc total spec emissions
     &          EANAM, NIPPA, MTMPNAME,            ! pol/act names
     &          CFDEV,                             ! control factor file
     &          SRCGRPFLAG, SGDEV, ISRCGRP,        ! source groups
     &          EMGGRD, EMGGRDSPC, EMGGRDSPCT,     ! emissions by source group
     &          SGINLNNAME, NSGOUTPUT              ! source group emissions output file

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: RPDFLAG, RPHFLAG, RPVFLAG, RPPFLAG, MOPTIMIZE,
     &          TVARNAME, METNAME,
     &          NREFSRCS, REFSRCS, NSRCCELLS, SRCCELLS, SRCCELLFRACS,
     &          EMPOLIDX,
     &          NEMTEMPS, EMTEMPS, EMXTEMPS, EMTEMPIDX,
     &          RPDEMFACS, RPHEMFACS, RPVEMFACS, RPPEMFACS,
     &          SPDFLAG, SPDPRO, MISCC,
     &          MSNAME_L, MSMATX_L, MNSMATV_L, GRDENV,
     &          MSNAME_S, MSMATX_S, MNSMATV_S,
     &          EANAMREP, CFPRO, CFFLAG,
     &          NMVSPOLS, MVSPOLNAMS,              ! MOVES lookup poll/spc names
     &          TEMPBIN, MTMP_OUT

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVCFIP, NINVSCC, INVSCC

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: MICNY, NCOUNTY, NSTATE

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NGRID

C.........  This module is used for reference county information
        USE MODMBSET, ONLY: NREFC, MCREFIDX,
     &                      NREFF, FMREFSORT, NFUELC, FMREFLIST

C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: SPEED, CSCC, VPOP, CIFIP, TZONES

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
        CHARACTER(2)    CRLF
        INTEGER         INDEX1
        INTEGER         FINDCFIRST
        INTEGER         FINDC
        INTEGER         WKDAY
        INTEGER         ENVINT
        INTEGER         STR2INT

        EXTERNAL    HHMMSS, INDEX1, FINDCFIRST, FINDC, WKDAY, 
     &              STR2INT, ENVINT, CRLF

C.........  LOCAL PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.6_Sep2018$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

C...........   Local arrays for per-source information
        INTEGER, ALLOCATABLE :: DAYBEGT( : )   ! daily start time for each source
        INTEGER, ALLOCATABLE :: DAYENDT( : )   ! daily end time for each source
        LOGICAL, ALLOCATABLE :: LDAYSAV( : )   ! true: src uses DST

C...........   Local arrays for hourly data
        REAL, ALLOCATABLE :: VMT( : )
        REAL, ALLOCATABLE :: HOTEL( : )
        REAL, ALLOCATABLE :: TEMPG( : )    ! Temp array for RPD and RPV
        REAL, ALLOCATABLE :: MAXTEMP( : )     ! max temp array for RPP
        REAL, ALLOCATABLE :: MINTEMP( : )     ! min temp array for RPP
        REAL, ALLOCATABLE :: MTMP_INVT( :,:,: )  ! tmp hourly emissoin for output temporal int output file
        REAL, ALLOCATABLE :: EMGRD( :,: )     ! emissions for each grid cell and species
        REAL, ALLOCATABLE :: TEMGRD( :,:,:)   ! emissions for each grid cell, species, and time steps
        REAL, ALLOCATABLE :: TMPEMGRD( :,: )  ! tmp emissions for each grid cell and species
        REAL, ALLOCATABLE :: TMPEMGGRD( : )   ! tmp emissions for output source groups

C...........   Logical names and unit numbers (not in MODMERGE)
        INTEGER         LDEV
     
C...........   Other local variables
    
        INTEGER          I, J, K, L1, L2, M, N, NG, V, NV, S, T ! counters and indices

        INTEGER          BIN1, BIN2    ! speed bins for current source
        INTEGER          CELL          ! current grid cell
        INTEGER          DAY           ! day-of-week index (monday=1)
        INTEGER          DAYMONTH      ! day-of-month
        INTEGER          DAYIDX        ! current day value index
        INTEGER          FMON, NMON    ! tmp fuel month
        INTEGER          FUELMONTH     ! current fuel month
        INTEGER          LFUELMONTH    ! last fuel month
        INTEGER          HOURIDX       ! current hour of the day
        INTEGER          IDX1, IDX2    ! temperature indexes for current cell
        INTEGER          IOS           ! tmp I/O status
        INTEGER          PDATE         ! Julian date (YYYYDDD) for RPP mode
	INTEGER          DDATE,DTIME   ! local date and time
        INTEGER          JDATE         ! Julian date (YYYYDDD)
        INTEGER          JTIME         ! time (HHMMSS)
        INTEGER          RDATE         ! last reporting Julian date (YYYYDDD)
        INTEGER          RTIME         ! last reporting time (HHMMSS)
        INTEGER       :: K1 = 0        ! tmp index
        INTEGER       :: K5 = 0        ! tmp index
        INTEGER          KM            ! tmp index to src-category species
        INTEGER          LDATE         ! Julian date from previous iteration
        INTEGER          MJDATE        ! mobile-source Julian date for by-day
        INTEGER          MONTH         ! current month
        INTEGER          LMONTH        ! month from previous iteration
        INTEGER       :: PDAY = 0      ! previous iteration day no.
        INTEGER          NPOLSPC       ! total number of pol/spc
        INTEGER          SIIDX         ! tmp pollutant idx
        INTEGER          SPIDX         ! tmp species idx
        INTEGER          POLIDX        ! current pollutant index
        INTEGER          SCCIDX        ! current SCC index
        INTEGER          SRC           ! current source number
        INTEGER          UUIDX, UOIDX, OUIDX, OOIDX  ! indexes for matching profiles
        INTEGER          USTART, UEND, OSTART, OEND
        INTEGER          MXWARN        !  maximum number of warnings
        INTEGER       :: NWARN = 0     !  current number of warnings
        INTEGER          GIDX          ! index to source group

        REAL             F1, F2        ! tmp conversion
        REAL             GFRAC         ! grid cell fraction
        REAL             SPEEDVAL      ! average speed value for current source
        REAL             TEMPVAL       ! temperature value for current grid cell
        REAL             SPDFAC        ! speed interpolation factor
        REAL             TEMPFAC       ! temperature interpolation factor
        REAL             MINTVAL, MAXTVAL  ! min and max temperature for current source
        REAL             UMIN, OMIN    ! bounding minimum temperature profile values
        REAL             UMAX, OMAX    ! bounding maximum temperature profile values
        REAL             MINFAC, MAXFAC    ! min and max temp interpolation factors
        REAL             PDIFF, TDIFF  ! temperature differences
        REAL             EFVAL1, EFVAL2, EFVALA, EFVALB, EFVAL   ! emission factor values
        REAL             EMVAL, EMFAC  ! emissions value
        REAL             EMVALSPC      ! speciated emissions value
        REAL          :: CFFAC = 1.0   ! control factor

        LOGICAL       :: NO_INTRPLT = .FALSE.   ! true: single interploation, false: bi-interpolation

        CHARACTER(300)     MESG    ! message buffer
        CHARACTER( 4 )     YEAR    ! modelin year
        CHARACTER( 7 )     TDATE   ! tmp julinan date
        CHARACTER(IOVLEN3) CPOL    ! tmp pollutant or emission type name
        CHARACTER(IOVLEN3) CSPC    ! tmp species name
        CHARACTER(IOVLEN3) LSPC    ! previous tmp species name
        CHARACTER(SCCLEN3) SCC     ! current source SCC
        CHARACTER(PLSLEN3) VBUF    ! pol to species or pol description buffer

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
        CALL RDSTCY( CDEV, NINVIFIP, INVCFIP )

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

C........ when not optimize memory 
        IF ( MOPTIMIZE ) THEN
            ALLOCATE( EMGRD( NGRID, NMSPC ), STAT=IOS )     ! gridded emissions
            CALL CHECKMEM( IOS, 'EMGRD', PROGNAME )
        ELSE
            ALLOCATE( TEMGRD( NGRID, NMSPC, NSTEPS ), STAT=IOS )     ! gridded emissions
            CALL CHECKMEM( IOS, 'TEMGRD', PROGNAME )
            TEMGRD = 0.  ! array
        END IF

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

        IF( LREPANY ) THEN        
            ALLOCATE( MEBSUM( NMSRC, NMSPC+NIPPA ), STAT=IOS )    ! source totals
            CALL CHECKMEM( IOS, 'MEBSUM', PROGNAME )
        END IF

        IF( MTMP_OUT ) THEN
            ALLOCATE( MTMP_INVT( NMSRC, NIPPA, NSTEPS ), STAT=IOS )    ! source totals
            CALL CHECKMEM( IOS, 'MTMP_INVT', PROGNAME )
            MTMP_INVT = 0.0
        END IF

        ALLOCATE( EMPOLIDX( NIPPA + NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMPOLIDX', PROGNAME )
        EMPOLIDX = 0   ! array

        IF( SMATCHK ) THEN
            ALLOCATE( MSMATX_L( NMSRC, MNSMATV_L ), STAT=IOS )    ! mole speciation matrix
            CALL CHECKMEM( IOS, 'MSMATX_L', PROGNAME )

            ALLOCATE( MSMATX_S( NMSRC, MNSMATV_S ), STAT=IOS )    ! mass speciation matrix
            CALL CHECKMEM( IOS, 'MSMATX_S', PROGNAME )
        END IF

        IF( RPDFLAG ) THEN
            ALLOCATE( VMT( NMSRC ), STAT=IOS )     ! hourly VMT
            CALL CHECKMEM( IOS, 'VMT', PROGNAME )
        END IF

        IF( RPHFLAG ) THEN
            ALLOCATE( HOTEL( NMSRC ), STAT=IOS )     ! hourly hotelling 
            CALL CHECKMEM( IOS, 'HOTEL', PROGNAME )
        END IF
        
        IF( RPDFLAG .OR. RPHFLAG .OR. RPVFLAG ) THEN
            ALLOCATE( TEMPG( NGRID ), STAT=IOS )    ! hourly temperatures
            CALL CHECKMEM( IOS, 'TEMPG', PROGNAME )
        ELSE
            ALLOCATE( MAXTEMP( NGRID ), STAT=IOS )    ! hourly temperatures
            CALL CHECKMEM( IOS, 'MAXTEMP', PROGNAME )
            ALLOCATE( MINTEMP( NGRID ), STAT=IOS )    ! hourly temperatures
            CALL CHECKMEM( IOS, 'MINTEMP', PROGNAME )
        END IF

C.........  Determine sources that observe DST
        CALL GETDYSAV( NMSRC, CIFIP, LDAYSAV )

C.........  Read gridding matrix
        CALL RDGMAT( MGNAME, NGRID, MNGMAT, MNGMAT,
     &               MGMATX( 1 ), MGMATX( NGRID + 1 ),
     &               MGMATX( NGRID + MNGMAT + 1 ) )

C.........  Build list of grid cells for each source
        CALL BLDSRCCELL( NMSRC, NGRID, MNGMAT, MGMATX( 1 ), 
     &                   MGMATX( NGRID + 1 ), MGMATX( NGRID + MNGMAT + 1 ) )

C.........  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )

C.........  Set up reference county information
        CALL SETREFCNTY

C.........  Build indicies for pollutant/species groups
        CALL BLDMRGIDX

C.........  Intialize state/county summed emissions to zero
        CALL INITSTCY

C.........  Read source group cross-reference file and assign sources to groups
        IF ( SRCGRPFLAG ) THEN
            CALL RDSRCGRPS( SGDEV, .TRUE., .NOT. MOPTIMIZE )

            IF ( MOPTIMIZE ) THEN
                ALLOCATE( TMPEMGGRD( NSGOUTPUT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'TMPEMGGRD', PROGNAME )
            END IF
        END IF

C.........  Open NetCDF output files, open ASCII report files, and write headers
        CALL OPENMRGOUT

C.........  Read control factor data
        IF( CFFLAG ) CALL RDCFPRO( CFDEV )

C.........  Read speciation matrices for current variable
        IF( SMATCHK ) THEN 
            DO V = 1, NSMATV

C.............  Extract name of variable
            VBUF = TSVDESC( V )

                CALL RDSMAT( MSNAME_L, VBUF, MSMATX_L( 1,V ) )
                CALL RDSMAT( MSNAME_S, VBUF, MSMATX_S( 1,V ) )

C.................  Switch SPC matrix (mole/mass) based on MRG_GRDOUT_UNIT, MRG_TOTOUT_UNIT
                IF( INDEX( GRDENV, 'mole' ) < 1 ) THEN
                    MSMATX_L( 1,V ) = MSMATX_S( 1,V )
                END IF

             END DO     
        END IF

C.........  Write out message with list of species
        CALL POLMESG( NMSPC, EMNAM )

C.........  Year of SDATE  
        WRITE( TDATE, '(I7)' ) SDATE
        YEAR = TDATE( 1:4 ) 

C.........  Loop over reference counties
        DO I = 1, NREFC

C.........  Skip if ref county is out of modeling domain for speed up
            IF( NREFSRCS( I ) < 1 ) CYCLE

C.............  Determine fuel month for current time step and reference county
            N = FINDCFIRST( MCREFIDX( I,1 ), NREFF, FMREFSORT( :,1 ) )
            M = FINDC( MCREFIDX( I,1 ), NFUELC, FMREFLIST( :,1 ) )

C.............  Determine month
            IF( N .LT. 0 .OR. M .LT. 0 ) THEN
                WRITE( MESG, 94010 )'ERROR: No fuel month data for ' //
     &            'reference county ' //  MCREFIDX( I,1 )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
                
C.............  Initializations before main time loop 
            JDATE   = SDATE
            JTIME   = STIME
            LDATE   = 0
            LMONTH  = 0
            DAY     = 1
            FUELMONTH = 0
            LFUELMONTH  = 0

C.............  Loop through output time steps
            DO T = 1, NSTEPS

                IF ( MOPTIMIZE ) THEN
                    EMGRD = 0.  ! array
                    IF ( SRCGRPFLAG ) THEN
                        EMGGRDSPC = 0.  ! array
                    END IF
                END IF
                TMPEMGRD = 0.  ! array

C.................  Determine month
                CALL DAYMON( JDATE, MONTH, DAYMONTH )

                IF ( MONTH .NE. LMONTH ) THEN
                    FUELMONTH = 0
                    NMON = STR2INT( FMREFLIST( M,2 ) )
                    DO J = N, N + NMON
                        FMON = STR2INT( FMREFSORT(J,3) )
                        IF( FMON == MONTH ) THEN
                            FUELMONTH = STR2INT( FMREFSORT( J,2 ) )
                            EXIT
                        END IF
                    END DO
                
                    IF( FUELMONTH == 0 ) THEN
                        WRITE( MESG, 94010 )'ERROR: Could not determine ' //
     &                    'fuel month for reference county ' //MCREFIDX( I,1 )//
     &                    ' and episode month', MONTH
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

C.....................  Read emission factors for reference county and month
C.....................  Update emission factors when there is a change of the fuel month
C                       at the last hour of the last day in fuel month to process the first day of next fuel month.

                    IF ( LFUELMONTH .NE. FUELMONTH ) THEN
                        WRITE( MESG,94010 ) 'Processing MOVES lookup ' //
     &                      'tables for reference county '// MCREFIDX( I,1 )//  
     &                      ' of fuel month:', FUELMONTH 
                        CALL M3MSG2( MESG )
 
                        IF( RPDFLAG ) THEN
                            CALL RDRPDEMFACS( I, FUELMONTH )
                        END IF

                        IF( RPHFLAG ) THEN
                            CALL RDRPHEMFACS( I, FUELMONTH )
                        END IF
             
                        IF( RPVFLAG ) THEN
                            CALL RDRPVEMFACS( I, FUELMONTH )
                        END IF
                
                        IF( RPPFLAG ) THEN
                            CALL RDRPPEMFACS( I, FUELMONTH )
                        END IF
                        LFUELMONTH = FUELMONTH
                    END IF
                    LMONTH = MONTH
                END IF

C.................  Write out message for new day.
C.................  Set start hour of day for all sources
                IF( JDATE .NE. LDATE ) THEN
                    CALL SETSRCDY( NMSRC, JDATE, TZONES, LDAYSAV, .TRUE.,
     &                             DAYBEGT, DAYENDT )
                END IF

C.................  Write out files that are being used for by-day treatment
                IF( RPDFLAG .OR. RPHFLAG .AND. DAY .NE. PDAY ) THEN
                    IF( MFLAG_BD ) THEN
                        MESG = '   with MTMP file ' // MTNAME( DAY )
                        CALL M3MSG2( MESG )
                    END IF
                    PDAY = DAY
                END IF

C.................  For new hour...
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

C.................  In RPH mode, read HOTELLING for current hour
                IF( RPHFLAG ) THEN
                    IF( .NOT. READSET( MTNAME( DAY ), 'HOTELLING', 1, ALLFILES,
     &                                 MJDATE, JTIME, HOTEL ) ) THEN
                        MESG = 'Could not read HOTELLING from MTMP file'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF
                END IF

C.................  In RPD and RPV modes, read temperatures for current hour
                IF( RPDFLAG .OR. RPHFLAG .OR. RPVFLAG ) THEN
                    IF( .NOT. READ3( METNAME, TVARNAME, 1, 
     &                               JDATE, JTIME, TEMPG ) ) THEN
                        MESG = 'Could not read ' // TRIM( TVARNAME ) //
     &                         ' from ' // TRIM( METNAME )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

                ELSE   ! RPP MODE
                    IF( .NOT. READ3( METNAME, 'MAXTEMP', 1, 
     &                               JDATE, 0, MAXTEMP ) ) THEN

C.........................  Look for a backup max temp from the first day of processing year
                        WRITE( TDATE, '(I7)' ) JDATE
                        PDATE = STR2INT( YEAR // TDATE( 5:7 ) )
                        IF( .NOT. READ3( METNAME, 'MAXTEMP', 1,
     &                                   PDATE, 0, MAXTEMP ) ) THEN
                            MESG = 'Could not read MAXTEMP'//
     &                             ' from ' // TRIM( METNAME )
                            CALL M3EXIT( PROGNAME, JDATE, 0, MESG, 2 )
                        END IF

                    END IF

                    IF( .NOT. READ3( METNAME, 'MINTEMP', 1, 
     &                               JDATE, 0, MINTEMP ) ) THEN

C.........................  Look for a backup min temp from the first day of processing year
                        IF( .NOT. READ3( METNAME, 'MINTEMP', 1,
     &                                   PDATE, 0, MINTEMP ) ) THEN
                            MESG = 'Could not read MINTEMP' //
     &                             ' from ' // TRIM( METNAME )
                            CALL M3EXIT( PROGNAME, JDATE, 0, MESG, 2 )
                        END IF
                    END IF

                END IF

C.................  Loop over sources in reference county
                DO S = 1, NREFSRCS( I )
                
                    SRC = REFSRCS( I,S )
                
C.....................  Determine hour index based on source's local time
                    HOURIDX = ( JTIME - DAYBEGT( SRC ) ) / 10000
                    IF( HOURIDX < 0 ) THEN
                        HOURIDX = HOURIDX + 24
                    END IF
                    HOURIDX = HOURIDX + 1  ! array index is 1 to 24

C.....................  Determine weekday index (Monday is 1)
                    DDATE = JDATE
                    DTIME = JTIME
                    CALL NEXTIME( DDATE, DTIME, -DAYBEGT( SRC ) ) 
                    DAY = WKDAY( DDATE )
                    IF( DAY .GT. 5 ) THEN
                        DAYIDX = 1   ! weekend (6-7)
                    ELSE
                        DAYIDX = 2   ! weekday (1-5)
                    END IF

C.....................  Determine SCC index for source
                    SCCIDX = MISCC( SRC )

C.....................  Determine speed bins for source
                    IF( RPDFLAG ) THEN
                        SPEEDVAL = BADVAL3
                        IF( SPDFLAG ) THEN
                            SPEEDVAL = SPDPRO( MCFIP( SRC ), SCCIDX, DAYIDX, HOURIDX )
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

C.....................  Loop over grid cells for this source
                    DO NG = 1, NSRCCELLS( SRC )

                        CELL = SRCCELLS( NG, SRC )
                        GFRAC = SRCCELLFRACS( NG, SRC )

                        IF( RPDFLAG ) EMFAC = VMT( SRC ) * GFRAC
                        IF( RPHFLAG ) EMFAC = HOTEL( SRC ) * GFRAC
                        IF( RPVFLAG .OR. RPPFLAG ) THEN
                            EMFAC = VPOP( SRC ) * GFRAC
                        END IF

C.............................  Determine temperature indexes for cell
                        IF( RPDFLAG .OR. RPHFLAG .OR. RPVFLAG ) THEN

                            TEMPVAL = ( TEMPG( CELL ) - CTOK ) * CTOF + 32.

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
                            IF( IDX1 == 0 ) THEN
                                IF( TEMPVAL >= (EMTEMPS( 1 ) - TEMPBIN) ) THEN
                                    IDX1 =1
                                    IDX2 =1
                                    IF( NWARN < MXWARN ) THEN
                                    WRITE( MESG, 94040 ) 'Grid cell ' //
     &                               'temperature', TEMPVAL, 'out of lowest limit ' //
     &                               'of emission factor data',  EMTEMPS( 1 ),
     &                               CRLF() // BLANK10 //'Lowest profile emission factor ' //
     &                               'is used based on temperature buffer bin', TEMPBIN 
                                    CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
                                    NWARN = NWARN + 1
                                    END IF

                                ELSE
                                    WRITE( MESG, 94040 ) 'ERROR: Grid cell ' //
     &                             'temperature', TEMPVAL, 'out of range ' //
     &                             'of minimum emission factor data', EMTEMPS( 1 )
                                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                               END IF
                            END IF

                            IF( IDX2 > NEMTEMPS ) THEN
                                IF( TEMPVAL <= (EMTEMPS( NEMTEMPS ) + TEMPBIN) ) THEN
                                    IDX1 = NEMTEMPS
                                    IDX2 = NEMTEMPS
                                    IF( NWARN < MXWARN ) THEN
                                    WRITE( MESG, 94040 ) 'Grid cell ' //
     &                               'temperature', TEMPVAL, 'out of highest limit ' //
     &                               'of emission factor data',  EMTEMPS( NEMTEMPS ),
     &                               CRLF()//BLANK10 //'Highest profile emission factor '//
     &                               'is used based on temperature buffer bin', TEMPBIN 
                                    CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
                                    NWARN = NWARN + 1
                                    ENDIF
                                ELSE
                                    WRITE( MESG, 94040 ) 'ERROR: Grid cell ' //
     &                             'temperature', TEMPVAL, 'out of range ' //
     &                             'of maximum emission factor data', EMTEMPS( NEMTEMPS )
                                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                               END IF

                            END IF
                            
                            IF( IDX1 .NE. IDX2 ) THEN
                                TEMPFAC = ( TEMPVAL - EMTEMPS( IDX1 ) ) /
     &                                    ( EMTEMPS( IDX2 ) - EMTEMPS( IDX1 ) )
                            ELSE
                                TEMPFAC = 1.
                            END IF

                        ELSE
C.........................  Determine profiles for current inventory county for RPP mode
C                           There will be 4 profiles used in total:
C                             UU - both min and max profile temps are under county temps
C                             UO - min profile temp is under county min, max profile temp is over county max
C                             OU - min profile temp is over county min, max profile temp is under county max
C                             OO - both min and max profile temps are over county temps
                            MINTVAL = MINTEMP( CELL )
                            MAXTVAL = MAXTEMP( CELL )

C.............................  MCFIP(SRC) and maxmum values within the index
                            IF ( (MINTVAL .LT. EMTEMPS( EMTEMPIDX( 1 ) ) )  
     &                          .AND. (MINTVAL .GE. (EMTEMPS( EMTEMPIDX( 1 ) ) - TEMPBIN )) ) THEN
                                IF( NWARN < MXWARN ) THEN
                                    WRITE( MESG, 94040 ) 'Lowest profile ' //
     &                               'minimum temperature ', EMTEMPS( EMTEMPIDX( 1 )) ,
     &                               ' is higher than county minimum temperature:', MINTVAL,
     &                               CRLF() // BLANK10 // 'minimum value will be set to ',
     &                               EMTEMPS( EMTEMPIDX( 1 )), 'based on temperature buffer bin ', TEMPBIN
                                    CALL M3WARN(PROGNAME, JDATE, JTIME, MESG )
                                    NWARN = NWARN + 1
                                END IF 
                                MINTVAL = EMTEMPS( EMTEMPIDX( 1 ) )
                            END IF
                            IF ( (MAXTVAL .GT. EMXTEMPS( EMTEMPIDX( NEMTEMPS ) ))
     &                          .AND. ( MAXTVAL .LE. (EMXTEMPS( EMTEMPIDX( NEMTEMPS ) ) + TEMPBIN) ) ) THEN
                                IF( NWARN < MXWARN ) THEN
                                     WRITE( MESG, 94040 ) 'Highest profile ' // 
     &                                  'maximum temperature ', EMXTEMPS( EMTEMPIDX( NEMTEMPS )) ,
     &                                  ' is lower than county maximum temperature:',  MAXTVAL,
     &                                  CRLF() // BLANK10 // 'maximum value will be set to ',
     &                                  EMTEMPS( EMTEMPIDX( NEMTEMPS )), 'based on temperature buffer bin ', TEMPBIN
                                     CALL M3WARN(PROGNAME, JDATE, JTIME, MESG )
                                     NWARN = NWARN + 1
                                END IF
                                MAXTVAL = EMXTEMPS( EMTEMPIDX( NEMTEMPS ))
                            END IF
                            
C.............................  Check that min and max county temps were found
                            IF( MINTVAL .LT. AMISS3 .OR. MAXTVAL .LT. AMISS3 ) THEN
                                WRITE( MESG, 94010 ) 'Could not find minimum ' //
     &                            'or maximum temperatures for county ' // CIFIP( SRC ) 
     &                            // ' and episode month', MONTH
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                            END IF

C.............................  Find indexes of bounding minimum temperatures
                            USTART = 0
                            OSTART = 0
                            PDIFF = 999
                            NO_INTRPLT = .FALSE.     ! no interpolation flag
                            DO K = 1, NEMTEMPS

                                TDIFF = MINTVAL - EMTEMPS( EMTEMPIDX( K ) )

C.................................  Once profile min temp is greater than county temp, this loop is done                            
                                IF( TDIFF .LT. 0 ) THEN
                                    OSTART = K
                                    OMIN = EMTEMPS( EMTEMPIDX( K ) )
                                    EXIT
                                END IF

C.................................  If current profile min temp is closer to county temp, store index                            
                                IF( TDIFF .LT. PDIFF ) THEN
                                    USTART = K
                                    UMIN = EMTEMPS( EMTEMPIDX( K ) )
                                END IF
                        
                                PDIFF = TDIFF
                            END DO

C.............................  Check that appropriate minimum temperatures were found
                            IF( USTART == 0 ) THEN
                                WRITE( MESG, 94040 ) 'ERROR: Lowest profile ' //
     &                            'minimum temperature', 
     &                            EMTEMPS( EMTEMPIDX( 1 ) ),
     &                            'is higher than county minimum temperature',
     &                            MINTVAL 
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                            END IF
                        
                            IF( OSTART == 0 ) THEN
                                IF ( MINTVAL .EQ. EMTEMPS( EMTEMPIDX( NEMTEMPS ) ) ) THEN
                                    OSTART = NEMTEMPS 
                                ELSE
                                    WRITE( MESG, 94040 ) 'ERROR: Highest ' //
     &                                 'profile minimum temperature', 
     &                                 EMTEMPS( EMTEMPIDX( NEMTEMPS ) ),
     &                                 'is lower than county minimum temperature',
     &                                 MINTVAL 
                                   CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                                END IF
                            END IF

C.............................  Find indexes of bounding maximum temperatures
                            UUIDX = 0
                            UOIDX = 0
                            UMAX = EMXTEMPS( EMTEMPIDX( USTART ) )
                            DO K = USTART, NEMTEMPS

C.................................  Check that profile minimum temperature hasn't changed
                                IF( EMTEMPS( EMTEMPIDX( K ) ) .NE. UMIN ) THEN
                                    EXIT
                                END IF

                                OMAX = EMXTEMPS( EMTEMPIDX( K ) )
                                IF( OMAX .EQ. MAXTVAL ) THEN
                                    UOIDX = EMTEMPIDX( K )
                                    UUIDX = EMTEMPIDX( K )
                                    EXIT
                                END IF

                                IF( OMAX > MAXTVAL ) THEN
                                    UOIDX = EMTEMPIDX( K )
                                    IF( K > USTART ) THEN
                                        UUIDX = EMTEMPIDX( K - 1 )
                                    END IF
                                    EXIT
                                END IF
                            END DO

C.............................  Check that appropriate maximum temperatures were found
                            IF( UUIDX .EQ. 0 .AND. UOIDX .NE. 0 ) THEN
                                WRITE( MESG, 94040 ) 'ERROR: Lowest profile of ' //
     &                            'max temperature', UMAX,'- min temperature',UMIN,
     &                            ' is higher than county maximum temperature:',
     &                            MAXTVAL
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                            END IF

                            IF ( UUIDX .EQ. 0 .AND. UOIDX .EQ. 0 ) THEN
                                WRITE( MESG, 94040 ) 'ERROR: Highest profile of ' //
     &                            'max temperature', OMAX,'- min temperature',UMIN, 
     &                            ' is lower than county maximum temperature',
     &                            MAXTVAL
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                            END IF

                            OUIDX = 0
                            OOIDX = 0
                            UMAX = EMXTEMPS( EMTEMPIDX( OSTART ) )
                            DO K = OSTART, NEMTEMPS

C.................................  Check that profile minimum temperature hasn't changed
                                IF( EMTEMPS( EMTEMPIDX( K ) ) .NE. OMIN ) THEN
                                    EXIT
                                END IF

                                OMAX = EMXTEMPS( EMTEMPIDX( K ) )
C...............................  If MAXTVAL equal to maximum value
                                IF( OMAX .EQ. MAXTVAL ) THEN
                                    OOIDX = EMTEMPIDX( K )
                                    OUIDX = EMTEMPIDX( K )
                                    EXIT
                                END IF

                                IF( OMAX > MAXTVAL ) THEN
                                    OOIDX = EMTEMPIDX( K )
                                    IF( K > OSTART ) THEN
                                        OUIDX = EMTEMPIDX( K - 1 )
                                    END IF
                                    EXIT
                                END IF
                            END DO

C.............................  Determine whether min/max temp within one temp bin
                            IF( OMAX .EQ. OMIN ) NO_INTRPLT = .TRUE.

C.............................  Check that appropriate maximum temperatures were found
                            IF( .NOT. NO_INTRPLT ) THEN
                              IF( OUIDX .EQ. 0 .AND. OOIDX .NE. 0 ) THEN
                                WRITE( MESG, 94040 ) 'ERROR: Lowest profile of ' //
     &                            'max temperature', UMAX,'- min temperature',OMIN,
     &                            ' is higher than county maximum temperature:',
     &                            MAXTVAL
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                              END IF

                              IF( OUIDX .EQ. 0 .AND. OOIDX .EQ. 0 ) THEN
                                WRITE( MESG, 94040 ) 'ERROR: Highest profile of ' //
     &                            'max temperature', OMAX,'- min temperature',OMIN,
     &                            ' is lower than county maximum temperature',
     &                            MAXTVAL
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                              END IF

C...............................  Check that maximum temperatures of profiles match
                              IF( EMXTEMPS( UUIDX ) .NE. EMXTEMPS( OUIDX ) .OR.
     &                          EMXTEMPS( UOIDX ) .NE. EMXTEMPS( OOIDX ) ) THEN
                                MESG = 'ERROR: Inconsistent temperature ' //
     &                            'profiles.'
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                              END IF

                              IF ( EMTEMPS( OUIDX ) .EQ. EMTEMPS( UUIDX ) ) THEN
                                  MAXFAC = 1.
                              ELSE
                                  MINFAC = ( MINTVAL - EMTEMPS( UUIDX ) ) /
     &                               ( EMTEMPS( OUIDX ) - EMTEMPS( UUIDX ) )
                              END IF
                                                 
                              IF ( EMXTEMPS( OOIDX ) .EQ. EMXTEMPS( OUIDX ) ) THEN
                                  MAXFAC = 1.
                              ELSE
                                  MAXFAC = ( MAXTVAL - EMXTEMPS( OUIDX ) ) /
     &                               ( EMXTEMPS( OOIDX ) - EMXTEMPS( OUIDX ) )
                              END IF

                          END IF

                        END IF

C.........................  Loop through pollutant-species combos
                        NPOLSPC = NIPPA + NMSPC
                        IF( SMATCHK ) NPOLSPC = NSMATV

                        DO V = 1, NPOLSPC

C.............................  Lookup poll/species index from MOVES lookup EF 
                            SIIDX = 0
                            SPIDX = 0
                            IF( SMATCHK ) THEN
                                SIIDX = SIINDEX( V,1 )
                                SPIDX = SPINDEX( V,1 )
                                POLIDX = EMPOLIDX( SIIDX )
                            ELSE
                                IF( V <= NIPPA ) THEN
                                    SIIDX = V
                                ELSE
                                    SPIDX = V - NIPPA
                                END IF
                                POLIDX = EMPOLIDX( V )
                            END IF

                            IF( CFFLAG ) THEN
                                CFFAC = CFPRO(MCFIP(SRC), SCCIDX, V, MONTH )
                                IF( SMATCHK ) THEN
                                    CFFAC = CFPRO(MCFIP(SRC), SCCIDX, SIIDX, MONTH )
                                END IF
                            END IF

C.............................  Check if emission factors exist for this process/pollutant
                            IF( POLIDX .EQ. 0 ) CYCLE

C.............................  Calculate interpolated emission factor if process/pollutant has changed
                            IF( RPDFLAG ) THEN
                                IF( BIN1 .NE. BIN2 ) THEN
                                    EFVAL1 = RPDEMFACS( SCCIDX, BIN1, IDX1, POLIDX )
                                    EFVAL2 = RPDEMFACS( SCCIDX, BIN2, IDX1, POLIDX )
                                    EFVALA = SPDFAC * (EFVAL2 - EFVAL1) + EFVAL1
        
                                    EFVAL1 = RPDEMFACS( SCCIDX, BIN1, IDX2, POLIDX )
                                    EFVAL2 = RPDEMFACS( SCCIDX, BIN2, IDX2, POLIDX )
                                    EFVALB = SPDFAC * (EFVAL2 - EFVAL1) + EFVAL1
                                ELSE
                                    EFVALA = RPDEMFACS( SCCIDX, BIN1, IDX1, POLIDX )
                                    EFVALB = RPDEMFACS( SCCIDX, BIN1, IDX2, POLIDX )
                                END IF
                                EFVAL = TEMPFAC * (EFVALB - EFVALA) + EFVALA
                            END IF

                            IF( RPHFLAG ) THEN
                                EFVALA = RPHEMFACS( SCCIDX, IDX1, POLIDX )
                                EFVALB = RPHEMFACS( SCCIDX, IDX2, POLIDX )
                                EFVAL = TEMPFAC * (EFVALB - EFVALA) + EFVALA
                            END IF
                                
                            IF( RPVFLAG ) THEN
                                EFVALA = RPVEMFACS( DAYIDX, SCCIDX, HOURIDX, IDX1, POLIDX )
                                EFVALB = RPVEMFACS( DAYIDX, SCCIDX, HOURIDX, IDX2, POLIDX )

                                EFVAL = TEMPFAC * (EFVALB - EFVALA) + EFVALA
                            END IF
                                
                            IF( RPPFLAG ) THEN
                                IF( NO_INTRPLT ) THEN
                                    EFVAL = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, UOIDX, POLIDX )
                                ELSE
                                    EFVAL1 = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, UUIDX, POLIDX )
                                    EFVAL2 = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, UOIDX, POLIDX )
                                    EFVALA = MAXFAC * (EFVAL2 - EFVAL1) + EFVAL1
                                        
                                    EFVAL1 = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, OUIDX, POLIDX )
                                    EFVAL2 = RPPEMFACS( DAYIDX, SCCIDX, HOURIDX, OOIDX, POLIDX )
                                    EFVALB = MAXFAC * (EFVAL2 - EFVAL1) + EFVAL1
    
                                    EFVAL = MINFAC * (EFVALB - EFVALA) + EFVALA
                                    
                                END IF
                            END IF

                            IF( CFFLAG ) EFVAL = EFVAL * CFFAC 

C.............................  Calculate gridded, hourly emissions (g/hr/cell)
                            EMVAL = EFVAL * EMFAC

C.............................  Set units conversion factor
                            IF( SPIDX > 0 ) THEN
                                F1 = GRDFAC( SPIDX )
                                
                                IF( SMATCHK ) THEN
                                    EMVALSPC = EMVAL * MSMATX_L( SRC,V ) * F1
                                ELSE
                                    EMVALSPC = EMVAL * F1
                                ENDIF

                                IF( MOPTIMIZE ) THEN
                                    EMGRD( CELL,SPIDX ) = 
     &                                  EMGRD( CELL,SPIDX ) + EMVALSPC
     
                                    IF( SRCGRPFLAG ) THEN
                                        GIDX = ISRCGRP( SRC )
                                        EMGGRDSPC( CELL,GIDX,SPIDX ) =
     &                                     EMGGRDSPC( CELL,GIDX,SPIDX ) + EMVALSPC
                                    END IF

C...................................  If not use memory optimize
                                ELSE                      
                                    TEMGRD( CELL,SPIDX,T ) =
     &                                TEMGRD( CELL,SPIDX,T ) + EMVALSPC
     
                                    IF ( SRCGRPFLAG ) THEN
                                        GIDX = ISRCGRP( SRC )
                                        EMGGRDSPCT( CELL,GIDX,SPIDX,T ) =
     &                                    EMGGRDSPCT( CELL,GIDX,SPIDX,T ) + EMVALSPC
                                    END IF

                                END IF

                            END IF

C.............................  Store Temporal intermediate hourly emissions
                            IF( MTMP_OUT ) THEN
                                IF( SMATCHK ) THEN
                                    IF( EANAMREP( V ) ) THEN
                                        MTMP_INVT( SRC,SIIDX,T ) =
     &                                       MTMP_INVT( SRC,SIIDX,T ) +
     &                                       EMVAL * GM2TON       ! g/hr-cell * ton/g = ton/hr-cell
                                    END IF
                                ELSE
                                    MTMP_INVT( SRC,SIIDX,T ) =
     &                                   MTMP_INVT( SRC,SIIDX,T ) +
     &                                   EMVAL * GM2TON       ! g/hr-cell * ton/g = ton/hr-cell
                                END IF
                            END IF

C.............................  Add this cell's emissions to source totals in unit of moles/hr (instead of moles/s)
C                               not by applying 1/3600 factor (hr to sec)
                            IF( LREPANY ) THEN
                                IF( .NOT. (T == NSTEPS .AND. JTIME == 0) ) THEN
                                    IF( SMATCHK ) THEN
                                        EMVALSPC = EMVAL * MSMATX_S( SRC,V ) * TOTFAC( SPIDX )
                                        MEBSUM( SRC,SPIDX ) =  MEBSUM( SRC,SPIDX ) + 
     &                                                         EMVALSPC * GM2TON
                                        IF( EANAMREP( V ) ) THEN
                                            EMVALSPC = EMVAL * TOTFAC( NMSPC+SIIDX )
                                            MEBSUM( SRC,NMSPC+SIIDX ) =
     &                                          MEBSUM( SRC,NMSPC+SIIDX ) + EMVALSPC
                                        END IF
                                    ELSE
                                        IF( SPIDX > 0 ) THEN
                                            EMVALSPC = EMVAL
                                            MEBSUM( SRC,SPIDX ) =  MEBSUM( SRC,SPIDX ) + 
     &                                                             EMVALSPC
                                        ELSE   ! sum of inv emission in unit of tons/hr 
                                            EMVALSPC = EMVAL * TOTFAC( NMSPC+SIIDX )
                                            MEBSUM( SRC,NMSPC+SIIDX ) =
     &                                          MEBSUM( SRC,NMSPC+SIIDX ) + EMVALSPC
                                        END IF
                                    END IF
                                    RDATE = JDATE     ! last reporting date
                                    RTIME = JTIME     ! last reporting hour
                                END IF
                            END IF

                        END DO    ! end loop over pollutant-species combos

                    END DO    ! end loop over grid cells for source
                
                END DO    ! end loop over sources in inv. county

C.................  Read out old data if not first county
                IF ( MOPTIMIZE ) THEN
                    DO V = 1, NMSPC 
                        CSPC = EMNAM( V )

C.........................  sum old county data with new county
                        IF ( I > 1 ) THEN
                            IF(.NOT. READSET( MONAME, CSPC, 1,  ALLFILES,
     &                                JDATE, JTIME, TMPEMGRD( 1,V ) ) )THEN
                                 MESG = 'Could not read "' // CSPC // '" ' //
     &                             'from file "' // MONAME // '"'
                                 CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                            END IF

                            EMGRD( :,V ) = EMGRD( :,V ) + TMPEMGRD( :,V )

                        END IF

                        IF( LGRDOUT ) THEN
                            IF( .NOT. WRITESET( MONAME, CSPC, ALLFILES,
     &                              JDATE, JTIME, EMGRD( 1,V ) ) ) THEN
                                MESG = 'Could not write "' // CSPC // '" ' //
     &                           'to file "' // MONAME // '"'
                                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                            END IF
                        END IF
                        
                        IF( SRCGRPFLAG ) THEN
                            EMGGRD( :,: ) = EMGGRDSPC( :,:,V )
                            IF( I > 1 ) THEN
                                TMPEMGGRD = 0.  ! array
                                IF( .NOT. READ3( SGINLNNAME, CSPC, 1, 
     &                                           JDATE, JTIME, TMPEMGGRD ) ) THEN
                                    MESG = 'Could not read "' // CSPC // '" ' //
     &                                     'from file "' // SGINLNNAME // '"'
                                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                                END IF
                                
                                CALL WRSRCGRPS( CSPC, JDATE, JTIME, .TRUE., TMPEMGGRD )
                            ELSE
                                CALL WRSRCGRPS( CSPC, JDATE, JTIME, .FALSE., 0 )
                            END IF
                        END IF
                    END DO
                END IF   ! end memory optimize
           
                LDATE = JDATE
                CALL NEXTIME( JDATE, JTIME, TSTEP )     !  update model clock

            END DO   ! End loop on time steps

        END DO   ! end loop over inventory counties

C.........  Write state, county, and SCC emissions (all that apply)
C.........  The subroutine will only write for certain hours and
C           will reinitialize the totals after output
        IF( LREPANY ) CALL WRMRGREP( RDATE, RTIME )

C.........  Output optional hourly emissions for inventory pollutants for Temporal program
        IF( MTMP_OUT ) THEN
            JDATE  = SDATE
            JTIME  = STIME
            DO T = 1, NSTEPS

                DO V = 1, NIPPA 

                    CPOL = EANAM( V )  ! Pollutant

                    IF( .NOT. WRITE3( MTMPNAME, CPOL, JDATE, JTIME,
     &                                MTMP_INVT( :,V,T ) )         ) THEN
                        MESG = 'Could not write '//  TRIM( CPOL ) //
     &                         ' from ' // TRIM( MTMPNAME )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

                END DO

                CALL NEXTIME( JDATE, JTIME, TSTEP )     !  update model clock

            END DO

        END IF

C.........  Output gridded houlry emissions 
        IF( .NOT. MOPTIMIZE ) THEN
            JDATE  = SDATE
            JTIME  = STIME
            DO T = 1, NSTEPS    

                DO V = 1, NMSPC

                    CSPC = EMNAM( V )

                    TMPEMGRD( :,V ) = TEMGRD( :,V,T )

C.....................  Write out gridded data
                    IF( LGRDOUT ) THEN
                        IF( .NOT. WRITESET( MONAME, CSPC, ALLFILES,
     &                            JDATE, JTIME, TMPEMGRD( 1,V ) ) ) THEN
                            MESG = 'Could not write "' // CSPC //'" '//
     &                      'to file "' // MONAME // '"'
                            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                        END IF
                    END IF
                    
                    IF( SRCGRPFLAG ) THEN
                        EMGGRD( :,: ) = EMGGRDSPCT( :,:,V,T )
                        CALL WRSRCGRPS( CSPC, JDATE, JTIME, .FALSE., 0 )
                    END IF

                END DO

                CALL NEXTIME( JDATE, JTIME, TSTEP )     !  update model clock

            END DO   ! End loop on time steps

        END IF

C.........  Close output file
        IF( .NOT. CLOSESET( MONAME ) ) THEN
            MESG = 'Could not close file:"'// MONAME // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        DEALLOCATE( DAYBEGT, DAYENDT, LDAYSAV )

        IF ( MOPTIMIZE ) THEN
            DEALLOCATE( EMGRD, TMPEMGRD, MGMATX, TMPEMGGRD )
        ELSE
            DEALLOCATE( TEMGRD, TMPEMGRD, MGMATX )
        END IF

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

