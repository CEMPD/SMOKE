
        SUBROUTINE GENPDOUT( FDEV, TZONE, SDATE, STIME, NSTEPS, INSTEP, 
     &                       OUTSTEP, NVAR, NVSP, MXPDSRC, TYPNAM, 
     &                       FNAME, EAIDX, SPIDX )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine reads and writes the day-specific or hour-specific
C      emissions.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
C
C*************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  EXTERNAL FUNCTIONS
        CHARACTER*2  CRLF
        LOGICAL      ENVYN


C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV      ! file unit no.
        INTEGER     , INTENT (IN) :: TZONE     ! output time zone
        INTEGER     , INTENT (IN) :: SDATE     ! Julian starting date in TZONE
        INTEGER     , INTENT (IN) :: STIME     ! start time of data in TZONE
        INTEGER     , INTENT (IN) :: NSTEPS    ! no. time steps
        INTEGER     , INTENT (IN) :: INSTEP    ! expected data time step HHMMSS
        INTEGER     , INTENT (IN) :: OUTSTEP   ! output time step HHMMSS
        INTEGER     , INTENT (IN) :: NVAR      ! no. period-specific variables
        INTEGER     , INTENT (IN) :: NVSP      ! no. period-spec special vars
        INTEGER     , INTENT (IN) :: MXPDSRC   ! maximum period-specific sources
        CHARACTER(*), INTENT (IN) :: TYPNAM    ! 'day' or 'hour'
        CHARACTER(*), INTENT (IN) :: FNAME     ! logical file name
        INTEGER     , INTENT (IN) :: EAIDX( NIPPA ) ! index to EANAM
        INTEGER     , INTENT (IN) :: SPIDX( MXSPDAT ) ! index to SPDATNAM

C.........  Local allocatable arrays
        LOGICAL, ALLOCATABLE :: EASTAT( : )    ! true: act/pol present in data

C.........  Local arrays
        LOGICAL         SPSTAT( MXSPDAT )     ! true: special data variable used

C...........   Other local variables
        INTEGER          S, T

        INTEGER          IOS                  ! i/o status
        INTEGER          JDATE                ! tmp Julian date
        INTEGER          JTIME                ! tmp time HHMMSS
        INTEGER          NPDSRC               ! number of day/hour-spec sources
        INTEGER          PDEMDIM              ! dim for PDEMOUT

        LOGICAL       :: DFLAG    = .FALSE.  ! true: day-specific
        LOGICAL       :: EFLAG    = .FALSE.  ! true: error found
        LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called
        LOGICAL, SAVE :: OFLAG    = .FALSE.  ! true: PFLAG & hourly
        LOGICAL, SAVE :: PFLAG    = .FALSE.  ! true: create hourly profiles
        LOGICAL, SAVE :: SFLAG    = .FALSE.  ! true: create daily totals

        CHARACTER*300 :: MESG = ' '          ! message buffer

        CHARACTER(LEN=NAMLEN3) ONAME         ! output file name
 
        CHARACTER*16 :: PROGNAME = 'GENPDOUT' !  program name

C***********************************************************************
C   begin body of program GENPDOUT

C.........  For the first time the routine is called...
        IF( FIRSTIME ) THEN

C.............  Get environment variable for creating daily data from the hourly
            MESG = 'Use daily totals only from hourly data file'
            SFLAG = ENVYN( 'HOURLY_TO_DAILY', MESG, .FALSE., IOS )

C.............  Get environment variable for creating hourly profiles from the
C               hourly data
            MESG = 'Create hourly profiles from hourly data'
            PFLAG = ENVYN( 'HOURLY_TO_PROFILE', MESG, .FALSE., IOS )

            IF( SFLAG .AND. PFLAG ) THEN
                MESG = 'WARNING: Ignoring HOURLY_TO_PROFILE "Y" ' //
     &                 'value because HOURLY_TO_DAILY is set to "Y"'
                CALL M3MSG2( MESG )
                PFLAG = .FALSE.

            END IF

            FIRSTIME = .FALSE.

        END IF

C.........  Perform case-specific settings
        OFLAG = .FALSE.
        SELECT CASE( TYPNAM )
        CASE( 'day' ) 
            DFLAG = .TRUE.

        CASE( 'hour' )
            IF( PFLAG ) OFLAG = .TRUE.
            DFLAG = .FALSE.

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: Do not know type ' // TYPNAM 
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.........  Allocate memory for logical status array for pol/act even though
C           it does not need to be set because EAIDX has already been 
C           determined.
        ALLOCATE( EASTAT( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EASTAT', PROGNAME )
        EASTAT = .FALSE.  ! array

C.........  Allocate memory for reading data
        ALLOCATE( MXPDPT( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MXPDPT', PROGNAME )
        ALLOCATE( NPDPT ( NSTEPS )        , STAT=IOS )
        CALL CHECKMEM( IOS, 'NPDPT', PROGNAME )
        ALLOCATE( CODEA ( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CODEA', PROGNAME )
        ALLOCATE( IDXSRC( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRC', PROGNAME )
        ALLOCATE( SPDIDA( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPDIDA', PROGNAME )
        ALLOCATE( EMISVA( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISVA', PROGNAME )
        ALLOCATE( DYTOTA( MXPDSRC,NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DYTOTA', PROGNAME )
        ALLOCATE( LPDSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LPDSRC', PROGNAME )

C.........  Initialize arrays
        MXPDPT = 0        ! array
        NPDPT  = 0        ! array
        CODEA  = 0        ! array
        IDXSRC = 0        ! array
        SPDIDA = 0        ! array
        EMISVA = BADVAL3  ! array
        DYTOTA = BADVAL3  ! array
        LPDSRC = .FALSE.  ! array

C.........  Message before reading the input file (list of files)
        MESG = 'Reading ' // TYPNAM // '-specific data...'
        CALL M3MSG2( MESG )

C.........  Loop through input files and actually read the data
        CALL RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, DFLAG, 
     &                 FNAME, SDATE, STIME, NSTEPS, EASTAT, SPSTAT )

C.........  Determine the actual number of day-specific or hour-specific sources
        NPDSRC = 0
        DO S = 1, NSRC
            IF( LPDSRC( S ) ) NPDSRC = NPDSRC + 1
        END DO

C.........  Make sure that that actual number of sources over all sources does
C           not exceed the maximum number of sources over all hours
        IF( NPDSRC .GT. MXPDSRC ) THEN

            WRITE( MESG,94010 ) 'INTERNAL ERROR: Actual number of ' //
     &             TYPNAM // 'sources, NPDSRC=', NPDSRC, CRLF() // 
     &             BLANK10 // 'dimensioned number, MXPDSRC =', MXPDSRC,
     &             '. Fix by ensuring all period-specific' // CRLF() //
     &             BLANK10 // 'sources in file for at the same day '//
     &             'or hour.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        ELSE IF( NPDSRC .EQ. 0 ) THEN

            MESG = 'No period-specific sources found in input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Allocate memory for daily or hourly output arrays.  Allocate 
C           memory as one block which will be separated into an integer section
C           and a real section when WRPDEMIS is called.  This permits
C           writing with a single WRITE3 statement.
        ALLOCATE( PDEMOUT( NPDSRC,NVAR+NVSP+1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PDEMOUT', PROGNAME )
        ALLOCATE( PDTOTL( NPDSRC,NVAR+NVSP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PDTOTL', PROGNAME )

C.........  Open day-specific or hour-specific output file, depending on value
C           of TYPNAM
        CALL OPENPDOUT( NPDSRC, NVAR, TZONE, SDATE, STIME, OUTSTEP, 
     &                  TYPNAM, OFLAG, EAIDX, SPSTAT, ONAME )

C.........  Loop through time steps and output emissions and other data

        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS

            CALL WRPDEMIS( JDATE, JTIME, T, NPDSRC, NVAR, NVSP, ONAME, 
     &                     OFLAG, EAIDX, SPIDX, PDEMOUT( 1,1 ), 
     &                     PDEMOUT( 1,2 ), EFLAG )

            CALL NEXTIME( JDATE, JTIME, OUTSTEP )

        END DO     !  End of loop over time steps

C.............  Abort if error found
	IF ( EFLAG ) THEN
            MESG = 'Problem with input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	END IF

C.........  Deallocate local memory
        DEALLOCATE( EASTAT )

C.........  Deallocate global memory
        DEALLOCATE( MXPDPT, NPDPT, CODEA, IDXSRC, SPDIDA, EMISVA, 
     &              DYTOTA, LPDSRC, PDEMOUT, PDTOTL )
 
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GENPDOUT
