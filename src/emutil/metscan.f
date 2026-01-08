
        PROGRAM METSCAN

C***********************************************************************
C  program body starts at line  
C
C  DESCRIPTION:
C       This program scans meteorology data for any time range within a year
C       and creates a first/last freeze date file.
C
C  PRECONDITIONS REQUIRED:
C       Postprocessed MM5 meteorology that contains temperature data
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       HRBIO, PREBMET
C
C  REVISION  HISTORY:
C       Created 10/2000 by M. Houyoux
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES:

C        INCLUDE 'PARMS3.EXT'      ! I/O API constants
C        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
C        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     ! emissions constants
        INCLUDE 'CONST3.EXT'      ! physical constants for I/O API

C...........   PARAMETERS and their descriptions:

        CHARACTER(50)  SCCSW          ! SCCS string with version number at end
 
        PARAMETER ( SCCSW   = '%W%'    )
     
C...........   EXTERNAL FUNCTIONS and their descriptions:

C       CHARACTER(2)    CRLF   
        INTEGER         GETFLINE
C       INTEGER         JULIAN
C       CHARACTER(16)   PROMPTMFILE
C       REAL            YR2DAY

C        EXTERNAL        CRLF, GETFLINE, JULIAN, PROMPTMFILE, YR2DAY
        EXTERNAL     GETFLINE

C.........  Local parameters
        REAL, PARAMETER :: FREEZEVAL = CTOK

C.........  Gridded temperature data
        REAL, ALLOCATABLE :: TMPR ( : )    !  air temperature (K)

C...........   Gridded freezing information
        INTEGER, ALLOCATABLE :: FREEZEGRP( : )  ! freeze group
        INTEGER, ALLOCATABLE :: SEASON   ( : )  ! 1=summer; 0=winter
        INTEGER, ALLOCATABLE :: FRZRANGE1( : )  ! freeze day first range value
        INTEGER, ALLOCATABLE :: FRZRANGE2( : )  ! freeze day second range value

        LOGICAL, ALLOCATABLE :: LFREEZE  ( :,: )! freeze status for all days

C...........   Logical names and unit numbers

        INTEGER         LDEV    ! unit number for log device

        CHARACTER(16)   MNAME   ! logical name for meteorology file
        CHARACTER(16)   ONAME   ! output file logical name

C...........   Other variables and their descriptions:

        INTEGER         C, D, L1, L2, T  !  indices and counters

        INTEGER         HR      ! current simulation hour

        INTEGER         IOS     ! temporay IO status
        INTEGER         DEC31   ! Julian day - December 31st
        INTEGER         EDATE   ! end Julian date (YYYYDDD)
        INTEGER         EDAY    ! end Julian day
        INTEGER         ETIME   ! end time (HHMMSS)
        INTEGER         JAN01   ! Julian day - January 1st
        INTEGER         JDATE   ! current simulation date (YYYYDDD)
        INTEGER         JUL31   ! Julian day - July 31st
        INTEGER         JTIME   ! current simulation time (HHMMSS)
        INTEGER         SDATE   ! start Julian date (YYYYDDD)
        INTEGER         STIME   ! start time (HHMMSS)
        INTEGER         NGRID   ! no. of grid cells
        INTEGER         NDAYS   ! number of days in met data
        INTEGER         NDAYYR  ! number of days in year of input file
        INTEGER         NSTEPS  ! run duration (hours) met or emis
        INTEGER         SDAY    ! start Julian day
        INTEGER         TSTEP   ! Met file time step
        INTEGER         YEAR    ! 4-digit year

        LOGICAL      :: EFLAG = .FALSE.   !  true: error found
        LOGICAL      :: FFLAG = .FALSE.   !  true: all freezing found
        LOGICAL      :: NFLAG = .TRUE.    !  true: all non-freezing found

        CHARACTER(5)    HEMISW  ! hemisphere switch (NORTH or SOUTH)
        CHARACTER(300)  MESG    ! message buffer for M3EXIT()

        CHARACTER(IOVLEN3) VARNAM   ! temperature variable name

C.......   Input met and grid variables:
        CHARACTER(16) :: PROGNAME = 'METSCAN'   !  program name

C***********************************************************************
C   begin body of program METSCAN

        LDEV = INIT3()
 
C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get environment variable settings...
C.........  Get name of temperature variable to use in analysis
        MESG = 'Variable name for temperature'
        CALL ENVSTR( 'TMPR_VAR', MESG, 'TA', VARNAM, IOS )
        CALL UPCASE( VARNAM )

        MESG = 'Northern or Southern Hemisphere switch'
        CALL ENVSTR( 'N_S_HEMI', MESG, 'NORTH', HEMISW, IOS )
        CALL UPCASE( HEMISW )

C.........  Open input meteorology file
        MESG = 'Enter logical name for the METEOROLOGY file'
        MNAME = 'MET_CRO_3D'
        MNAME = PROMPTMFILE( MESG, FSREAD3, MNAME, PROGNAME )

C.........  Get header from file and store necessary info
        IF ( .NOT. DESC3( MNAME ) ) THEN
            MESG = 'Could not get description of file ' // MNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE
            NSTEPS = MXREC3D
            SDATE  = SDATE3D
            STIME  = STIME3D
            TSTEP  = TSTEP3D
            NGRID  = NCOLS3D * NROWS3D
            YEAR   = SDATE / 1000

C.............  Set the ending date
            EDATE = SDATE
            ETIME = STIME 
            CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

C.............  Error if year is not the same
            IF( EDATE .LT. SDATE .OR. 
     &          EDATE/1000 .NE. SDATE/1000 ) THEN
                MESG = 'Cannot process multi-year data'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Determine number of days in input file (round up) and year
            NDAYS  = EDATE - SDATE + 1
            NDAYYR = INT( 1. / YR2DAY( YEAR ) )

C.............  Set day range
            SDAY = SDATE - ( SDATE / 1000 ) * 1000    ! integer math
            EDAY = EDATE - ( EDATE / 1000 ) * 1000    ! integer math

        END IF

C.........  Abort of hemisphere switch is not set
        IF( HEMISW .EQ. ' '  ) THEN

            MESG = 'Program not yet configured for multiple ' //
     &             'hemispheres. Must set N_S_HEMI.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Allocate memory
        ALLOCATE( TMPR( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPR', PROGNAME )
        ALLOCATE( FREEZEGRP( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FREEZEGRP', PROGNAME )
        ALLOCATE( FRZRANGE1( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FRZRANGE1', PROGNAME )
        ALLOCATE( FRZRANGE2( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FRZRANGE2', PROGNAME )
        ALLOCATE( SEASON( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SEASON', PROGNAME )
        ALLOCATE( LFREEZE( NDAYYR, NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LFREEZE', PROGNAME )

C.........  Initialize freeze status a never having freeze cells in the whole 
C           year
        LFREEZE = .FALSE.   ! array

C.........  Initialize freeze range depending on hemisphere specified
        IF( HEMISW .EQ. 'NORTH' ) THEN
            FRZRANGE1 = 0    ! array
            FRZRANGE2 = 400  ! array

        ELSE IF( HEMISW .EQ. 'SOUTH' ) THEN
            FRZRANGE1 = 400  ! array
            FRZRANGE2 = 0    ! array

        END IF

C.........  Loop through time steps of Met data
        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS

C.............  Set day index
            D = JDATE - ( JDATE / 1000 ) * 1000    ! integer math            

C.............  Read in temperature data
            IF( .NOT. READ3( MNAME, VARNAM, 1,
     &                       JDATE, JTIME, TMPR      ) ) THEN

                L1 = LEN_TRIM( VARNAM )
                L2 = LEN_TRIM( MNAME )
                MESG = 'ERROR: Could not read variable "' //
     &                 VARNAM( 1:L1 ) // '" from file  "' // 
     &                 MNAME( 1:L2 )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

C.............  Flag cells/dates with freezing temperatures
C.............  Set global all cell-time freeze and no-freeze flags
            DO C = 1, NGRID

                IF( TMPR( C ) .LE. FREEZEVAL ) THEN
                    LFREEZE( D,C ) = .TRUE.
                    NFLAG = .FALSE.
                ELSE
                    FFLAG = .FALSE.
                END IF

            END DO

            CALL NEXTIME( JDATE, JTIME, TSTEP )
            
        END DO

C.........  Assign cells into four groups:
C           Group 1: no freeze days
C           Group 2: Northern Hemisphere winter pattern
C           Group 3: Southern Hemisphere winter pattern
C           Group 4: all freeze days

        DO C = 1, NGRID

C.............  All cells are non-freezing at all times
            IF ( NFLAG ) THEN
                FREEZEGRP( C ) = 1
                CYCLE

C.............  All cells are freezing at all times
            ELSE IF( FFLAG ) THEN
                FREEZEGRP( C ) = 4
                CYCLE

            ENDIF

C.............  If hemisphere has been specified, then set group accordingly
            IF( HEMISW .EQ. 'NORTH' ) THEN
                FREEZEGRP( C ) = 2

            ELSE IF( HEMISW .EQ. 'SOUTH' ) THEN
                FREEZEGRP( C ) = 3

            END IF

        END DO

C.........  Determine first/last Julian day for July 
        JAN01 = 1
        JUL31 = JULIAN( YEAR, 7 , 31 )
        DEC31 = JULIAN( YEAR, 12, 31 )

C.........  For groups 2 and 3, determine date ranges for freeze periods
        DO C = 1, NGRID

            SELECT CASE( FREEZEGRP( C ) )

C.............  Northern hemisphere
            CASE( 2 )

C.................  Set last freeze day of spring
                DO D = JUL31, JAN01, -1

                    IF( D .GE. SDAY .AND. D .LE. EDAY ) THEN
                        IF( LFREEZE( D,C ) ) THEN
                            FRZRANGE1( C ) = D
                            EXIT
                        END IF
                    END IF

                END DO

C.................  Set first freeze date of fall
                DO D = JUL31, DEC31

                    IF( D .GE. SDAY .AND. D .LE. EDAY ) THEN
                        IF( LFREEZE( D,C ) ) THEN
                            FRZRANGE2( C ) = D
                            EXIT
                        END IF
                    END IF

                END DO

C.............  Southern hemisphere
            CASE( 3 )

C.................  Set first freeze date of year
                DO D = JAN01, JUL31

                    IF( D .GE. SDAY .AND. D .LE. EDAY ) THEN
                        IF( LFREEZE( D,C ) ) THEN
                            FRZRANGE1( C ) = D
                            EXIT
                        END IF
                    END IF

                END DO
                
C.................  Set last freeze date of year
                DO D = DEC31, JUL31, -1

                    IF( D .GE. SDAY .AND. D .LE. EDAY ) THEN
                        IF( LFREEZE( D,C ) ) THEN
                            FRZRANGE2( C ) = D
                            EXIT
                        END IF
                    END IF

                END DO

            END SELECT

        END DO

C.........  Open output file
        NVARS3D = 1
        NLAYS3D = 1
        STIME3D = 0
        TSTEP3D = 240000
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3

        VNAME3D( 1 ) = 'SEASON'
        UNITS3D( 1 ) = 'n/a'
        VTYPE3D( 1 ) = M3INT
        VDESC3D( 1 ) = '1=summer, 0=winter'

        MESG = 'Enter logical name for the BIOGENIC SEASONS file'
        ONAME = 'BIOSEASON'
        ONAME = PROMPTMFILE( MESG, FSUNKN3, ONAME, PROGNAME )
 
C.........  Loop through days and cells and output freeze date status for
C           each day and grid cell

        JDATE = SDATE
        JTIME = 0
        DO D = SDAY, EDAY

C.............  Set season status of grid cells based on freeze information
            DO C = 1, NGRID

                SELECT CASE( FREEZEGRP( C ) )
                CASE( 1 )
                    SEASON( C ) = 1

                CASE( 2 )
                    IF( D .LE. FRZRANGE1( C ) .OR.
     &                  D .GE. FRZRANGE2( C )      ) THEN
                        SEASON( C ) = 0
                    ELSE
                        SEASON( C ) = 1
                    END IF

                CASE( 3 )
                    IF( D .GE. FRZRANGE1( C ) .OR.
     &                  D .LE. FRZRANGE2( C )      ) THEN
                        SEASON( C ) = 0
                    ELSE
                        SEASON( C ) = 1
                    END IF

                CASE( 4 )
                    SEASON( C ) = 0

                END SELECT

            END DO

C.............  Write gridded season data
            IF( .NOT. WRITE3( 
     &          ONAME, 'SEASON', JDATE, JTIME, SEASON ) ) THEN

                MESG = 'Problem writing "SEASON" to file "' //
     &                 ONAME( 1:LEN_TRIM( ONAME ) ) // '."'

                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

C.............  Get next output time step
            CALL NEXTIME( JDATE, JTIME, TSTEP3D )

        END DO

C.........  End of program:

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

        END PROGRAM METSCAN  

