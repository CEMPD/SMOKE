
      PROGRAM GEOFAC

C***********************************************************************
C
C  DESCRIPTION: Takes gridded SMOKE emissions file and applies
C               a user-supplied factor for each individual species
C               and applies it for a certain geographical region
C               obtained from input mask file..
C               The resulting emissions are output to a netCDF file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C             10/00 : Prototype by JMV
C             04/24 : Improvised by Huy Tran (UNC-IE)
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
C*************************************************************************

      USE MODFILESET

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     ! Emissions constants
        INCLUDE 'SETDECL.EXT'     ! FileSetAPI function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         GETFLINE
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        INTEGER         TRIMLEN

        EXTERNAL  GETFLINE, PROMPTFFILE, PROMPTMFILE,
     &            TRIMLEN
        LOGICAL, EXTERNAL :: ENVYN          ! from IOAPI
           
C...........   PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

C...........  LOCAL VARIABLES

C       REAL  FEMIS                              ! temporary value for emissions
        REAL, ALLOCATABLE :: EMIS( :, : , : )    ! emissions 
        REAL, ALLOCATABLE :: SPCFAC( : )         ! species factors

        INTEGER  TSTEP                           ! time step
        INTEGER  I, J, K , L, M, N               ! counters
        INTEGER, ALLOCATABLE ::  IMASK ( :, : )  ! mask values
        INTEGER  LDEV                            ! log file unit number
        INTEGER  RDEV                            ! species factors unit number
        INTEGER  NSPECS                          ! number of species
        INTEGER  HR                              ! hour loop counter
        INTEGER  NLINES                          ! number of species factors 
        INTEGER  NSTEPS                          ! number of time steps
        INTEGER  IOS                             ! iostat
        INTEGER  SDATE                           ! start date
        INTEGER  STIME                           ! start time
        INTEGER  NLAYS                           ! number of layers in emis file

        CHARACTER(16), ALLOCATABLE :: SPCNAM ( : ) ! species names with facs
        CHARACTER(16)  ENAME                       ! logical name for gridded emis input file
        CHARACTER(16)  ONAME                       ! logical name for output file
        CHARACTER(16)  MNAME                       ! logical name for mask file
        CHARACTER(300) MESG                        ! message buffer for M3EXIT()

        CHARACTER(16) :: PROGNAME = 'GEOFAC'       !  program name

        LOGICAL :: ZEROYN                          ! flag to zero out emissions for grid-cell that are not interested
        LOGICAL :: MASKFRAC                        ! flag to use fractional mask value; if false then use unity mask value
        CHARACTER(16) :: MASKVAR                   ! Variable name in mask file indicating masked region
        REAL :: MAXMSK

C***********************************************************************
C   begin body of program

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Prompt for name of NetCDF input file

        ENAME = PROMPTSET(
     &       'Enter logical name for SMOKE gridded input (NetCDF) file',
     &        FSREAD3, 'INFILE', PROGNAME )

C Huy Tran 04/24: ! Call DESCSET with index -1 will cause it to use 
C /NUMBER OF VARIABLES/ in FILEDESC attribute which may differ from the actual 
C number of variable in SMOKE gridded input file
C Modified to use index 1 to force using number of variable in the header of 
C SMOKE gridded input file
C       IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN
        IF ( .NOT. DESCSET( ENAME, 1 ) ) THEN

            MESG = 'Could not get description of file "' //
     &             ENAME( 1:TRIMLEN( ENAME ) ) // '"'
            CALL M3EXIT( PROGNAME , 0, 0, MESG, 2 )

        ENDIF

C.........  Assign local variables

        TSTEP  = TSTEP3D
        SDATE  = SDATE3D
        STIME  = STIME3D
        NSPECS = NVARS3D
        NSTEPS = MXREC3D
        NLAYS  = NLAYS3D

C.......   Get the mask file name

        MNAME = PROMPTMFILE(
     &       'Enter logical name for mask input (NetCDF) file',
     &        FSREAD3, 'GEOMASK', PROGNAME )

C.......   Get variable name that represent mask area in mask file
        CALL ENVSTR( 'MASKVAR', 'Variable indicating masked region',
     &               'MASK', MASKVAR, IOS )

C.......   Check if to use fractional mask value or unity mask value
        MASKFRAC = ENVYN( 'MASKFRAC', 'Whether or not using fractional mask value',
     &               .TRUE., IOS )

C.......   Get the species fac file name

        RDEV = PROMPTFFILE(
     &           'Enter logical name for factors file',
     &           .TRUE., .TRUE., 'SPECFACS', PROGNAME )

C.......   Get whether or not to zero out emission for grid-cell not = 1 in mask file
        ZEROYN = ENVYN( 'ZEROEMIS', 'Whether or not zeroing out emiss',
     &                       .FALSE., IOS )

C.........  Read mask file header      

        ALLOCATE( IMASK ( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMASK', PROGNAME )
     
        IF ( .NOT. READ3( MNAME, MASKVAR, 1,
     &                            0, 0, IMASK ) ) THEN
                  CALL M3ERR( PROGNAME , 0, 0,
     &                       'Error reading masked region from file '
     &                       // MNAME , .TRUE. )
        ENDIF
        MAXMSK = maxval(IMASK)
 
C.........  Read species factors

        NLINES = GETFLINE( RDEV, 'Species factors file' )

        ALLOCATE( SPCNAM( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCNAM', PROGNAME )

        ALLOCATE( SPCFAC( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCFAC', PROGNAME )

C.........  Create file descriptions for output file

        FDESC3D  = '  ' 
        FDESC3D( 1 ) = '/' // PROGNAME // '/'

C........   Read in factors for species

        DO N = 1, NLINES
         READ( RDEV, 93000 ) SPCNAM ( N ) , SPCFAC ( N )
         M = N + 1 
         WRITE ( FDESC3D( M ) , 94000 ) SPCNAM( N ), SPCFAC( N )  
        ENDDO

C.........  Open output file

        ONAME = PROMPTSET(
     &          'Enter logical name for OUTPUT gridded netCDF file '
     &          , FSUNKN3, 'OUTFILE',  PROGNAME      )

C.........  Allocate memory for emissions

        ALLOCATE( EMIS( NCOLS3D, NROWS3D, NLAYS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMIS', PROGNAME )

C.........  Loop over time steps

        DO HR = 1, NSTEPS

C...........  Write to screen because WRITE3 only writes to LDEV
C...........  Loop through cells and apply factor to cells in mask

          DO  L = 1, NSPECS
C...........  Read input file for time and species of interest

              IF ( .NOT. READSET( ENAME, VNAMESET( L ), ALLAYS3,
     &              ALLFILES, SDATE, STIME, EMIS   ) ) THEN
                 CALL M3ERR( PROGNAME , 0, 0,
     &               'Error reading ' // VNAMESET( L ) //
     &               ' from file ' // ENAME , .TRUE. )
              ENDIF

              DO I = 1, NCOLS3D 
              DO J = 1, NROWS3D

C...........    If grid cell is NOT in geographical region of interest
                IF (IMASK ( I, J ) .EQ. 0 ) GOTO 100
                IF ( .NOT. MASKFRAC .AND.
     &                IMASK ( I, J ) .NE. MAXMSK ) GOTO 100
                GOTO 200
100             IF (ZEROYN) EMIS ( I, J, : ) = 0.0 ! set emission to zero if needed
                CYCLE                              ! skip to next grid cell

200             CONTINUE

                DO M = 1, NLINES

C...........       Find out if factor for this species exists
                   IF ( VNAMESET ( L ) .EQ. SPCNAM ( M )  ) THEN

                        EMIS (I,J,:) = EMIS ( I, J, :  ) * SPCFAC( M ) 

                   ENDIF

                ENDDO
    
              ENDDO ! Finished NROWS
              ENDDO ! Finished NCOLS

C............  Write out new emissions
              
C           IF ( .NOT. WRITESET( ONAME, VNAMESET( L ), -1, SDATE,
            IF ( .NOT. WRITESET( ONAME, VNAMESET( L ),  1, SDATE,
     &                           STIME, EMIS ) ) THEN

              CALL M3EXIT( PROGNAME , SDATE, STIME,
     &                          'Could not write "' //
     &                           VNAMESET( L ) //
     &                          '" to ' // ONAME, 2 )

            END IF         !  if writeset() failed

           ENDDO           ! Finished for this species

C............  Next time step

           CALL NEXTIME( SDATE, STIME, TSTEP)

        ENDDO 

C.........   End of program:

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

C...........   Informational (LOG) message formats... 92xxx

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A, F6.3 ) 

C...........   Internal buffering formats............ 94xxx

94000   FORMAT ( '/SPECIES FAC/ ', A, F6.3 )
     
        END PROGRAM  GEOFAC
