
      PROGRAM SMK2EMIS

C***********************************************************************
C  program body starts at line 152
C
C  DESCRIPTION:
C       Use SMOKE NetCDF gridded input to produce UAM EMISSIONS input 
C       files.  Outputs 2-d files only.
C
C  PRECONDITIONS REQUIRED:
C       M3IO single layer gridded data
C       Tested for lat-lon and UTM projections only
C       Tested for CAMx and UAM-V compatiblity
C       Not fully tested for UAM-IV compatiblity
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       GETYN, PROMPTFFILE, PROMPTMFILE, NEXTIME
C
C  REVISION  HISTORY:
C       Prototype  4/00 by JMV
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

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CRL

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   PARAMETERS and their descriptions:

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.8_Oct2020$' ! CVS release tag

C.........  EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)    CRLF
        INTEGER         ENVINT
        REAL            ENVREAL
        LOGICAL         ENVYN
        INTEGER         GETFLINE
        INTEGER         GETNUM
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         TRIMLEN
 
        EXTERNAL    CRLF, ENVINT, ENVREAL, ENVYN, GETFLINE, GETNUM, 
     &              INDEX1, PROMPTFFILE, PROMPTMFILE, SECSDIFF, TRIMLEN

C.........  Allocatable arrays
        CHARACTER(4), ALLOCATABLE :: SPCNAM ( :,: )  ! output species names
        REAL,         ALLOCATABLE :: EMIS   ( :,: )  ! input and output variable

        CHARACTER(IOVLEN3), ALLOCATABLE :: OUTNAMS( : ) ! output spec names

C.........  Static arrays
        CHARACTER(4)    CFTYPE( 10 )  ! UAM filename
        CHARACTER(4)    CFNOTE( 60 )  ! UAM note

        CHARACTER(IOVLEN3) SEGMENT( 2 )

C.........  Unit numbers and logical file names
        INTEGER :: LDEV          ! unit number for log file
        INTEGER :: MDEV          ! unit number for mapping file
        INTEGER :: ODEV          ! UAM file desigator for output file

        CHARACTER(16)  ENAME     !  NetCDF logical input file name

C.........  Other local variables

        INTEGER    I, J, K, L1, L2, HR   ! indices and counters

        INTEGER    BTIME         ! beginning time of output file (YYYYDDD)
        INTEGER    ETIME         ! ending time of output file (HHMMSS)
        INTEGER    IBD           ! start date for each output time period (YYDDD)
        INTEGER    IED           ! end date for output time period (YYDDD)
        INTEGER    IOS           ! i/o status
        INTEGER    IREC          ! line count
        INTEGER    IUTMZONE      ! UTM zone
        INTEGER :: IXY = 0       ! segment origin always 0
        INTEGER    JDATE         ! date variable (YYYYDDD)
        INTEGER    JTIME         ! time variable (HHMMSS)
        INTEGER :: NLAYS = 1     ! number of vertical layers in UAM file
        INTEGER    NLINE         ! no. lines input ascii file
        INTEGER :: NSEG = 1      ! number of segments always 1
        INTEGER    NSTEPS        ! number of time steps
        INTEGER    NSPECS        ! number of species
        INTEGER    NZLOWR, NZUPPR ! UAM-IV diffbreak layer variables
        INTEGER    NDAYS, NHRS   
        INTEGER    SDATE         ! UAM file starting date
        INTEGER    SDATEP1       ! date variable plus one time step (YYYYDDD)
        INTEGER    STIME         ! UAM file starting hour
        INTEGER    STIMEP1       ! time variable plus one time step (HHMMSS)
        INTEGER    TSTEP         ! time step (HHMMSS)
        INTEGER    TZONE         ! time zone of output

        LOGICAL :: EFLAG   = .FALSE. ! true: error found
        LOGICAL :: MAPFLAG = .FALSE. ! true: use variable mapping input file

        REAL       BT            ! Beginning time for a time step period
        REAL       ET            ! Ending time for a time step period
        REAL       HTSUR, HTLOWR, HTUPPR  ! UAM-IV layer ht variables 
        REAL       UTMX, UTMY    ! 
        REAL       XORIG, YORIG  ! Grid origin
        REAL       XCELL, YCELL  ! Grid resolution

        CHARACTER      TEMPCH  !  For integer conversion
        CHARACTER(10)  HDRKEY  !  UAM file name label for output file
        CHARACTER(44)  NOTEDEF !  UAM file note default
        CHARACTER(44)  UNOTE   !  UAM file note from env variable
        CHARACTER(60)  FNOTE   !  UAM file description
        CHARACTER(80)  LINE    !  line buffer
        CHARACTER(256) MESG    !  Temporary message array

        CHARACTER(16) :: PROGNAME = 'SMK2EMIS'   !  program name

C***********************************************************************
C   begin body of program SMK2EMIS

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get environment variables that control this program
        MAPFLAG = ENVYN( 'SMK2EMIS_VMAP_YN', 'Use names mapping file', 
     &                   .FALSE., IOS )

C.........  Set source category based on environment variable setting
        CALL GETCTGRY

C.........  Prompt for name of NetCDF input file

        ENAME = PROMPTMFILE(
     &       'Enter logical name for SMOKE gridded input (NetCDF) file',
     &        FSREAD3, CRL // 'GTS_L', PROGNAME )

C.........  Prompt for variable mapping file, if needed
        IF( MAPFLAG ) THEN

            MDEV = PROMPTFFILE( 
     &             'Enter logical name for VARIABLE NAME MAP file',
     &             .TRUE., .TRUE., 'VNAMMAP', PROGNAME )

        END IF

C.........  Prompt for keyword name of output file     
C.........  Only EMISSIONS supported at this point

        MESG = 'UAM file label for output file'
        CALL ENVSTR( 'FLABEL', MESG, 'EMISSIONS', HDRKEY, IOS )

C.........  Convert UAM keyword name to character format (as in UAM)
        CFTYPE = ' '
        DO J = 1, 10
            CFTYPE( J ) = HDRKEY( J:J )
        END DO

C.........  Prompt for name of output file

        ODEV = PROMPTFFILE(
     &         'Enter logical name for ' // HDRKEY // 
     &         ' UAM output file ',
     &         .FALSE., .FALSE., 'UAM_' // CRL // 'GTS', PROGNAME )

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........   Get setting from environment variables
C.........   for UTM and UAM4 variables

        IUTMZONE = ENVINT( 'UTM_ZONE', 'UTM zone ',
     &                   -9 , IOS )
      
        NZLOWR = ENVINT( 'UAM4_LAYBELOW', 'Layers below diffbreak',
     &                   0 , IOS )
        NZUPPR = ENVINT( 'UAM4_LAYABOVE', 'Layers above diffbreak',
     &                   0 , IOS )

        MESG = 'Height of surface layer [meters]'
        HTSUR = ENVREAL( 'UAM4_HTSFC', MESG, 0., IOS )

        MESG = 'Min. ht. of cells between sfc and diffbreak [meters]'
        HTLOWR = ENVREAL( 'UAM4_HTLOWR', MESG, 0., IOS )

        MESG = 'Min. ht. of cells between diffbreak and top [meters]'
        HTUPPR = ENVREAL( 'UAM4_HTUPPR', MESG, 0., IOS )

C.........  Create FNOTE from grid name and UAM_NOTE env. variable
C.........  Convert FNOTE to character format (as in UAM)

        NOTEDEF = ' UAM gridded emissions from ' // PROGNAME

        MESG = 'UAM file note for output file'
        CALL ENVSTR( 'UAM_NOTE', MESG, NOTEDEF , UNOTE, IOS )

        FNOTE  = GDNAM3D // UNOTE ( 1 : 44 )  

        CFNOTE = ' '
        DO J = 1, 60
            CFNOTE( J ) = FNOTE( J:J )
        END DO

C.........  Read NetCDF header information:
C.........  coordinate info, date, time, timestep , and variables

        IF ( .NOT. DESC3( ENAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME , 0, 0, MESG, 2 )

        END IF

        NSPECS = NVARS3D
        NSTEPS = MXREC3D
        TSTEP  = TSTEP3D
        JDATE  = SDATE3D
        JTIME  = STIME3D
        XORIG  = XORIG3D
        YORIG  = YORIG3D
        XCELL  = XCELL3D
        YCELL  = YCELL3D

C...  If time step other than hourly then exit 
C...  Assumes UAM/CAMx wants hourly emissions for now

        IF ( TSTEP .NE. 10000 ) THEN

           WRITE( MESG,94010 ) 'ERROR: Gridded SMOKE emissions '//
     &             'time-step not hourly :  TSTEP =  ', TSTEP   
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.............  Prompt for starting date of UAM output file

        SDATE  = GETNUM( JDATE, 9999999, JDATE,
     &                  'Enter starting date (YYYYDDD)' )

        IF( SDATE .NE. JDATE ) JTIME = 0  ! Reset because now could be 0 -24

C.............  Prompt for starting time of UAM output file

        STIME  = GETNUM( JTIME, 235959, JTIME,
     &                 'Enter starting time (HHMMSS)' )
        BTIME  = INT( STIME / TSTEP )

C.............  Prompt for time number of time steps 
C.............  Number of hours to output to UAM file   

        I = SECSDIFF( JDATE, JTIME, SDATE, STIME ) / 3600
        NSTEPS = NSTEPS - I
        MESG   = 'Enter number of time steps'
        NSTEPS = GETNUM( 1, NSTEPS, NSTEPS, MESG )

C.........  Calculate ending date and time

        NDAYS = NSTEPS / 24 
        NHRS  = MOD( NSTEPS , 24 ) 

        IF ( SDATE .LT. 2000000 ) THEN
          IBD = SDATE  - 1900000
        ELSE
          IBD = SDATE  - 2000000
        END IF
       
        IED = IBD + NDAYS
        ETIME = MOD( ( BTIME + NHRS - 1 ), 24 )

C.........  Allocate memory for and read remapping information
        ALLOCATE( OUTNAMS( NSPECS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OUTNAMS', PROGNAME )

C.........  Initialize names map as input names
        OUTNAMS( 1:NSPECS ) = VNAME3D( 1:NSPECS )

        IF ( MAPFLAG ) THEN

            NLINE = GETFLINE( MDEV, 'Variable mapping file' )

            CALL M3MSG2( 'Reading variable mapping file...' )

            IREC = 0
            DO I = 1, NLINE

                READ( MDEV, 93000, END=999, IOSTAT=IOS ) LINE
                IREC = IREC + 1

                IF( LINE .EQ. ' ' ) CYCLE

                IF ( IOS .NE. 0 ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                  'I/O error', IOS, 
     &                  'reading variable mapping file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE

                END IF

                CALL PARSLINE( LINE, 2, SEGMENT )

C.................  Find variable name in list of species
                K = INDEX1( SEGMENT( 1 ), NSPECS, VNAME3D ) 

C.................  If not found, write message
                IF( K .LE. 0 ) THEN

                    L1 = LEN_TRIM( SEGMENT( 1 ) )
                    MESG = 'WARNING: Ignoring "' // 
     &                     SEGMENT( 1 )( 1:L1 ) //
     &                     '" in variable mapping file'
                    CALL M3MSG2( MESG )

C.................  If found, write message and reassign
                ELSE

                    L1 = LEN_TRIM( SEGMENT( 2 ) )
                    L2 = LEN_TRIM( OUTNAMS( K ) )
                    MESG = 'NOTE: Renaming variable "' // 
     &                     OUTNAMS( K )( 1:L2 ) // '" as "' //
     &                     SEGMENT( 2 )( 1:L1 ) // '"'
                    CALL M3MSG2( MESG )

                    OUTNAMS( K ) = SEGMENT( 2 )

                END IF

            END DO

        END IF

        IF ( EFLAG ) THEN

            MESG = 'Problem reading mapping names file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Allocate for species name in character format
        ALLOCATE( SPCNAM( NSPECS, 10  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPCNAM', PROGNAME )

C.........  Convert output species names to character format
        DO I = 1, NSPECS
          DO J = 1, 10
              SPCNAM( I,J ) = OUTNAMS( I )( J:J )
          END DO
        END DO

C.........  Allocate memory for emissions 
        ALLOCATE( EMIS( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMIS', PROGNAME )

C.........  Write UAM Binary file  .........................................

C....   Write UAM EMISSIONS header

        WRITE(ODEV,ERR=7010) CFTYPE, CFNOTE, NSEG, NSPECS, IBD,
     &                       REAL( BTIME ), IED , REAL( ETIME )

        WRITE(ODEV,ERR=7010) UTMX, UTMY, IUTMZONE,
     &                       XORIG, YORIG, XCELL,
     &                       YCELL, NCOLS3D, NROWS3D, NLAYS,
     &                       NZLOWR, NZUPPR, HTSUR, HTLOWR, HTUPPR 


        WRITE(ODEV,ERR=7010) IXY, IXY, NCOLS3D, NROWS3D

C.........  Write out character form of output variable names

        WRITE(ODEV,ERR=7010) ((SPCNAM(J,I),I=1,10),J=1,NSPECS)

C.........  Initialize time and date for time step after STIME and SDATE

        SDATEP1 = SDATE
        STIMEP1 = STIME
        CALL NEXTIME (SDATEP1, STIMEP1, TSTEP)

C.........  Loop over time steps

        DO HR = 1, NSTEPS 

C.........  Initialize output file special time parameters
C.........  Note 20th and 21st centuries only are supported 

            BT  = REAL( STIME/10000 )
            ET  = BT + 1.
            IF ( SDATE .LT. 2000000 ) THEN
              IBD = SDATE   - 1900000
              IED = SDATEP1 - 1900000 
            ELSE
              IBD = SDATE   - 2000000
              IED = SDATEP1 - 2000000
            END IF

            WRITE (*,93200) SDATE, STIME

C.............  Set UAM formatted times and dates

            WRITE(ODEV,ERR=7011) IBD, BT, IED, ET

C.............  Loop over output variables

            DO K = 1, NSPECS     
       
                EMIS = 0.0   !  array
 
C.....................  Read input file for time and species of interest

                IF ( .NOT. READ3( ENAME, VNAME3D( K ), 1, 
     &                            SDATE, STIME, EMIS )) THEN
                  CALL M3ERR( PROGNAME , 0, 0,
     &                       'Error reading ' // VNAME3D( K ) // 
     &                       ' from file ' // ENAME , .TRUE. )
                END IF 

C.............  Write UAM formatted output data for all variables

                WRITE(ODEV,ERR=7011) NSEG, (SPCNAM( K,J ),J=1,10),
     &                     (( EMIS( I,J ),I = 1,NCOLS3D ), J=1,NROWS3D )

            END DO       ! Output variables loop

C.............  Update times and dates for next pass through loop

            SDATE = SDATEP1
            STIME = STIMEP1
            CALL NEXTIME( SDATEP1, STIMEP1, TSTEP) 

        END DO ! Time step loop 


C.........  End of program

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of variable ' // CRLF() // BLANK5 //
     &         'mapping file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************    ERROR MESSAGES     ******************************

7010  WRITE (*,*) 'ERROR:  WRITING UAM HEADER OF OUTPUT'
      STOP
7011  WRITE (*,*) 'ERROR:  WRITING UAM RECORD OUTPUT'
      STOP

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000 FORMAT ( A )

93100 FORMAT ( A1 )

93200 FORMAT (/, 5X, 'Writing out for date ',i7,' and time ',i6,'.')

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )


      END PROGRAM SMK2EMIS

