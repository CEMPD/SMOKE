
      PROGRAM SMK2EMIS

C***********************************************************************
C  program body starts at line 127
C
C  DESCRIPTION:
C       Use SMOKE NetCDF gridded input to produce UAM EMISSIONS input 
C       files.  Originally developed for EMISSIONS input, but potential 
C       use for other UAM input files.
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
C*************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !

C...........   PARAMETERS and their descriptions:

        CHARACTER*50, PARAMETER :: SCCSW = '@(#)$Id$'

C.........  EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         ENVINT
        REAL            ENVREAL
        INTEGER         GETNUM
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         SECSDIFF
        INTEGER         TRIMLEN
 
        EXTERNAL    ENVINT, ENVREAL, GETNUM, PROMPTFFILE, 
     &              PROMPTMFILE, SECSDIFF, TRIMLEN

C.........  Local variables

        INTEGER   BTIME         ! Beginning time of output file (YYYYDDD)
        INTEGER   ETIME         ! Ending time of output file (HHMMSS)
        INTEGER   IBD           ! Start date for each output time period (YYDDD)
        INTEGER   IED           ! End date for output time period (YYDDD)
        INTEGER   IFTYPE( 10 )  ! UAM filename in integer format
        INTEGER   IFNOTE( 60 )  ! UAM note in integer format
        INTEGER   SDATE         ! UAM file starting date
        INTEGER   STIME         ! UAM file starting hour
        INTEGER   JDATE         ! Date variable (YYYYDDD)
        INTEGER   SDATEP1       ! Date variable plus one time step (YYYYDDD)
        INTEGER   JTIME         ! Time variable (HHMMSS)
        INTEGER   STIMEP1       ! Time variable plus one time step (HHMMSS)
        INTEGER   :: NSEG = 1   ! Number of segments always 1
        INTEGER   ODEV          ! UAM file desigator for output file
        INTEGER   TSTEP         ! Time step (HHMMSS)
        INTEGER   IUTMZONE      ! UTM zone
        INTEGER   TZONE         ! time zone of output
        INTEGER   NSTEPS        ! number of time steps
        INTEGER   NSPECS        ! number of species
        INTEGER   NZLOWR, NZUPPR ! UAM-IV diffbreak layer variables
        INTEGER   NLAYS         ! number of vertical layers in UAM file
        INTEGER   :: IXY = 0    ! segment origin always 0
        INTEGER   LDEV          !  unit number for log file
        INTEGER   IOS
        INTEGER   I, J, K, HR
        INTEGER   NDAYS, NHRS   

        INTEGER, ALLOCATABLE ::  INTNAM ( : , : )  ! species names

        REAL    BT            ! Beginning time for a time step period
        REAL    ET            ! Ending time for a time step period
        REAL    UTMX, UTMY    ! 
        REAL    HTSUR, HTLOWR, HTUPPR  ! UAM-IV layer ht variables 
        REAL    XORIG, YORIG  ! Grid origin
        REAL    XCELL, YCELL  ! Grid resolution

        REAL,ALLOCATABLE ::  EMIS( :, : ) ! Input and output variable
  
        CHARACTER*16   ENAME   !  NetCDF logical input file name
        CHARACTER*10   HDRKEY  !  UAM file name label for output file
        CHARACTER*1    TEMPCH  !  For integer conversion
        CHARACTER*256  MESG    !  Temporary message array
        CHARACTER*60   FNOTE   !  UAM file description
        CHARACTER*44   NOTEDEF ! UAM file note default
        CHARACTER*44   UNOTE   ! UAM file note from env variable

        CHARACTER*16 :: PROGNAME = 'SMK2EMIS'   !  program name

C***********************************************************************
C   begin body of program SMK2EMIS

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Prompt for name of NetCDF input file

        ENAME = PROMPTMFILE(
     &       'Enter logical name for SMOKE gridded input (NetCDF) file',
     &        FSREAD3, 'AGTS', PROGNAME )

C.........  Prompt for keyword name of output file     
C.........  Only EMISSIONS supported at this point

        MESG = 'UAM file label for output file'
        CALL ENVSTR( 'FLABEL', MESG, 'EMISSIONS', HDRKEY, IOS )

C..... Convert UAM keyword name to integer format (as in UAM)

        DO J= 1, 10
            TEMPCH = HDRKEY( J:J )
            READ(TEMPCH, 91000) IFTYPE( J )
        ENDDO

C.........  Prompt for name of output file

        ODEV = PROMPTFFILE(
     &         'Enter logical name for ' // HDRKEY // 
     &         ' UAM output file ',
     &         .FALSE., .FALSE., 'UAMGTS', PROGNAME )

C.........  Read NetCDF header information:
C.........  coordinate info, date, time, timestep , and variables

        IF ( .NOT. DESC3( ENAME ) ) THEN

            MESG = 'Could not get description of file "' //
     &             ENAME( 1:TRIMLEN( ENAME ) ) // '"'
            CALL M3EXIT( PROGNAME , 0, 0, MESG, 2 )

        ENDIF

        NSPECS = NVARS3D

C.........  Get the time zone for output of the emissions
        TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )

C.........   Get setting from environment variables
C.........   for UTM and UAM4 variables

        IUTMZONE = ENVINT( 'UTM_ZONE', 'UTM zone ',
     &                   -9 , IOS )
      
        NLAYS = ENVINT( 'UAM_LAYERS', 'Number of vertical layers',
     &                   NLAYS3D , IOS )

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

C..... Create FNOTE from grid name and UAM_NOTE env. variable
C..... Convert FNOTE to integer format (as in UAM)

        NOTEDEF = ' UAM gridded emissions from ' // PROGNAME

        MESG = 'UAM file note for output file'
        CALL ENVSTR( 'UAM_NOTE', MESG, NOTEDEF , UNOTE, IOS )

        FNOTE  = GDNAM3D // UNOTE ( 1 : 44 )  

        DO J = 1,60
            TEMPCH = FNOTE( J: J )
            READ(TEMPCH, 91000) IFNOTE( J )
        ENDDO

C..... Allocate for species name in interger format

        ALLOCATE( INTNAM( NSPECS, 10  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INTNAM', PROGNAME )

C..... Convert output species names to integer format

        DO I = 1, NSPECS
          DO J = 1, 10
              TEMPCH = VNAME3D( I )( J:J )
              READ( TEMPCH , 91000) INTNAM ( I, J ) 
          ENDDO
        ENDDO

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

        ENDIF

C.............  Prompt for starting date of UAM output file

        SDATE  = GETNUM( JDATE, 9999999, JDATE,
     &                  'Enter starting date (YYYYDDD)' )

        IF( SDATE .NE. JDATE ) JTIME = 0  ! Reset because now could be 0 -24

C.............  Prompt for starting time of UAM output file

        STIME  = GETNUM( JTIME, 235959, JTIME,
     &                 'Enter starting time (HHMMSS)' )

C.............  Prompt for time number of time steps 
C.............  Number of hours to output to UAM file   

        I = SECSDIFF( JDATE, JTIME, SDATE, STIME ) / 3600
        NSTEPS = NSTEPS - I
        MESG   = 'Enter number of time steps'
        NSTEPS = GETNUM( 1, NSTEPS, NSTEPS, MESG )

C...   Calculate ending date and time

        NDAYS = NSTEPS / 24 
        NHRS  = MOD( NSTEPS , 24 ) 

        IF ( SDATE .LT. 2000000 ) THEN
          IBD = SDATE  - 1900000
        ELSE
          IBD = SDATE  - 2000000
        ENDIF
       
        IED = IBD + NDAYS
        ETIME = STIME + ( NHRS * TSTEP ) 

C....   Write UAM EMISSIONS header

        WRITE(ODEV,ERR=7010) IFTYPE, IFNOTE, NSEG, NSPECS, IBD,
     &                       BTIME, IED, ETIME

        WRITE(ODEV,ERR=7010) UTMX, UTMY, IUTMZONE,
     &                       XORIG, YORIG, XCELL,
     &                       YCELL, NCOLS3D, NROWS3D, NLAYS,
     &                       NZLOWR, NZUPPR, HTSUR, HTLOWR, HTUPPR 


        WRITE(ODEV,ERR=7010) IXY, IXY, NCOLS3D, NROWS3D

C.........  Write out integer form of output variable names

        WRITE(ODEV,ERR=7010) ((INTNAM(J,I),I=1,10),J=1,NSPECS)

C.........  Allocate memory for emissions 

        ALLOCATE( EMIS( NCOLS3D, NROWS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMIS', PROGNAME )

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
            ENDIF

            WRITE (*,93000) SDATE, STIME

C.............  Set UAM formatted times and dates

            WRITE(ODEV,ERR=7011) IBD, BT, IED, ET

C.............  Loop over output variables

            DO K = 1, NSPECS     
       
                EMIS = 0.0   !  array
 
C.....................  Read input file for time and species of interest

                IF ( .NOT. READ3( ENAME, VNAME3D( K ), 1, 
     &                            SDATE, STIME, EMIS ) ) THEN
                  CALL M3ERR( PROGNAME , 0, 0,
     &                       'Error reading ' // VNAME3D( K ) // 
     &                       ' from file ' // ENAME , .TRUE. )
                ENDIF 

C.............  Write UAM formatted output data for all variables

                WRITE(ODEV,ERR=7011) NSEG, (INTNAM( K,J ),J=1,10),
     &                     (( EMIS( I,J ),I = 1,NCOLS3D ), J=1,NROWS3D )

            ENDDO       ! Output variables loop

C.............  Update times and dates for next pass through loop

            SDATE = SDATEP1
            STIME = STIMEP1
            CALL NEXTIME( SDATEP1, STIMEP1, TSTEP) 

        ENDDO ! Time step loop 


C.........  End of program

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


C******************    ERROR MESSAGES     ******************************

7010  WRITE (*,*) 'ERROR:  WRITING UAM HEADER OF OUTPUT'
      STOP
7011  WRITE (*,*) 'ERROR:  WRITING UAM RECORD OUTPUT'
      STOP

C******************  FORMAT  STATEMENTS   ******************************

91000 FORMAT (A1)
93000 FORMAT (/, 5X, 'Writing out for date ',i7,' and time ',i6,'.')

C...........   Internal buffering formats............ 94xxx
94010   FORMAT( 10 ( A, :, I5, :, 2X ) )


      END PROGRAM SMK2EMIS

