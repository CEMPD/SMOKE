
        SUBROUTINE CHKEMFAC( GNAME, TDEV, EDEV, NTP, PFLAG, MODELNAM,
     &                       TZONE, TDMAX )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine is for perform steps needed for using activities and emission factors.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 9/2005 by B. Baek 
C
C**************************************************************************
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
C**************************************************************************

C.........  Modules for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC, ONLY: NEFS, INPUTHC, OUTPUTHC, EMTNAM, EMTPOL,
     &                      NEPOL, NETYPE, EFDAYS, EFIDX, EFLIST,
     &                      EFLOGS, EFTYPE, USETIME

C.........  This module contains the information about the source category
        USE MODINFO,  ONLY: NSRC, NIACT 

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVDVTS

C.........  This module is used for MOBILE6 setup information 
        USE MODMBSET, ONLY: DAILY, WEEKLY, MONTHLY, EPISLEN

C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: ITDATE, STTIME, RUNLEN

C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: TZONES

        IMPLICIT NONE
 
C.........  INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions


C.........  EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         CHKINT
        CHARACTER(2)    CRLF
        INTEGER         GETDATE
        INTEGER         GETFLINE
        CHARACTER(14)   MMDDYY
        INTEGER         PROMPTFFILE
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        LOGICAL         SETENVVAR

        EXTERNAL    CHKINT, CRLF,GETDATE, GETFLINE, PROMPTFFILE, 
     &              SECSDIFF, STR2INT, SETENVVAR

C........  Subroutine Arguments

        CHARACTER(*) , INTENT( IN ) :: GNAME           ! ungridding matrix
        INTEGER      , INTENT( IN ) :: TDEV            ! unit number for emission processes file
        INTEGER      , INTENT( IN ) :: EDEV            ! unit number for ef file list
        INTEGER      , INTENT( IN ) :: NTP             ! starting Julian date
        LOGICAL      , INTENT( IN ) :: PFLAG           ! number of output time steps
        CHARACTER(*) , INTENT( IN ) :: MODELNAM        ! emission factor model name
        INTEGER      , INTENT( IN ) :: TZONE           ! output-file time zone
        INTEGER      , INTENT( IN ) :: TDMAX           ! max no. days in episode

C........   Emission factor arrays        
        INTEGER       , ALLOCATABLE :: SRCS( : )       ! temporary array for sources in each ef file
        INTEGER       , ALLOCATABLE :: UMAT( : )       ! contiguous ungridding matrix

C........   Other local variables

        INTEGER         I, II, J, K, L, L2, N, S
        INTEGER         AVERTYPE            ! time period averaging type
        INTEGER         EARLYDATE           ! earliest starting date based on time zones
        INTEGER         EARLYTIME           ! earliest starting time based on time zones
        INTEGER         EDATE, ETIME        ! ending Julian date and time
        INTEGER         EFSDATE, EFEDATE    ! start and end date of current ef file
        INTEGER         ENDPOS              ! ending position in ef day array
        INTEGER         IOS                 ! i/o status
        INTEGER         LATEDATE            ! latest ending date based on time zones
        INTEGER         LATETIME            ! latest ending time
        INTEGER         NDAYS               ! no. days in episode
        INTEGER         NMATX               ! size of ungridding matrix
        INTEGER         NLINES              ! no. lines in ef list file
        INTEGER         NSTEPS              ! number of output time steps
        INTEGER         SDATE, STIME        ! starting Julian date and time
        INTEGER         STPOS               ! starting position in ef day array
        INTEGER         TSTEP               ! output time step
        INTEGER         VDEV                ! unit no. for inventory data table
        INTEGER         TZMIN               ! minimum time zone in inventory      
        INTEGER         TZMAX               ! maximum time zone in inventory

        LOGICAL      :: FNDOUTPUT = .FALSE. ! true: found output hydrocarbon
        LOGICAL      :: ENDFLAG             ! true: couldn't find file end date
        LOGICAL      :: EFLAG = .FALSE.  !  error-flag

        CHARACTER(256)       CURFNM    ! current emission factor file name
        CHARACTER(16)        CURLNM    ! current ef logical file name
        CHARACTER(3)         INTBUF    ! buffer for integer
        CHARACTER(300)       MESG      ! buffer for M3EXIT() messages
        CHARACTER(20)        SEARCHSTR ! string used in search
        CHARACTER(MXDLEN3)   TEMPLINE  ! line from file description
        CHARACTER(IOVLEN3)   VOLNAM    ! volatile pollutant name

        CHARACTER(16) :: PROGNAME = 'CHKEMFAC' ! program name

C***********************************************************************
C   begin body of program CHKEMFAC

C.....  Perform steps needed for using activities and emission factors

C.....  Allocate memory for emission factor arrays
        ALLOCATE( EFDAYS( NTP,TDMAX,4 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFDAYS', PROGNAME )
        ALLOCATE( EFTYPE( NTP,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFFILE', PROGNAME )
        ALLOCATE( EFIDX( NTP,NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFIDX', PROGNAME )

        EFDAYS = 0
        EFTYPE = ' '
        EFIDX  = 0

C.....  Allocate memory for USETIME( 4 ) (true : time period is used)
        ALLOCATE( USETIME( NTP,4 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'USETIME', PROGNAME )
        USETIME = .FALSE.                    ! true: time period is used

C.....  Read list of emission factor files
        NLINES = GETFLINE( EDEV, 'Emission factor file list' )

        ALLOCATE( EFLIST( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFLIST', PROGNAME )
        ALLOCATE( EFLOGS( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EFLOGS', PROGNAME )
      
        EFLIST = ' '
        EFLOGS = ' '

C.........  Read header of ungridding matrix
        IF( .NOT. DESC3( GNAME ) ) THEN
            MESG = 'Could not get description for file ' // GNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store number of ungridding factors
        NMATX = NCOLS3D

C.........  Allocate memory for ungridding matrix
    	IF( ALLOCATED( UMAT ) ) DEALLOCATE( UMAT )
        ALLOCATE( UMAT( NSRC + 2*NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UMAT', PROGNAME )

C.........  Read ungridding matrix
        CALL RDUMAT( GNAME, NSRC, NMATX, NMATX, UMAT( 1 ),
     &               UMAT( NSRC+1 ), UMAT( NSRC+NMATX+1 )  )

C.........  Read emission processes file.  Populate array in MODEMFAC.
        CALL RDEPROC( TDEV )

C.........  Store sources that are outside the grid
        DO S = 1, NSRC
            IF( UMAT( S ) == 0 ) THEN
                EFIDX( II,S ) = -9
            END IF
        END DO

C.....  Perform steps needed for using activities and emission factors
        CALL RDLINES( EDEV, 'Emission factor file list', NLINES, 
     &                EFLIST )

C.........  Loop through activities and...
C.........  NOTE - this is not fully implemented for multiple activities. 
C           To do this, the data structures and RDEFACS will need to be 
C           updated. Also, the variable names in the emission factor file
C           are not truly supporting 16-character pollutant and 
C           emission process names, because it is only set up for MOBILE5
        DO I = 1, NIACT

C.............  Skip activities that do not have emissions types
            IF( NETYPE( I ) .LE. 0 ) CYCLE            
C.............  Set up emission process variable names
            CALL EFSETUP( 'NONE', MODELNAM, NEFS, VOLNAM )

        END DO

C.....  Read inventory table
        VDEV = PROMPTFFILE( 
     &       'Enter logical name for INVENTORY DATA TABLE file',
     &       .TRUE., .TRUE., 'INVTABLE', PROGNAME )
        CALL RDCODNAM( VDEV )
        
C.........  Check if processing NONHAP values

C.........  Set input and output hydrocarbon names
        INPUTHC = TRIM( VOLNAM )
        OUTPUTHC = 'NONHAP' // TRIM( INPUTHC )

        FNDOUTPUT = .FALSE.
        K = 0

C.........  Loop through all pollutants        
        DO I = 1, MXIDAT
        
            IF( INVDNAM( I ) == OUTPUTHC ) THEN
                FNDOUTPUT = .TRUE.
                CYCLE
            END IF

C.............  If requested hydrocarbon is not TOG or VOC, skip rest of loop
            IF( INPUTHC /= 'TOG' .AND. INPUTHC /= 'VOC' ) EXIT
      
            IF( INVDVTS( I ) /= 'N' ) THEN
        
C.................  Check that pollutant is generated by MOBILE6   
                DO J = 1, NEPOL
                    IF( INVDNAM( I ) == EMTPOL( J ) ) THEN
                        IF( INVDVTS( I ) == 'V' ) THEN
                            K = K + 1
                        ELSE IF( INPUTHC == 'TOG' ) THEN
                            K = K + 1
                        END IF
                        EXIT
                    END IF
                END DO
            END IF
        END DO

C.........  If output was not found, set name to blank        
        IF( .NOT. FNDOUTPUT .OR. K == 0 ) THEN
            OUTPUTHC = ' '             
        END IF

C.........  Rename emission factors if necessary
        IF( OUTPUTHC /= ' ' ) THEN
            DO I = 1, SIZE( EMTNAM,1 )
                L = INDEX( EMTNAM( I,1 ), ETJOIN )
                L2 = LEN_TRIM( ETJOIN )
                
                IF( EMTNAM( I,1 )( L+L2:IOVLEN3 ) == INPUTHC ) THEN
                    EMTNAM( I,1 )( L+L2:IOVLEN3 ) = OUTPUTHC
                    CYCLE
                END IF
            END DO
        END IF

C.........  Loop through EF files
        MESG = 'Checking emission factor files...'
        CALL M3MSG2( MESG )

C.........  Define the minimum and maximum time zones in the inventory
        TZMIN = MINVAL( TZONES )
        TZMAX = MAXVAL( TZONES )

C.........  Adjust TZMIN for possibility of daylight savings
        TZMIN = MAX( TZMIN - 1, 0 )

C.........  Loop over entire episode periods
        DO II = 1, NTP  

C.............  Determine number of days in episode
           IF( PFLAG ) THEN
               SDATE = ITDATE( II )
               STIME = STTIME( II )
               TSTEP  = 10000  ! Only 1-hour time steps supported
               NSTEPS = RUNLEN ( II ) / TSTEP
           ELSE

C.........  Get episode settings from the Models-3 environment variables
C           when $GE_DAT/procdates.txt is not available for episode time periods
               SDATE  = 0
               STIME  = 0
               NSTEPS = 1
               CALL GETM3EPI( TZONE, SDATE, STIME, TSTEP, NSTEPS )
               TSTEP  = 10000  ! Only 1-hour time steps supported
            END IF

C.............  Earliest day is start time in maximum time zone
            EARLYDATE = SDATE
            EARLYTIME = STIME
            CALL NEXTIME( EARLYDATE, EARLYTIME, 
     &                   -( TZMAX - TZONE )*10000 )
            
C.............  If time is before 6 am, need previous day also
            IF( EARLYTIME < 60000 ) EARLYDATE = EARLYDATE - 1
            
C.............  Latest day is end time in minimum time zone
            EDATE = SDATE
            ETIME = STIME
            CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

            LATEDATE = EDATE
            LATETIME = ETIME
            CALL NEXTIME( LATEDATE, LATETIME, 
     &                   -( TZMIN - TZONE )*10000 )
C.............  If time is before 6 am, don't need last day
            IF( LATETIME < 60000 ) LATEDATE = LATEDATE - 1

            NDAYS = SECSDIFF( EARLYDATE, 0, LATEDATE, 0 ) / ( 24*3600 )
            NDAYS = NDAYS + 1

            DO N = 1, NLINES

                CURFNM = EFLIST( N )

C.................  Skip any blank lines
                IF( CURFNM == ' ' ) CYCLE

C.................  Determine file type
                IF( INDEX( CURFNM, 'daily' ) > 0 ) THEN
                    AVERTYPE = DAILY
                ELSE IF( INDEX( CURFNM, 'weekly' ) > 0 ) THEN
                    AVERTYPE = WEEKLY
                ELSE IF( INDEX( CURFNM, 'monthly' ) > 0 ) THEN
                    AVERTYPE = MONTHLY
                ELSE IF( INDEX( CURFNM, 'episode' ) > 0 ) THEN
                    AVERTYPE = EPISLEN
                ELSE
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not determine time period ' //
     &                     'of file ' // TRIM( CURFNM )
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
C.................  Assign and store logical file name
                WRITE( INTBUF,94030 ) N
                CURLNM = 'EMISFAC_' // ADJUSTL( INTBUF )
                EFLOGS( N ) = CURLNM

C.................  Set logical file name
                IF( .NOT. SETENVVAR( CURLNM, CURFNM ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not set logical file name ' //
     &                     'for file ' // CRLF() // BLANK10 // '"' //
     &                     TRIM( CURFNM ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

                USETIME( II,AVERTYPE ) = .TRUE.

C.................  Try to open file   
                IF( .NOT. OPENSET( CURLNM, FSREAD3, PROGNAME ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not open emission factors ' //
     &                     'file ' // CRLF() // BLANK10 // '"' //
     &                     TRIM( CURFNM ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Read file description
                IF( .NOT. DESCSET( CURLNM, ALLFILES ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not get description for ' // 
     &                     'file ' // CRLF() // BLANK10 // '"' //
     &                     TRIM( CURFNM ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
                EFSDATE = SDATE3D
                
C.................  Find end date in file description
                SEARCHSTR = '/END DATE/ '
                L = LEN_TRIM( SEARCHSTR ) + 1
                ENDFLAG = .FALSE.
                
                DO I = 1, MXDESC3
                   IF( INDEX( FDESC3D( I ), 
     &                        SEARCHSTR( 1:L ) ) > 0 ) THEN
                       TEMPLINE = FDESC3D( I )
                       IF( CHKINT( TEMPLINE( L+1:L+8 ) ) ) THEN
                           EFEDATE = STR2INT( TEMPLINE( L+1:L+8 ) )
                           EXIT
                       ELSE
                           ENDFLAG = .TRUE.
                           EXIT
                       END IF
                   END IF
                   
                   IF( I == MXDESC3 ) THEN
                       ENDFLAG = .TRUE.
                   END IF
                END DO

                IF( ENDFLAG ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not get ending date of ' //
     &                     'file ' // CRLF() // BLANK10 // '"' //
     &                     TRIM( CURFNM ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Determine starting and ending positions in array
                STPOS  = SECSDIFF( EARLYDATE, 0, EFSDATE, 0 )/(24*3600)
                STPOS  = STPOS + 1
                ENDPOS = SECSDIFF( EARLYDATE, 0, EFEDATE, 0 )/(24*3600)
                ENDPOS = ENDPOS + 1

C.................  Make sure starting and ending positions are valid
                IF( STPOS < 1 ) THEN
                    IF( ENDPOS > 0 ) THEN
                        STPOS = 1
                    ELSE
                        CYCLE
                    END IF
                END IF 
                
                IF( ENDPOS > NDAYS ) THEN
                    IF( STPOS <= NDAYS ) THEN
                        ENDPOS = NDAYS
                    ELSE
                        CYCLE
                    END IF
                END IF

C.................  Store day info
                DO I = STPOS, ENDPOS

                    EFDAYS( II, I, AVERTYPE ) = N
                    
                END DO

C.................  Allocate memory for temporary source info
                IF( ALLOCATED( SRCS ) ) DEALLOCATE( SRCS )
                ALLOCATE( SRCS( NROWS3D ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SRCS', PROGNAME )
            
                SRCS = 0
            
C.................  Read source information
                IF( .NOT. READSET( CURLNM, 'SOURCES', ALLAYS3, 
     &                             ALLFILES, SDATE3D, STIME3D, 
     &                             SRCS ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not read SOURCES ' // 
     &                     'from file ' // CRLF() // BLANK10 // '"' //
     &                     TRIM( CURFNM ) // '".'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

C.................  Store source information
                DO S = 1, NROWS3D
                                    
                    IF( SRCS( S ) /= 0 ) THEN

C.........................  Make sure source number is valid
                        IF( SRCS( S ) < 1 .OR. SRCS( S ) > NSRC ) CYCLE
                        
C.........................  Skip sources that are outside the grid
                        IF( EFIDX( II,SRCS( S ) ) == -9 ) CYCLE
                    
                        EFIDX( II,SRCS( S ) ) = S
                        WRITE( EFTYPE( II,SRCS( S ) ),'(I1)' ) AVERTYPE

                    END IF

                END DO
                
C.................  Close current file
                IF( .NOT. CLOSESET( CURLNM ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not close file ' // 
     &                     TRIM( CURFNM )
                    CALL M3MESG( MESG )
                    CYCLE
                END IF
                
            END DO    ! end of NLINES
 
C.............  Exit if there was a problem with the emission factor files
            IF( EFLAG ) THEN
                MESG = 'Problem checking emission factor files'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
C.............  Make sure all days are covered
            DO I = DAILY, EPISLEN
                IF( USETIME( II,I ) .EQV. .TRUE. ) THEN 
                    IF( MINVAL( EFDAYS( II,1:NDAYS,I ) ) == 0 ) THEN
                        MESG = 'ERROR: Emission factor files do not ' //
     &                         'cover requested time period.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                END IF
            END DO

C............  Print warning for sources that don't have emission factors
           DO I = 1, NSRC
               IF( EFIDX( II,I ) == 0 ) THEN
!                   WRITE( MESG,94070 ) 'WARNING: No VMT or emission ' //
!     &                    'factors available for' // CRLF() // 
!     &                    BLANK10 // 'Region: ', IFIP( I ),
!     &                    ' SCC: ' // CSCC( I )
!                   CALL M3MESG( MESG )
                   EFIDX( II,I ) = -1
               END IF
           END DO

        END DO    ! end of loop for entire episode periods



C.........  Successful completion
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94030   FORMAT( I3 )

94070   FORMAT( A, I5, A )

        END SUBROUTINE CHKEMFAC
