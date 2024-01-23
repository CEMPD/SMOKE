        SUBROUTINE SETOUTDATE( G_SDATE, G_STIME, G_NSTEPS,
     &                         NFILE, SDATE, STIME, NSTEP,
     &                         NAMES, MRGDIFF, USEFIRST )
     
C***********************************************************************
C  subroutine SETOUTDATE body starts at line 90
C
C  DESCRIPTION:
C      This subroutine determines the output date, time, and number
C      of time steps for the programs Mrggrid and Mrgelev.
C
C  PRECONDITIONS REQUIRED:
C      Arrays containing start date, time, and time steps for each file
C      If merging different days, already retrieved environment variables
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      M3MESG, M3MSG2, M3EXIT
C      NEXTIME
C      SECSDIFF
C      SEC2TIME
C
C  REVISION  HISTORY:
C      Created 5/2005 by C. Seppanen based on code from mrggrid.f
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

        IMPLICIT NONE

C.........  INCLUDES
        INCLUDE 'EMCNST3.EXT'   ! emissions constant parameters
        
C.........  EXTERNAL FUNCTIONS
        CHARACTER(2) CRLF
        INTEGER      SECSDIFF
        INTEGER      SEC2TIME
        
        EXTERNAL     CRLF, SECSDIFF, SEC2TIME
        
C.........  SUBROUTINE ARGUMENTS
        INTEGER,        INTENT(INOUT) :: G_SDATE           ! start date of output
        INTEGER,        INTENT(INOUT) :: G_STIME           ! start time of output
        INTEGER,        INTENT(INOUT) :: G_NSTEPS          ! number of time steps in output
        INTEGER,        INTENT(IN)    :: NFILE             ! number of input files
        INTEGER,        INTENT(IN)    :: SDATE( NFILE )    ! start dates of files
        INTEGER,        INTENT(IN)    :: STIME( NFILE )    ! start times of files
        INTEGER,        INTENT(IN)    :: NSTEP( NFILE )    ! no. of time steps in files
        CHARACTER(256), INTENT(IN)    :: NAMES( NFILE )    ! file names (may be logical or physical)
        LOGICAL,        INTENT(IN)    :: MRGDIFF           ! true: allow different days to be merged
        LOGICAL,        INTENT(OUT)   :: USEFIRST( NFILE ) ! true: use first time step of file
        
C.........  LOCAL VARIABLES
        INTEGER          I                        ! counter
        INTEGER          SECS                     ! number of seconds
        INTEGER          SECSMAX                  ! max. number of seconds
        INTEGER          SECSMIN                  ! min. number of seconds
        INTEGER          TMPDATE                  ! temporary date
        INTEGER          TMPTIME                  ! temporary time
        INTEGER          TMPSTEP                  ! temporary no. of time steps
        
        LOGICAL          ISLOGNAM                 ! true: file names are logical names
        LOGICAL       :: EFLAG = .FALSE.          ! true: an error has happened

        CHARACTER(512)   MESG                     ! message buffer
        CHARACTER(256)   NAM                      ! current file name
        
        CHARACTER(16) :: PROGNAME = 'SETOUTDATE'  ! program name

C***********************************************************************
C   begin body of subroutine SETOUTDATE

C.........  Check if we have logical or physical file names; this will
C           determine how error messages are printed
        IF( LEN_TRIM( NAMES( 1 ) ) <= 16 ) THEN
            ISLOGNAM = .TRUE.
        ELSE
            ISLOGNAM = .FALSE.
        END IF

C.........  If different days can be merged, then check file consistency
C           and confirm that environment settings for start date, start time,
C           and number of time steps match files
        IF( MRGDIFF ) THEN
        
C.............  Check that all files have the same start time
            TMPTIME = STIME( 1 )
            
            DO I = 2, NFILE
                NAM = NAMES( I )
                
                IF( STIME( I ) /= TMPTIME ) THEN
                    EFLAG = .TRUE.

                    WRITE( MESG,94010 ) 'ERROR: Start time ', 
     &                  STIME( I ), 'in file'

                    IF( ISLOGNAM ) THEN
                        MESG = TRIM( MESG ) // '"' // TRIM( NAM ) // '"'
                    ELSE
                        MESG = TRIM( MESG ) // CRLF() // BLANK10 // 
     &                         TRIM( NAM )
                    END IF

                    WRITE( MESG,94010 ) TRIM( MESG ) // CRLF() // 
     &                  BLANK10 // 'is inconsistent with the first ' //
     &                  'file value of ', TMPTIME

                    CALL M3MSG2( MESG )
                END IF
            END DO

            IF( EFLAG ) THEN
                MESG = 'Inconsistent start time in input files'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Check that environment settings are consistent with files
            IF( TMPTIME /= G_STIME ) THEN
                WRITE( MESG,94010 ) 'ERROR: Value for G_STTIME ',
     &              G_STIME, 'is inconsistent with the start' //
     &              CRLF() // BLANK10 // 'time of the files', TMPTIME
                CALL M3MSG2( MESG )
                
                MESG = 'Inconsistent environment settings'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            DO I = 1, NFILE
                NAM = NAMES( I )
                
C.................  Check that all files are at least long enough to
C                   to cover requested output duration
                IF( NSTEP( I ) < G_NSTEPS ) THEN
                    EFLAG = .TRUE.
                    
                    WRITE( MESG,94010 ) 'ERROR: Number of time steps ',
     &                  NSTEP( I ), 'in file'
                    
                    IF( ISLOGNAM ) THEN
                        MESG = TRIM( MESG ) // '"' // TRIM( NAM ) // '"'
                    ELSE
                        MESG = TRIM( MESG ) // CRLF() // BLANK10 //
     &                         TRIM( NAM )
                    END IF
                    
                    WRITE( MESG,94010 ) TRIM( MESG ) // CRLF() //
     &                  BLANK10 // 'is insufficient to cover the ' //
     &                  'requested number of output' // CRLF() //
     &                  BLANK10 // 'time steps ', G_NSTEPS
                    
                    CALL M3MSG2( MESG )
                END IF

C.................  Check if file contains output start date and 
C                   enough data to cover output duration
                IF( .NOT. EFLAG ) THEN
                    SECS = SECSDIFF( SDATE( I ), STIME( I ),
     &                               G_SDATE, G_STIME )
                    TMPSTEP = SECS / 3600
                    
                    IF( TMPSTEP > 0 .AND. TMPSTEP < NSTEP( I ) ) THEN
                        IF( NSTEP( I ) - TMPSTEP >= G_NSTEPS ) THEN
                            USEFIRST( I ) = .FALSE.
                        ELSE
                            USEFIRST( I ) = .TRUE.
                            
                            MESG = 'WARNING: Input file'
                            
                            IF( ISLOGNAM ) THEN
                                MESG = TRIM( MESG ) // '"' //
     &                                 TRIM( NAM ) // '"'
                            ELSE
                                MESG = TRIM( MESG ) // CRLF() //
     &                                 BLANK10 // TRIM( NAM )
                            END IF
                            
                            WRITE( MESG,94010 ) TRIM( MESG ) //
     &                        CRLF() // BLANK10 // 'contains the ' //
     &                        'requested output start date ', G_SDATE,
     &                        'and' // 
     &                        CRLF() // BLANK10 // 'start time ', 
     &                        G_STIME, 'but does not contain ' //
     &                        'enough data to ' //
     &                        CRLF() // BLANK10 // 'cover the ' //
     &                        'requested number of time steps',
     &                        G_NSTEPS, '.'
                            CALL M3MSG2( MESG )
                            
                            WRITE( MESG,94010 ) BLANK5 // 'Data ' //
     &                        'from the start date of the file ',
     &                        SDATE( I ), 'will be used instead.'
                            CALL M3MSG2( MESG )
                        END IF
                    ELSE

C.........................  File doesn't contain output start date and time
                        USEFIRST( I ) = .TRUE.
                    END IF
                END IF
            END DO
            
            IF( EFLAG ) THEN
                MESG = 'Problem with duration of input files'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.........  Otherwise, only merge input data from same day; set start date
C           and time as latest start date and time in input files and number
C           of time steps based on shortest input file
        ELSE

C.............  Generate reference date to compare against files (first day
C               of year)
            TMPDATE = ( SDATE( 1 )/1000 ) * 1000 + 1
            
C.............  Find maximum difference between files and reference date
            SECSMAX = 0
            DO I = 1, NFILE
                SECS = SECSDIFF( TMPDATE, STIME( 1 ), 
     &                           SDATE( I ), STIME( I ) )
                SECSMAX = MAX( SECSMAX, SECS )
            END DO
            
C.............  Calculate latest start date and time of all files
            TMPTIME = SEC2TIME( SECSMAX )
            
            G_SDATE = TMPDATE
            G_STIME = STIME( 1 )
            CALL NEXTIME( G_SDATE, G_STIME, TMPTIME )

C.............  Initialize number of output time steps with longest
C               possible duration
            TMPDATE = G_SDATE
            TMPTIME = G_STIME
            CALL NEXTIME( TMPDATE, TMPTIME, MAXVAL( NSTEP ) * 10000 )
            SECSMIN = SECSDIFF( G_SDATE, G_STIME, TMPDATE, TMPTIME )

C.............  Find minimum duration based on already determined start
C               date and time
            DO I = 1, NFILE
                TMPDATE = SDATE( I )
                TMPTIME = STIME( I )
                CALL NEXTIME( TMPDATE, TMPTIME, NSTEP( I ) * 10000 )
                SECS = SECSDIFF( G_SDATE, G_STIME, TMPDATE, TMPTIME )
                SECSMIN = MIN( SECSMIN, SECS )
            END DO

C.............  Calculate shortest number of time steps            
            TMPTIME = SEC2TIME( SECSMIN )
            
            G_NSTEPS = TMPTIME / 10000

            IF( G_NSTEPS <= 0 ) THEN
                MESG = 'ERROR: At least two input files contain ' //
     &                 'data for time periods that do' //
     &                 CRLF() // BLANK10 // 'not overlap. No output ' //
     &                 'file can be created.'
                CALL M3MSG2( MESG )
                
                MESG = 'Problem with input file time periods'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

        WRITE( MESG,94010 ) 'NOTE: Output file will start on ', 
     &      G_SDATE, 'at time ', G_STIME, 'and will contain' //
     &      CRLF() // BLANK10 // 'data for ', G_NSTEPS, 'time steps.'
        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE SETOUTDATE
