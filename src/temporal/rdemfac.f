
        SUBROUTINE RDEMFAC( II, N, SRC, NGPS, NGST, EARLYDATE, JDATE,
     &                      JTIME, NDAYS, TMAX, TMIN, TZONE, EMFAC, 
     &                      NAMOUT ) 

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
C
C Pathname: $Source: /afs/isis/depts/cep/em
C smoke@unc.educ/apps/archive/smoke/smoke/src/temporal/rdemfac.f,v $
C Last updated: $Date$ 
C
C**************************************************************************

C.........  Modules for public variables

C.........  This module contains emission factor tables and related
        USE MODEMFAC, ONLY: EFDAYS, EFIDX, EFLIST, EFLOGS, EFTYPE,
     &                      TEMPEF, USETIME

C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: TZONES

C.........  This module is used for MOBILE6 setup information 
        USE MODMBSET, ONLY: DAILY, WEEKLY, MONTHLY, EPISLEN


        IMPLICIT NONE
 
C.........  INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions


C..........  EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(2)    CRLF
        INTEGER         SECSDIFF
        INTEGER         STR2INT
        LOGICAL         SETENVVAR

        EXTERNAL    CRLF, SECSDIFF, STR2INT, SETENVVAR

C.........  Subroutine Arguments        

        INTEGER     , INTENT( IN ) :: II             ! idex of episode time periods
        INTEGER     , INTENT( IN ) :: N              ! current no. of pol/emis-types groups
        INTEGER     , INTENT( IN ) :: SRC            ! no. of sources
        INTEGER     , INTENT( IN ) :: NGPS           ! no. of pol/emis-types groups
        INTEGER     , INTENT( IN ) :: NGST           ! no. of pols/emis-types per group 
        INTEGER     , INTENT( IN ) :: EARLYDATE      ! earliest starting date based on time zones
        INTEGER     , INTENT( IN ) :: JDATE          ! ending Julian date and time
        INTEGER     , INTENT( IN ) :: JTIME          ! ending Julian date and time
        INTEGER     , INTENT( IN ) :: NDAYS          ! no. days in episode
        INTEGER     , INTENT( IN ) :: TMAX           ! maximum time zone in inventory      
        INTEGER     , INTENT( IN ) :: TMIN           ! minimum time zone in inventory      
        INTEGER     , INTENT( IN ) :: TZONE          ! output-file time zone
        REAL        , INTENT( IN OUT ) :: EMFAC ( SRC , NGST ) ! emission factors
        CHARACTER(*), INTENT( IN OUT ) :: NAMOUT( NGST, NGPS ) ! inv pol names

C.........  Reshaped input variables and output variables

C...........   Other local variables

        INTEGER         I, J, K, L, L1, S

        INTEGER         IOS                 ! i/o status
        INTEGER         FDATE, FTIME        ! emission factor date and time
        INTEGER         STPOS               ! starting position in ef day array
        REAL            RTMP                ! tmp float

        LOGICAL      :: EFLAG = .FALSE.  !  error-flag

        CHARACTER(IOVLEN3)   CBUF      ! pollutant name temporary buffer 
        CHARACTER(256)       CURFNM    ! current emission factor file name
        CHARACTER(16)        CURLNM    ! current ef logical file name
        CHARACTER(300)       MESG      ! buffer for M3EXIT() messages

        CHARACTER(16) :: PROGNAME = 'RDEMFAC' ! program name

C***********************************************************************
C   begin body of program RDEMFAC

C........  Loop through pollutants/emission-types in this group
       DO I = 1, NGST

           CBUF = NAMOUT( I,N )
           L1   = LEN_TRIM( CBUF )

C............  Skip blanks that can occur when NGRP > 1
           IF ( CBUF .EQ. ' ' ) CYCLE

C............  Check that this pollutant uses emission factors
C              Look for double underscore in pollutant name
           K = INDEX( CBUF, ETJOIN )

           IF( K == 0 ) THEN
               EMFAC( :,I ) = -1
               CYCLE
           END IF 

C............  Loop through time zones
           DO J = TMIN, TMAX
C................  Adjust time zone J based on output time zone and account
C                  for 6 AM starting time in files
               K = J - TZONE + 6

               FDATE = JDATE
               FTIME = JTIME

               CALL NEXTIME( FDATE, FTIME, -K * 10000 )

C................  Use date and time to find appropriate ef file
               STPOS = SECSDIFF( EARLYDATE, 0, FDATE, 0 )
               STPOS = STPOS / ( 24*3600 )
               STPOS = STPOS + 1

               DO L = DAILY, EPISLEN
               
                   IF( USETIME( II,L ) .EQV. .FALSE. ) CYCLE

                   IF( STPOS <= 0 .OR. STPOS > NDAYS ) THEN
                       MESG = 'ERROR: Invalid position'
                       CALL M3EXIT( PROGNAME, FDATE, FTIME, MESG, 2 )
                   END IF

                   CURFNM = EFLIST( EFDAYS( II,STPOS,L ) )
                   CURLNM = EFLOGS( EFDAYS( II,STPOS,L ) )

C....................  Set logical file name
                   IF( .NOT. SETENVVAR( CURLNM, CURFNM ) ) THEN
                       EFLAG = .TRUE.
                       MESG = 'ERROR: Could not set ' //
     &                        'logical file name for ' //
     &                        'file ' // CRLF() // BLANK10
     &                        // '"' // TRIM( CURFNM ) // 
     &                        '".'
                       CALL M3EXIT( PROGNAME, FDATE, FTIME, MESG, 2 )
                   END IF

C....................  Open current file
                   IF( .NOT. OPENSET( CURLNM, FSREAD3, 
     &                                PROGNAME ) ) THEN
                       EFLAG = .TRUE.
                       MESG = 'ERROR: Could not open ' //
     &                        'emission factors file ' //
     &                        CRLF() // BLANK10 // '"' //
     &                        TRIM( CURFNM ) // '".'
                       CALL M3EXIT( PROGNAME, FDATE, FTIME, MESG, 2 )
                   END IF

C....................  Read file description
                   IF( .NOT. DESCSET( CURLNM, 
     &                                ALLFILES ) ) THEN
                       MESG = 'ERROR: Could not get ' //
     &                        'description for file ' //
     &                        CRLF() // BLANK10 // '"' // 
     &                        TRIM( CURFNM ) // '".'
                       CALL M3EXIT( PROGNAME, FDATE, FTIME, MESG, 2 )
                   END IF

C....................  Read emission factors from current file
                   IF( .NOT. READSET( CURLNM, CBUF, ALLAYS3, 
     &                              ALLFILES, SDATE3D, 
     &                              FTIME, TEMPEF ) ) THEN
                       EFLAG = .TRUE.
                       MESG = 'Error reading "'// 
     &                        CBUF(1:L1) //
     &                        '" from file ' // 
     &                        CRLF() // BLANK10 // '"' // 
     &                        TRIM( CURFNM ) // '."'
                       CALL M3EXIT( PROGNAME, FDATE, FTIME, MESG, 2 )
                   END IF

C....................  Store emission factors by source                            
                   DO S = 1, SRC

C........................  Skip sources that are outside the grid or
C                          don't use emission factors
                       IF( EFIDX( II,S ) == -9 .OR. 
     &                     EFIDX( II,S ) == -1 ) THEN
                           EMFAC( S,I ) = 0
                           CYCLE
                       END IF

                       IF( TZONES( S ) == J .AND. 
     &                     STR2INT(EFTYPE( II,S )) == L ) THEN
                           EMFAC( S,I ) = TEMPEF( EFIDX( II,S ) )
                       END IF
                   END DO   ! end source loop

               END DO   ! end time period loop

           END DO   ! end time zone loop

C............  If there are any missing values in the data, give an
C              error to avoid problems in genhemis routine
           RTMP = MINVAL( EMFAC( 1:SRC,I ) )
           IF( RTMP == IMISS3 ) THEN
               EFLAG = .TRUE.
               MESG = 'ERROR: Missing emission ' //
     &            'factors(s) for "'// CBUF( 1:L1 ) // '".'
               CALL M3MSG2( MESG )
           END IF

       END DO  ! End loop on pollutants/emission-types I in this group

C........  Abort if error found
       IF( EFLAG ) THEN
           MESG = 'Problem with emission factors.'
           CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
       END IF

C.........  Successful completion
       RETURN


       END SUBROUTINE RDEMFAC
