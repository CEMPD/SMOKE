
        SUBROUTINE WRIDAPOL( CATEGORY, VNAM, NSRC, VCNT, POLALL, FDEV, 
     &                       STATUS )

C***********************************************************************
C  subroutine body starts at line 76
C
C  DESCRIPTION:
C      Write inventory pollutant-specific data for variables listed in 
C      subroutine arguments to a temporary file. This subroutine opens the
C      temporary files as well.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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
C***************************************************************************

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS:
        INTEGER         GETEFILE
        EXTERNAL        GETEFILE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*), INTENT (IN) :: CATEGORY            ! src category
        CHARACTER*(*), INTENT (IN) :: VNAM                ! variable name
        INTEGER      , INTENT (IN) :: NSRC                ! no. of sources
        INTEGER      , INTENT (IN) :: VCNT                ! no. vars to write
        REAL         , INTENT (IN) :: POLALL( NSRC,VCNT ) ! var-specific data
        INTEGER      , INTENT(OUT) :: FDEV                ! unit number
        INTEGER      , INTENT(OUT) :: STATUS              ! exit status

C...........   Other local variables

        INTEGER          I, S     ! counters and indices

        INTEGER          IOS      ! i/o status
        INTEGER, SAVE :: LP       ! length of PATHNM

        LOGICAL       :: FIRSTIME = .TRUE.  ! true: first time routine is called

        CHARACTER*300          MESG           ! tmp message buffer
        CHARACTER*300, SAVE :: PATHNM = ' '   ! tmp path name
        CHARACTER*340          FILENM         ! tmp file name (with path)

        CHARACTER*16 :: PROGNAME = 'WRIDAPOL' ! program name

C***********************************************************************
C   begin body of subroutine WRIDAPOL

C.........  First time routine is called...
        IF( FIRSTIME ) THEN

C.............  Get directory for temporary files
            MESG = 'Path for temporary files will be written'
            CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )
            LP = LEN_TRIM( PATHNM )

            FIRSTIME = .FALSE.

        END IF

C.........  Initialize exit status
        STATUS = 0

C.........  Open temprary file
        FILENM = PATHNM( 1:LP ) // '/grwinven_tmp_' // VNAM
        FDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )
        
C.........  Append current pollutant to output buffer, w/ correct format
        SELECT CASE ( CATEGORY )
        CASE( 'AREA' )

            DO S = 1, NSRC

                WRITE( FDEV, 93200, ERR=999 ) 
     &               ( POLALL( S,I ), I = 1, VCNT )

            END DO ! End loop over sources

        CASE( 'MOBILE' )

            DO S = 1, NSRC

                WRITE( FDEV, 93210, ERR=999 ) 
     &               ( POLALL( S,I ), I = 1, VCNT )

            END DO ! End loop over sources

        CASE( 'POINT' ) 

            DO S = 1, NSRC

                WRITE( FDEV, 93220, ERR=999 ) 
     &               ( POLALL( S,I ), I = 1, VCNT )

            END DO ! End loop over sources

        END SELECT

C.........  Normal completion of subroutine

        RETURN

C.........  Exit with errors
999     MESG = 'ERROR: Problem writing to temporary file ' //
     &         'in subroutine ' // PROGNAME
        CALL M3MSG2( MESG )

        STATUS = 1

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   File read/write formats............ 93xxx

C......... NOTE - this three formats need to be consistent with wridaout.f
93200   FORMAT( 1X, E13.6, 1X, E13.6, 1X, E13.6, 1X, E13.6,          ! area
     &          1X, F4.0, 1X, F4.0 )

93210   FORMAT( 1X, E20.13, 1X, E20.13 )                             ! mobile

93220   FORMAT( 1X, E13.6, 1X, E13.6, 1X, E13.6, 1X, F4.0, 1X, E13.6, ! point
     &          1X, E10.2, 1X, E10.3 )

        END

