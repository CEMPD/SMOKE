
        SUBROUTINE OPENCTMP( PKTTYP, IDEV )

C***********************************************************************
C  subroutine body starts at line 75
C
C  DESCRIPTION:
C      This subroutine opens the temporary files that will contain the
C      indices to the control packet data tables
C
C      
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C
C***************************************************************************
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
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER(2)  CRLF
        INTEGER       GETEFILE
        INTEGER       JUNIT

        EXTERNAL      CRLF, GETEFILE, JUNIT

C...........   SUBROUTINE ARGUMENTS:

        CHARACTER(*), INTENT (IN) :: PKTTYP   ! packet type
        INTEGER     , INTENT(OUT) :: IDEV     ! logical file number


C...........   Other local variables

        INTEGER          IOS                   ! i/o status
        INTEGER, SAVE :: LP                    ! length of tmp file path

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time subroutine called

        CHARACTER(200), SAVE :: PATHNM         ! path name for tmp file
        CHARACTER(220)          FILENM                ! file name
        CHARACTER(256)          MESG                  ! message buffer
        CHARACTER(16)        :: PROGNAME = 'OPENCTMP' ! program name

C***********************************************************************
C   Begin body of subroutine OPENCTMP

        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

            MESG = 'Path where temporary control files will be written'
            CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )
            LP = LEN_TRIM( PATHNM )

            IF( IOS .NE. 0 ) THEN
                MESG = 'WARNING: Large temporary files being placed '//
     &                 'executable directory because ' // CRLF() //
     &                 BLANK10 // 'environment variable SMK_TMPPATH '//
     &                 'is not set properly'
                CALL M3MSG2( MESG )
            END IF

        END IF

        SELECT CASE( PKTTYP )

           CASE ( 'CTG' )

              FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_ctg'
              IDEV = JUNIT()
              OPEN( IDEV, ERR=1006, FILE=FILENM )

           CASE ( 'CONTROL', 'EMS_CONTROL' )

              FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_ctl'
              IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

           CASE ( 'ALLOWABLE' )

              FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_alw'
              IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

           CASE ( 'PROJECTION' )

              FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_proj'
              IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )
              
           CASE ( 'MACT' )
           
              FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_mact'
              IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

        END SELECT

        RETURN

C.........  Open file error
1006    MESG = 'Could not open temporary ASCII file for control info.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENCTMP
