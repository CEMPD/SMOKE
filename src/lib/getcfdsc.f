        CHARACTER(*) FUNCTION GETCFDSC( FILEINFO, KEY, REQUIRED )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION: 
C     Retreives a character string from the FDESC array
C
C  PRECONDITIONS REQUIRED:
C     
C
C  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
C
C  REVISION  HISTORY:
C       prototype 1/99 by M Houyoux
C
C***********************************************************************
C
C Project Title: EDSS Tools Library
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
C****************************************************************************
        USE M3UTILIO

        IMPLICIT NONE

C...........   Include files
        INCLUDE 'IOCNST3.EXT'   !  I/O API constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   ARGUMENTS and their descriptions:

        CHARACTER(*), INTENT (IN) :: FILEINFO( * ) ! Array of file information
        CHARACTER(*), INTENT (IN) :: KEY           ! Key to find in FILEINFO
        LOGICAL     , INTENT (IN) :: REQUIRED      ! true: key must be found

C...........   EXTERNAL FUNCTIONS:
C       CHARACTER(2) CRLF
C       INTEGER      STR2INT

C        EXTERNAL     CRLF, STR2INT

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       L1, L2       ! length of strings
        INTEGER       I            ! loop index
        INTEGER       K            ! description string position of key
        INTEGER       LENGTH       ! length of function

        CHARACTER(300) CVAL        ! temporary output buffer
        CHARACTER(300) BUFFER      ! key buffer
        CHARACTER(300) MESG        ! message buffer

        CHARACTER(16) :: PROGNAME = 'GETCFDSC'    ! Program name

C***********************************************************************
C   begin body of subroutine GETCFDSC

C.........  Ensure left-justified keyword
        L1 = LEN_TRIM( KEY )
        BUFFER = ADJUSTL( KEY( 1:L1 ) )

        L1 = LEN_TRIM( BUFFER )

C.........  Get length of function
        LENGTH = LEN( GETCFDSC )

        DO I = 1, MXDESC3

            K = INDEX( FILEINFO( I ), BUFFER( 1:L1 ) )

            IF( K .GT. 0 ) THEN
                
                L2 = MAX( K+L1, LEN_TRIM( FILEINFO( I ) ) )
                CVAL = FILEINFO( I )( K+L1:L2 )

                L2 = LEN_TRIM( CVAL )

                IF( L2 .GT. LENGTH ) THEN
                    WRITE( MESG,94010 )
     &                     'INTERNAL WARNING: Length of string used ' //
     &                     'to call function "' //
     &                     PROGNAME( 1:LEN_TRIM( PROGNAME ) ) // '"' //
     &                     CRLF() // BLANK16 //
     &                     'is not long enough. Needs to be at least',
     &                     LENGTH, '. Trimming FDESC3D entry.'
                    CALL M3MSG2( MESG )

                END IF

                GETCFDSC = ADJUSTL( CVAL( 1:LENGTH ) )
                RETURN

            END IF

        END DO

C.........  If we get here, then key was not found in FDESC, so if it was
C           required, then abort.

        IF( REQUIRED ) THEN
            MESG = 'FDESC3D packet "' // KEY( 1:L1 ) // 
     &             '" was not found in NetCDF file!'
            CALL M3EXIT( MESG, 0, 0, PROGNAME, 2 )

        ELSE
            GETCFDSC = ' '
            RETURN

        END IF
    
C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END
