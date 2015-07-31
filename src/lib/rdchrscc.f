
        SUBROUTINE RDCHRSCC( FDEV, NSCC, SCCLIST )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C       Reads character SCC codes
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C       Written  1/99 by M. Houyoux
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

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'

C...........   ARGUMENTS and their descriptions: actually-occurring ASC table

        INTEGER      FDEV             !  unit number for actual-SCC file
        INTEGER      NSCC             !  number of lines in actual-SCC file
        CHARACTER(*) SCCLIST( NSCC )  !  list of SCCs

C...........   SCRATCH LOCAL VARIABLES and their descriptions:
  
        INTEGER         I                 !  counters and indices
        INTEGER         IOS               !  I/O Status
        INTEGER         IREC              !  input line (record) number

        LOGICAL      :: EFLAG = .FALSE.   !  error flag

        CHARACTER(7)    FMTSCC        !  read format for SCCs
        CHARACTER(300)  MESG          !  message buffer
        CHARACTER(SCCLEN3) SBUFFA     !  read SCC buffer

        CHARACTER(16) :: PROGNAME = 'RDCHRSCC' ! program name

C***********************************************************************
C   begin body of subroutine  RDCHRSCC

C.........  Write SCC read format
        WRITE( FMTSCC, '(A,I2.2,A)' ) '(A', SCCLEN3, ')'

        IREC = 0
        DO I = 1, NSCC    !  head of the FDEV-read loop

            READ( FDEV, FMTSCC, END=999, IOSTAT=IOS ) SBUFFA
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading SCC CODES file at line', IREC
                CALL M3MESG( MESG )
                CYCLE   !  to end of loop

            END IF      !  if i/o error; else if out-of-order

            SCCLIST( I ) = SBUFFA           

        ENDDO       !  end of the FDEV-read loop

        IF( EFLAG ) THEN

            MESG = 'Problem reading actual SCC file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

        REWIND( FDEV )

        RETURN

C.........  Special exit from loop for end of file

999     WRITE( MESG,94010 ) 'Unexpected end-of-file at line', IREC,
     &                      'of actual SCC file.'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( I7, I3 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

        END SUBROUTINE RDCHRSCC
