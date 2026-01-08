C copied by: mhouyoux
C origin: rdascc.F 3.2

        SUBROUTINE  RDASCC( ADEV, NDIM, NASC, ASCA7, ASCA3 )

C***********************************************************************
C  subroutine body starts at line  72
C
C  DESCRIPTION:
C       Reads formatted actual-ASC file.
C
C  PRECONDITIONS REQUIRED:
C       Actual-ASC file opened on unit ADEV.
C       Actual-ASC file is sorted, and formatted (I7,I3)
C
C  REVISION  HISTORY:
C       Prototype  12/96 by CJC for area-source submodel in SMOKE
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C****************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES:

C        INCLUDE 'PARMS3.EXT'      ! I/O API constants
C        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
C        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations


C...........   ARGUMENTS and their descriptions: actually-occurring ASC table

        INTEGER     ADEV          !  unit number for actual-ASC file
        INTEGER     NDIM          !  max dimensioned number of ASCs
        INTEGER     NASC          !  actual number of ASCs returned
        INTEGER     ASCA7( NDIM ) !  leading-7 digits
        INTEGER     ASCA3( NDIM ) !  trailing-3 digits


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         IREC            !  input line (record) number
        INTEGER         IOS             !  I/O Status
        INTEGER         I               !  loop counter
        INTEGER         ID7,  ID3
        INTEGER         LID7, LID3
        LOGICAL         EFLAG   !  input error flag
        CHARACTER(300)  MESG    !  message buffer for M3MESG() and M3EXIT()

        CHARACTER(16) :: PROGNAME = 'RDASCC' ! program name


C***********************************************************************
C   begin body of subroutine  RDASCC

        CALL M3MSG2( 'Reading ACTUAL ASCs file...' )

        IREC  =  0
        I     =  0
        NASC  =  0
        LID3  = -1
        LID7  = -1
        EFLAG = .FALSE.

11      CONTINUE                        !  head of the ADEV-read loop

            READ( ADEV, 93020, END=99, IOSTAT=IOS ) ID7, ID3

            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading ACTUAL ASC file at line', IREC
                CALL M3MESG( MESG )
                GO TO  11                   !  to head of loop

            ELSE IF ( ( LID7 .GT. ID7 ) .OR.
     &                ( LID7 .EQ. ID7 .AND. LID3 .GE. ID3 ) ) THEN 

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'ASC table out of order at line', IREC
                CALL M3MESG( MESG )
                GO TO  11

            END IF              !  if i/o error; else if out-of-order

            I = I + 1
            IF ( I .LE. NDIM ) THEN 

                ASCA7( I ) = ID7
                ASCA3( I ) = ID3

            END IF              !  if I in bounds
     
            LID7 = ID7
            LID3 = ID3
            GO TO  11                   !  to head of loop

99      CONTINUE                        !  end of the ADEV-read loop

        WRITE( MESG,94010 )
     &      'Dimensioned ASC TABLE size', NDIM, 
     &      'Actual size', I
        CALL M3MSG2( MESG )

        IF ( I .GT. NDIM ) THEN 
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'ACTUAL ASC table overflow', 2 )
        ELSE IF ( EFLAG ) THEN
            CALL M3EXIT( 'SPCAMAT', 0, 0, 
     &                   'Error reading ACTUAL ASC file.', 2 )
        END IF

        NASC = I

      RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93020   FORMAT( I7, I3 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )


        END

