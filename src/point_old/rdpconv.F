
C Version "@(#)$Id$ $Source$ $Date$ 

        SUBROUTINE RDCONV( KDEV, FIP1, SCC1, FAC1, N1,
     &                      SCC2, FAC2, N2 )

C***********************************************************************
C  function body starts at line  77
C
C  FUNCTION:
C       Reads the indicated volatile organic conversion file,
C       sorts it, and returns the number of profile entries found.
C       In case of overflow, returns negative of number of profile
C       entries found.
C
C  PRECONDITIONS REQUIRED:
C       volatile conversion file already opened on unit KDEV
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       prototype 6/96 by MRH
C
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'CHDIMS3.EXT'     ! emis chem info (both inventory and model)
      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C...........   ARGUMENTS and their descriptions:

        INTEGER     KDEV           !  unit number for speciation profile file
        INTEGER     N1             !  number of entries in group 1
        INTEGER     N2             !  number of entries in group 2
        INTEGER     FIP1( MXSCNV ) ! sorted state codes in group 1
        INTEGER     SCC1( MXSCNV ) ! sorted SCC codes in group 1
        INTEGER     SCC2( MXSCNV ) ! sorted SCC codes in group 2
        REAL        FAC1( MXSCNV ) ! sorted conversion factors in group 1
        REAL        FAC2( MXSCNV ) ! sorted conversion factors in group 2

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     STR2INT
        REAL        STR2REAL

        EXTERNAL    STR2INT, STR2REAL

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         FIP              ! Temporary FIP code
        INTEGER         I, I1, I2, J     ! Counters
        INTEGER         IOS              ! I/O Status code
        INTEGER         IREC             ! Record number of input file
        INTEGER         NTOT             ! Total items read in
        INTEGER         SCC              ! Temporary SCC code
        INTEGER         SID              ! Temporary state ID

        INTEGER         FIPA( MXSCNV )   ! unsorted FIPS
        INTEGER         INDX( MXSCNV )   ! index for sorting
        INTEGER         SCCA( MXSCNV )   ! unsorted SCCS

        REAL            FAC              ! Temporary conversion factor
        REAL            FACA( MXSCNV )   ! unsorted conversion factors

        LOGICAL         EFLAG            ! Error flag

        CHARACTER*128   LINE
        CHARACTER*256   MESG

C***********************************************************************
C   begin body of function  RDCONV

C.........  Read volatile pollutants conversion factors file
        IREC  =  0
        I     =  0
        EFLAG = .FALSE.
77      CONTINUE        !  head of the KDEV-read loop

            IREC = IREC + 1
            READ( KDEV, 93010, END=99, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading VOLATILE CONVERSION line', IREC
                CALL M3MESG( MESG )

            ELSE

                IF( LINE( 1:2 ) .EQ. ' .' ) THEN
                    SID = 0
                ELSE
                    SID = STR2INT( LINE( 1:2 ) )
                    IF( SID .LT. 0 ) SID = 0
                ENDIF

                SCC = STR2INT( LINE( 3:17 ) )
                IF( SCC .LT. 0 ) SCC = 0

                FAC = STR2REAL( LINE( 19:23 ) )
                IF( FAC .LT. 0. ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94020 )
     &                'Bad factor at line', IREC,
     &                'of VOLATILE CONVERSION FACTOR file =',
     &                FAC
                    GO TO 77  ! to head of read loop
                ENDIF

                I = I + 1
                IF ( I .LE. MXSCNV ) THEN

                    INDX( I ) = I
                    FIPA( I ) = SID * 1000
                    SCCA( I ) = SCC
                    FACA( I ) = FAC

                ELSE

                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                  'Max CONVERSION table size', MXSCNV, 
     &                  'exceeded at line', IREC
                    CALL M3MESG( MESG )

                END IF

            END IF

            GO TO  77   !  to head of loop

99      CONTINUE        !  end of the KDEV-read loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RDCONV', 0, 0, 
     &           'Error reading VOLATILE CONVERSION file', 2 )

        ENDIF

        NTOT = I
                
        CALL SORTI2( NTOT, INDX, FIPA, SCCA )

        I1 =  0
        I2 =  0
        DO 101 I = 1, NTOT

            J   = INDX( I )
            FIP = FIPA( J )
            SCC = SCCA( J )
            FAC = FACA( J )

            IF( FIP .GT. 0 ) THEN 

                I1 = I1 + 1
                FIP1( I1 ) = FIP
                SCC1( I1 ) = SCC
                FAC1( I1 ) = FAC

            ELSE

                I2 = I2 + 1
                SCC2( I2 ) = SCC
                FAC2( I2 ) = FAC

            END IF

101     CONTINUE

        N1 = I1
        N2 = I2
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xx

93010   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

94020   FORMAT( A, I8, 2X, 10 ( A, :, E10.3, :, 1X ) )


        END

