
        SUBROUTINE RDMVADJ( FDEV, MAX1, MAX2, NCOLS, NROWS, N1, N2,
     &                      OUTFIP, OUTRDT, OUTCEL, FAC1, FAC2  )

C***********************************************************************
C  subroutine body starts at line 127
C
C  DESCRIPTION:
C       Reads the mobile sources adjustments file that contains the information
C       on how to adjust the gridding matrix (not ungridding) entries so
C       that the user can adjust mobile emissions by grid cell, county, state,
C       road class, or urban/rural road class designation.
C
C  PRECONDITIONS REQUIRED:
C       File FDEV opened. Maximum county-based and cell-based arrays 
C       previously set. Memory allocated for output arrays.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       I/O API functions and subroutines
C
C  REVISION  HISTORY:
C       Written 4/98 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'MBDIMS3.EXT'     ! mobile-specific information

C...........   FUNCTION ARGUMENTS and their descriptions:

        INTEGER     FDEV           !  unit number for mobile adjustments file
        INTEGER     MAX1           !  max for FIP/RoadClass entries
        INTEGER     MAX2           !  max for cell-based entries
        INTEGER     NCOLS          !  number of columns
        INTEGER     NROWS          !  number of rows
        INTEGER     N1             !  actual number of FIP/RoadClass entries
        INTEGER     N2             !  actual number of cell-based entries
        INTEGER     OUTFIP( MAX1 ) !  sorted FIP codes 
        INTEGER     OUTRDT( MAX1 ) !  sorted road class codes 
        INTEGER     OUTCEL( MAX2 ) !  sorted cell numbers 
        REAL        FAC1  ( MAX1 ) !  factors OUTFIP and OUTRDT
        REAL        FAC2  ( MAX2 ) !  factors for OUTCEL

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2 CRLF
        INTEGER     FIND1
        INTEGER     STR2INT
        INTEGER     TRIMLEN

        EXTERNAL    CRLF, FIND1, STR2INT, TRIMLEN

C...........   LOCAL PARAMETERS:

        CHARACTER*5  BLANK5
        CHARACTER*10 BLANK10
        PARAMETER  ( BLANK5  = ' '
     &             , BLANK10 = ' ' )

C...........   Unsorted tables and indices for FIP/RoadClass-based factors
        INTEGER         INDX1( MXADJ1 )
        INTEGER         IFIPA( MXADJ1 )
        INTEGER         IRDTA( MXADJ1 )

        REAL            FAC1A( MXADJ1 )

C...........   Unsorted tables and indices for cell-based factors
        INTEGER         INDX2( MXADJ2 )
        INTEGER         ICELA( MXADJ2 )

        REAL            FAC2A( MXADJ2 )

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I,  J            ! Counters and indices

        INTEGER         FIP              ! Temporary FIP code
        INTEGER         IOS              ! I/O Status code
        INTEGER         IREC             ! Record number of input file
        INTEGER         ISTART, IEND     ! start/end index of MROADS3 for expndg
        INTEGER         IX, IY           ! tmp x-cell and y-cell
        INTEGER         JX, JY           ! position of '-' in cell range string
        INTEGER         LX, LY           ! length of cell range string
        INTEGER         XBEG, XEND       ! starting and ending X-cells in range
        INTEGER         YBEG, YEND       ! starting and ending Y-cells in range

        INTEGER         RDT              ! Temporary road class

        REAL            FAC              ! Temporary adjustment factor

        LOGICAL         CFLAG            ! True: records is cell-type
        LOGICAL         DFLAG            ! True: default record
        LOGICAL         EFLAG            ! Error flag
        LOGICAL         WFLAG            ! True: default FIP/RDT entry found

        CHARACTER*1     URI              ! urban or rural indicator
        CHARACTER*25    CELGRPX          ! x-cell range
        CHARACTER*25    CELGRPY          ! y-cell range
        CHARACTER*256   MESG

C***********************************************************************
C   begin body of subroutine RDMVADJ

        EFLAG = .FALSE.
        WFLAG = .FALSE.

C.........  Make sure that calling maxima are the same as compiled
        IF( MAX1 .GT. MXADJ1 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Calling maximum ' //
     &             'for table 1 =', MAX1, 'but compiled (MXADJ1)=',
     &             MXADJ1
            CALL M3MSG2( MESG )

        END IF

        IF( MAX2 .GT. MXADJ2 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Calling maximum ' //
     &             'for table 2 =', MAX2, 'but compiled (MXADJ2)=',
     &             MXADJ2
            CALL M3MSG2( MESG )

        END IF

        IF( EFLAG ) THEN
            MESG = 'Problem calling routine RDMVADJ'
            CALL M3EXIT( 'RDMVADJ', 0, 0, MESG, 2 )
        END IF

C.........  Initialize variables before start of loop
        IREC  =  0
        N1    =  0
        N2    =  0

C.........  Loop through lines of mobile adjustments factors file
99      CONTINUE        !  head of the FDEV-read loop

            READ( FDEV, *, END=111, IOSTAT=IOS ) 
     &            FIP, RDT, URI, CELGRPX, CELGRPY, FAC

            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )
     &              'I/O error', IOS, 
     &              'reading MOBILE ADJUSTMENTS file at line', IREC,
     &              CRLF() // BLANK10 // 'Remember to use quotes ' //
     &              'around character strings!'
                CALL M3MESG( MESG )

            END IF

C.............  Filter out -9 ot 0 for cell group
            IF( CELGRPX .EQ. '-9' .OR. CELGRPX .EQ. '0' ) CELGRPX = ' '
            IF( CELGRPY .EQ. '-9' .OR. CELGRPY .EQ. '0' ) CELGRPY = ' '

C.............  Set default to FIP-based record
            CFLAG = .FALSE.

C.............  Check if record is a cell-based record
            IF( CELGRPX .NE. ' ' ) CFLAG = .TRUE.

C.............  Process cell-based record
            IF( CFLAG ) THEN

                JX = INDEX( CELGRPX, '-' )
                JY = INDEX( CELGRPY, '-' )
                LX = TRIMLEN( CELGRPX )
                LY = TRIMLEN( CELGRPY )

C.................  Convert range strings to integer for loop indices
                IF( JX .GT. 0 ) THEN
                    XBEG = STR2INT( CELGRPX(    1:JX-1 ) )
                    XEND = STR2INT( CELGRPX( JX+1:LX   ) )
                ELSE
                    XBEG = STR2INT( CELGRPX( 1:LX ) )
                    XEND = XBEG
                END IF
                
                IF( JY .GT. 0 ) THEN
                    YBEG = STR2INT( CELGRPY(    1:JY-1 ) )
                    YEND = STR2INT( CELGRPY( JY+1:LY   ) )
                ELSE
                    YBEG = STR2INT( CELGRPY( 1:LY ) )
                    YEND = YBEG
                END IF

C.................  If any of these are negative, then there was a problem
                IF( XBEG .LT. 0 .OR. XEND .LT. 0 .OR.
     &              YBEG .LT. 0 .OR. YEND. LT. 0      ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Could not read ' //
     &                     'coordinate ranges at line ', IREC
                    CALL M3MESG( MESG )
                    GO TO 99   ! to head of main read loop

                END IF

C.................  Ensure that the "ends" are > the "begs"
                IF( XBEG .GT. XEND ) THEN
                    I = XEND
                    XEND = XBEG
                    XBEG = I
                END IF

                IF( YBEG .GT. YEND ) THEN
                    I = YEND
                    YEND = YBEG
                    YBEG = I
                END IF

C.................  Loop through cells of range
               DO 88 IX = XBEG, XEND

                    DO 77 IY = YBEG, YEND

                        IF( IX .GT. NCOLS ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: X-cell out ' //
     &                       'of range ( 1,', NCOLS, ') at line', IREC
                            CALL M3MESG( MESG )
                            GO TO 99   ! to head of main read loop

                        ELSE IF( IY .GT. NROWS ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Y-cell out ' //
     &                       'of range ( 1,', NROWS, ') at line', IREC
                            CALL M3MESG( MESG )
                            GO TO 99   ! to head of main read loop

                        ELSE IF( IY .LE. 0 ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: X-cell ' //
     &                       'indicates record is for cell-based ' //
     &                       'factor, but' // CRLF() // BLANK10 // 
     &                       'Y-cell is less than 0 at line', IREC
                            CALL M3MESG( MESG )
                            GO TO 77    ! to end of Y loop

                        END IF 

                        N2 = N2 + 1

                        IF( N2 .LE. MXADJ2 ) THEN
                            INDX2( N2 ) = N2
                            ICELA( N2 ) = ( IY - 1 ) * NCOLS + IX
                            FAC2A( N2 ) = FAC
                        END IF

77                  END DO  ! end of Y loop
88              END DO      ! end of X loop

C.............  Otherwise, process FIP-based record
            ELSE

                DFLAG = .FALSE.   ! Initialize as non-default

                IF( FIP .LE. 0   .AND. RDT .LE. 0   .AND.
     &              URI .NE. 'U' .AND. URI .NE. 'R'      ) THEN

                    IF( WFLAG ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Multiple entries '// 
     &                     'with FIPS code and Road Class both ' //
     &                     'default at line', IREC
                        CALL M3MESG( MESG )
                        GO TO 99   ! to head of loop

                    ELSE
                        WFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'WARNING: X-cell ' // 
     &                     'indicates record is for FIP/RoadClass-' //
     &                     'based factor, but' // CRLF() // BLANK10 // 
     &                     'FIPS code and Road Class are defaults ' //
     &                     'at line', IREC
                        CALL M3MSG2( MESG )
                    END IF

                END IF

C.................  Set up loop for expanding road class records for urban or
C                   rural (if indicator is defined)
                IF( URI .EQ. 'U' ) THEN
                    ISTART = 7
                    IEND   = 12

                ELSE IF( URI .EQ. 'R' ) THEN
                    ISTART = 1
                    IEND   = 6

                ELSE IF( URI .NE. '-' ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Invalid value for ' //
     &                     'urban/rural indicator at line', IREC, 
     &                     CRLF() // BLANK10 // 'Valid values are ' //
     &                     '"U", "R", and "-" for none.' 
                    CALL M3MESG( MESG )
                    GO TO 99   ! to head of loop

C................. If not set, then set index by finding value in list
                ELSE

                    ISTART = FIND1( RDT, NRCLAS, MROADS3 )

                    IF( RDT .EQ. 0  ) THEN
                        DFLAG = .TRUE.
                        ISTART = 1
                        IEND   = 1

                    ELSE IF( ISTART .GT. 0 ) THEN
                        IEND   = ISTART

                    ELSE 
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Invalid road ' //
     &                         'class code value', RDT, 'at line', IREC 
                        CALL M3MESG( MESG )
                        GO TO 99   ! to head of loop

                    END IF

                END IF

                DO 81 J = ISTART, IEND

                    N1 = N1 + 1

                    IF( N1 .LE. MXADJ1 ) THEN
                        INDX1( N1 ) = N1
                        IFIPA( N1 ) = FIP
                        FAC1A( N1 ) = FAC
                        IF( DFLAG ) THEN
                            IRDTA( N1 ) = 0
                        ELSE
                            IRDTA( N1 ) = MROADS3( J )
                        END IF
                    END IF

81              CONTINUE  

            END IF

101         GO TO  99   !  to head of loop

111     CONTINUE        !  end of the FDEV-read loop

        IF( N1 .GT. 0 .AND. N1 .GT. MXADJ1 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 
     &             'Max mobile adjustments FIP/RoadClass table ' //
     &             'size (MXADJ1)=', MXADJ1, CRLF() // BLANK5 //
     &             'Actual number = ', N1
            CALL M3MSG2( MESG )
        END IF

        IF( N2 .GT. 0 .AND. N2 .GT. MXADJ2 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 
     &             'Max mobile adjustments cell table ' //
     &             'size (MXADJ2)=', MXADJ2, CRLF() // BLANK5 //
     &             'Actual number = ', N2
            CALL M3MSG2( MESG )
        END IF

        IF ( EFLAG ) THEN
            MESG = 'Problem reading mobile adjustments file.'
            CALL M3EXIT( 'RDMVADJ', 0, 0, MESG, 2 )
        ENDIF


C.........  Sort FIP/RoadClass-based table
        CALL SORTI2( N1, INDX1, IFIPA, IRDTA )

C.........  Sort Cell-based table
        CALL SORTI1( N2, INDX2, ICELA )

C.........  Store FIP/RoadClass-based sorted table for output
        DO 151 I = 1, N1

            J = INDX1( I )

            OUTFIP( I ) = IFIPA( J )
            OUTRDT( I ) = IRDTA( J )
            FAC1  ( I ) = FAC1A( J )

151     CONTINUE

C.........  Store Cell-based sorted table for output
        DO 181 I = 1, N2

            J = INDX2( I )

            OUTCEL( I ) = ICELA( J )
            FAC2  ( I ) = FAC2A( J )

181     CONTINUE

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

        END

