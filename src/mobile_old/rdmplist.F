
        INTEGER FUNCTION RDMPLIST( FDEV, MPFIP, MPRDT, MPLNK, 
     &                             MPREF, REFCNT, REFLIST ) 

C***********************************************************************
C  function body starts at line  108
C
C  DESCRIPTION:  
C       read and sort file MPLIST for TMPMOBIL
C
C  PRECONDITIONS REQUIRED:  MPLIST opened on unit FDEV, and positioned
C       at the start of the file
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       SORTI3
C
C  REVISION  HISTORY:
C       Prototype  1/96 by CJC
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

        INCLUDE 'CHDGEN3.EXT'     ! generalized speciation parms
        INCLUDE 'MBDIMS3.EXT'     ! mobile-source parameters
        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations


C...........   ARGUMENTS and their descriptions:

        INTEGER   FDEV
        INTEGER   MPFIP( MXPLST )
        INTEGER   MPRDT( MXPLST )
        INTEGER   MPLNK( MXPLST )
        INTEGER   MPREF( 0:23, MXPLST )    !  hour&source-specific xrefs
        INTEGER   REFCNT                   !  count of total unique PSIs
        INTEGER   REFLIST( MXPSI )         !  list of unique PSIs


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER   FIND1
        INTEGER   STR2INT
        INTEGER   TRIMLEN

        EXTERNAL  FIND1, STR2INT, TRIMLEN

C...........   PARAMETERS and their descriptions

        INTEGER    MXFLDS         ! Maximum number of columns in MPLIST

        PARAMETER( MXFLDS = 27 )

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER   I, J, K, L, P, L1, S
        INTEGER   IOS
        INTEGER	  FIP, RDT, LNK, REF( 0:23 )
        INTEGER   IREC
        INTEGER   LENGTH    ! Length of LINE, not including trailing blanks
        INTEGER   MPCNT
        INTEGER   MPDEXA( MXPLST )
        INTEGER   MPFIPA( MXPLST )
        INTEGER   MPRDTA( MXPLST )
        INTEGER   MPLNKA( MXPLST )
        INTEGER   MPREFA( 0:23, MXPLST )   !  hour&source-specific xrefs
        INTEGER   NREF      ! Number of PSIs in a segment of LINE
        INTEGER   NCOL      ! Number of column in LINE
        INTEGER   SBPP
        INTEGER   SB( MXFLDS )  ! Beginning position in LINE of each segment
        INTEGER   SE( MXFLDS )  ! Ending    position in LINE of each segment

        LOGICAL	  EFLAG

        CHARACTER*1     THISCHAR  ! Current character in LINE
        CHARACTER*1     PREVCHAR  ! Previous character in LINE
        CHARACTER*256	LINE      ! Line of input from MPLIST
        CHARACTER*256	MESG

C***********************************************************************
C   begin body of subroutine  RDMPLIST

        MPCNT  = 0
        REFCNT = 0
        IREC   = 0
        I      = 0
        EFLAG  = .FALSE.

11      CONTINUE	!  head of loop reading FDEV

            READ( FDEV, 93000, END = 99, IOSTAT=IOS ) LINE

            IREC = IREC + 1

C.............  Check for line read errors

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 ) 
     &               'I/O error in MPLIST at line', IREC,
     &               'IOSTAT=', IOS

                CALL M3EXIT( 'RDMPLIST', 0, 0, MESG, 2 )

            END IF

C.............  Parse first segment of LINE to ensure input using ' '

            LENGTH   = TRIMLEN( LINE ) 
            PREVCHAR = ' '
            S        = 1
            L1       = 1

            DO 15 P = 1, LENGTH

                THISCHAR = LINE( P:P )

C.................  Check for error in '*' based columns
                IF( ( THISCHAR .EQ. '*' .AND. PREVCHAR .EQ. ' ' ) .OR.
     &              ( THISCHAR .EQ. ' ' .AND. PREVCHAR .EQ. '*' ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                   'An illegal " " is adjoins a "*" at line', 
     &                   IREC, 'of the MPLIST file'
                    CALL M3MESG( MESG )
                    GO TO 11          ! To head of read loop

C.................  Set begining of column S
                ELSEIF( THISCHAR .NE. ' ' .AND. PREVCHAR .EQ. ' ' ) THEN
                   L1 = P

C.................  Store segment beginning and end & increment segment counter
                ELSEIF( THISCHAR .EQ. ' ' .AND. PREVCHAR .NE. ' ' ) THEN

                    IF( S .LE. MXFLDS ) THEN

                        SB( S ) = L1
                        SE( S ) = P - 1

                    ELSE

                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                   'More than', MXFLDS, 'columns at line', IREC,
     &                   'of MPLIST file'
                        CALL M3MESG( MESG )
                        GO TO 11          ! To head of read loop

                    ENDIF

                    S = S + 1 
    
                ENDIF

                PREVCHAR = THISCHAR

15          CONTINUE

C.............  Get last column 
            SB( S ) = L1
            SE( S ) = LENGTH

C.............  Check for valid number of columns

            IF( S .EQ. 1 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &               'Columns are not separated by a " " at line', IREC,
     &               'of MPLIST file'
                CALL M3MESG( MESG )
                GO TO 11          ! To head of read loop

            ELSEIF( S .LT. 3 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &               'Incomplete record at line', IREC,
     &               'of MPLIST file'
                CALL M3MESG( MESG )
                GO TO 11          ! To head of read loop

            ELSE

                NCOL = S

            ENDIF

C.............  Process segments of LINE into known parts

            FIP = STR2INT( LINE( SB( 1 ): SE( 1 ) ) )
            RDT = STR2INT( LINE( SB( 2 ): SE( 2 ) ) )
            LNK = STR2INT( LINE( SB( 3 ): SE( 3 ) ) )

            J = -1             !  J starts at 0
            DO 33 S = 4, NCOL  !  head of segment processing loop

                P = INDEX( LINE( SB( S ): SE( S ) ), '*' )

C.................  For only one PSI given
                IF( P .LE. 0 ) THEN 
             
                    J        = J + 1
                    REF( J ) = STR2INT( LINE( SB( S ): SE( S ) ) )
 
C.................  For multiple PSIs given
                ELSE

                    SBPP = SB( S ) + P
                    NREF = STR2INT( LINE( SB( S ): SBPP - 2 ) )

                    DO 31 L = 1, NREF
                        J        = J + 1
                        REF( J ) = STR2INT( LINE( SBPP: SE( S ) ) )
31                  CONTINUE

                ENDIF

33          CONTINUE

C.............  Check for too many or too few records in line
            IF( J .NE. 23 ) THEN  

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &               'Need PSIs for 24 hrs at line', IREC,
     &               'of MPLIST (instead of', J, ').'
                CALL M3MESG( MESG )
                GO TO 11         ! To head of file read loop

            ENDIF

            I = I + 1

            IF ( I .LE. MXPLST ) THEN

                MPDEXA( I ) = I
                MPFIPA( I ) = FIP
                MPRDTA( I ) = RDT
                MPLNKA( I ) = LNK

                DO  43  J = 0, 23

                    MPREFA( J,I ) = REF( J )
   
C.....................  If haven't encounterd this PSI before, insert in
C.....................  sorted list 
                    IF( FIND1( REF( J ), REFCNT, REFLIST ) .LE. 0 ) THEN

                        REFCNT = REFCNT + 1
                        IF( REFCNT.EQ.1 ) THEN
                            REFLIST( 1 ) = REF( J )
                        ELSE
                            DO  41  K = REFCNT, 2, -1
                                IF ( REFLIST( K-1 ) .GT. REF( J ) ) THEN
                                     REFLIST( K ) = REFLIST( K-1 )
                                ELSE
                                     REFLIST( K ) = REF( J )
                                     GO TO  42
                                END IF
41                          CONTINUE
                            REFLIST( 1 ) = REF( J )
                        END IF
42                      CONTINUE !  exit from insertion

                    END IF  	 !  if new hour index

43              CONTINUE	 !  end loop on hour indexes J

            END IF               !  end dimensioning check

        GO TO 11                 !  to head of reading FDEV

99      CONTINUE	         !  end of loop reading FDEV

        IF( I .GT. MXPLST ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &           'MPLIST entries -- dimensioned (MXPLST):', MXPLST,
     &           ', actual:', I
            CALL M3MSG2( MESG )
        ENDIF

        IF ( EFLAG ) THEN
            RDMPLIST = -I
            RETURN
        ELSE
            RDMPLIST = I
        END IF


C.......   Now sort the entries:

        CALL SORTI3( I, MPDEXA, MPFIPA, MPRDTA, MPLNKA )

        DO  122  J = 1, I
            K = MPDEXA( J )
            MPFIP( J ) = MPFIPA( K )
            MPRDT( J ) = MPRDTA( K )
            MPLNK( J ) = MPLNKA( K )
            DO  111  L = 0, 23
                MPREF( L,J ) = MPREFA( L,K )
111         CONTINUE
122     CONTINUE

        RETURN

C******************  FORMAT  STATEMENTS   ******************************
C...........   I/O formats........................... 93xxx

93000	FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010	FORMAT( 10( A, :, I7, :, 2X ) )

        END

