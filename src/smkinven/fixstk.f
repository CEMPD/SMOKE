
        SUBROUTINE  FIXSTK( FDEV, NSRC, IFIP, ISCC, CSOURC,
     &                      STKHT, STKDM, STKTK, STKVE )

C***********************************************************************
C  subroutine body starts at line 157
C
C  DESCRIPTION:
C	Read and use replacement stack parameters from file PSTK to fill in
C	stack parameters which are "missing" (i.e., negative). Also use
C       ultimate defaults (set as local parameters) when no other stack 
C       parameters are available.
C
C  PRECONDITIONS REQUIRED:
C	Opened file with unit FDEV
C       Memory of arrays already allocated and with NSRC
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C	Subroutines: I/O API subroutines, CHECKMEM, FMTCSRC
C       Function: I/O API functions, GETFLINE
C
C  REVISION  HISTORY:
C	prototype 12/95 by CJC
C       copied by: mhouyoux
C       origin: fixstk.F 4.3
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
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        INTEGER		FIND1, FIND2
        INTEGER		GETFLINE

        EXTERNAL	CRLF, FIND1, FIND2, GETFLINE

C...........   ARGUMENTS and their descriptions:

        INTEGER       FDEV             ! unit number for stack parameter file PSTK
        INTEGER       NSRC             ! actual number of sources
        INTEGER       IFIP  ( NSRC )   ! FIP codes
        INTEGER       ISCC  ( NSRC )   ! SCC codes
        CHARACTER*(*) CSOURC( NSRC )   ! concat source chars
        REAL	      STKHT ( NSRC )   ! stack height (m)
        REAL	      STKDM ( NSRC )   ! stack diameter (m)
        REAL	      STKTK ( NSRC )   ! stack exhaust temperature (K)
        REAL	      STKVE ( NSRC )   ! stack exhaust velocity (m/s)

C...........   PARAMETERS and their descriptions:

        REAL       MINHT        ! Mininum stack height (m)
        REAL       MINDM        ! Mininum stack diameter (m)
        REAL       MINTK        ! Mininum stack exit temperature (K)
        REAL       MINVE        ! Mininum stack exit velocity (m/s)
        REAL       MAXHT        ! Maximum stack height (m)
        REAL       MAXDM        ! Maximum stack diameter (m)
        REAL       MAXTK        ! Maximum stack exit temperature (K)
        REAL       MAXVE        ! Maximum stack exit velocity (m/s)

        PARAMETER( MINHT = 0.5,
     &             MINDM = 0.01,
     &             MINTK = 260.,
     &             MINVE = 0.0001,
     &             MAXHT = 2100.,
     &             MAXDM = 100.,
     &             MAXTK = 2000.,
     &             MAXVE = 500.    )

C...........    LOCAL VARIABLES and their descriptions:

        INTEGER, ALLOCATABLE:: INDXA( : ) !  Sorting index
        INTEGER, ALLOCATABLE:: SFIPA( : ) !  Unsorted FIP state code from PSTK
        INTEGER, ALLOCATABLE:: SSCCA( : ) !  Unsorted SCC code from PSTK
        REAL   , ALLOCATABLE:: SHTA ( : ) !  Unsorted height from PSTK
        REAL   , ALLOCATABLE:: SDMA ( : ) !  Unsorted diameter from PSTK
        REAL   , ALLOCATABLE:: STKA ( : ) !  Unsorted temperature from PSTK
        REAL   , ALLOCATABLE:: SVEA ( : ) !  Unsorted velocity from PSTK

        REAL                   HT0      !  ultimate fallback height
        REAL                   DM0      !  ultimate fallback diameter
        REAL                   TK0      !  ultimate fallback temperature
        REAL                   VE0      !  ultimate fallback velocity

        INTEGER                NR1      !  size of SCC-only table
        INTEGER, ALLOCATABLE:: SC1( : ) !  SCC code
        INTEGER, ALLOCATABLE:: ID1( : ) !  Index to unsorted arrays from PSTK

        INTEGER                NR2      !  size of SCC-state table
        INTEGER, ALLOCATABLE:: FP2( : ) !  FIP state code
        INTEGER, ALLOCATABLE:: SC2( : ) !  SCC code
        INTEGER, ALLOCATABLE:: ID2( : ) !  Index to unsorted arrays from PSTK

        INTEGER                NR3      !  size of FIP-SCC table
        INTEGER, ALLOCATABLE:: FP3( : ) !  FIP code
        INTEGER, ALLOCATABLE:: SC3( : ) !  SCC code
        INTEGER, ALLOCATABLE:: ID3( : ) !  Index to unsorted arrays from PSTK
        
        REAL		HT	!  temporary height
        REAL		DM	!  temporary diameter
        REAL		TK	!  temporary exit temperature
        REAL		VE	!  temporary velocity

        INTEGER		FIP	!  temporary FIPs code
        INTEGER		I, J, S, K	!  source subscript
        INTEGER		IOS	!  I/O error status
        INTEGER		IREC	!  current record number
        INTEGER		L2	!  buffer length
        INTEGER		LDEV	!  log file unit number
        INTEGER		NLINE	!  Number of lines
        INTEGER		NPSTK	!  Number of PSTK entries
        INTEGER		SCC	!  temporary SCC code
        INTEGER		SID	!  temporary state ID

        LOGICAL		EFLAG   !  error flag
        LOGICAL		DFLAG( NSRC ) ! true if source getting default parms 
        DATA            EFLAG / .FALSE. /

        CHARACTER*300	BUFFER  !  temporary buffer
        CHARACTER*300	MESG    !  message buffer
        CHARACTER*16 :: PROGNAME = 'FIXSTK'  ! program name
     
C***********************************************************************
C   begin body of subroutine FIXSTK

C.........   Get LOG file unit, so can write to directly (using M3MESG would
C            add too many spaces for some messages)
        LDEV = INIT3()

        CALL M3MSG2( 'Reading default stack parameters...' )

C.........  Get dimensions of input file
        NLINE = GETFLINE( FDEV, 'Stack replacement file')

C.........  Allocate memory for arrays.  Since this file is not likely to be
C           large, allocate all arrays based on number of lines in the file.
        ALLOCATE( INDXA( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXA', PROGNAME )
        ALLOCATE( SFIPA( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SFIPA', PROGNAME )
        ALLOCATE( SSCCA( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SSCCA', PROGNAME )
        ALLOCATE( SHTA( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SHTA', PROGNAME )
        ALLOCATE( SDMA( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SDMA', PROGNAME )
        ALLOCATE( STKA( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'STKA', PROGNAME )
        ALLOCATE( SVEA( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SVEA', PROGNAME )
        ALLOCATE( SC1( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SC1', PROGNAME )
        ALLOCATE( ID1( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ID1', PROGNAME )
        ALLOCATE( FP2( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FP2', PROGNAME )
        ALLOCATE( SC2( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SC2', PROGNAME )
        ALLOCATE( ID2( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ID2', PROGNAME )
        ALLOCATE( FP3( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FP3', PROGNAME )
        ALLOCATE( SC3( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SC3', PROGNAME )
        ALLOCATE( ID3( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ID3', PROGNAME )

C.........  Read the PSTK file until hit the end of the file

        IREC  = 0
        I = 0
        DO        !  head of input loop

            READ( FDEV, *, END=22, IOSTAT=IOS ) FIP, SCC, HT, DM, TK, VE

            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN	!  I/O error

                WRITE( MESG,94010 ) 'Error', IOS, 
     &                              'reading PSTK at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE

                I = I + 1
                INDXA( I ) = I
                SFIPA( I ) = FIP
                SSCCA( I ) = SCC
                SHTA ( I ) = HT
                SDMA ( I ) = DM
                STKA ( I ) = TK
                SVEA ( I ) = VE

            ENDIF

        ENDDO

22      CONTINUE        !  end of input loop

        NPSTK = I

C.........  Sort PSTK data 
        CALL SORTI2( NPSTK, INDXA, SFIPA, SSCCA )

C.........  Disaggregate PSTK data into 4 categories
        NR1   = 0
        NR2   = 0
        NR3   = 0
        DO I = 1, NPSTK

            J   = INDXA( I )
            FIP = SFIPA( J )
            SCC = SSCCA( J ) 

            IF( FIP .EQ. 0 .AND. SCC .EQ. 0 ) THEN  ! fallback default
                HT0 = SHTA ( J )
                DM0 = SDMA ( J )
                TK0 = STKA ( J )
                VE0 = SVEA ( J )

            ELSEIF( FIP .EQ. 0 ) THEN               !  SCC only
                NR1 = NR1 + 1
                SC1( NR1 ) = SCC
                ID1( NR1 ) = J

            ELSE IF( MOD( FIP, 1000 ) .EQ. 0 ) THEN !  state and SCC
                NR2 = NR2 + 1
                FP2( NR2 ) = FIP / 1000
                SC2( NR2 ) = SCC
                ID2( NR2 ) = J

            ELSE                                    !  FIP and SCC
                NR3 = NR3 + 1
                FP3( NR3 ) = FIP
                SC3( NR3 ) = SCC
                ID3( NR3 ) = J

            END IF

        ENDDO ! End loop on NPSTK
 
C.........  Bound stack parameters to minima and maxima values
C.........  This is in a separate loop to permit better reporting
C.........  Watch out for negative numbers or zeroes, because these are 
C.........  the missing stack parameters, which should get defaults.

        CALL M3MSG2( 'Bounding MIN and MAX stack parameters...' )

        DO S = 1, NSRC

            HT = STKHT( S )
            DM = STKDM( S )
            TK = STKTK( S )
            VE = STKVE( S )

            IF ( HT .GT. MAXHT .OR.
     &         ( HT .LT. MINHT .AND. HT .GT. 0 ) .OR.
     &           DM .GT. MAXDM .OR.
     &         ( DM .LT. MINDM .AND. DM .GT. 0 ) .OR.
     &           TK .GT. MAXTK .OR.
     &         ( TK .LT. MINTK .AND. TK .GT. 0 ) .OR.
     &           VE .GT. MAXVE .OR.
     &         ( VE .LT. MINVE .AND. VE .GT. 0 ) ) THEN

                CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2 )
                WRITE( MESG,94010 ) BUFFER( 1:L2 ) // ' SCC: ', ISCC(S)
                CALL M3MESG( MESG )

            ENDIF

            IF ( HT .GT. MAXHT ) THEN
                WRITE( LDEV,94030 ) 'Height', HT, MAXHT
                HT = MAXHT

            ELSEIF( HT .LT. MINHT .AND. HT .GT. 0 ) THEN
                WRITE( LDEV,94040 ) 'Height', HT, MINHT
                HT = MINHT

            END IF

            IF ( DM .GT. MAXDM ) THEN
                WRITE( LDEV,94030 ) '  Diam', DM, MAXDM
                DM = MAXDM

            ELSEIF( DM .LT. MINDM .AND. DM .GT. 0 ) THEN
                WRITE( LDEV,94040 ) '  Diam', DM, MINDM
                DM = MINDM

            END IF

            IF ( TK .GT. MAXTK )THEN
                WRITE( LDEV,94030 ) '  Temp', TK, MAXTK
                TK = MAXTK

            ELSEIF( TK .LT. MINTK .AND. TK .GT. 0 ) THEN 
                WRITE( LDEV,94040 ) '  Temp', TK, MINTK
                TK = MINTK

            END IF

            IF ( VE .GT. MAXVE )THEN
                WRITE( LDEV,94030 ) ' Veloc', VE, MAXVE
                VE = MAXVE

            ELSEIF( VE .LT. MINVE .AND. VE .GT. 0 ) THEN
                WRITE( LDEV,94040 ) ' Veloc', VE, MINVE
                VE = MINVE
 
            END IF

            STKHT( S ) = HT
            STKDM( S ) = DM
            STKTK( S ) = TK
            STKVE( S ) = VE

        ENDDO ! Loop on sources

        CALL M3MSG2( 'Fixing MISSING stack parameters...' )

C...........   Now do replacements of MISSING stack parameters:
C...........   4 passes -- ht, dm, tk, ve
C...........   Treat parameters equal to 0 as missing

        DO S = 1, NSRC

            K = 0                ! Initialize K to test if replacements made
            DFLAG( S ) = .FALSE. ! Initialize DFLAG to test if defaults used

            IF ( STKHT( S ) .LE. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF( K .LE. 0 ) K = FIND2( FIP, 0, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    J  = ID3 ( K )
                    HT = SHTA( J )
                    DM = SDMA( J )
                    TK = STKA( J )
                    VE = SVEA( J )
                ELSE
                    SID = FIP/1000
                    K   = FIND2( SID, SCC, NR2, FP2, SC2 )

                    IF( K .LE. 0 ) K = FIND2( SID, 0, NR2, FP2, SC2 )

                    IF ( K .GT. 0 ) THEN
                        J  = ID2 ( K )
                        HT = SHTA( J )
                        DM = SDMA( J )
                        TK = STKA( J )
                        VE = SVEA( J )
                    ELSE

                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            J  = ID1 ( K )
                            HT = SHTA( J )
                            DM = SDMA( J )
                            TK = STKA( J )
                            VE = SVEA( J )
                        ELSE
                            DFLAG( S ) = .TRUE.

                        END IF 
                    END IF 
                END IF 

                IF( .NOT. DFLAG( S ) ) THEN

                    CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                    WRITE( MESG,94010 ) 
     &                     BUFFER( 1:L2 ) // ' SCC: ', ISCC( S ),
     &                     CRLF() // BLANK5 // 
     &                     '             Old        New'
                    CALL M3MESG( MESG )

                    WRITE( LDEV,94020 )     'Height', STKHT( S ), HT
                    STKHT( S ) = HT

                    IF ( STKDM( S ) .LE. 0 ) THEN
                        WRITE( LDEV,94020 ) '  Diam', STKDM( S ), DM
                        STKDM( S ) = DM
                    ENDIF

                    IF ( STKTK( S ) .LE. 0 ) THEN
                        WRITE( LDEV,94020 ) '  Temp', STKTK( S ), TK
                        STKTK( S ) = TK
                    ENDIF 

                    IF ( STKVE( S ) .LE. 0 ) THEN
                        WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                        STKVE( S ) = VE
                    ENDIF
                ENDIF

            END IF	!  if stack height bad

            IF ( STKDM( S ) .LE. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF( K .LE. 0 ) K = FIND2( FIP, 0, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    J  = ID3 ( K )
                    DM = SDMA( J )
                    TK = STKA( J )
                    VE = SVEA( J )
                ELSE
                    SID = FIP/1000
                    K   = FIND2( SID, SCC, NR2, FP2, SC2 )

                    IF( K .LE. 0 ) K = FIND2( SID, 0, NR2, FP2, SC2 )

                    IF ( K .GT. 0 ) THEN
                        J  = ID2 ( K )
                        DM = SDMA( J )
                        TK = STKA( J )
                        VE = SVEA( J )
                    ELSE
                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            J  = ID1 ( K )
                            DM = SDMA( J )
                            TK = STKA( J )
                            VE = SVEA( J )
                        ELSE
                            DFLAG( S ) = .TRUE.

                        END IF 
                    END IF 
                END IF 

                IF( .NOT. DFLAG( S ) ) THEN

                    CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                    WRITE( MESG,94010 ) 
     &                     BUFFER( 1:L2 ) // ' SCC: ', ISCC( S ),
     &                     CRLF() // BLANK5 // 
     &                     '             Old        New'
                    CALL M3MESG( MESG )

                    WRITE( LDEV,94020 ) '  Diam', STKDM( S ), DM
                    STKDM( S ) = DM

                    IF ( STKTK( S ) .LE. 0 ) THEN
                        WRITE( LDEV,94020 ) '  Temp', STKTK( S ), TK
                        STKTK( S ) = TK
                    ENDIF 

                    IF ( STKVE( S ) .LE. 0 ) THEN
                        WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                        STKVE( S ) = VE
                    ENDIF
                ENDIF

            END IF  	!  if stack diameter bad

            IF ( STKTK( S ) .LE. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF( K .LE. 0 ) K = FIND2( FIP, 0, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    J  = ID3 ( K )
                    TK = STKA( J )
                    VE = SVEA( J )
                ELSE
                    SID = FIP/1000
                    K   = FIND2( SID, SCC, NR2, FP2, SC2 )

                    IF( K .LE. 0 ) K = FIND2( SID, 0, NR2, FP2, SC2 )

                    IF ( K .GT. 0 ) THEN
                        J  = ID2 ( K )
                        TK = STKA( J )
                        VE = SVEA( J )
                    ELSE
                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            J  = ID1 ( K )
                            TK = STKA( J )
                            VE = SVEA( J )
                        ELSE
                            DFLAG( S ) = .TRUE.

                        END IF 
                    END IF 
                END IF 

                IF( .NOT. DFLAG( S ) ) THEN

                    CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                    WRITE( MESG,94010 ) 
     &                     BUFFER( 1:L2 ) // ' SCC: ', ISCC( S ),
     &                     CRLF() // BLANK5 // 
     &                     '             Old        New'
                    CALL M3MESG( MESG )

                    WRITE( LDEV,94020 ) '  Temp', STKTK( S ), TK
                    STKTK( S ) = TK

                    IF ( STKVE( S ) .LE. 0 ) THEN
                        WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                        STKVE( S ) = VE
                    ENDIF
                ENDIF

            END IF	!  if stack exhaust temperature bad

            IF ( STKVE( S ) .LE. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF( K .LE. 0 ) K = FIND2( FIP, 0, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    J  = ID3 ( K )
                    VE = SVEA( J )
                ELSE
                    SID = FIP/1000
                    K   = FIND2( FIP/1000, SCC, NR2, FP2, SC2 )

                    IF( K .LE. 0 ) K = FIND2( SID, 0, NR2, FP2, SC2 )

                    IF ( K .GT. 0 ) THEN
                        J  = ID2 ( K )
                        VE = SVEA( J )
                    ELSE
                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            J  = ID1 ( K )
                            VE = SVEA( J )
                        ELSE
                            DFLAG( S ) = .TRUE.

                        END IF 
                    END IF 
                END IF 

                IF( .NOT. DFLAG( S ) ) THEN

                    CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                    WRITE( MESG,94010 ) 
     &                     BUFFER( 1:L2 ) // ' SCC: ', ISCC( S ),
     &                     CRLF() // BLANK5 // 
     &                     '             Old        New'
                    CALL M3MESG( MESG )

                    WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                    STKVE( S ) = VE
                ENDIF

            END IF	!  if stack exhaust velocity bad

        ENDDO  !  end loop on sources for fixing missing stack parameters

C.........  Apply ultimate fallback parameters, and write report
C.........  This is in a separate loop to permit better reporting
     
        CALL M3MESG( 'Ultimate fallback stack parameters report:' )

        DO S = 1, NSRC

            IF( DFLAG( S ) ) THEN

                CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                WRITE( MESG,94010 ) 
     &                 BUFFER( 1:L2 ) // ' SCC: ', ISCC( S ),
     &                 CRLF() // BLANK5 // 
     &                 '             Old        New'
                CALL M3MESG( MESG )

                IF ( STKHT( S ) .LE. 0 ) THEN
                    WRITE( LDEV,94020 ) 'Height', STKHT( S ), HT0
                    STKHT( S ) = HT0
                ENDIF

                IF ( STKDM( S ) .LE. 0 ) THEN
                    WRITE( LDEV,94020 ) '  Diam', STKDM( S ), DM0
                    STKDM( S ) = DM0
                ENDIF

                IF ( STKTK( S ) .LE. 0 ) THEN
                    WRITE( LDEV,94020 ) '  Temp', STKTK( S ), TK0
                    STKTK( S ) = TK0
                ENDIF 

                IF ( STKVE( S ) .LE. 0 ) THEN
                    WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE0
                    STKVE( S ) = VE0
                ENDIF

            ENDIF

        ENDDO  ! Loop through sources for applying ultimate fallbacks

        DEALLOCATE( INDXA )
        DEALLOCATE( SFIPA )
        DEALLOCATE( SSCCA )
        DEALLOCATE( SHTA )
        DEALLOCATE( SDMA )
        DEALLOCATE( STKA )
        DEALLOCATE( SVEA )
        DEALLOCATE( SC1 )
        DEALLOCATE( ID1 )
        DEALLOCATE( FP2 )
        DEALLOCATE( SC2 )
        DEALLOCATE( ID2 )
        DEALLOCATE( FP3 )
        DEALLOCATE( SC3 )
        DEALLOCATE( ID3 )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010	FORMAT( 10 ( A, :, I8, :, 1X ) )

94020   FORMAT( 7X, A, 2X, E10.3, 1X, E10.3 )

94030   FORMAT( 7X, A6, 1X, '> max.  Change from ', 
     &          E10.3, ' to ', E10.3 )

94040   FORMAT( 7X, A6, 1X, '< min.  Change from ', 
     &          E10.3, ' to ', E10.3 )

        END

