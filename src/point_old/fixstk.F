
        SUBROUTINE  FIXSTK( FDEV, NSRC, IFIP, ISCC, IPLT, ISTK,
     &                      STKHT, STKDM, STKTK, STKVE, NPOL, EMISV )

C***********************************************************************
C  subroutine body starts at line 159
C
C  FUNCTION:
C	Use replacement stack parameters from file PSTK to fill inn
C	stack parameters which are "missing" (i.e., negative).
C
C  PRECONDITIONS REQUIRED:
C	Correctly set logical name for the PSTK file
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C	PROMPTFFILE
C
C  REVISION  HISTORY:
C	prototype 12/95 by CJC
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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

        INCLUDE 'PTDIMS3.EXT'      ! point source dimensioning constants
        INCLUDE 'PARMS3.EXT'       ! 
        INCLUDE 'FDESC3.EXT'       ! 
        INCLUDE 'IODECL3.EXT'      ! 


C...........   ARGUMENTS and their descriptions:

        INTEGER FDEV	        !  unit number for stack parameter file PSTK
        INTEGER NSRC            !  actual number of sources
        INTEGER IFIP ( NSRC )	!  FIP codes
        INTEGER ISCC ( NSRC )	!  SCC codes
        INTEGER IPLT ( NSRC )	!  Plant ID codes (for writing message only)
        INTEGER ISTK ( NSRC )	!  Stack ID codes (for writing message only)
        INTEGER NPOL            !  actual number of pollutants
        REAL	STKHT( NSRC )	!  stack height (m)
        REAL	STKDM( NSRC )	!  stack diameter (m)
        REAL	STKTK( NSRC )	!  stack exhaust temperature (K)
        REAL	STKVE( NSRC )	!  stack exhaust velocity (m/s)
        REAL    EMISV( NSRC, NPOL )  ! emissions


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

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        INTEGER		FIND1, FIND2
        EXTERNAL	CRLF, FIND1, FIND2


C...........   LOCAL PARAMETERS and their descriptions:
        CHARACTER*5     BLANK5
        INTEGER         MXFPSC
        INTEGER         MXSTSC

        PARAMETER(      BLANK5 = ' ',
     &                  MXFPSC = NPFIP * NPSCC,
     &                  MXSTSC = NPSID * NPSCC  ) 

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        REAL		HT0	!  ultimate fallback height
        REAL		DM0	!  ultimate fallback diameter
        REAL		TK0	!  ultimate fallback temperature
        REAL		VE0	!  ultimate fallback velocity

        INTEGER		NR1		!  size of SCC-only table
        INTEGER		SC1( NPSCC )	!  SCC code
        REAL		HT1( NPSCC )	!  SCC-only height
        REAL		DM1( NPSCC )	!  SCC-only diameter
        REAL		TK1( NPSCC )	!  SCC-only temperature
        REAL		VE1( NPSCC )	!  SCC-only velocity

        INTEGER		NR2	        !  size of SCC-state table
        INTEGER		FP2( MXSTSC )   !  FIP state code
        INTEGER		SC2( MXSTSC )   !  SCC code
        REAL		HT2( MXSTSC )   !  SCC-state height
        REAL		DM2( MXSTSC )   !  SCC-state diameter
        REAL		TK2( MXSTSC )   !  SCC-state temperature
        REAL		VE2( MXSTSC )   !  SCC-state velocity

        INTEGER		NR3		!  size of FIP-SCC table
        INTEGER		FP3( MXFPSC )	!  FIP code
        INTEGER		SC3( MXFPSC )	!  SCC code
        REAL		HT3( MXFPSC )	!  FIP-SCC height
        REAL		DM3( MXFPSC )	!  FIP-SCC diameter
        REAL		TK3( MXFPSC )	!  FIP-SCC temperature
        REAL		VE3( MXFPSC )	!  FIP-SCC velocity
        
        INTEGER		I, S, K	!  source subscript
        INTEGER		IOS	!  I/O error status
        INTEGER		LINE	!  current line number
        
        REAL		HT	!  temporary height
        REAL		DM	!  temporary diameter
        REAL		TK	!  temporary exit temperature
        REAL		VE	!  temporary velocity

        INTEGER		LDEV	!  log file unit number
        INTEGER		LFIP	!  previous FIPs code
        INTEGER		LSCC	!  previous SCC code
        INTEGER		FIP	!  temporary FIPs code
        INTEGER		SCC	!  temporary SCC code
        INTEGER		SID	!  temporary state ID

        LOGICAL		EFLAG   !  error flag
        LOGICAL		DFLAG( NPSRC ) ! true if source getting default parms 
        CHARACTER*256	MESG	!  error-message buffer

C***********************************************************************
C   begin body of subroutine  FIXSTK

        LDEV = INIT3()  ! Need for message writing

        CALL M3MSG2( 'Reading default stack parameters...' )

C.......   First, read the (first and therefore) ultimate fallback record:

        READ( FDEV,*, IOSTAT=IOS )  FIP, SCC, HT0, DM0, TK0, VE0

        IF ( IOS .NE. 0 ) THEN
             CALL M3EXIT( 'FIXSTK', 0, 0,
     &              'Error reading PSTK at line 1', 1 )
        ELSE IF ( FIP .NE. 0  .OR.  SCC .NE. 0 ) THEN
             CALL M3EXIT( 'FIXSTK', 0, 0,
     &              'No fallback record in PSTK', 2 )
        END IF

C.......   Now read the rest of the file:

        EFLAG = .FALSE.
        LINE  = 1
        NR1   = 0
        NR2   = 0
        NR3   = 0

11      CONTINUE        !  head of input loop

            LFIP = FIP
            LSCC = SCC
            LINE = LINE + 1

            READ( FDEV,*, END=22, IOSTAT=IOS )  FIP, SCC, HT, DM, TK, VE

            IF ( IOS .NE. 0 ) THEN	!  I/O error

                WRITE( MESG,94010 ) 'Error reading PSTK at line', LINE
                CALL M3EXIT( 'FIXSTK', 0, 0, MESG, 1 )

            ELSE IF ( LFIP .GT. FIP  .OR.
     &                ( LFIP .EQ. FIP  .AND.  
     &                  LSCC .GT. SCC ) ) THEN	!  out of order

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'PSTK out of order at line', LINE
                CALL M3MESG( MESG )

            ELSE IF ( FIP .EQ. 0 ) THEN		!  SCC only

                NR1 = NR1 + 1
                IF ( NR1 .LE. NPSCC ) THEN
                    SC1( NR1 ) = FIP
                    HT1( NR1 ) = HT
                    DM1( NR1 ) = DM
                    TK1( NR1 ) = TK
                    VE1( NR1 ) = VE
                END IF

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN	!  state and SCC

                NR2 = NR2 + 1
                IF ( NR2 .LE. MXSTSC ) THEN
                    FP2( NR2 ) = FIP / 1000
                    SC2( NR2 ) = SCC
                    HT2( NR2 ) = HT
                    DM2( NR2 ) = DM
                    TK2( NR2 ) = TK
                    VE2( NR2 ) = VE
                END IF

            ELSE					!  FIP and SCC

                NR3 = NR3 + 1
                IF ( NR3 .LE. MXFPSC ) THEN
                    FP3( NR3 ) = FIP
                    SC3( NR3 ) = SCC
                    HT3( NR3 ) = HT
                    DM3( NR3 ) = DM
                    TK3( NR3 ) = TK
                    VE3( NR3 ) = VE
                END IF

            END IF	!  if I/O error, or out of order or ...

            GO TO  11

22      CONTINUE        !  end of input loop

C...........   Report dimensions of segmented stack parms

        WRITE( MESG,94010 )
     &      'Number of STACK FIX PARMS by SCC entries--' //
     &      CRLF() // BLANK5 // '   dimensioned (NPSCC):', NPSCC,
     &      'actual:', NR1
        CALL M3MSG2( MESG )
 
        WRITE( MESG,94010 )
     &      'Number of STACK FIX PARMS by STATE/SCC entries--' //
     &      CRLF() // BLANK5 // '   dimensioned (MXSTSC):', MXSTSC,
     &      'actual:', NR2
        CALL M3MSG2( MESG )
 
        WRITE( MESG,94010 )
     &      'Number of STACK FIX PARMS by FIP/SCC entries--' //
     &      CRLF() // BLANK5 // '   dimensioned (MXFPSC):', MXFPSC,
     &      'actual:', NR3
        CALL M3MSG2( MESG )
 
C...........   If there is an overflow, abort
        IF( NR1 .GT. NPSCC  .OR.
     &      NR2 .GT. MXSTSC .OR.
     &      NR3 .GT. MXFPSC      ) THEN
 
            CALL M3EXIT( 'FIXSTK', 0, 0,
     &                   'STACK FIX PARAMETERS table overflow', 2 )
 
        END IF

C.........  Bound stack parameters to minima and maxima values
C.........  This is in a separate loop to permit better reporting
C.........  Watch out for negative numbers or zeroes, because these are 
C.........  the missing stack parameters, which should get defaults.

        CALL M3MSG2( 'Bounding MIN and MAX stack parameters...' )

        DO 29 S = 1, NSRC

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
                WRITE( MESG,94050 ) 
     &                 'FIPS code:', IFIP( S ), 'SCC:', ISCC( S ),
     &                 'PLT:', IPLT( S ), 'STK:', ISTK( S ),
     &                 'EMIS:', ( EMISV( S,I ), I = 1, NPOL )

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

29      CONTINUE


        CALL M3MSG2( 'Fixing MISSING stack parameters...' )

C...........   Get LOG file unit, so can write to directly (using M3MESG would
C..........    add too many spaces
        LDEV = INIT3()

C...........   Now do replacements of MISSING stack parameters:
C...........   4 passes -- ht, dm, tk, ve
C...........   Treat parameters equal to 0 as missing

        DO  33  S = 1, NSRC

            K = 0                ! Initialize K to test if replacements made
            DFLAG( S ) = .FALSE. ! Initialize DFLAG to test if defaults used

            IF ( STKHT( S ) .LE. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF( K .LE. 0 ) K = FIND2( FIP, 0, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    HT = HT3( K )
                    DM = DM3( K )
                    TK = TK3( K )
                    VE = VE3( K )
                ELSE
                    SID = FIP/1000
                    K   = FIND2( SID, SCC, NR2, FP2, SC2 )

                    IF( K .LE. 0 ) K = FIND2( SID, 0, NR2, FP2, SC2 )

                    IF ( K .GT. 0 ) THEN
                        HT = HT2( K )
                        DM = DM2( K )
                        TK = TK2( K )
                        VE = VE2( K )
                    ELSE

                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            HT = HT1( K )
                            DM = DM1( K )
                            TK = TK1( K )
                            VE = VE1( K )
                        ELSE
                            DFLAG( S ) = .TRUE.

                        END IF 
                    END IF 
                END IF 

                IF( .NOT. DFLAG( S ) ) THEN

                    WRITE( MESG,94050 ) 
     &                'FIPS code:', FIP, 'SCC:', SCC,
     &                'PLT:', IPLT( S ), 'STK:', ISTK( S ),
     &                'EMIS:', ( EMISV( S,I ), I = 1, NPOL ), 
     &                CRLF() // BLANK5 //  '             Old        New'
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
                    DM = DM3( K )
                    TK = TK3( K )
                    VE = VE3( K )
                ELSE
                    SID = FIP/1000
                    K   = FIND2( SID, SCC, NR2, FP2, SC2 )

                    IF( K .LE. 0 ) K = FIND2( SID, 0, NR2, FP2, SC2 )

                    IF ( K .GT. 0 ) THEN
                        DM = DM2( K )
                        TK = TK2( K )
                        VE = VE2( K )
                    ELSE
                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            DM = DM1( K )
                            TK = TK1( K )
                            VE = VE1( K )
                        ELSE
                            DFLAG( S ) = .TRUE.

                        END IF 
                    END IF 
                END IF 

                IF( .NOT. DFLAG( S ) ) THEN

                    WRITE( MESG,94050 ) 
     &                'FIPS code:', FIP, 'SCC:', SCC,
     &                'PLT:', IPLT( S ), 'STK:', ISTK( S ),
     &                'EMIS:', ( EMISV( S,I ), I = 1, NPOL ), 
     &                CRLF() // BLANK5 //  '             Old        New'
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
                    TK = TK3( K )
                    VE = VE3( K )
                ELSE
                    SID = FIP/1000
                    K   = FIND2( SID, SCC, NR2, FP2, SC2 )

                    IF( K .LE. 0 ) K = FIND2( SID, 0, NR2, FP2, SC2 )

                    IF ( K .GT. 0 ) THEN
                        TK = TK2( K )
                        VE = VE2( K )
                    ELSE
                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            TK = TK1( K )
                            VE = VE1( K )
                        ELSE
                            DFLAG( S ) = .TRUE.

                        END IF 
                    END IF 
                END IF 

                IF( .NOT. DFLAG( S ) ) THEN

                    WRITE( MESG,94050 ) 
     &                'FIPS code:', FIP, 'SCC:', SCC,
     &                'PLT:', IPLT( S ), 'STK:', ISTK( S ),
     &                'EMIS:', ( EMISV( S,I ), I = 1, NPOL ), 
     &                CRLF() // BLANK5 //  '             Old        New'
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
                    VE = VE3( K )
                ELSE
                    SID = FIP/1000
                    K   = FIND2( FIP/1000, SCC, NR2, FP2, SC2 )

                    IF( K .LE. 0 ) K = FIND2( SID, 0, NR2, FP2, SC2 )

                    IF ( K .GT. 0 ) THEN
                        VE = VE2( K )
                    ELSE
                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            VE = VE1( K )
                        ELSE
                            DFLAG( S ) = .TRUE.

                        END IF 
                    END IF 
                END IF 

                IF( .NOT. DFLAG( S ) ) THEN

                    WRITE( MESG,94050 ) 
     &                'FIPS code:', FIP, 'SCC:', SCC,
     &                'PLT:', IPLT( S ), 'STK:', ISTK( S ),
     &                'EMIS:', ( EMISV( S,I ), I = 1, NPOL ), 
     &                CRLF() // BLANK5 //  '             Old        New'
                    CALL M3MESG( MESG )

                    WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                    STKVE( S ) = VE
                ENDIF

            END IF	!  if stack exhaust velocity bad

33      CONTINUE        !  end loop fixing missing stack parameters

C.........  Apply ultimate default parameters, and write report
C.........  This is in a separate loop to permit better reporting
     
        CALL M3MESG( 'Ultimate fallback stack parameters report:' )
        DO 44 S = 1, NSRC

            FIP = IFIP( S )
            SCC = ISCC( S )

            IF( DFLAG( S ) ) THEN

                WRITE( MESG,94050 ) 
     &                'FIPS code:', FIP, 'SCC:', SCC,
     &                'PLT:', IPLT( S ), 'STK:', ISTK( S ),
     &                'EMIS:', ( EMISV( S,I ), I = 1, NPOL ), 
     &                CRLF() // BLANK5 //  '             Old        New'
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

44      CONTINUE


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010	FORMAT( 10 ( A, :, I8, :, 1X ) )

94020   FORMAT( 7X, A, 2X, E10.3, 1X, E10.3 )

94030   FORMAT( 7X, A6, 1X, '> max.  Change from ', 
     &          E10.3, ' to ', E10.3 )

94040   FORMAT( 7X, A6, 1X, '< min.  Change from ', 
     &          E10.3, ' to ', E10.3 )

94050   FORMAT( A, I5.5, 1X, A, I8.8, 1X, A, I5, 1X, A, I5, : , 1X,
     &          A, <NPOL>( F9.2, 1X ), :, A )

        END

