C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C ****************************************************************************/

        SUBROUTINE  FIXSTK( FDEV, NSRC, IFIP, ISCC, 
     &                      STKHT, STKDM, STKTK, STKVE )

C***********************************************************************
C  subroutine body starts at line 120
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

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PTDIMS3.EXT'      ! point source dimensioning constants


C...........   ARGUMENTS and their descriptions:

        INTEGER FDEV	        !  unit number for stack parameter file PSTK
        INTEGER NSRC            !  actual number of sources
        INTEGER IFIP ( NSRC )	!  FIP codes
        INTEGER ISCC ( NSRC )	!  SCC codes
        REAL	STKHT( NSRC )	!  stack height (m)
        REAL	STKDM( NSRC )	!  stack diameter (m)
        REAL	STKTK( NSRC )	!  stack exhaust temperature (K)
        REAL	STKVE( NSRC )	!  stack exhaust velocity (m/s)


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER		FIND1, FIND2
        EXTERNAL	FIND1, FIND2


C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER		S, K	!  source subscript
        INTEGER		IOS	!  I/O error status
        INTEGER		LINE	!  current line number
        
        INTEGER		FIP	!  contents of input record
        INTEGER		SCC	!  contents of input record
        INTEGER		LFIP	!  contents of input record
        INTEGER		LSCC	!  contents of input record
        REAL		HT	!  contents of input record
        REAL		DM	!  contents of input record
        REAL		TK	!  contents of input record
        REAL		VE	!  contents of input record

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
        INTEGER		FP2( NPSRC )	!  FIP state code
        INTEGER		SC2( NPSRC )	!  SCC code
        REAL		HT2( NPSRC )	!  SCC-state height
        REAL		DM2( NPSRC )	!  SCC-state diameter
        REAL		TK2( NPSRC )	!  SCC-state temperature
        REAL		VE2( NPSRC )	!  SCC-state velocity

        INTEGER		NR3		!  size of FIP-SCC table
        INTEGER		FP3( NPSRC )	!  FIP code
        INTEGER		SC3( NPSRC )	!  SCC code
        REAL		HT3( NPSRC )	!  FIP-SCC height
        REAL		DM3( NPSRC )	!  FIP-SCC diameter
        REAL		TK3( NPSRC )	!  FIP-SCC temperature
        REAL		VE3( NPSRC )	!  FIP-SCC velocity
        
        LOGICAL		EFLAG   !  error flag
        CHARACTER*160	MESG	!  error-message buffer

C***********************************************************************
C   begin body of subroutine  FIXSTK

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
                ELSE
                    WRITE( MESG,94010 ) 
     &              'SCC-only table overflow at line', LINE
                    CALL M3MESG( MESG )
                END IF

            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN	!  state and SCC

                NR2 = NR2 + 1
                IF ( NR2 .LE. NPSRC ) THEN
                    FP2( NR2 ) = FIP / 1000
                    SC2( NR2 ) = SCC
                    HT2( NR2 ) = HT
                    DM2( NR2 ) = DM
                    TK2( NR2 ) = TK
                    VE2( NR2 ) = VE
                ELSE
                    WRITE( MESG,94010 ) 
     &              'SCC-state table overflow at line', LINE
                    CALL M3MESG( MESG )
                END IF

            ELSE					!  FIP and SCC

                NR3 = NR3 + 1
                IF ( NR3 .LE. NPSRC ) THEN
                    FP3( NR3 ) = FIP
                    SC3( NR3 ) = SCC
                    HT3( NR3 ) = HT
                    DM3( NR3 ) = DM
                    TK3( NR3 ) = TK
                    VE3( NR3 ) = VE
                ELSE
                    WRITE( MESG,94010 ) 
     &              'FIP-SCC table overflow at line', LINE
                    CALL M3MESG( MESG )
                END IF

            END IF	!  if I/O error, or out of order or ...

            GO TO  11

22      CONTINUE        !  end of input loop


C...........   Now do replacements of missing stack parameters:
C...........   4 passes -- ht, dm, tk, ve

        DO  33  S = 1, NSRC

            IF ( STKHT( S ) .LT. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    HT = HT3( K )
                    DM = DM3( K )
                    TK = TK3( K )
                    VE = VE3( K )
                ELSE
                    K = FIND2( FIP/1000, SCC, NR2, FP2, SC2 )
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
                            HT = HT0
                            DM = DM0
                            TK = TK0
                            VE = VE0
                        END IF 
                    END IF 
                END IF 

                STKHT( S ) = HT
                IF ( STKDM( S ) .LT. 0 ) STKDM( S ) = DM
                IF ( STKTK( S ) .LT. 0 ) STKTK( S ) = TK
                IF ( STKVE( S ) .LT. 0 ) STKVE( S ) = VE

            END IF	!  if stack height bad

            IF ( STKDM( S ) .LT. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    DM = DM3( K )
                    TK = TK3( K )
                    VE = VE3( K )
                ELSE
                    K = FIND2( FIP/1000, SCC, NR2, FP2, SC2 )
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
                            DM = DM0
                            TK = TK0
                            VE = VE0
                        END IF 
                    END IF 
                END IF 

                STKDM( S ) = DM
                IF ( STKTK( S ) .LT. 0 ) STKTK( S ) = TK
                IF ( STKVE( S ) .LT. 0 ) STKVE( S ) = VE

            END IF  	!  if stack diameter bad

            IF ( STKTK( S ) .LT. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    TK = TK3( K )
                    VE = VE3( K )
                ELSE
                    K = FIND2( FIP/1000, SCC, NR2, FP2, SC2 )
                    IF ( K .GT. 0 ) THEN
                        TK = TK2( K )
                        VE = VE2( K )
                    ELSE
                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            TK = TK1( K )
                            VE = VE1( K )
                        ELSE
                            TK = TK0
                            VE = VE0
                        END IF 
                    END IF 
                END IF 

                STKTK( S ) = TK
                IF ( STKVE( S ) .LT. 0 ) STKVE( S ) = VE

            END IF	!  if stack exhaust temperature bad

            IF ( STKVE( S ) .LT. 0.0 ) THEN
                FIP = IFIP( S )
                SCC = ISCC( S )
                K = FIND2( FIP, SCC, NR3, FP3, SC3 )

                IF ( K .GT. 0 ) THEN
                    VE = VE3( K )
                ELSE
                    K = FIND2( FIP/1000, SCC, NR2, FP2, SC2 )
                    IF ( K .GT. 0 ) THEN
                        VE = VE2( K )
                    ELSE
                        K = FIND1( SCC, NR1, SC1 )
                        IF ( K .GT. 0 ) THEN
                            VE = VE1( K )
                        ELSE
                            VE = VE0
                        END IF 
                    END IF 
                END IF 

                STKVE( S ) = VE

            END IF	!  if stack exhaust velocity bad

33      CONTINUE        !  end loop fixing stack parameters


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010	FORMAT( 10 ( A, :, I5, :, 1X ) )


        END

