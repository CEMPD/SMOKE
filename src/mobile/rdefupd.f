
        SUBROUTINE RDEFUPD( FDEV, NPSI, VALIDPSI, MXUPDAT, 
     &                      NUPDAT, UPDATE )
   
C***********************************************************************
C  subroutine RDEFUPD body starts at line < >
C
C  DESCRIPTION:
C      Read targeted update file. This file contains a list of PSIs that can
C      be used to force updating of all emission factors for that PSI.  The
C      routine screens out PSIs that are not in the master list.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C 
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C 
C See file COPYRIGHT for conditions of use.
C 
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'    !  emissions constant parameters

C.........  EXTERNAL functions
        CHARACTER*2 CRLF
        INTEGER     FIND1
        INTEGER     STR2INT

        EXTERNAL    CRLF, FIND1, STR2INT

C.........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN)  :: FDEV              ! input update file
        INTEGER, INTENT (IN)  :: NPSI              ! no. of valid PSIs
        INTEGER, INTENT (IN)  :: VALIDPSI( NPSI )  ! valid PSIs
        INTEGER, INTENT (IN)  :: MXUPDAT           ! max no. update PSIs
        INTEGER, INTENT (OUT) :: NUPDAT            ! no. update 
        INTEGER, INTENT (OUT) :: UPDATE( MXUPDAT ) ! array of PSIs

C.........  Variables dimensioned by subroutine arguments
        INTEGER         INDEXA( MXUPDAT )  ! sorting index
        INTEGER         UPDATA( MXUPDAT )  ! unsorted update PSI list

C.........  Local variables
        INTEGER         I, J, K, KK, L, N

        INTEGER         IREC              ! line counter
        INTEGER         IOS               ! I/O status
        INTEGER         PSI               ! tmp PSI

        LOGICAL      :: EFLAG = .FALSE.   ! true: error flag

        CHARACTER*300   MESG    ! message buffer
        CHARACTER*300   LINE    ! line read buffer

        CHARACTER*16 :: PROGNAME = 'RDEFUPD' ! program name

C***********************************************************************
C   begin body of subroutine RDEFUPD

        IREC = 0
        I    = 0
        J    = 0
        DO

            READ( FDEV, 93000, END=101, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            L = INDEX( LINE, '-' ) 

            IF( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'Error ', IOS,
     &                 'reading TARGETED EF UPDATE file at line', IREC
                CALL M3MESG( MESG )
                CYCLE           ! To head of read loop

            ELSE IF( L .LE. 0 ) THEN
                K  = STR2INT( LINE( 1:LEN_TRIM( LINE ) ) )
                KK = K 

            ELSE
                K  = STR2INT( LINE( 1  :L-1 ) )
                KK = STR2INT( LINE( L+1:LEN_TRIM( LINE ) ) )

            END IF

            DO PSI = K, KK

                N = FIND1( PSI, NPSI, VALIDPSI )

                IF( N .LE. 0 ) THEN
                    WRITE( MESG,94010 ) 'WARNING: PSI', PSI, 
     &                     'in EF update file dropped because PSI is '//
     &                     CRLF() // BLANK10 // 
     &                     'not in EF cross-reference file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

                J = J + 1

                IF( J .LE. MXUPDAT ) THEN
                    INDEXA( J ) = J
                    UPDATA( J ) = PSI
                ENDIF

            END DO

        END DO     ! To head of read loop

101     CONTINUE   ! exit from read loop

        NUPDAT = J

C.........  Make sure that number of update arrays
        IF( NUPDAT .GT. MXUPDAT ) THEN

            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'storing EF update information was', MXUPDAT,
     &             CRLF() // BLANK10 // 'but actually needed', NUPDAT
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

C.........  Sort list that we have so far
        CALL SORTI1( NUPDAT, INDEXA, UPDATA )

        DO I = 1, NUPDAT
            J           = INDEXA( I )
            UPDATE( I ) = UPDATA( J )
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

        END SUBROUTINE RDEFUPD
