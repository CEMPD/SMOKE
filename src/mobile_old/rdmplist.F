
C Version "@(#)$Id$ $Source$ $Date$ 

        INTEGER FUNCTION RDMPLIST( FDEV, MPFIP, MPRDT, MPLNK, 
     &                             MPREF, REFCNT, REFLIST ) 

C***********************************************************************
C  subroutine body starts at line  71
C
C  FUNCTION:  read and sort file MPLIST for TMPMOBIL
C
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

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
      INCLUDE 'CHDIMS3.EXT'     ! chemistry parameters
      INCLUDE 'MBDIMS3.EXT'     ! mobile-source parameters


C...........   ARGUMENTS and their descriptions:

        INTEGER   FDEV
        INTEGER   MPFIP( NLIST )
        INTEGER   MPRDT( NLIST )
        INTEGER   MPLNK( NLIST )
        INTEGER   MPREF( 0:23, NLIST )    !  hour&source-specific xrefs
        INTEGER   REFCNT
        INTEGER   REFLIST( NLIST )


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER   FIND1
        EXTERNAL  FIND1

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER   IREC
        INTEGER   MPCNT
        INTEGER   MPDEXA( NLIST )
        INTEGER   MPFIPA( NLIST )
        INTEGER   MPRDTA( NLIST )
        INTEGER   MPLNKA( NLIST )
        INTEGER   MPREFA( 0:23, NLIST )   !  hour&source-specific xrefs
        INTEGER   ISTAT
        INTEGER   I, J, K, L
        INTEGER	  FIP, RDT, LNK, REF( 0:23 )

        LOGICAL	  EFLAG

        CHARACTER*256	MESG

C***********************************************************************
C   begin body of subroutine  RDMPLIST

        MPCNT  = 0
        REFCNT = 0
        IREC   = 0
        I      = 0
        EFLAG  = .FALSE.

11      CONTINUE	!  head of loop reading FDEV

            IREC = IREC + 1
            READ( FDEV, *, END = 99, IOSTAT=ISTAT )
     &              FIP,
     &              RDT,
     &              LNK,
     &            ( REF( J ), J = 0, 23 )

            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG,94010 ) 
     &          'I/O error in MPLIST at line', IREC,
     &          'IOSTAT=', ISTAT
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  11
            END IF

            I = I + 1

            IF ( I .GT. NLIST ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Emission Factor XREF table overflow at line', IREC
                CALL M3MESG( MESG )
                GO TO  11
            END IF

            MPDEXA( I ) = I
            MPFIPA( I ) = FIP
            MPRDTA( I ) = RDT
            MPLNKA( I ) = LNK

            DO  23  J = 0, 23

                MPREFA( J,I ) = REF( J )

                IF( FIND1( REF( J ), REFCNT, REFLIST ) .LE. 0 ) THEN

                    REFCNT = REFCNT + 1
                    DO  21  K = REFCNT, 1, -1
                        IF ( REFLIST( K-1 ) .GT. REF( J ) ) THEN
                            REFLIST( K ) = REFLIST( K-1 )
                        ELSE
                            REFLIST( K ) = REF( J )
                            GO TO  22
                        END IF
21                  CONTINUE
22                  CONTINUE	!  exit from insertion

                END IF  	!  if new hour index

23          CONTINUE	!  end loop on hour indexes J

        GO TO 11        !  to head of reading FDEV

99      CONTINUE	!  end of loop reading FDEV

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
C...........   Internal buffering formats............ 94xxx

94010	FORMAT( 10( A, :, I7, :, 2X ) )

        END

