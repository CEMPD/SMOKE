
        SUBROUTINE SETEFTMPR( PIDX, STARTOVR, TI, TMMI, TF, TMIN, TMAX )
   
C***********************************************************************
C  subroutine SETEFTMPR body starts at line < >
C
C  DESCRIPTION:
C     Set the ambient, minimum, and maximum temperature to use for computing
C     emission factors.
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
 
C...........   MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS
c        INTEGER       FIND1

c        EXTERNAL      FIND1

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT    (IN) :: PIDX      ! position of this PSI in list
        LOGICAL     , INTENT(IN OUT) :: STARTOVR  ! true: caller ends loop
        INTEGER     , INTENT   (OUT) :: TI        ! index to VLDTMPR
        INTEGER     , INTENT   (OUT) :: TMMI      ! index to VLDTMIN
        REAL        , INTENT   (OUT) :: TF        ! ambient tmpr
        REAL        , INTENT   (OUT) :: TMIN      ! minimum diurnal tmpr
        REAL        , INTENT   (OUT) :: TMAX      ! maximum diurnal tmpr

C...........   Local allocatable arrays. VLDTMPR is in MODMET.
        INTEGER, ALLOCATABLE, SAVE :: TFEFIDX( : )! ptr to VLDTMPR per PSI/tmpr
        INTEGER, ALLOCATABLE, SAVE :: TFEFPTR( : )! ptr to TFEFIDX per PSI
        INTEGER, ALLOCATABLE, SAVE :: TFEFEND( : )! end position in TFEFIDX per PSI

C...........   Local variables
        INTEGER         I, J, K

        INTEGER         IMM                 ! tmp index for min/max tmpr
        INTEGER         IOS                 ! i/o status
        INTEGER         ITF                 ! tmp index for ambient tmpr
        INTEGER         NTFEF               ! no. ambient tmprs for EFs

        REAL            TMAXCHK             ! tmp max temperature
        REAL            TMINCHK             ! tmp min temperature

        LOGICAL      :: EFLAG    = .FALSE.  ! true: error found
        LOGICAL      :: FIRSTIME = .TRUE.   ! true: first time called

        CHARACTER*300   MESG                ! message buffer

        CHARACTER*16 :: PROGNAME = 'SETEFTMPR' ! program name

C***********************************************************************
C   begin body of subroutine SETEFTMPR

C.........  Allocate local memory for ambient temperatures for each PSI.
C.........  In future, these arrays could be created based on met
C           preprocessing, as with the min/max temperature combinations
        IF( FIRSTIME ) THEN

            NTFEF = NTMPR * NPSIALL
            ALLOCATE( TFEFIDX( NTFEF ),STAT=IOS )
            CALL CHECKMEM( IOS, 'TFEFIDX', PROGNAME )
            ALLOCATE( TFEFPTR( NPSIALL ),STAT=IOS )
            CALL CHECKMEM( IOS, 'TFEFPTR', PROGNAME )
            ALLOCATE( TFEFEND( NPSIALL ),STAT=IOS )
            CALL CHECKMEM( IOS, 'TFEFEND', PROGNAME )

C.............  Initialize arrays so that all PSIs are using all ambient
C               temperatures in the list from MODMET
            DO I = 1, NPSIALL 
                DO J = 1, NTMPR

                    K = ( I-1 ) * NTMPR + J
                    TFEFIDX( K ) = J                

                END DO

                TFEFPTR( I ) = ( I-1 ) * NTMPR + 1
                TFEFEND( I ) = TFEFPTR( I ) + NTMPR - 1

            END DO

            FIRSTIME = .FALSE.

        END IF

C.........  If ambient temperatures for current PSI have not been used up, 
C           then increment to next temperature and set min/max temperatures
C           based on it
        
C.........  Retrieve index for ambient temperature and increment pointer
        ITF = 0
        IF( TFEFPTR( PIDX ) .LE. TFEFEND( PIDX ) ) THEN
            ITF = TFEFPTR( PIDX )
            TFEFPTR( PIDX ) = TFEFPTR( PIDX ) + 1
        END IF

C.........  Retrieve index for min/max temperature combination (no increment)
        IMM = 0
        IF( MMTEFPTR( PIDX ) .LE. MMTEFEND( PIDX ) ) THEN
            IMM = MMTEFPTR( PIDX )
        END IF

C.........  Set missing min/max temperatures
        TMIN = AMISS3
        TMAX = AMISS3

C.........  Case in which ambient temperatures are still being computed
C           takes precendece
        IF( ITF .GT. 0 ) THEN

            TI = TFEFIDX( ITF )
            TF = VLDTMPR( TI )

C.............  Retrieve next available min/max temperature combination
            IF( IMM .GT. 0 ) THEN
                TMMI = MMTEFIDX( IMM )
                TMIN = VLDTMIN ( TMMI )
                TMAX = VLDTMAX ( TMMI )
            END IF

C.............  Evaluate criteria min/max temperatures
            TMINCHK = MAX( TF - 2.*TMMINVL, MINT_MIN )
            TMINCHK = MIN( TMINCHK, MINT_MAX )
            TMAXCHK = MIN( TF + 2.*TMMINVL, MAXT_MAX )
            TMAXCHK = MAX( TMAXCHK, MAXT_MIN )

C.............  Check to see if min/max just happens to be mean the 2*TMMINVL
C               criteria.  If so, increment min/max point and use set flag
C               to indicate diurnal emissions should be stored
            IF( TMIN .EQ. TMINCHK .AND. TMAX .EQ. TMAXCHK ) THEN

                MMTEFPTR( PIDX ) = MMTEFPTR( PIDX ) + 1

C.............  Set min/max so 2*TINTV less/greater than ambient
            ELSE

                TMMI = 0
                TMIN = TMINCHK
                TMAX = TMAXCHK

            END IF

C.............  Reset restart flag
            STARTOVR = .FALSE.

C.........  Now process case in which no ambient temperatures are left for the
C           PSI
        ELSE IF( IMM .GT. 0 ) THEN

            TI   = 0
            TMMI = MMTEFIDX( IMM )
            TMIN = VLDTMIN ( TMMI )
            TMAX = VLDTMAX ( TMMI )

            TF = TMIN + ( TMAX - TMIN ) / 2.0

            MMTEFPTR( PIDX ) = MMTEFPTR( PIDX ) + 1

C.............  Reset restart flag
            STARTOVR = .FALSE.

C.........  If both ambient and min/max temperatures are expended, then
C           set restart flag and return.
        ELSE

            STARTOVR = .TRUE.

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE SETEFTMPR

