
        SUBROUTINE TMPRINFO( ENVFLAG, CONTROL )

C***********************************************************************
C  subroutine TMPRINFO body starts at line 83
C
C  DESCRIPTION:
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
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

C...........   MODULES for public variables
C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants

C...........   EXTERNAL FUNCTIONS
        REAL      ENVREAL
        REAL      GETRFDSC

        EXTERNAL  ENVREAL, GETRFDSC

C.........  SUBROUTINE ARGUMENTS
        LOGICAL     , INTENT(IN) :: ENVFLAG  ! true: initialize tmpr info from
                                             !       environment variables
                                             ! false: get values from FDESC3D
        CHARACTER(*), INTENT(IN) :: CONTROL  ! 'NOMINMAX' | 'MINMAX' | 'BOTH'

C...........   Other local variables
        INTEGER     I, J, K     ! counters and indices

        INTEGER     IOS         ! ENVREAL return status
        INTEGER     NMAXALL     ! no. of all possible maximum values
        INTEGER     NMINALL     ! no. of all possible minimum values
        INTEGER     T1, T2      ! tmp counters for min/max tmpr setting

        DOUBLE PRECISION :: MNTMN       ! tmp min tmpr min
        DOUBLE PRECISION :: MNTMX       ! tmp min tmpr max
        DOUBLE PRECISION :: MXTMN       ! tmp max tmpr min
        DOUBLE PRECISION :: MXTMX       ! tmp max tmpr max
        DOUBLE PRECISION :: TMXIT       ! tmp max interval
        DOUBLE PRECISION :: TINTV       ! tmp tmpr interval
        DOUBLE PRECISION :: TDIFF       ! tmp temperature difference
        DOUBLE PRECISION :: TMAX        ! tmp maximum temperature value
        DOUBLE PRECISION :: TMIN        ! tmp mimimum temperature value

        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'TMPRINFO' ! program name

C***********************************************************************
C   begin body of subroutine TMPRINFO

C.........  Get the values from the environment...

        IF( ENVFLAG ) THEN

C.............  Get ranges for minimum and maximum daily temperatures...
C.............  Defaults are from MOBILE5a/b
            MESG = 'Minimum allowed minimum daily temperature [deg F]'
            MINT_MIN = ENVREAL( 'SMK_MINT_MIN', MESG, 0., IOS )

            MESG = 'Maximum allowed minimum daily temperature [deg F]'
            MINT_MAX = ENVREAL( 'SMK_MINT_MAX', MESG, 100., IOS )

            MESG = 'Minimum allowed maximum daily temperature [deg F]'
            MAXT_MIN = ENVREAL( 'SMK_MAXT_MIN', MESG, 10., IOS )

            MESG = 'Maximum allowed maximum daily temperature [deg F]'
            MAXT_MAX = ENVREAL( 'SMK_MAXT_MAX', MESG, 120., IOS )

C.............  Get tmpr intervals to use for computing emission factors
            MESG = 'Temperature interval [deg F]'
            TMMINVL = ENVREAL( 'SMK_TF_INTERVAL', MESG, 2., IOS )

            MESG = 'Maximum permitted temperature interval [deg F]'
            TMXINVL = ENVREAL( 'SMK_TF_MAXINTVL', MESG, 40., IOS )


C.............  Convert temperature parameters to proper units
C.............  Kelvin for now - future can be dependant on met input units

            MINT_MIN = FTOC * ( MINT_MIN - 32. ) + CTOK
            MINT_MAX = FTOC * ( MINT_MAX - 32. ) + CTOK
            MAXT_MIN = FTOC * ( MAXT_MIN - 32. ) + CTOK
            MAXT_MAX = FTOC * ( MAXT_MAX - 32. ) + CTOK
            TMMINVL  = FTOC * TMMINVL
            TMXINVL  = FTOC * TMXINVL

C.........  Get the values from the FDESC3D array
        ELSE

            MINT_MIN = GETRFDSC( FDESC3D, '/MINT_MIN/', .TRUE. )
            MAXT_MAX = GETRFDSC( FDESC3D, '/MAXT_MAX/', .TRUE. )
            TMMINVL  = GETRFDSC( FDESC3D, '/T_INTERVAL/', .TRUE. )

            IF( CONTROL .NE. 'NOMINMAX' ) THEN
                MINT_MAX = GETRFDSC( FDESC3D, '/MINT_MAX/', .TRUE. )
                MAXT_MIN = GETRFDSC( FDESC3D, '/MAXT_MIN/', .TRUE. )
                TMXINVL  = GETRFDSC( FDESC3D, '/T_MAXINTVL/', .TRUE. )
            END IF

        END IF

C.........  Store temperatures as double precision variables
        MNTMN = DBLE( MINT_MIN )
        MNTMX = DBLE( MINT_MAX )
        MXTMN = DBLE( MAXT_MIN )
        MXTMX = DBLE( MAXT_MAX )
        TINTV = DBLE(  TMMINVL )
        TMXIT = DBLE(  TMXINVL )

C.........  For non min-max type call, create list of value temperatures, and
C           return to calling routine
        IF ( CONTROL .NE. 'MINMAX' ) THEN

C.............  Determine the number of valid temperatures
            NTMPR = ( MAXT_MAX - MINT_MIN ) / TMMINVL + 1

C.............  Allocate memory for valid temperatures and set them
            IF( .NOT. ALLOCATED( VLDTMPR ) ) THEN
                ALLOCATE( VLDTMPR( NTMPR ), STAT=IOS )
                CALL CHECKMEM( IOS, 'VLDTMPR', PROGNAME )

                DO I = 1, NTMPR
                    VLDTMPR( I ) = REAL( MNTMN + ( I - 1 ) * TINTV )
                END DO
            END IF

        END IF

C.............  If not doing minimum/maximum array, then return
        IF ( CONTROL .EQ. 'NOMINMAX' ) RETURN

C.........  Allocate maximum amount of memory for arrays of min/max
C           temperatures
        NMINALL = ( MINT_MAX - MINT_MIN ) / TMMINVL + 1
        NMAXALL = ( MAXT_MAX - MAXT_MIN ) / TMMINVL + 1 

C.........  Allocate memory for valid temperatures and set them
        IF( .NOT. ALLOCATED( VLDTMIN ) ) THEN
            ALLOCATE( VLDTMIN( NMINALL * NMAXALL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VLDTMIN', PROGNAME )
            ALLOCATE( VLDTMAX( NMINALL * NMAXALL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VLDTMAX', PROGNAME ) 

C.............  Initialize
            VLDTMIN = 0.
            VLDTMAX = 0.

C.............  Build arrays of min/max temperatures
C.............  Start TMAX computation from MNTMN to make sure all values use
C               the same temporal basis.
            K = 0        
            T1 = 0
            TMIN = MNTMN
            DO WHILE( TMIN .LE. MNTMX )

        	T2 = 0
        	TMAX = MNTMN
        	DO 
                    IF( TMAX                .GT. MXTMN .OR. 
     &                  ABS( TMAX - MXTMN ) .LE. 0.01        ) THEN


                        TDIFF = TMAX - TMIN

C.........................  Increment table and store for valid tmpr combos
                        IF( TMAX  .GT. TMIN  .AND.
     &                      TDIFF .LE. TMXIT      ) THEN

                	    K = K + 1
                	    VLDTMIN( K ) = REAL( TMIN )
                	    VLDTMAX( K ) = REAL( TMAX )

                	END IF
                    END IF

                    IF( TMAX                .GT. MXTMX .OR. 
     &                  ABS( TMAX - MXTMX ) .LE. 0.01        ) EXIT

C.....................  Increase maximum
                    T2 = T2 + 1
                    TMAX = MNTMN + T2 * TINTV

                END DO
 
C.................  Increase minimum
        	T1 = T1 + 1
        	TMIN = MNTMN + T1 * TINTV

            END DO

            NVLDTMM = K

        END IF  ! already been set or not
        
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )
 
        END SUBROUTINE TMPRINFO
