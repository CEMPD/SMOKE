
        SUBROUTINE GETPDINFO( FDEV, TZONE, TSTEP, TYPNAM, SDATE, STIME, 
     &                        NSTEPS, NPDVAR, MXPDSRC, EAIDX )

C***************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine gets the vital information from the day-specific or 
C      hour-specific input files so that memory can be allocated to read in
C      the data.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 12/99 by M. Houyoux
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR

        IMPLICIT NONE

C.........  INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN):: FDEV          ! file unit no.
        INTEGER     , INTENT (IN):: TZONE         ! output time zone
        INTEGER     , INTENT (IN):: TSTEP         ! time step HHMMSS
        CHARACTER(*), INTENT (IN):: TYPNAM        ! name of processing type
        INTEGER     , INTENT(OUT):: SDATE         ! Julian start date in TZONE
        INTEGER     , INTENT(OUT):: STIME         ! start time of data in TZONE
        INTEGER     , INTENT(OUT):: NSTEPS        ! no. time steps
        INTEGER     , INTENT(OUT):: NPDVAR        ! no. pol/act variables
        INTEGER     , INTENT(OUT):: MXPDSRC       ! max. no. srcs over all times
        INTEGER     , INTENT(OUT):: EAIDX( NIPPA )! index to EANAM

C.........  Local allocatable arrays
        LOGICAL, ALLOCATABLE :: EASTAT( : )    ! true: act/pol present in data

C.........  Other local variables
        INTEGER         N, V             ! counters and indices
        INTEGER         IOS              ! i/o status

        LOGICAL       :: DFLAG    = .FALSE.  ! true: day-specific

        CHARACTER*300   MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'GETPDINFO' !  program name

C***********************************************************************
C   begin body of program GETPDINFO

C.........  Allocate memory for logical status array for pol/act
        ALLOCATE( EASTAT( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EASTAT', PROGNAME )
        EASTAT = .FALSE.  ! array

        MESG = 'Determining number of time steps for ' // TYPNAM //
     &         '-specific files...'
        CALL M3MSG2( MESG )
        
C.........  Perform case-specific settings
        SELECT CASE( TYPNAM )
        CASE( 'day' ) 
            DFLAG = .TRUE.

        CASE( 'hour' )
            DFLAG = .FALSE.

        CASE DEFAULT
            MESG = 'INTERNAL ERROR: Do not know type ' // TYPNAM 
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.........  Get the dates (in the output time zone) from the files and
C           flag the pollutants of interest
        CALL RDLOOPPD( FDEV, TZONE, TSTEP, MXPDSRC, DFLAG, SDATE, STIME,
     &                 NSTEPS, EASTAT )

C.........  Allocate memory and initialize for the maximum number of 
C           records per time step
        ALLOCATE( MXPDPT( NSTEPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MXPDPT', PROGNAME )
        MXPDPT = 0  ! array

        MESG = 'Determining number of sources for ' // TYPNAM //
     &         '-specific files...'
        CALL M3MSG2( MESG )

C.........  Get the maximum number of records per time step - i.e., populate
C           MXSRCPD
        CALL RDLOOPPD( FDEV, TZONE, TSTEP, MXPDSRC, DFLAG, SDATE, STIME, 
     &                 NSTEPS, EASTAT )

C.........  Create index to pollutant/activity names for current data files
        N = 0
        DO V = 1, NIPPA

            IF( EASTAT( V ) ) THEN
                N = N + 1
                EAIDX( N ) = V
            END IF

        END DO
        NPDVAR = N

C.........  Compute the maximum number of sources per time step
C.........  NOTE - MXPDPT is in the MODDAYHR module
        MXPDSRC = MAXVAL( MXPDPT )

C.........  Deallocate local memory
        DEALLOCATE( EASTAT )

C.........  Deallocate global memory that is no longer needed
        DEALLOCATE( MXPDPT )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GETPDINFO


