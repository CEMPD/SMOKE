
        SUBROUTINE GETPDINFO( FDEV, TZONE, INSTEP, OUTSTEP, TYPNAM, 
     &                        FNAME, SDATE, STIME, NSTEPS, NPDVAR, 
     &                        NPDVSP, MXPDSRC, EAIDX, SPIDX )

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
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NSPDAT

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: MXPDPT

        IMPLICIT NONE

C.........  INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER(2)    CRLF
        INTEGER         GETFLINE
        INTEGER         GETFORMT

        EXTERNAL        CRLF, GETFLINE, GETFORMT

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN):: FDEV          ! file unit no.
        INTEGER     , INTENT (IN):: TZONE         ! output time zone
        INTEGER     , INTENT (IN):: INSTEP        ! expected data time step HHMMSS
        INTEGER     , INTENT (IN):: OUTSTEP       ! output time step HHMMSS
        CHARACTER(*), INTENT (IN):: TYPNAM        ! name of processing type
        CHARACTER(*), INTENT (IN):: FNAME         ! logical file name
        INTEGER     , INTENT(OUT):: SDATE         ! Julian start date in TZONE
        INTEGER     , INTENT(OUT):: STIME         ! start time of data in TZONE
        INTEGER     , INTENT(OUT):: NSTEPS        ! no. time steps
        INTEGER     , INTENT(OUT):: NPDVAR        ! no. pol/act variables
        INTEGER     , INTENT(OUT):: NPDVSP        ! no. pol/act/special data
        INTEGER     , INTENT(OUT):: MXPDSRC       ! max. no. srcs over all times
        INTEGER     , INTENT(OUT):: EAIDX( NIPPA )! index to EANAM
        INTEGER     , INTENT(OUT):: SPIDX( MXSPDAT )! index to SPDATNAM

C.........  Local allocatable arrays...
        LOGICAL, ALLOCATABLE :: EASTAT( : )   ! true: act/pol present in data

C.........  Local arrats
        LOGICAL         SPSTAT( MXSPDAT )     ! true: special data variable used

C.........  Other local variables
        INTEGER         N, V             ! counters and indices
        INTEGER         FILFMT           ! format code of files in list
        INTEGER         INVFMT           ! inventory format code
        INTEGER         IOS              ! i/o status
        INTEGER         NLINE            ! number of lines

        LOGICAL       :: DFLAG    = .FALSE.  ! true: day-specific

        CHARACTER(300)  MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'GETPDINFO' !  program name

C***********************************************************************
C   begin body of program GETPDINFO

C.........  Allocate memory for logical status array for pol/act
        ALLOCATE( EASTAT( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EASTAT', PROGNAME )

        EASTAT = .FALSE.  ! array
        SPSTAT = .FALSE.  ! array

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

C.........  Ensure that input file is a list-formatted file
        INVFMT = GETFORMT( FDEV, -1 )

        IF( INVFMT .NE. LSTFMT ) THEN
            MESG = TYPNAM// '-specific input file is not provided by '//
     &             'a list of files OR ' // CRLF() // BLANK10 // 
     &             'files in list provided could not be found.'
            CALL M3MSG2( MESG )

            MESG = 'Problem reading inventory file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
         
C.........  Get the dates (in the output time zone) from the files, 
C           flag the pollutants of interest, and flag the special variables
C           contained in the file.
        CALL RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, DFLAG, 
     &                 FNAME, SDATE, STIME, NSTEPS, FILFMT, 
     &                 EASTAT, SPSTAT )

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
        CALL RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, DFLAG, 
     &                 FNAME, SDATE, STIME, NSTEPS, FILFMT, 
     &                 EASTAT, SPSTAT )

C.........  Create index to pollutant/activity names for current data files
        N = 0
        DO V = 1, NIPPA

            IF( EASTAT( V ) ) THEN
                N = N + 1
                EAIDX( N ) = V
            END IF

        END DO
        NPDVAR = N

C.........  Create index to special data variable names for current data files
C.........  The idex serves a different purpose from EAIDX and is constructed
C           differently intentionally.
        N = 0
        DO V = 1, MXSPDAT

            IF( SPSTAT( V ) ) THEN
                N = N + 1
                SPIDX( V ) = N
            END IF

        END DO
        NSPDAT = N
        NPDVSP = NPDVAR + NSPDAT

C.........  Compute the maximum number of sources per time step
C.........  NOTE - MXPDPT is in the MODDAYHR module
        MXPDSRC = MAXVAL( MXPDPT )

C.........  If no sources matched then error
        IF ( MXPDSRC .EQ. 0 ) THEN

            MESG = 'No ' // TYPNAM //'-specific sources matched ' //
     &             'the inventory for the time period processed.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Deallocate local memory
        DEALLOCATE( EASTAT )

C.........  Deallocate global memory that is no longer needed
        DEALLOCATE( MXPDPT )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GETPDINFO


