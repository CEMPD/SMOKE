
        SUBROUTINE RDINVPOL( FILNAM, NSRC, VCNT, JDATE, JTIME, VNAMES, 
     &                       POLDAT, STATUS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Reads inventory pollutant-specific data for variables listed in VNAMES
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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

C...........   INCLUDE FILES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER     INDEX1
        EXTERNAL    INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*), INTENT (IN) :: FILNAM           ! Logical file name	
        INTEGER      , INTENT (IN) :: NSRC             ! Number of sources
        INTEGER      , INTENT (IN) :: VCNT             ! Number of variables
        INTEGER      , INTENT (IN) :: JDATE            ! Julian date
        INTEGER      , INTENT (IN) :: JTIME            ! time (HHMMSS)
        CHARACTER*(*), INTENT (IN) :: VNAMES( VCNT )   ! Variable names
        REAL         , INTENT(OUT) :: POLDAT( NSRC,VCNT ) ! Data
        INTEGER      , INTENT(OUT) :: STATUS           ! Exit status

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: IREAD ( : )  ! integer read array

C...........   Other local variables

        INTEGER         K, LV, V  ! counters and indices
        INTEGER         IOS       ! i/o status

        CHARACTER(LEN=IOVLEN3)   VARBUF
        CHARACTER*300   MESG 

        CHARACTER*16 :: PROGNAME = 'RDINVPOL' ! program name

C***********************************************************************
C   begin body of subroutine RDINVPOL

        STATUS = 0

C.........  Get description of file header so we can get variable types
        IF ( .NOT. DESC3( FILNAM ) ) THEN
            MESG = 'Could not get description of file "'
     &             // FILNAM( 1:LEN_TRIM( FILNAM ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
       

C.........  Read variables by type, and store as REAL in POLDAT
        DO V = 1, VCNT

            VARBUF = VNAMES( V )
            LV = LEN_TRIM( VARBUF )

C.............  Find variable name in list to get type
            K = INDEX1( VARBUF, NVARS3D, VNAME3D )

            IF( VTYPE3D( K ) .EQ. M3INT ) THEN

C.................  If memory is not allocated for integer read array, then
C                   allocate it
                IF( .NOT. ALLOCATED( IREAD ) ) THEN
                    ALLOCATE( IREAD( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IREAD', PROGNAME )
                END IF

                IF( .NOT. READ3( FILNAM, VARBUF, ALLAYS3,
     &                           JDATE, JTIME, IREAD ) ) THEN
                    STATUS = 1
                    MESG = 'ERROR: Could not read "' //
     &                      VARBUF( 1:LV ) // '" from file.'
                    CALL M3MSG2( MESG )

                ELSE

                    POLDAT( :,V ) = REAL( IREAD )   ! Array

                END IF

            ELSE IF( .NOT. READ3( FILNAM, VARBUF, ALLAYS3,
     &                       JDATE, JTIME, POLDAT( 1,V ) ) ) THEN
                STATUS = 1
                MESG = 'ERROR: Could not read "' //
     &                 VARBUF( 1:LV ) // '" from file.'
                CALL M3MSG2( MESG )

            END IF

        END DO

C.........  Deallocate local memory
        IF( ALLOCATED( IREAD ) ) DEALLOCATE( IREAD )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

        END

