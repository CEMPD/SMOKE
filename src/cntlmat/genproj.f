
        SUBROUTINE GENPROJ( NSRC, NIPPA, BYEAR, PYEAR, ENAME, 
     &                      USEPOL, EANAM )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine processes the projection data and writes out
C      the projection matrix.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     
C
C************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         PROMPTFFILE

        EXTERNAL   PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC   ! no. sources
        INTEGER     , INTENT (IN) :: NIPPA  ! no. inventory pollutants
        INTEGER     , INTENT (IN) :: BYEAR  ! base year
        INTEGER     , INTENT (IN) :: PYEAR  ! projection year
        CHARACTER(*), INTENT (IN) :: ENAME  ! emission inventory file name
        LOGICAL     , INTENT (IN) :: USEPOL( NIPPA )  ! true: pol used in pkt
        CHARACTER(*), INTENT (IN) :: EANAM ( NIPPA )  ! all pollutant names        

C...........   Local arrays allocated by subroutine arguments
        INTEGER          ISPRJ ( NSRC ) ! projection control data table index
        REAL             PRJFAC( NSRC ) ! projection factor

C...........   Logical names

        CHARACTER*16     PNAME      ! logical name for projection matrix

C...........   Other local variables

        INTEGER          K, L, S    ! counters and indices
        INTEGER          IDUM       ! dummy integer
        INTEGER          IOS        ! i/o error status

        LOGICAL       :: EFLAG    = .FALSE.   ! true: error has occurred

        CHARACTER*16     VARNAM 
        CHARACTER*300    MESG                 ! message buffer

        CHARACTER*16  :: PROGNAME = 'GENPROJ' ! program name

C***********************************************************************
C   begin body of subroutine GENPROJ


C.........  Determine which projection packet goes to each source.
C           Since the projection packet is not pollutant specific,
C           simply send the first pollutant in the pollutant list
C           to ASGNCNTL in order to avoid looping through all pollutants.

         CALL ASGNCNTL( NSRC, 1, 'PROJECTION', USEPOL, EANAM(1), 
     &                  IDUM, ISPRJ )

C.........  Allocate memory for the projection matrix.  Use data
C           structures for point sources, but this routine can be used for area
C           sources or mobile sources as well. 

C        ALLOCATE( PRJFAC( NSRC ), STAT=IOS )
C        CALL CHECKMEM( IOS, 'PRJFAC', PROGNAME )

C.........  Loop through all sources and store projection information for
C           those that have it. Otherwise, set projection factor=1.

c        ISPRJ = 1  ! array

        DO S = 1, NSRC

            K = ISPRJ( S )       ! index to projection data tables

            IF( K .GT. 0 ) THEN

C.................  Store projection factor

                PRJFAC( S ) = PRJFC( K )

            ELSE
C.................  If source does not have projection info., set to 1.

                PRJFAC( S ) = 1.0

            END IF

        END DO

C.........  Set up and open output projection matrices

        CALL OPENPMAT( ENAME, BYEAR, PYEAR, PNAME )

C.........  Write the projection matrix

C.........  Initialize message to use in case there is an error

        MESG = 'Problem writing to output file "' //
     &         PNAME( 1:LEN_TRIM( PNAME ) ) // '"'

        L = LEN_TRIM( MESG )

C.........  Write the I/O API variables for the non-speciation data

        IF( .NOT. WRITE3( PNAME, 'PFAC', 0, 0, PRJFAC ) ) THEN
            MESG = MESG( 1:L ) // ' for variable "PRJFAC"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF


c note: these deallocations cause an unexplained error on SGI.  There is
c    n: a memory problem somewhere that is probably causing this, but
c    n: I was not able to find it.

c        DEALLOCATE( PRJFC )


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        END SUBROUTINE GENPROJ
