
        SUBROUTINE RD3MASK( FNAME, JDATE, JTIME, NDIM, NVLIST, VNAMES, 
     &                      VINDX, OUTVAR )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads in certain variables from a list of variables, 
C      depending on the index, which indicates if those variables are expected
C      to be present in a file.  The index is a pointer which says which
C      part of the output array to read the variable into. The variables are
C      assumed to be reals.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         INDEX1

        EXTERNAL        CRLF, INDEX1

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME     ! i/o api file name
        INTEGER     , INTENT (IN) :: JDATE     ! Julian date (YYYYDDD)
        INTEGER     , INTENT (IN) :: JTIME     ! time (HHMMSS)
        INTEGER     , INTENT (IN) :: NDIM      ! dimension for output array
        INTEGER     , INTENT (IN) :: NVLIST    ! variable description to read
        CHARACTER(*), INTENT (IN) :: VNAMES( NVLIST ) ! variable names
        INTEGER     , INTENT (IN) :: VINDX ( NVLIST ) ! var index to OUTVAR
        REAL        , INTENT(OUT) :: OUTVAR( NDIM,* ) ! coeffs for sources

C.........  Other local variables
        INTEGER         J, L, L1, V       !  counters and indices

        CHARACTER(LEN=IOVLEN3 ) VBUF    !  variable name buffer
        CHARACTER*300           MESG    !  message buffer

        CHARACTER*16 :: PROGNAME = 'RD3MASK' ! program name

C***********************************************************************
C   begin body of subroutine RD3MASK

C.........  Loop through variables that are possibilities for reading
        DO V = 1, NVLIST

C.............  Check variable index to see if this variable is to be read
            J = VINDX( V ) 
            IF( J .EQ. 0 ) CYCLE  ! to go bottom of loop

C.............  Read variable and print nice error message if cannot
            VBUF = VNAMES( J )
            IF ( .NOT. READ3( 
     &           FNAME, VBUF, 1, JDATE, JTIME, OUTVAR( 1,J ) ) ) THEN

                L  = LEN_TRIM( FNAME )
                L1 = LEN_TRIM( VBUF )
                MESG = 'Could not read variable "' // VBUF( 1:L1 ) //
     &                 '" from file "' // FNAME( 1:L )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF    !  if read3() failed for file

        END DO        !  End loop on possible variables

        RETURN

        END SUBROUTINE RD3MASK
