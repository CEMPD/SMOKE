        SUBROUTINE RDTMPU( UINAME, FNAME, JDATE, JTIME, NDIM, NVLIST,
     &                     VNAMES, VINDX, OUTVAR )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine calls RD3MASK twice to read both certain and uncerain 
C      variables into their proper place from a list 
C      of variables, 
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
C.........  This module contains the major data structure and control flags
        USE MODMERGE

C.........  This module contains the global variables for uncertainty
        USE MODUNCERT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file desc. data structures
c        INCLUDE 'IOSTRG3.EXT'   !  I/O API string lengths

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         INDEX1
        CHARACTER*16    PROMPTMFILE  

        EXTERNAL        CRLF, INDEX1, PROMPTMFILE

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: UINAME      ! invetory file name
        CHARACTER(*), INTENT (IN) :: FNAME       ! i/o api file name
        INTEGER     , INTENT (IN) :: JDATE       ! Julian date (YYYYDDD)
        INTEGER     , INTENT (IN) :: JTIME       ! time (HHMMSS)
        INTEGER     , INTENT (IN) :: NDIM        ! dimension for output array
        INTEGER     , INTENT (IN) :: NVLIST      ! variable description to read
        CHARACTER(*), INTENT (IN) :: VNAMES( NVLIST ) ! variable names
        INTEGER     , INTENT (IN) :: VINDX ( NVLIST ) ! var index to OUTVAR
        REAL        , INTENT(OUT) :: OUTVAR( NDIM,* ) ! coeffs for sources

C.........  Other local variables
        INTEGER         J, L, L1, LASTJ, N, S, V ! counters and indices
        INTEGER         IOS                      ! I/O status variables
        INTEGER         LASTJ                    ! save J value

        REAL,  ALLOCATABLE       ::  INTMPU( :,: )    ! uncert emis facts from temporal
        REAL,  ALLOCATABLE       ::  INITV( : )       ! used to initialize vector to 0

        CHARACTER*16                 NAMBUF      ! temp name buffer
        CHARACTER*16                 VBUF        ! variable buffer
        CHARACTER*300                MESG        ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDTMPU' ! program name

C***********************************************************************
C   begin body of subroutine RDTMPU

C.........  Prompt for and open I/O API temporal uncertainty file(s)...
        MESG = 'Reading temporal uncertainty input file...'
        NAMBUF = PROMPTMFILE( MESG, FSREAD3, FNAME, PROGNAME )

C.........  Get header description of non-diurnal file to get date for reading
        IF( .NOT. DESC3( NAMBUF ) ) THEN
            MESG = 'Could not get description of file "' 
     &              // NAMBUF // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for emission factors and initialize to zero
        ALLOCATE( INTMPU( UNSRC, NIPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INTMPU', PROGNAME )

        CALL RDU3MASK( NAMBUF, JDATE, JTIME, UNSRC, NVLIST,
     &                 VNAMES, UI_EXIST, INTMPU )

C.........  Loop through variables that are possibilities for reading
        LASTJ = 0
        DO N = 1, NVLIST

            J = VINDX( N ) 
            IF( J .EQ. 0 .OR. LASTJ .EQ. J) CYCLE
            
            LASTJ = J

            DO V = 1, NDIM

                OUTVAR( V,J ) = 0.0 

                IF( USTAT( V ) ) THEN

                    IF( INSRC( SRCNUM( V ) ) .NE. V ) THEN 
                        WRITE( MESG,94010) 'debug: discrepancy for ', V
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                    OUTVAR( V,J ) = INTMPU( SRCNUM( V ),J ) 

                END IF
            END DO

        END DO

        IF( AFLAG ) THEN

            A_EXIST = UI_EXIST

            IF( SFLAG ) THEN
                DO V = 1, NVLIST
                    IF( VINDX( V ) .EQ. 0 ) AS_EXIST( V, 1) = 0
                END DO
            END IF

        END IF

        DEALLOCATE( INTMPU )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

        END SUBROUTINE RDTMPU
