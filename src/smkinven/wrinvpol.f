
        SUBROUTINE WRINVPOL( ENAME, CATEGORY, NSRC, IPCNT, NPPOL, 
     &                       POLNAM, POLBUF )

C***********************************************************************
C  subroutine body starts at line 86
C
C  DESCRIPTION:
C      This subroutine writes the inventory pollutant-based variables
C      and/or VMT to the I/O API files
C
C  PRECONDITIONS REQUIRED:
C      Logical file name ENAME opened
C      Number of sources NSRC defined correctly
C      Pollutant count IPCNT and output names POLNAM defined correctly
C      Output array POLVAL populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine, BLDENAMS
C      Functions: I/O API subroutine
C
C  REVISION  HISTORY:
C      Created 11/98 by M. Houyoux
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

      IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  Subroutine arguments and their descriptions:
        CHARACTER(*), INTENT (IN) :: ENAME            !  I/O API file name
        CHARACTER(*), INTENT (IN) :: CATEGORY         !  source category
        INTEGER     , INTENT (IN) :: NSRC             !  no. sources
        INTEGER     , INTENT (IN) :: IPCNT            !  pollutant count
        INTEGER     , INTENT (IN) :: NPPOL            !  no. vars per pol
        CHARACTER(*), INTENT (IN) :: POLNAM( IPCNT )  !  names of output pols
        REAL        , INTENT (IN) :: POLBUF( NSRC,IPCNT*NPPOL )
                                            !  pol-based emissions & other data
C...........   Output variable information
        INTEGER                 EOTYPES( IPCNT,NPPOL ) ! types output vars
        CHARACTER(LEN=IOVLEN3)  EONAMES( IPCNT,NPPOL ) ! names for pol-spec
        CHARACTER(LEN=IODLEN3)  CDUM   ( IPCNT,NPPOL ) ! char dummy

C...........   Temporary integer array for output of integer variables
        INTEGER, ALLOCATABLE :: INTBUF ( : )

C...........   Other local variables
        INTEGER                 I, J, L, V
        INTEGER                 IOS     !  i/o status

        CHARACTER*300           MESG    !  message buffer
        CHARACTER(LEN=IOVLEN3 ) VAR     !  tmp variable name
 
        CHARACTER*16 :: PROGNAME = 'WRINVPOL' !  program name

C***********************************************************************
C   begin body of program WRINVPOL

C.........  Create message to use in case there is an error
        MESG = 'Error writing output file "' //
     &         ENAME( 1:LEN_TRIM( ENAME ) ) // '"'

C.........  Get the list of variable names per pollutant
        CALL BLDENAMS( CATEGORY, IPCNT, NPPOL, POLNAM, EONAMES,
     &                 CDUM, EOTYPES, CDUM )

C.........  Write the I/O API file, one variable at a time
        DO V = 1, IPCNT

            DO I = 1, NPPOL

                J = ( V - 1 ) * NPPOL + I
                VAR = EONAMES( V,I )
                L   = LEN_TRIM( VAR )

C.................  Convert real variable to an integer, if needed
                IF( EOTYPES( V,I ) .EQ. M3INT ) THEN

                    IF( .NOT. ALLOCATED( INTBUF ) ) THEN
                        ALLOCATE( INTBUF( NSRC ), STAT=IOS )
                        CALL CHECKMEM( IOS, 'INTBUF', PROGNAME )
                    END IF

                    INTBUF = INT( POLBUF( 1,J ) )  ! arrays

                    IF( .NOT. WRITE3( ENAME,VAR(1:L),0,0,INTBUF )) THEN
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

C.................  Simply write out real values directly
                ELSE
                
                    IF( .NOT. WRITE3( ENAME,VAR(1:L),0,0,
     &                                POLBUF(1,J)         ) ) THEN
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END IF

            END DO

        END DO

        DEALLOCATE( INTBUF )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93100   FORMAT( I2 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2 ) )

        END SUBROUTINE WRINVPOL
