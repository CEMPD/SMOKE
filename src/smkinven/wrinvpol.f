
        SUBROUTINE WPNTSPOL( ENAME, NPSRC, IPCNT, POLNAM, POLBUF )

C***********************************************************************
C  subroutine body starts at line 89
C
C  DESCRIPTION:
C      This subroutine writes the point source pollutant-based variables
C      to the I/O API files
C
C  PRECONDITIONS REQUIRED:
C      Logical file name ENAME opened
C      Number of sources NPSRC defined correctly
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

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         TRIMLEN

        EXTERNAL        TRIMLEN

C.........  Subroutine arguments and their descriptions:
        CHARACTER*(*) ENAME                 !  I/O API file name
        INTEGER       NPSRC                 !  actual source count
        INTEGER       IPCNT                 !  pollutant count
        CHARACTER*(*) POLNAM( IPCNT )       !  names of pollutants for output
        REAL          RULPEN( NPSRC,IPCNT ) !  rule penetration   fraction
        REAL          RULEFF( NPSRC,IPCNT ) !  rule effectiveness fraction
        REAL          CTLEFF( NPSRC,IPCNT ) !  control efficiency fraction
                                            !  pol-based emissions & other data:
        REAL          POLBUF( NPSRC,IPCNT*NPTPPOL3 )

C...........   Other local variables
        INTEGER                 I, J, V

        INTEGER                 IDUM   ( IPCNT,NPTPPOL3 ) ! Int dummy

        CHARACTER(LEN=IOVLEN3 ) VAR  !  tmp variable names
        CHARACTER(LEN=IOVLEN3)  ECNAMES( IPCNT,NPTPPOL3 ) ! Names for pol-spec
        CHARACTER(LEN=IODLEN3)  CDUM   ( IPCNT,NPTPPOL3 ) ! Char dummy
 
        CHARACTER*300 MESG             !  message buffer

        CHARACTER*16 :: PROGNAME = 'WPNTSPOL' !  program name

C***********************************************************************
C   begin body of program WPNTSPOL

C.........  Create message to use in case there is an error
        MESG = 'Error writing output file "' //
     &         ENAME( 1:TRIMLEN( ENAME ) ) // '"'

C.........  Get the list of variable names per pollutant
        CALL BLDENAMS( 'POINT', IPCNT, NPTPPOL3, POLNAM, ECNAMES,
     &                  CDUM, IDUM, CDUM )

C.........  Write the I/O API file, one variable at a time
        DO V = 1, IPCNT

            DO I = 1, NPTPPOL3

                J = ( V - 1 ) * NPTPPOL3 + I
                VAR = ECNAMES( V,I )

                IF( .NOT. WRITE3( ENAME, VAR, 0, 0, POLBUF(1,J) ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            ENDDO

        ENDDO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93100   FORMAT( I2 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2 ) )

        END
