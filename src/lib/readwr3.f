
        SUBROUTINE READWR3( INFILE, OUTFILE, VNAME, LAYSVAL, 
     &                      JDATE, JTIME, VTYPE, NDIM, STATUS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Reads variable from INFILE and write variable to OUTFILE
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

C...........   INCLUDE FILES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS:
        INTEGER         TRIMLEN

        EXTERNAL        TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*)   INFILE           ! Name of file being read	
        CHARACTER*(*)   OUTFILE          ! Name of file being written	
        CHARACTER*(*)   VNAME            ! Variable name being read/written
        INTEGER         LAYSVAL          ! layer number or value of ALLAYS3
        INTEGER         JDATE            ! Julian date being read/written
        INTEGER         JTIME            ! Julian time being read/written
        INTEGER         VTYPE            ! Integer code for variable type
        INTEGER         NDIM             ! Dimension of var being read/written
        INTEGER         STATUS           ! Exit status

C...........   Allocatable arrays
        INTEGER         INTVAL ( NDIM )  !  Integer value
        REAL            REALVAL( NDIM )  !  Real value

C...........   Other local variables
        CHARACTER*300   MESG 

        CHARACTER*16 :: PROGNAME = 'READWR3' ! program name

C***********************************************************************
C   begin body of subroutine READWR3

        STATUS = 0

        IF( VTYPE .EQ. M3INT ) THEN
            IF( .NOT. READ3( INFILE, VNAME, LAYSVAL,
     &          JDATE, JTIME, INTVAL ) ) THEN
                STATUS = 1
                MESG = 'ERROR: Could not read "' //
     &                 VNAME( 1:TRIMLEN( VNAME ) ) // '" from file.'
                CALL M3MESG( MESG )
            END IF

            IF( .NOT. WRITE3( OUTFILE, VNAME, LAYSVAL,
     &          JDATE, JTIME, INTVAL ) ) THEN
                STATUS = 1
                MESG = 'ERROR: Could not write "' //
     &                 VNAME( 1:TRIMLEN( VNAME ) ) // '" from file.'
                CALL M3MESG( MESG )
            END IF

        ELSEIF( VTYPE .EQ. M3REAL ) THEN
            IF( .NOT. READ3( INFILE, VNAME, LAYSVAL,
     &          JDATE, JTIME, REALVAL ) ) THEN
                STATUS = 1
                MESG = 'ERROR: Could not read "' //
     &                 VNAME( 1:TRIMLEN( VNAME ) ) // '" from file.'
                CALL M3MESG( MESG )
            END IF

            IF( .NOT. WRITE3( OUTFILE, VNAME, LAYSVAL,
     &          JDATE, JTIME, REALVAL ) ) THEN
                STATUS = 1
                MESG = 'ERROR: Could not write "' //
     &                 VNAME( 1:TRIMLEN( VNAME ) ) // '" from file.'
                CALL M3MESG( MESG )
            END IF

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

