
        SUBROUTINE WCNTLREP( ODEV, ADEV, CDEV, GDEV, LDEV )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine write the report file for emission controls by source
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: %W%
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
C Pathname: %P%
C Last updated: %G% %U%
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'FLTERR.EXT'    !  functions for comparing two numbers

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER         GETEFILE

        EXTERNAL   GETEFILE

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT (IN) :: ODEV   ! file unit no. for output report
        INTEGER     , INTENT (IN) :: ADEV   ! file unit no. for tmp ADD file
        INTEGER     , INTENT (IN) :: CDEV   ! file unit no. for tmp CTL file 
        INTEGER     , INTENT (IN) :: GDEV   ! file unit no. for tmp CTG file
        INTEGER     , INTENT (IN) :: LDEV   ! file unit no. for tmp ALW file

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: ALWINDX ( :,: ) ! indices to ALW controls table
        INTEGER, ALLOCATABLE :: CTGINDX ( :,: ) ! indices to CTG controls table
        INTEGER, ALLOCATABLE :: CTLINDX ( :,: ) ! indices to CTL controls table
        INTEGER, ALLOCATABLE :: PLTINDX ( : )   ! index from sources to plants

        REAL   , ALLOCATABLE :: BACKOUT ( : )   ! factor used to account for pol
                                                ! specific control info that is
                                                ! already in the inventory
        REAL   , ALLOCATABLE :: CTLEFF  ( : )   ! control efficiency
        REAL   , ALLOCATABLE :: EMIS    ( : )   ! base inventory emissions
        REAL   , ALLOCATABLE :: FACTOR  ( : )   ! multiplicative controls
        REAL   , ALLOCATABLE :: RULEFF  ( : )   ! rule effectiveness
        REAL   , ALLOCATABLE :: RULPEN  ( : )   ! rule penetration
        REAL   , ALLOCATABLE :: PLTINEM ( :,: ) ! initial emissions
        REAL   , ALLOCATABLE :: PLTOUTEM( :,: ) ! controlled emissions

        LOGICAL, ALLOCATABLE :: PLTFLAG ( : )   ! true: plant controlled

C.........  Local arrays
        INTEGER                 OUTTYPES( NVCMULT,6 ) ! var type:int/real

        CHARACTER(LEN=IOVLEN3)  OUTNAMES( NVCMULT,6 ) ! var names
        CHARACTER(LEN=IOULEN3)  OUTUNITS( NVCMULT,6 ) ! var units
        CHARACTER(LEN=IODLEN3)  OUTDESCS( NVCMULT,6 ) ! var descriptions

C...........   Other local variables
        INTEGER          S, V  ! counters and indices

        INTEGER          CIDX     ! control plant index

        REAL             E_IN   ! emissions before controls
        REAL             E_OUT  ! emissions after controls
        REAL             FAC    ! control factor

        CHARACTER*300          MESG       ! message buffer
        CHARACTER(LEN=IOVLEN3) PNAM       ! tmp pollutant name

        CHARACTER*16  :: PROGNAME = 'WCNTLREP' ! program name

C***********************************************************************
C   begin body of subroutine WCNTLREP

C.........  Rewind temporary files
        IF( ADEV .GT. 0 ) REWIND( ADEV )
        IF( CDEV .GT. 0 ) REWIND( CDEV )
        IF( GDEV .GT. 0 ) REWIND( GDEV )
        IF( LDEV .GT. 0 ) REWIND( LDEV )

C.........  For each pollutant that receives controls, obtain variable
C             names for control efficiency, rule effectiveness, and, in the
C             case of AREA sources, rule penetration. These variable names
C             will be used in reading the inventory file.

c note: updated for all pollutants that get controls

        CALL BLDENAMS( CATEGORY, NVCMULT, 6, PNAMMULT, OUTNAMES,
     &                 OUTUNITS, OUTTYPES, OUTDESCS )

C...........  Read in indices from temporary files. No error checking is
C             performed because it is assumed that the program has already
C             successfully written the temporary files.

C.........  Loop through pollutants
C note: must change NVCMULT to other variable for global program pollutants.
        DO V = 1, NVCMULT

C.............  Loop through sources and output 
            DO S = 1, NSRC

C.................  If ADDITIVE packet applies for this pollutant
c note: must be added

C.................  If CONTROL or EMS CONTROL packet applies for this pollutant
                IF( PCTLFLAG( V, 1 ) ) THEN

                    READ( CDEV,* ) CIDX, PNAM, E_IN, E_OUT, FAC
        
                    CALL CONTROL_MESG( ODEV, 'CONTROL', S, CIDX, PNAM, 
     &                                 E_IN, E_OUT, FAC )
            
                END IF

C.................  If CTG packet applies for this pollutant
                IF( PCTLFLAG( V, 2 ) ) THEN
                    READ( GDEV,* ) CIDX, PNAM, E_IN, E_OUT, FAC

                    CALL CONTROL_MESG( ODEV, 'CTG', S, CIDX, PNAM, 
     &                                 E_IN, E_OUT, FAC )

                END IF

C.................  If ALLOWABLE packet applies for this pollutant
                IF( PCTLFLAG( V, 3 ) ) THEN

                    READ( LDEV,* ) CIDX, PNAM, E_IN, E_OUT, FAC

C.....................  Interpret control information and write message
                    CALL CONTROL_MESG( ODEV, 'ALLOWABLE', S, CIDX, PNAM, 
     &                                 E_IN, E_OUT, FAC )

                END IF

C.................  If REACTIVITY packet applies for this pollutant
c note: must be added

            END DO
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram 
            SUBROUTINE CONTROL_MESG( FDEV, CTYPE, SMKID, CIDX, PNAM,
     &                               EMIS_IN, EMIS_OUT, FACTOR )

C.............  Subprogram arguments
            INTEGER     , INTENT (IN) :: FDEV     ! output unit number
            CHARACTER(*), INTENT (IN) :: CTYPE    ! control type
            INTEGER     , INTENT (IN) :: SMKID    ! smoke ID
            INTEGER     , INTENT (IN) :: CIDX     ! control index
            CHARACTER(*), INTENT (IN) :: PNAM     ! pollutant name
            REAL        , INTENT (IN) :: EMIS_IN  ! emissions before control
            REAL        , INTENT (IN) :: EMIS_OUT ! emissions after control
            REAL        , INTENT (IN) :: FACTOR   ! EMIS_IN*FACTOR = EMIS_OUT

C.............  Local variables
            INTEGER           L, L2

            INTEGER, SAVE ::  PSMKID = 0   ! previous smoke ID

            CHARACTER*300     BUFFER                 ! message buffer
            CHARACTER*300     MESG                   ! message buffer

            CHARACTER(LEN=IOVLEN3), SAVE :: LNAM  ! previous pollutant name

C----------------------------------------------------------------------

C.............  Skip records that are not controlled
            IF( CIDX .EQ. 0 ) RETURN

C.............  Write message for pollutant for each new pollutant
            IF( PNAM .NE. LNAM ) THEN

                MESG = 'Controls for pollutant ' // PNAM
                L = LEN_TRIM( MESG )
                WRITE( FDEV,93000 ) ' '
                WRITE( FDEV,93000 ) MESG( 1:L )
                WRITE( FDEV,93000 ) ' '
                LNAM = PNAM

            END IF

C.............  Only write out header line if source is different from previous
C
            IF( SMKID .NE. PSMKID ) THEN

                CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                WRITE( FDEV,93000 ) BLANK5 // BUFFER( 1:L2 )
                PSMKID = SMKID

            END IF

C.............  Write warning if source is controled and factor is 1.0
            IF( FACTOR .EQ. 1. ) THEN
                WRITE( MESG,93380 ) CTYPE, EMIS_IN, EMIS_OUT
                L = LEN_TRIM( MESG )
                WRITE( FDEV,93000 ) MESG( 1:L )

C.............  Otherwise, write standard control message
            ELSE
                WRITE( MESG,93400 ) CTYPE, EMIS_IN, EMIS_OUT, FACTOR
                L = LEN_TRIM( MESG )
                WRITE( FDEV,93000 ) MESG( 1:L )

            END IF

            RETURN

C*************  SUBPROGRAM FORMAT  STATEMENTS   **************************

C...........   Formatted file I/O formats............ 93xxx

93000       FORMAT( A )

93380       FORMAT( 10X, A, ' Packet. Before: ', F10.3, ' After: ', 
     &              F10.3, ' [tons/yr]. WARNING: Control factor of 1.' )

93400       FORMAT( 10X, A, ' Packet. Before: ', F10.3, ' After: ',  
     &              F10.3, ' [tons/yr]. Factor:', F5.2 )

            END SUBROUTINE CONTROL_MESG

C----------------------------------------------------------------------

        END SUBROUTINE WCNTLREP
