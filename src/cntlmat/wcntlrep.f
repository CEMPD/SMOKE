
        SUBROUTINE WCNTLREP( CDEV, GDEV, LDEV, MDEV )

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
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: %W%
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
C Pathname: %P%
C Last updated: %G% %U%
C
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CSOURC

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: NVCMULT, PNAMMULT, RPTDEV, PCTLFLAG

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CRL, NSRC, NCHARS

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       INTEGER         GETEFILE
C       INTEGER         PROMPTFFILE

C        EXTERNAL   GETEFILE, PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT (IN) :: CDEV   ! file unit no. for tmp CTL file 
        INTEGER     , INTENT (IN) :: GDEV   ! file unit no. for tmp CTG file
        INTEGER     , INTENT (IN) :: LDEV   ! file unit no. for tmp ALW file
        INTEGER     , INTENT (IN) :: MDEV   ! file unit no. for tmp MACT file

C.........  Local arrays
        INTEGER             OUTTYPES( NVCMULT,6 ) ! var type:int/real

        CHARACTER(IOVLEN3)  OUTNAMES( NVCMULT,6 ) ! var names
        CHARACTER(IOULEN3)  OUTUNITS( NVCMULT,6 ) ! var units
        CHARACTER(IODLEN3)  OUTDESCS( NVCMULT,6 ) ! var descriptions

C...........   Other local variables
        INTEGER          S, V  ! counters and indices

        INTEGER          CIDX   ! control plant index
        INTEGER          ODEV   ! file unit no. for output report

        REAL             E_IN   ! emissions before controls
        REAL             E_OUT  ! emissions after controls
        REAL             FAC    ! control factor

        CHARACTER(256)     MESG       ! message buffer
        CHARACTER(IOVLEN3) PNAM       ! tmp pollutant name

        CHARACTER(16) :: PROGNAME = 'WCNTLREP' ! program name

C***********************************************************************
C   begin body of subroutine WCNTLREP

C.........  Rewind temporary files
        IF( CDEV .GT. 0 ) REWIND( CDEV )
        IF( GDEV .GT. 0 ) REWIND( GDEV )
        IF( LDEV .GT. 0 ) REWIND( LDEV )
        IF( MDEV .GT. 0 ) REWIND( MDEV )

C.........  Open reports file
        IF( MAX( CDEV, GDEV, LDEV, MDEV ) .GT. 0 ) THEN
            RPTDEV( 2 ) = PROMPTFFILE( 
     &                'Enter logical name for SUMMARY ' //
     &                'CONTROLS REPORT',
     &                .FALSE., .TRUE., CRL // 'CSUMREP', PROGNAME )
            ODEV = RPTDEV( 2 )
        END IF

C.........  For each pollutant that receives controls, obtain variable
C             names for control efficiency, rule effectiveness, and, in the
C             case of AREA sources, rule penetration. These variable names
C             will be used in reading the inventory file.
        
C.........  Check that NVCMULT does not equal 0, otherwise some systems will get confused
        IF( NVCMULT == 0 ) RETURN

        CALL BLDENAMS( CATEGORY, NVCMULT, 6, PNAMMULT, OUTNAMES,
     &                 OUTUNITS, OUTTYPES, OUTDESCS )

C...........  Read in indices from temporary files. No error checking is
C             performed because it is assumed that the program has already
C             successfully written the temporary files.

C.........  Loop through pollutants
        DO V = 1, NVCMULT

C.............  Loop through sources and output 
            DO S = 1, NSRC

C.................  If MACT packet applies for this pollutant
                IF( PCTLFLAG( V, 4 ) ) THEN
                
                    READ( MDEV,* ) CIDX, PNAM, E_IN, E_OUT, FAC
                    
                    CALL CONTROL_MESG( ODEV, 'MACT', S, CIDX, PNAM,
     &                                 E_IN, E_OUT, FAC )
     
                END IF

C.................  If CONTROL packet applies for this pollutant
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

            CHARACTER(300)    BUFFER                 ! message buffer
            CHARACTER(300)    MESG                   ! message buffer

            CHARACTER(IOVLEN3), SAVE :: LNAM  ! previous pollutant name

C----------------------------------------------------------------------

C.............  Skip records that are not controlled
            IF( CIDX .EQ. 0 ) RETURN

C.............  Write message for pollutant for each new pollutant
            IF( PNAM .NE. LNAM ) THEN

                MESG = 'Source controls for pollutant ' // PNAM
                L = LEN_TRIM( MESG )
                WRITE( FDEV,93000 ) ' '
                WRITE( FDEV, 93000 ) REPEAT( '-', 80 )
                WRITE( FDEV,93000 ) MESG( 1:L )
                WRITE( FDEV,93000 ) ' '

            END IF

C.............  Only write out header line if source is different from previous
C               or if the message for a new pollutant was written
            IF( SMKID .NE. PSMKID .OR. PNAM .NE. LNAM ) THEN

                CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                WRITE( FDEV,93000 ) BLANK5 // BUFFER( 1:L2 )
                PSMKID = SMKID

            END IF
            
            LNAM = PNAM

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

93380       FORMAT( 10X, A, ' Packet. Before: ', E11.4, ' After: ',
     &              E11.4, 
     &              ' [tons/day]. WARNING: Control factor of 1.' )

93400       FORMAT( 10X, A, ' Packet. Before: ', E11.4, ' After: ',  
     &              E11.4, ' [tons/day]. Factor:', E9.3 )

            END SUBROUTINE CONTROL_MESG

C----------------------------------------------------------------------

        END SUBROUTINE WCNTLREP
