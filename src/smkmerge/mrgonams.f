
        SUBROUTINE MRGONAMS

C***********************************************************************
C  subroutine MRGONAMS body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to set the default names of the
C      output files. It sets them regardless of whether they will be used
C      or not.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C
C***********************************************************************
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
C****************************************************************************

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: ARFLAG, MRFLAG, PRFLAG,
     &                      AUFLAG, MUFLAG, PUFLAG,
     &                      TFLAG, SFLAG, LFLAG,
     &                      LREPSPC, LREPCTL, SUBOUTNAME,
     &                      AREPNAME, BREPNAME, MREPNAME, PREPNAME, 
     &                      TREPNAME, NGRPS, IGRPIDX, SUBSECFLAG,
     &                      AONAME, BONAME, MONAME, PONAME, TONAME,
     &                      PINGNAME, INLINENAME, PELVNAME, 
     &                      SGINLNNAME, NUNITS, GRDUNIT

        IMPLICIT NONE

C.........  EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        LOGICAL            ENVYN
        
        EXTERNAL        ENVYN

        INTEGER         I, J, N, K    ! indices and counters

        INTEGER         IOS           ! tmp I/O status
 
        LOGICAL      :: KFLAG = .FALSE.          ! true: simplied output file names
        LOGICAL      :: MOLEFLAG = .FALSE.       ! true: outputting moles

        CHARACTER(300)  MESG    ! message buffer

        CHARACTER(16) :: PROGNAME = 'MRGONAMS' ! program name

C***********************************************************************
C   begin body of subroutine MRGONAMS

C.........  Retrieve variable to indicate whether to use annual or average day data
        MESG = 'Use customized SMKMERGE output file names'
        KFLAG = ENVYN( 'SMKMERGE_CUSTOM_OUTPUT', MESG, .FALSE., IOS )

C.........  Set default output file name(s) depending on inputs.  Set both
C           report file(s) and gridded file(s)...

C.........  Initialize - everything will be gridded
        N = 1
        IF( SUBSECFLAG ) N = NGRPS
        ALLOCATE( SUBOUTNAME( N ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUBOUTNAME', PROGNAME )
        SUBOUTNAME = 'SUBOUT'

        DO K = 1, N

          IF( KFLAG ) THEN

            AREPNAME = 'REPAG'
            BREPNAME = 'REPBG'
            MREPNAME = 'REPMG'
            PREPNAME = 'REPPG'
            TREPNAME = 'REPEG'
    
            AONAME = 'AOUT'
            BONAME = 'BOUT'
            MONAME = 'MOUT'
            PONAME = 'POUT'
            TONAME = 'EOUT'
    
            PINGNAME = 'PING'
            PELVNAME = 'ELEV'
            INLINENAME = 'INLN'
            SGINLNNAME = 'SGINLN'
            IF( SUBSECFLAG ) WRITE( SUBOUTNAME( K ), '(A,I2.2)' ) 'SUBOUT', IGRPIDX( K )

        ELSE

            AREPNAME = 'REPAG'
            BREPNAME = 'REPBG'
            MREPNAME = 'REPMG'
            PREPNAME = 'REPPG'
            TREPNAME = 'REPEG'
    
            AONAME = 'AG'
            BONAME = 'BGTS'
            MONAME = 'MG'
            PONAME = 'PG'
            TONAME = 'EG'
    
            PINGNAME = 'PING'
            INLINENAME = 'INLN'
            PELVNAME = 'ELEV'
            SGINLNNAME = 'SGINLN'
            IF( SUBSECFLAG ) WRITE( SUBOUTNAME( K ), '(A,I2.2)' ) 'SUBOUT', IGRPIDX( K )
 
            IF( TFLAG ) THEN
    
                CALL TRIM_AND_CONCAT( AREPNAME, 'T' )
                CALL TRIM_AND_CONCAT( BREPNAME, 'T' )
                CALL TRIM_AND_CONCAT( MREPNAME, 'T' )
                CALL TRIM_AND_CONCAT( PREPNAME, 'T' )
                CALL TRIM_AND_CONCAT( TREPNAME, 'T' )
    
                CALL TRIM_AND_CONCAT( AONAME, 'T' )
                CALL TRIM_AND_CONCAT( MONAME, 'T' )
                CALL TRIM_AND_CONCAT( PONAME, 'T' )
                CALL TRIM_AND_CONCAT( TONAME, 'T' )
                CALL TRIM_AND_CONCAT( PINGNAME, 'T' )
                CALL TRIM_AND_CONCAT( INLINENAME, 'T' )
                CALL TRIM_AND_CONCAT( PELVNAME, 'T' )
                CALL TRIM_AND_CONCAT( SGINLNNAME, 'T' )
                CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), 'T' )
    
            END IF
    
            IF( SFLAG ) THEN
    
                IF( LREPSPC ) THEN
                    CALL TRIM_AND_CONCAT( AREPNAME, 'S' )
                    CALL TRIM_AND_CONCAT( BREPNAME, 'S' )
                    CALL TRIM_AND_CONCAT( MREPNAME, 'S' )
                    CALL TRIM_AND_CONCAT( PREPNAME, 'S' )
                    CALL TRIM_AND_CONCAT( TREPNAME, 'S' )
                ENDIF
    
                CALL TRIM_AND_CONCAT( AONAME, 'S' )
                CALL TRIM_AND_CONCAT( MONAME, 'S' )
                CALL TRIM_AND_CONCAT( PONAME, 'S' )
                CALL TRIM_AND_CONCAT( TONAME, 'S' )
                CALL TRIM_AND_CONCAT( PINGNAME, 'S' )
                CALL TRIM_AND_CONCAT( INLINENAME, 'S' )                
                CALL TRIM_AND_CONCAT( PELVNAME, 'S' )
                CALL TRIM_AND_CONCAT( SGINLNNAME, 'S' )
                CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), 'S' )
    
            END IF
    
            IF( AUFLAG .OR. ARFLAG ) THEN
    
                IF( LREPCTL ) THEN
                    CALL TRIM_AND_CONCAT( AREPNAME, 'C' )
                END IF
    
                CALL TRIM_AND_CONCAT( AONAME, 'C' )
    
            END IF
    
            IF( MUFLAG .OR. MRFLAG ) THEN
    
                IF( LREPCTL ) THEN
                    CALL TRIM_AND_CONCAT( MREPNAME, 'C' )
                END IF
    
                CALL TRIM_AND_CONCAT( MONAME, 'C' )
    
            END IF
    
            IF( PUFLAG .OR. PRFLAG ) THEN
    
                IF( LREPCTL ) THEN
                    CALL TRIM_AND_CONCAT( PREPNAME, 'C' )
                END IF
    
                CALL TRIM_AND_CONCAT( PONAME, 'C' )
                CALL TRIM_AND_CONCAT( PINGNAME, 'C' )
                CALL TRIM_AND_CONCAT( INLINENAME, 'C' )
                CALL TRIM_AND_CONCAT( PELVNAME, 'C' )
                CALL TRIM_AND_CONCAT( SGINLNNAME, 'C' )
                CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), 'C' )
    
            END IF
    
            IF( AUFLAG .OR. ARFLAG .OR.
     &          MUFLAG .OR. MRFLAG .OR.
     &          PUFLAG .OR. PRFLAG     ) THEN
                IF( LREPCTL ) THEN
                    CALL TRIM_AND_CONCAT( TREPNAME, 'C' )
                END IF
    
                CALL TRIM_AND_CONCAT( TONAME, 'C' )
    
            END IF
    
            IF( LFLAG ) THEN
    
                CALL TRIM_AND_CONCAT( PONAME, '3D' )
                CALL TRIM_AND_CONCAT( TONAME, '3D' )
    
            END IF
    
C.............  Now append mass or mole for I/O API files, depending on which
C               inputs were used
            IF( SFLAG ) THEN
    
C.................  Set flag if any of the output species are mole-based
                DO I = 1, NUNITS
                
                    J = INDEX( GRDUNIT( I ), 'mole' )
                    IF( J .GT. 0 ) THEN
                        MOLEFLAG = .TRUE.
                        EXIT
                    END IF
    
                END DO
    
C.................  Get output file names depending on if there are moles in units
                IF( MOLEFLAG ) THEN 
    
                    CALL TRIM_AND_CONCAT( AREPNAME, '_L' )
                    CALL TRIM_AND_CONCAT( BREPNAME, '_L' )
                    CALL TRIM_AND_CONCAT( MREPNAME, '_L' )
                    CALL TRIM_AND_CONCAT( PREPNAME, '_L' )
                    CALL TRIM_AND_CONCAT( TREPNAME, '_L' )
    
                    CALL TRIM_AND_CONCAT( AONAME, '_L' )
                    CALL TRIM_AND_CONCAT( BONAME, '_L_O' )
                    CALL TRIM_AND_CONCAT( MONAME, '_L' )
                    CALL TRIM_AND_CONCAT( PONAME, '_L' )
                    CALL TRIM_AND_CONCAT( TONAME, '_L' )
                    CALL TRIM_AND_CONCAT( PINGNAME, '_L' )
                    CALL TRIM_AND_CONCAT( INLINENAME, '_L' )
                    CALL TRIM_AND_CONCAT( PELVNAME, '_L' )
                    CALL TRIM_AND_CONCAT( SGINLNNAME, '_L' )
                    CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), '_L' )
    
                ELSE 
    
                    CALL TRIM_AND_CONCAT( AREPNAME, '_S' )
                    CALL TRIM_AND_CONCAT( BREPNAME, '_S' )
                    CALL TRIM_AND_CONCAT( MREPNAME, '_S' )
                    CALL TRIM_AND_CONCAT( PREPNAME, '_S' )
                    CALL TRIM_AND_CONCAT( TREPNAME, '_S' )
    
                    CALL TRIM_AND_CONCAT( AONAME, '_S' )
                    CALL TRIM_AND_CONCAT( BONAME, '_S_O' )
                    CALL TRIM_AND_CONCAT( MONAME, '_S' )
                    CALL TRIM_AND_CONCAT( PONAME, '_S' )
                    CALL TRIM_AND_CONCAT( TONAME, '_S' )
                    CALL TRIM_AND_CONCAT( PINGNAME, '_S' )
                    CALL TRIM_AND_CONCAT( INLINENAME, '_S' )
                    CALL TRIM_AND_CONCAT( PELVNAME, '_S' )
                    CALL TRIM_AND_CONCAT( SGINLNNAME, '_S' )
                    CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), '_S' )
    
                END IF
    
            END IF

          END IF

        END DO
    
        RETURN

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal subprogram trims and concatonates the first
C               string with the second string
            SUBROUTINE TRIM_AND_CONCAT( PART1, PART2 )

C.............  Subprogram arguments
            CHARACTER(*)  PART1
            CHARACTER(*)  PART2

C.............  Local variables
            INTEGER       L

C----------------------------------------------------------------------

            L = LEN_TRIM( PART1 )
            PART1 = PART1( 1:L ) // PART2
 
            END SUBROUTINE TRIM_AND_CONCAT

        END SUBROUTINE MRGONAMS
