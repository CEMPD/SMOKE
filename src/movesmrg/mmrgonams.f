
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
        USE MODMERGE, ONLY: MREPNAME, MONAME, SGINLNNAME, SUBOUTNAME,
     &                      NUNITS, GRDUNIT, NGRPS, IUGRPNUM, SUBSECFLAG

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
        MESG = 'Use customized MOVESMRG output file names'
        KFLAG = ENVYN( 'MOVESMRG_CUSTOM_OUTPUT', MESG, .FALSE., IOS )

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

            MREPNAME = 'REPMG'
            MONAME = 'MOUT'
            SGINLNNAME = 'SGINLN'
            IF( SUBSECFLAG ) WRITE( SUBOUTNAME( K ), '(A,I2.2)' ) 'SUBOUT', IUGRPNUM( K )

          ELSE

            MREPNAME = 'REPMG'
            MONAME = 'MG'
            SGINLNNAME = 'SGINLN'
            IF( SUBSECFLAG ) WRITE( SUBOUTNAME( K ), '(A,I2.2)' ) 'SUBOUT', IUGRPNUM( K )
    
            CALL TRIM_AND_CONCAT( MREPNAME, 'T' )
            CALL TRIM_AND_CONCAT( MONAME, 'T' )
            CALL TRIM_AND_CONCAT( SGINLNNAME, 'T' )
            CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), 'T' )
    
            CALL TRIM_AND_CONCAT( MREPNAME, 'S' )
            CALL TRIM_AND_CONCAT( MONAME, 'S' )
            CALL TRIM_AND_CONCAT( SGINLNNAME, 'S' )
            CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), 'S' )
    
C.............  Now append mass or mole for I/O API files, depending on which
C               inputs were used
    
C.............  Set flag if any of the output species are mole-based
            DO I = 1, NUNITS
            
                J = INDEX( GRDUNIT( I ), 'mole' )
                IF( J .GT. 0 ) THEN
                    MOLEFLAG = .TRUE.
                    EXIT
                END IF
    
            END DO
    
C.............  Get output file names depending on if there are moles in units
            IF( MOLEFLAG ) THEN 
    
                CALL TRIM_AND_CONCAT( MREPNAME, '_L' )
                CALL TRIM_AND_CONCAT( MONAME, '_L' )
                CALL TRIM_AND_CONCAT( SGINLNNAME, '_L' )
                CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), '_L' )
    
            ELSE 
    
                CALL TRIM_AND_CONCAT( MREPNAME, '_S' )
                CALL TRIM_AND_CONCAT( MONAME, '_S' )
                CALL TRIM_AND_CONCAT( SGINLNNAME, '_S' )
                CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), '_S' )
    
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
