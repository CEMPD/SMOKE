
        SUBROUTINE ASGNGRPS( SRCARRAY )

C***********************************************************************
C  subroutine body starts at line 78
C
C  DESCRIPTION:
C       Assigns temporal averaging type to each source
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
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
C***********************************************************************

C...........   MODULES for public variables
        USE MODMBSET, ONLY: DAILY, WEEKLY, MONTHLY, EPISLEN
        
        USE MODMET, ONLY: DYCODES, WKCODES, MNCODES, EPCODES
                
        USE MODINFO, ONLY: NSRC
        
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER, EXTERNAL :: STR2INT
        
C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: SRCARRAY( NSRC,2 )  ! array to hold county codes

C...........   Local arrays
        INTEGER IDX( NSRC )               ! index to sort SRCARRAY
        INTEGER ARRAYPOS( 4 )             ! position in array for each averaging type
        
C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices                     
        
        INTEGER CURRAVER                  ! current averaging type
        INTEGER CURRCNTY                  ! current county
        INTEGER PREVCNTY                  ! previous county

        LOGICAL :: EFLAG      = .FALSE.   ! true: error found
        
        CHARACTER(LEN=100)     LINE     !  line buffer
        CHARACTER(LEN=300)     MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'ASGNGRPS'   ! program name
        
C***********************************************************************
C   begin body of subroutine ASGNGRPS

C.........  Sort county-by-src array by averaging type and county
        DO I = 1, NSRC
            IDX( I ) = I
        END DO
            
        CALL SORTI1( NSRC, IDX, SRCARRAY( :,1 ) )

C.........  Initialize variables
        ARRAYPOS = 0   ! array
        PREVCNTY = 0

        DO I = 1, NSRC
     
C.............  Skip unused sources
            IF( SRCARRAY( IDX( I ),2 ) == 0 ) CYCLE
        
            CURRAVER = SRCARRAY( IDX( I ),2 )            
            CURRCNTY = SRCARRAY( IDX( I ),1 )
            
            IF( CURRCNTY == PREVCNTY ) CYCLE
            
            ARRAYPOS( CURRAVER ) = ARRAYPOS( CURRAVER ) + 1
            
            SELECT CASE ( CURRAVER )
            CASE( DAILY )
                DYCODES( ARRAYPOS( CURRAVER ) ) = SRCARRAY( IDX( I ),1 )
            CASE( WEEKLY )
                WKCODES( ARRAYPOS( CURRAVER ) ) = SRCARRAY( IDX( I ),1 )
            CASE( MONTHLY )
                MNCODES( ARRAYPOS( CURRAVER ) ) = SRCARRAY( IDX( I ),1 )
            CASE( EPISLEN )
                EPCODES( ARRAYPOS( CURRAVER ) ) = SRCARRAY( IDX( I ),1 )
            CASE DEFAULT
                MESG = 'Unrecognized temporal averaging value'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END SELECT 
            
            PREVCNTY = CURRCNTY
            
        END DO

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( I6, 1X, I6, 1X, I1 ) 

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )     
        
        END SUBROUTINE ASGNGRPS
        
