
        SUBROUTINE PROCEMSPT( FDEV, FNAME, CFLAG, WFLAG )

C***********************************************************************
C  subroutine body starts at line 232
C
C  DESCRIPTION:
C      This subroutine reads the EMS-95 point files and processes the
C      source characteristics.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Copied from rdemspt.f by C. Seppanen (2/03)
C
C****************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: LSTSTR

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         GETINVYR
        
        EXTERNAL        CRLF, GETINVYR

C...........   SUBROUTINE ARGUMENTS
        INTEGER,          INTENT (IN) :: FDEV   ! unit no. of inv file
        CHARACTER(LEN=*), INTENT (IN) :: FNAME  ! logical name of inv file
        LOGICAL,          INTENT (IN) :: CFLAG  ! true: recalc vel w/ flow and diam
        LOGICAL,          INTENT (IN) :: WFLAG  ! true: convert lat-lons to western hemisphere

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: UTMZONE( : )    ! UTM zone from facility file

C...........   Other local variables
        INTEGER              I                      ! counter
        INTEGER              IOS                    ! I/O status

        CHARACTER(LEN=8)     PFILTYPE               ! previous file type
        CHARACTER(LEN=300)   INFILE                 ! inventory file name
        CHARACTER(LEN=300)   MESG                   ! message buffer

        CHARACTER(LEN=16) :: PROGNAME = 'PROCEMSPT' ! program name

C***********************************************************************
C   begin body of subroutine PROCEMSPT

C.........  Allocate local arrays
        ALLOCATE( UTMZONE( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UTMZONE', PROGNAME )
        UTMZONE = 0

C.........  Set previous file type to stack (as if we've read a complete set of files)
        PFILTYPE = 'STACK'

C.........  Loop over all files in inventory list file
        DO I = 1, SIZE( LSTSTR )

C.............  Don't read emission file since it will be read in RDINVDATA
            IF( PFILTYPE == 'DEVICE' ) THEN
                PFILTYPE = 'EMISSION'
                CYCLE
            END IF

            INFILE = LSTSTR( I )
            
C.............  Check for INVYEAR packet
            IF( GETINVYR( INFILE ) > 0 ) CYCLE

C.............  Open current file
            OPEN( FDEV, FILE=INFILE, STATUS='OLD', IOSTAT=IOS )

            IF( IOS /= 0 ) THEN
                    
                WRITE( MESG,94010 ) 'Problem at line ', I, 'of ' // 
     &                 TRIM( FNAME ) // '.' // ' Could not open ' //
     &                 'file:' // CRLF() // BLANK5 // TRIM( INFILE ) 
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        
            ELSE
                WRITE( MESG,94010 ) 'Successful OPEN for ' //
     &                 'inventory file:' // CRLF() // BLANK5 // 
     &                 TRIM( INFILE )
                CALL M3MSG2( MESG ) 
        
            END IF

C.............  Call correct reader based on previous file type
            SELECT CASE( PFILTYPE )
            
            CASE( 'EMISSION' )
                PFILTYPE = 'FACILITY'
                CALL RDFACEMSPT( FDEV, UTMZONE )
                
            CASE( 'FACILITY' )
                PFILTYPE = 'PROCESS'
                CALL RDPROCEMSPT( FDEV )
            
            CASE( 'PROCESS' )
                PFILTYPE = 'STACK'
                CALL RDSTKEMSPT( FDEV, CFLAG, WFLAG, UTMZONE )
            
            CASE( 'STACK' )
                PFILTYPE = 'DEVICE'
                CALL RDDEVEMSPT( FDEV )
            
            END SELECT
            
            CLOSE( FDEV )          
            
        END DO  ! loop over all files in inventory list

C.........  Deallocate local memory
        DEALLOCATE( UTMZONE )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE PROCEMSPT
