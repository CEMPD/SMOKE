
        PROGRAM TMPBEIS3

C***********************************************************************
C  program body starts at line  187
C
C  DESCRIPTION:
C       Computes hourly time stepped gridded biogenic emissions using 
C       normalized gridded emissions from Normbeis3 and postprocessed MM5
C       meteorology.
C
C  PRECONDITIONS REQUIRED:
C       Postprocessed MM5 meteorology that contains temperature, 
C       solar radiation, and pressure data. 
C       Normalized gridded emissions B3GRD from Normbeis3 
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       HRBIO, PREBMET
C
C  REVISION  HISTORY:
C      3/01: Prototype by Jeff Vukovich
C            Tested only on 36km Lambert domain 
C            Summer/winter switch file option not tested
C      8/04: Incorporated BEIS v3.12 by C. Seppanen
C      4/06: changed to BEIS3.13 by G. Pouliot
C      3/08: changed to BEIS3.14 by G. Pouliot
C      ?/14: changed to BEIS3.60 by G. Pouliot
C      7/15  Changed from BEIS3.60 to BEIS3.61 by Baek
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

        IMPLICIT NONE

C.........  Includes
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

C.........  Local parameters
        CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv5.2.1_Sep2025$' ! CVS release tag
        
C.........  External functions
        LOGICAL, EXTERNAL :: ENVYN

C.........  Logical names and unit numbers
        INTEGER         LDEV    !  unit number for log device

C.........  Other local variables
        INTEGER         IOS     !  I/O status

        CHARACTER(16)   BEISVER !  version of BEIS3 to use
        CHARACTER(300)  MESG    !  message buffer for M3EXIT()

        CHARACTER(16) :: PROGNAME = 'TMPBEIS3'   !  program name

C***********************************************************************
C   begin body of program TMPBEIS3

        LDEV = INIT3()
 
C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get the BEIS3 model version to use
        MESG = 'Version of BEIS3 to use'
        CALL ENVSTR( 'BEIS_VERSION', MESG, '3.7', BEISVER, IOS )
        
        SELECT CASE( BEISVER )
        CASE( '3.7' )
            CALL TMPBEIS360( CVSW )
        CASE( '3.61' )
            CALL TMPBEIS360( CVSW )
        CASE( '3.14' )
            CALL TMPBEIS314( CVSW )
        CASE DEFAULT
            MESG = 'ERROR: Unrecognized BEIS_VERSION setting; valid ' //
     &             'settings are 3.14, 3.61 and 3.7'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END SELECT

C.........   End of program
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

        END PROGRAM TMPBEIS3  

