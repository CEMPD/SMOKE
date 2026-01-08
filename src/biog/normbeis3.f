
        PROGRAM NORMBEIS3

C***********************************************************************
C
C  DESCRIPTION:  Produces normalized biogenic emissions for use with
C                SMOKE-BEIS versions 3.14 and 3.6
C
C  SUBROUTINES AND FUNCTIONS CALLED: Calls Normbeis314 or Normbeis360
C
C  REVISION  HISTORY: 3/00 Prototype, Jeff Vukovich
C                     8/04 Integrated v3.12, C. Seppanen
C                     4/06 Changed Beis3.12 to BEIS3.13 G. Pouliot
C                     3/08 Changed Beis3.13 to BEIS3.14 G. Pouliot
C                     ?/14 Changed Beis3.14 to BEIS3.60 G. Pouliot
C                     7/15 Changed BEIS3.60 to BEIS3.61 by Baek
C                     8/20 Added BEIS3.70 with BELD5 and BEISFACS
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C.........  Includes
C        INCLUDE 'IODECL3.EXT'   ! I/O API function declarations
        INCLUDE 'IOSTRG3.EXT'     !
        
C.........  Local parameters
C       CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv5.2.1_Sep2025$' ! CVS release tag
  
C.........  External functions
C       LOGICAL, EXTERNAL :: ENVYN

C.........  Logical names and unit numbers
        INTEGER         LDEV    !  unit number for log device

C.........  Other local variables
        INTEGER         IOS     !  I/O status
       
        CHARACTER(16)   BEISVER !  version of BEIS3 to use
        CHARACTER(300)  MESG    !  message buffer
        
        CHARACTER(16) :: PROGNAME = 'NORMBEIS3'   !  program name

C***********************************************************************
C   begin body of program NORMBEIS3

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get the BEIS3 model version to use
        MESG = 'Version of BEIS3 to use'
        CALL ENVSTR( 'BEIS_VERSION', MESG, '3.7', BEISVER, IOS )
        
        SELECT CASE( BEISVER )
        CASE( '3.7' )
            CALL NORMBEIS370(  )
        CASE( '3.61' )
            CALL NORMBEIS360(  )
        CASE( '3.14' )
            CALL NORMBEIS312(  )
        CASE DEFAULT
            MESG = 'ERROR: Unrecognized BEIS_VERSION setting; valid ' //
     &             'settings are 3.14, 3.61, and 3.7'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END SELECT

C.........  End of program
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

        END PROGRAM NORMBEIS3 

