
        SUBROUTINE WREPINVEN( ADEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine writes out the REPINVEN file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 11/2002 by A. Holland
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel gsions (SMOKE) Modeling
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

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS
        
C.........  This module contains the arrays for the area-to-point x-form
        USE MODAR2PT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        
C...........   EXTERNAL FUNCTIONS and their descriptions:

	INTEGER		INDEX1
        
        EXTERNAL	INDEX1

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT (IN) :: ADEV  ! file unit no. for REPINVEN file

C...........   Local variables

	CHARACTER*1	KEEP
        CHARACTER(LEN=DDSLEN3)	DESC
        
        INTEGER		I, K

C...........   Other local variables

	CHARACTER*300	MESG

        CHARACTER*16  :: PROGNAME = 'WREPINVEN' ! program name

C***********************************************************************
C   begin body of subroutine WREPINVEN


C.........  Write out first report to REPINVEN file

        WRITE( ADEV, 93010 ) 'CAS Code', 'Keep', 'Nrecs', 
     &         'Emissions', 'CAS description'
     
        WRITE( ADEV, 93020 ) '[tons/year]'
        
        WRITE( ADEV, 93000 ) REPEAT( '-', 85 )
          
        DO I = 1, NUNIQCAS
          IF( UCASNPOL( I ) .EQ. UCASNKEP( I ) ) THEN
            KEEP = 'Y'
          ELSE IF( UCASNPOL( I ) .NE. UCASNKEP( I ) .AND.
     &             UCASNKEP( I ) .NE. 0 ) THEN
            KEEP = 'P'
          ELSE IF( UCASNKEP( I ) .EQ. 0 ) THEN
            KEEP = 'N'
          END IF
            
          K = INDEX1( UNIQCAS( I ), NINVTBL, ITCASA )
          IF( K .GT. 0 ) THEN
            DESC = ITCASDSCA( K )
          ELSE
            MESG = 'CAS code not found in raw list'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 1 )
          END IF
          
          WRITE( ADEV, 93030 ) UNIQCAS( I ), KEEP, RECSBYCAS( I ),
     &             EMISBYCAS( I ), DESC

	END DO
          
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( 1X, A8, 4X, A4, 4X, A5, 6X, A9, 9X, A15 )

93020   FORMAT( 32X, A11 )

93030   FORMAT( 1X, A8, 4X, A1, 7X, I5, 4X, F16.10, 4X, A40 )



        END SUBROUTINE WREPINVEN
