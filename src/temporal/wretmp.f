
        SUBROUTINE WRETMP( OUTFILE, NPELV, NPSRC, TMPEIDX, EMIST, 
     &                     NEOUT, INDXE, EMISE )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reallocates the elevated emissions storage to the 
C      types needed for the sparse output, copies the emissions from the
C      main array, and writes the output emissions.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C      Functions: I/O API functions
C
C  REVISION  HISTORY:
C      Created by M. Houyoux 11/98
C
C****************************************************************************/
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

        IMPLICIT NONE

C...........   INCLUDES

c        INCLUDE 'EMCNST3.EXT'   !  emissions parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
c        INTEGER         LBLANK
c        INTEGER         TRIMLEN

c        EXTERNAL LBLANK, TRIMLEN

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*(*)   OUTFILE          ! logical name of elevated output file
        INTEGER         NPELV            ! exact number of elevated entries
        INTEGER         NPSRC            ! number of point sources
        INTEGER         TMPEIDX( NPELV ) ! 
        REAL            EMIST  ( NPSRC )
        INTEGER         NEOUT
        INTEGER         INDXE  ( NPELV )
        REAL            EMISE  ( NPELV )

        CHARACTER*16 :: PROGNAME = 'WRETMP' !  program name

C***********************************************************************
C   begin body of subroutine WRETMP


        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93100   FORMAT( I2, ', "', A, '"' )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2.2 ) )

        END
