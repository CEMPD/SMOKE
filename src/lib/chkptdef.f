
        SUBROUTINE CHKPTDEF( CHKNCHAR, CHKJSCC )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine compares the source definition that is input through
C      the subroutine arguments to the defition from the inventory file 
C      that is stored in MODINFO.  This is used to check to ensure the
C      cross-reference files and the inventory use the same definition of
C      point sources.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Create 3/99 by M. Houyoux
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C.........  INCLUDE FILES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........  EXTERNAL FUNCTIONS
        CHARACTER*2   CRLF
        EXTERNAL      CRLF

C.........  SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: CHKNCHAR     ! number of chars in source def
        INTEGER, INTENT (IN) :: CHKJSCC      ! position of SCC in source def

C.........  Other local variables
        LOGICAL      :: EFLAG = .FALSE.      ! true: error found

        CHARACTER*300   MESG 

        CHARACTER*16 :: PROGNAME = 'CHKPTDEF' ! program name

C***********************************************************************
C   begin body of subroutine CHKPTDEF

C.........  Compare the input total number of sources with the inventory
        IF( CHKNCHAR .NE. NCHARS ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'ERROR: Source definition mismatch. Number of ' //
     &             'characteristics in file: ', CHKNCHAR, ', ' // 
     &             CRLF() // BLANK10 // 'but in ' // CATDESC // 
     &             ' source inventory file: ', NCHARS
            CALL M3MSG2( MESG )

        END IF

C.........  Compare the input position of SCCs with the inventory
        IF( CHKJSCC .NE. JSCC ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )
     &             'ERROR: Source definition mismatch. Position of ' //
     &             'SCC in characteristics in file: ', CHKJSCC, ', ' //
     &             CRLF() // BLANK10 // 'but in ' // CATDESC // 
     &             ' source inventory file: ', JSCC
            CALL M3MSG2( MESG )

        END IF 

C.........  If there is an error, exit
        IF( EFLAG ) THEN

            MESG = 'Inconsistent input files.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE CHKPTDEF

