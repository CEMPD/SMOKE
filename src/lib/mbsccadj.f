
        SUBROUTINE MBSCCADJ( IREC, TSCC, CRWT, CVID, TSCCINTL, EFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      Separates the parts of the mobile source SCC and creates an internal
C      SCC used by SMOKE for applying cross-reference data.
C
C  PRECONDITIONS REQUIRED:
C      RDMVINFO subroutine called previously
C      TSCC is defined properly
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 10/99 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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
C.........  This module is for mobile-specific data
        USE MODMOBIL

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FIND1
        INTEGER         STR2INT

        EXTERNAL    CRLF, FIND1, STR2INT

C.........   SUBROUTINE ARGUMENTS
        INTEGER               , INTENT (IN) :: IREC     ! line number
        CHARACTER(LEN=SCCLEN3), INTENT (IN) :: TSCC     ! external SCC
        CHARACTER(LEN=RWTLEN3), INTENT(OUT) :: CRWT     ! roadway type no.
        CHARACTER(LEN=VIDLEN3), INTENT(OUT) :: CVID     ! vehicle type ID no.
        CHARACTER(LEN=SCCLEN3), INTENT(OUT) :: TSCCINTL ! internal SCC
        LOGICAL               , INTENT(OUT) :: EFLAG    ! true: error found

C.........  Local variables
        INTEGER         K               !  find1 index
        INTEGER         RCL             !  tmp road class

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called

        CHARACTER*10 , SAVE :: RWTFMT   !  frmt to write roadway type to string
        CHARACTER*300          MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'MBSCCADJ' ! program name

C***********************************************************************
C   begin body of subroutine MBSCCADJ

C.........  Set up roadway type format
        IF( FIRSTIME ) THEN
            WRITE( RWTFMT, '("(I",I2.2,".",I2.2,")")' ) RWTLEN3, RWTLEN3
        END IF

C.........  Set road class from the TSCC
        RCL = STR2INT( TSCC( 8:10 ) )
    
C.........  Find road class in list of valid ones, to get index
        K   = FIND1( RCL, NRCLAS, AMSRDCLS )
        IF( RCL .NE. 0 .AND. K .LE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &             'ERROR: Road class ', RCL, 'at line', IREC, 
     &             'of surrogates cross-reference is' //
     &             CRLF() // BLANK10 // 'not in the list of '//
     &             'valid road class codes.'
            CALL M3MESG( MESG )
            RETURN

C.........  Use index to get valid roadway type and write to buffer
        ELSE IF( RCL .NE. 0 ) THEN

            WRITE( CRWT,RWTFMT ) RDWAYTYP( K )

        ELSE
            CRWT = REPEAT( '0', RWTLEN3 )

        END IF

C.........  Set vehicle type code from the TSCC
        CVID = TSCC( 3:6 )

C.........  Rearrange the SCC to reflect the desired heirarchy for
C           left-right SCC
        TSCCINTL = CRWT // CVID
        CALL PADZERO( TSCCINTL )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C.........  Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE MBSCCADJ
