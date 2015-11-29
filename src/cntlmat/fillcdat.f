
        SUBROUTINE FILLCDAT( PKTTYP, JT, PKTINFO )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine populates the control data table depending on the 
C      packet type being processed when the routine is called (PKTTYP)
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Started 3/99 by M. Houyoux
C
C***************************************************************************
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: FACCTG, CUTCTG, FACMACT, FACRACT, 
     &                      ICTLEQUP, CCTLSIC, FACCEFF, FACREFF, 
     &                      FACRLPN, CALWSIC, FACALW, EMCAPALW,
     &                      EMREPALW, EMREPREA, PRJFCREA, MKTPNREA,
     &                      CSCCREA, CSPFREA, CPRJSIC, PRJFC, CTLRPLC,
     &                      MACEXEFF, MACNWEFF, MACNWFRC, CMACSRCTYP, 
     &                      CTGCOMT, CTLCOMT, ALWCOMT,
     &                      REACOMT, PRJCOMT, MACCOMT

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'CPKTDAT.EXT'   !  control packet contents        

C...........   EXTERNAL FUNCTIONS:
        INTEGER       STR2INT
        EXTERNAL      STR2INT

C...........   SUBROUTINE ARGUMENTS:

        CHARACTER(*), INTENT (IN) :: PKTTYP    ! packet type 
        INTEGER     , INTENT (IN) :: JT        ! index to control data tables
        TYPE(CPACKET),INTENT (IN) :: PKTINFO   ! packet information

        CHARACTER(16) :: PROGNAME = 'FILLCDAT' ! program name

C***********************************************************************
C   Begin body of subroutine FILLCDAT

        SELECT CASE( PKTTYP )

        CASE( 'CTG' )
            FACCTG ( JT ) = PKTINFO%FAC1
            CUTCTG ( JT ) = PKTINFO%FAC2
            FACMACT( JT ) = PKTINFO%FAC3
            FACRACT( JT ) = PKTINFO%FAC4
            CTGCOMT( JT ) = PKTINFO%COMMENT

        CASE( 'CONTROL' )
            ICTLEQUP( JT ) = INT    ( PKTINFO%FAC1 )
            CCTLSIC ( JT ) = PKTINFO%CSIC
            FACCEFF ( JT ) = PKTINFO%FAC2
            FACREFF ( JT ) = PKTINFO%FAC3
            FACRLPN ( JT ) = PKTINFO%FAC4
            IF( PKTINFO%REPFLAG == 'R' .OR.
     &          PKTINFO%REPFLAG == 'Y'      ) THEN
                CTLRPLC( JT ) = .TRUE.
            END IF
            CTLCOMT( JT ) = PKTINFO%COMMENT

        CASE( 'ALLOWABLE' )
            CALWSIC ( JT ) = PKTINFO%CSIC
            FACALW  ( JT ) = PKTINFO%FAC1
            EMCAPALW( JT ) = PKTINFO%FAC2
            EMREPALW( JT ) = PKTINFO%FAC3
            ALWCOMT( JT ) = PKTINFO%COMMENT

        CASE( 'REACTIVITY' )
            EMREPREA( JT ) = PKTINFO%FAC1
            PRJFCREA( JT ) = PKTINFO%FAC2
            MKTPNREA( JT ) = PKTINFO%FAC3
            CSCCREA ( JT ) = PKTINFO%NSCC
            CSPFREA ( JT ) = PKTINFO%TMPPRF
            REACOMT( JT ) = PKTINFO%COMMENT

        CASE( 'PROJECTION' )
            CPRJSIC( JT ) = PKTINFO%CSIC
            PRJFC  ( JT ) = PKTINFO%FAC1
            PRJCOMT( JT ) = PKTINFO%COMMENT

        CASE( 'MACT' )
            MACEXEFF( JT ) = PKTINFO%FAC1
            MACNWEFF( JT ) = PKTINFO%FAC2
            MACNWFRC( JT ) = PKTINFO%FAC3
            
            CMACSRCTYP( JT ) = PKTINFO%CSTYP
            CALL PADZERO( CMACSRCTYP( JT ) )
            MACCOMT( JT ) = PKTINFO%COMMENT
            
C.............  Make sure src type is only 00, 01, or 02            
            IF( CMACSRCTYP( JT ) /= '01' .AND.
     &          CMACSRCTYP( JT ) /= '02'       ) THEN
                CMACSRCTYP( JT ) = '00'
            END IF

        END SELECT

        RETURN

        END SUBROUTINE FILLCDAT
