
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
     &                      ICTLEQUP, ICTLSIC, FACCEFF, FACREFF, 
     &                      FACRLPN, IALWSIC, FACALW, EMCAPALW,
     &                      EMREPALW, EMREPREA, PRJFCREA, MKTPNREA,
     &                      CSCCREA, CSPFREA, IPRJSIC, PRJFC, IEMSSIC,
     &                      BASCEFF, BASREFF, BASRLPN, EMSCEFF, EMSREFF,
     &                      EMSRLPN, EMSPTCF, EMSTOTL, CTLRPLC,
     &                      MACEXEFF, MACNWEFF, MACNWFRC, CMACSRCTYP

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

	CHARACTER(300)   MESG                  ! message buffer
        CHARACTER(16) :: PROGNAME = 'FILLCDAT' ! program name

C***********************************************************************
C   Begin body of subroutine FILLCDAT

        SELECT CASE( PKTTYP )

        CASE( 'CTG' )
            FACCTG ( JT ) = PKTINFO%FAC1
            CUTCTG ( JT ) = PKTINFO%FAC2
            FACMACT( JT ) = PKTINFO%FAC3
            FACRACT( JT ) = PKTINFO%FAC4

        CASE( 'CONTROL' )
            ICTLEQUP( JT ) = INT    ( PKTINFO%FAC1 )
            ICTLSIC ( JT ) = STR2INT( PKTINFO%CSIC )
            FACCEFF ( JT ) = PKTINFO%FAC2
            FACREFF ( JT ) = PKTINFO%FAC3
            FACRLPN ( JT ) = PKTINFO%FAC4
            IF( PKTINFO%REPFLAG == 'R' .OR.
     &          PKTINFO%REPFLAG == 'Y'      ) THEN
                CTLRPLC = .TRUE.
            END IF

        CASE( 'ALLOWABLE' )
            IALWSIC ( JT ) = STR2INT( PKTINFO%CSIC )
            FACALW  ( JT ) = PKTINFO%FAC1
            EMCAPALW( JT ) = PKTINFO%FAC2
            EMREPALW( JT ) = PKTINFO%FAC3

        CASE( 'REACTIVITY' )
            EMREPREA( JT ) = PKTINFO%FAC1
            PRJFCREA( JT ) = PKTINFO%FAC2
            MKTPNREA( JT ) = PKTINFO%FAC3
            CSCCREA ( JT ) = PKTINFO%NSCC
            CSPFREA ( JT ) = PKTINFO%TMPPRF

        CASE( 'PROJECTION' )
            IPRJSIC( JT ) = STR2INT( PKTINFO%CSIC )
            PRJFC  ( JT ) = PKTINFO%FAC1

        CASE( 'EMS_CONTROL' )
            IEMSSIC ( JT ) = STR2INT( PKTINFO%CSIC )
            BASCEFF ( JT ) = PKTINFO%FAC1
            BASREFF ( JT ) = PKTINFO%FAC2
            BASRLPN ( JT ) = PKTINFO%FAC3
            EMSCEFF ( JT ) = PKTINFO%FAC4
            EMSREFF ( JT ) = PKTINFO%FAC5
            EMSRLPN ( JT ) = PKTINFO%FAC6
            IF( PKTINFO%FAC7 .GT. 0. ) 
     &          EMSPTCF ( JT ) = PKTINFO%FAC7
            IF( PKTINFO%FAC8 .GT. 0. ) 
     &          EMSTOTL ( JT ) = PKTINFO%FAC8

        CASE( 'MACT' )
            MACEXEFF( JT ) = PKTINFO%FAC1
            MACNWEFF( JT ) = PKTINFO%FAC2
            MACNWFRC( JT ) = PKTINFO%FAC3
            
            CMACSRCTYP( JT ) = PKTINFO%CSTYP
            CALL PADZERO( CMACSRCTYP( JT ) )
            
C.............  Make sure src type is only 00, 01, 02, 03 or 04            
            IF( CMACSRCTYP( JT ) /= '01' .AND.
     &          CMACSRCTYP( JT ) /= '02' .AND.
     &          CMACSRCTYP( JT ) /= '03' .AND. 
     &          CMACSRCTYP( JT ) /= '04'      ) THEN
                CMACSRCTYP( JT ) = '00'
		WRITE( MESG,94010 ) 'WARNING: Source type '//
     &            'code at entry ', JT, ' is invalid. Source '//
     &            'type code will be set to 00.'
		CALL M3MESG( MESG )
            END IF

        END SELECT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

94010	FORMAT( 10( A, :, I4, :, 1X ) )

        END SUBROUTINE FILLCDAT
