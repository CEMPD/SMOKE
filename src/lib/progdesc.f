
        SUBROUTINE PROGDESC( LDEV, INPROGNM )

***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C       This program is used to write out a description of a program
C
C***********************************************************************
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
C****************************************************************************

        IMPLICIT NONE

C.........  Subroutine arguments
        INTEGER       LDEV      ! Log file unit number
        CHARACTER*(*) INPROGNM  ! Name of program to write description for

C.........  Local variables
        CHARACTER*16    NAME      ! Uppercase name of program to write descript

        CHARACTER*16 :: PROGNAME = 'PROGDESC' 

C***********************************************************************
C   begin body of subroutine PROGDESC

        NAME = INPROGNM
        CALL UPCASE( NAME )

        IF( NAME .EQ. 'RAWPOINT' ) THEN
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program RAWPOINT to take ASCII point source',
     &  'files in IDA, EPS2, EMS-95, or SMOKE list format',
     &  'and produce the SMOKE point source inventory file.',
     &      ' '

        ELSEIF( NAME .EQ. 'RAWAREA' ) THEN
        ELSEIF( NAME .EQ. 'RAWMOBIL' ) THEN
        ELSEIF( NAME .EQ. 'RAWBIO' ) THEN
        ELSEIF( NAME .EQ. 'GRDBIO' ) THEN
        ELSEIF( NAME .EQ. 'SPCPMAT' ) THEN
        ELSEIF( NAME .EQ. 'SPCAMAT' ) THEN
        ELSEIF( NAME .EQ. 'SPCMMAT' ) THEN
        ELSEIF( NAME .EQ. 'GRDPMAT' ) THEN
        ELSEIF( NAME .EQ. 'GRDAMAT' ) THEN
        ELSEIF( NAME .EQ. 'GRDMMAT' ) THEN
        ELSEIF( NAME .EQ. 'TMPPOINT' ) THEN
        ELSEIF( NAME .EQ. 'TMPAREA' ) THEN
        ELSEIF( NAME .EQ. 'TMPMOBIL' ) THEN
        ELSEIF( NAME .EQ. 'TMPBIO' ) THEN
        ELSEIF( NAME .EQ. 'CTLPMAT' ) THEN
        ELSEIF( NAME .EQ. 'CTLAMAT' ) THEN
        ELSEIF( NAME .EQ. 'CTLMMAT' ) THEN
        ELSEIF( NAME .EQ. 'ELEVPOINT' ) THEN
        ELSEIF( NAME .EQ. 'LAYPOINT' ) THEN
        ELSEIF( NAME .EQ. 'EMISFAC' ) THEN
        ELSEIF( NAME .EQ. 'PREDIUR' ) THEN
        ELSEIF( NAME .EQ. 'GRDPOINT' ) THEN
        ELSEIF( NAME .EQ. 'GRDAREA' ) THEN
        ELSEIF( NAME .EQ. 'GRDMOBIL' ) THEN
        ELSEIF( NAME .EQ. 'CSGPOINT' ) THEN
        ELSEIF( NAME .EQ. 'CSGAREA' ) THEN
        ELSEIF( NAME .EQ. 'CSGMOBIL' ) THEN
        ELSE
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'No program description is available for ' // NAME

        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

        END
