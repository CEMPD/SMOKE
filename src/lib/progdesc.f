
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
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program SPCPMAT to take a SMOKE point source inventory file,',
     &  'an actual SCCs file, a speciation profiles file, a',
     &  'speciation cross-reference file, an optional pollutant-to-,',
     &  'pollutant conversion file, and produce a mass-based and/or',
     &  'mole-based SMOKE speciation matrices for all inventory',
     &  'pollutants using run-time defined combinations of pollutants',
     &  'and model species. The output species are defined at run time',
     &  'by the speciation profiles file, permitting support of any',
     &  'chemical mechanism.',
     &      ' '

        ELSEIF( NAME .EQ. 'SPCAMAT' ) THEN
        ELSEIF( NAME .EQ. 'SPCMMAT' ) THEN
        ELSEIF( NAME .EQ. 'GRDMAT' ) THEN
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program GRDMAT to take a SMOKE area, mobile, or point source',
     &  'inventory file, gridding surrogates,  surrogate cross-',
     &  'reference, and an optional link definitions file, and produce',
     &  'a SMOKE gridding matrix for a grid defined at run time. For',
     &  'mobile sources, an "ungridding" matrix is also created to',
     &  'allow the use of gridded temperature data in assigning',
     &  'factors to mobile sources.',
     &      ' '

        ELSEIF( NAME .EQ. 'TMPPOINT' ) THEN
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program TMPPOINT to take a SMOKE point source inventory file,',
     &  'an actual SCCs file, a temporal profiles file, a temporal ',
     &  'cross-reference file, an optional SMOKE day-specific file, ',
     &  'and an optional SMOKE hour-specific point source file, and',
     &  'produce hourly low-level and optionally hourly elevated point',
     &  'source emissions for the requested episode.',
     &      ' '

        ELSEIF( NAME .EQ. 'TMPAREA' ) THEN
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program TMPAREA to take a SMOKE area source inventory file,',
     &  'an actual SCCs file, a temporal profiles file, and a temporal',
     &  'cross-reference file, and produce hourly area source',
     &  'emissions for the requested episode.',
     &      ' '

        ELSEIF( NAME .EQ. 'TMPMOBIL' ) THEN
        ELSEIF( NAME .EQ. 'TMPBIO' ) THEN
        ELSEIF( NAME .EQ. 'CNTLMAT' ) THEN
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program CNTLMAT to take a SMOKE inventory file, a control',
     &  'packets file, an optional list of sources to track for the',
     &  'report, and produce one or more of the following: ',
     &  '     1) multiplicative control matrix ',
     &  '     2) additive control matrix ',
     &  '     3) reactivity control matrix ',
     &  '     4) projection matrix ',
     &  '     5) controls report ',
     &      ' '

        ELSEIF( NAME .EQ. 'ELEVPOINT' ) THEN
        ELSEIF( NAME .EQ. 'LAYPOINT' ) THEN
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program LAYPOINT to take a SMOKE point source inventory file',
     &  'and dot- and cross-point meteorology files and construct a',
     &  'point source layer fractions matrix for all selected hours',
     &  'and an optional report of plume exceeding a user-defined',
     &  'layer number. The program uses a Briggs method that has been',
     &  'adapted for multiple layers. ',
     &      ' '

        ELSEIF( NAME .EQ. 'EMISFAC' ) THEN
        ELSEIF( NAME .EQ. 'PREDIUR' ) THEN
        ELSEIF( NAME .EQ. 'SMKMERGE' ) THEN
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program SMKMERGE to merge the inventory or hourly emission',
     &  'files with gridding matrices and with optionally any',
     &  'combination of speciation matrices, multiplicative control',
     &  'matrices, additive control matrices, or reactivity control',
     &  'matrices. The program can operate on one to four source ',
     &  'categories (area, biogenic, mobile, or point sources), or any',
     &  'combination of these.  Gridded and/or state reports and/or',
     &  'county reports can be written from this program. If a layer-',
     &  'fractions file is input, then the total emissions output file',
     &  'is three-dimensional.'

        ELSEIF( NAME .EQ. 'GETRECS' ) THEN
            WRITE( LDEV,92000 ) 
     &  ' ',
     &  'Program GETRECS searches for a specific source, for all',
     &  'sources in a specific cell, or for combinations of source',
     &  'keys.  It creates an ASCII file which lists all details about',
     &  'the source including source number, grid cell, if found in',
     &  'gridding matrix, temporalization factors, control factors,',
     &  'inventory pollutant emissions, and model species emissions.'

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
