
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

C.........  SMOKE programs - listed in alphabetical order
        
        SELECT CASE( NAME )
        CASE( 'CNTLMAT' )       

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

        CASE( 'ELEVPOINT' ) 


        CASE( 'EMISFAC' )
            WRITE( *,92000 )
     &      ' ',
     &  'Program EMISFAC drives the MOBILE5a/b program by supplying',
     &  'a range of ambient temperatures and a scenario-specific',
     &  'MOBILE5 parameter file MPREF. Using multiple calls to',
     &  'MOBILE5, EMISFAC creates a diurnal and non-diurnal emission',
     &  'factors table for the specified temperatures and input',
     &  'parameter combinations.',
     &      ' '
 
        CASE( 'GETRECS' )
            WRITE( LDEV,92000 ) 
     &  ' ',
     &  'Program GETRECS searches for a specific source, for all',
     &  'sources in a specific cell, or for combinations of source',
     &  'keys.  It creates an ASCII file which lists all details about',
     &  'the source including source number, grid cell, if found in',
     &  'gridding matrix, temporalization factors, control factors,',
     &  'inventory pollutant emissions, and model species emissions.'

        CASE( 'GRDMAT' )
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

        CASE( 'LAYPOINT' ) 
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program LAYPOINT to take a SMOKE point source inventory file',
     &  'and dot- and cross-point meteorology files and construct a',
     &  'point source layer fractions matrix for all selected hours',
     &  'and an optional report of plume exceeding a user-defined',
     &  'layer number. The program uses a Briggs method that has been',
     &  'adapted for multiple layers. ',
     &      ' '

        CASE( 'METSCAN' )
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program METSCAN scans meteorology temperature data for any',
     &  'time range to determine the dates between which freezing',
     &  'occurs (sometimes called the first and last freeze dates.',
     &      ' '

        CASE( 'MRGGRID' )
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program MRGGRID reads 2-D area, biogenic, mobile, and 3-D',
     &  'point source emissions and merges into a single 3-D file.',
     &  'The time period merged is adjusted based on the latest',
     &  'starting file and earliest ending file.  All variables are',
     &  'merged, even if different variables are in each file.',
     &      ' '

        CASE( 'MVCONDNS' ) 
            WRITE( LDEV,92000 ) 
     &  ' ',
     &  'Program MVCONDNS condenses a mobile source inventory file to',
     &  'preprocess for preparing the MPLIST and MPREF files.'

        CASE( 'PREMOBL' )
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program PREMOBL to input gridded, time-dependent temperature',
     &  'data, an emission factors cross-reference, and an ungridding',
     &  'matrix. Program determines the minimum and maximum ',
     &  'temperatures per day for each source and for each emission',
     &  'factor.',
     &      ' '

        CASE( 'RAWBIO' ) 
            WRITE( LDEV,92000 )
     &      ' ',
     &  'Program RAWBIO to take the county level biomass,',
     &  'the emissions factors, and the surrogate factors,',
     &  'and produce gridded normalized biogenic emissions.',
     &  ' '

        CASE( 'SMK2EMIS' )
            WRITE( LDEV,92000 )
     &  ' ',
     &  'Program SMK2EMIS converts a NetCDF gridded, hourly emissions ',
     &  'file into a UAM ready gridded emissions file.  Program has ', 
     &  'been tested for lat-lon and UTM projections only.'

        CASE( 'SMKINVEN' ) 
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program SMKINVEN to take ASCII area or point source files',
     &  'in IDA, EPS2, EMS-95, or SMOKE list format, or mobile files',
     &  'in IDA format, and produce the I/O API and ASCII SMOKE',  
     &  'inventory files and list of unique SCCs in the inventory.',
     &      ' '

        CASE( 'SMKMERGE' )
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
     &  'is three-dimensional.',
     &      ' '

        CASE( 'SMKREPORT' )
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program SMKREPORT to take any combination of SMOKE intermed-',
     &  'iate files (but at least an inventory file) and generate',
     &  'one or more user-defined reports in any combination of output',
     &  'files.',
     &      ' '

        CASE( 'SPCMAT' )
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program SPCMAT to take a SMOKE area, mobile, or point source',
     &  'inventory file, a speciation profiles file, a speciation',
     &  'cross-reference file, an optional pollutant-to-pollutant,',
     &  'conversion file, and produce mass-based and/or mole-based',
     &  'SMOKE speciation matrices for all inventory pollutants',
     &  'using run-time defined combinations of pollutants and model',
     &  'species. The output species are defined at run time by the',
     &  'speciation profiles file, permitting support of any chemical',
     &  'mechanism.',
     &      ' '

        CASE( 'TEMPORAL' ) 
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program TEMPORAL to take a SMOKE area, mobile, or point',
     &  'source inventory file, a temporal profiles file, a temporal',
     &  'cross-reference file, an optional SMOKE day-specific file, ',
     &  'and an optional SMOKE hour-specific point source file, and',
     &  'produce hourly low-level and optionally hourly elevated point',
     &  'source emissions for the requested episode.',
     &      ' '

        CASE( 'TMPBIO' ) 
            WRITE( LDEV,92000 )
     &      ' ',
     &  'Program TMPBIO takes postprocessed MM5 meteorology and ',
     &  'normalized gridded emissions from RAWBIO or GRDBIO, and ', 
     &  'produces time stepped gridded speciated biogenic ',
     &  'emissions.',
     &      ' '

        CASE( 'UAM2NCF' ) 
            WRITE( LDEV,92000 )
     &  ' ',
     &  'Program UAM2NCF converts a UAM gridded, hourly emissions ',
     &  'into a NetCDF gridded, hourly emissions file.  Program has ',
     &  'been tested for the lat-lon projection only.'

        CASE DEFAULT
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'No program description is available for ' // NAME

        END SELECT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

        END
