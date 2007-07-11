
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
C****************************************************************************

        IMPLICIT NONE

C.........  Subroutine arguments
        INTEGER      LDEV      ! Log file unit number
        CHARACTER(*) INPROGNM  ! Name of program to write description for

C.........  Local variables
        CHARACTER(16)   NAME      ! Uppercase name of program to write descript

        CHARACTER(16) :: PROGNAME = 'PROGDESC' 

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
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program ELEVPOINT specifies elevated point sources in the',
     &  'inventory by doing one of the following:',
     &  '     1) Reads a list of preset "major" point sources',
     &  '     2) Uses a "cutoff" height and an analytical plume',
     &  '        rise equation to determine elevated sources for UAM-',
     &  '        style point source processing.',
     &  'It also writes one of the plume-in-grid source files for',
     &  'the CMAQ model.',
     &      ' '

        CASE( 'EMISFAC' )
            WRITE( *,92000 )
     &      ' ',
     &  'Program EMISFAC drives the MOBILE6 model using custom',
     &  'input files created from county-based hourly temperature',
     &  'profiles and MOBILE6 input scenarios. Separate runs of',
     &  'EMISFAC are needed for each temperature averaging type.',
     &  'After running MOBILE6, EMISFAC stores the source-based',
     &  'emission factors.',
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

        CASE( 'GRWINVEN' )
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program GRWINVEN reads a SMOKE area, mobile, or point source',
     &  'inventory file and one or more control and/or projection',
     &  'files, applies the controls and projections to the inventory,',
     &  'and writes out a new SMOKE inventory file and/or an IDA-',
     &  'formatted inventory file.',
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

        CASE( 'MBSETUP' )
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program MBSETUP prepares intermediate files needed to run',
     &  'PREMOBL and EMISFAC using the county cross-reference and',
     &  'settings files. It checks that each county in the inventory',
     &  'and within the grid has been assigned a reference county',
     &  'and creates the speed summary file grouping sources by county',
     &  'and speed.',
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
     &      ' ',
     &  'Program MVCONDNS condenses a mobile source inventory file to',
     &  'preprocess for preparing the MPLIST and MPREF files.',
     &      ' '

        CASE( 'MVSETUP' ) 
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program MVSETUP uses the SMOKE mobile source file to create',
     &  'a condensed list of all sources in the inventory by FIPS',
     &  'code, roadtype, vehicle type, and including the speed',
     &  'from the inventory (if any).',
     &      ' '

        CASE( 'PREMOBL' )
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program PREMOBL uses gridded, time-dependent temperature',
     &  'data and an ungridding matrix to create county-based',
     &  '24-hour temperature profiles.',
     &      ' '

        CASE( 'RAWBIO' ) 
            WRITE( LDEV,92000 )
     &      ' ',
     &  'Program RAWBIO uses the biogenic emission factors, and',
     &  'either:',
     &  '   1) county-level biogenic land use data and gridding',
     &  '      surrogate factors, or',
     &  '   2) gridded biogenic land use data',
     &  'to produce gridded normalized biogenic emissions.',
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
     &  'in IDA, EMS-95, or SMOKE list format, or mobile files',
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

         CASE( 'SURGTOOL' )       
            WRITE( LDEV,92000 ) 
     &      ' ',
     &  'Program CROSSGRID to take a SMOKE gridding surrogate file ',
     &  'for a "fine" grid, a grid definition for a "coarse" grid, ',
     &  'and produce the approximate "coarse" gridding surrogate file.',
     &  ' ',
     &  'NOTES:',
     &  '   (1)  Current version is for Lat-Lon, Lambert, and UTM',
     &  '        projections only. Can perform Lambert-to-Lambert',
     &  '        and UTM zone-to-zone transformations)',
     &  ' ',
     &  '   (4)  Inputs and outputs only SMOKE formatted files.',
     &  ' '

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
     &  'emissions. It optionally takes is a gridded winter-summer',
     &  'designation to apply winter and summer emission factors by',
     &  'grid cell.',
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
