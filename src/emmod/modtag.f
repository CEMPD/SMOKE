        MODULE MODTAG

!***********************************************************************
!  Module body starts at line 40
!
!  DESCRIPTION:
!     This module contains the public allocatable arrays for tagging 
!     tables
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 2/2009 by M. Houyoux
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! Pathname: $Source$
! Last updated: $Date$ 
!
!****************************************************************************

        IMPLICIT NONE

        INCLUDE 'EMPRVT3.EXT'   !  private emissions string widths parameters

        INTEGER, PUBLIC :: MXTAG = 0  !  maximum number of tags for any species

!.........  Tag numbers by species,pollutant
        INTEGER, ALLOCATABLE, PUBLIC :: TAGNUM( :,: ) ! number of tags for each species/pol combo

!.........  Tag names by tag number, species, pollutant
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGNAME( :,:,: )  ! tag labels for all tags for each species/pol combo

!.........  Cross-reference table sizes for tagging tables
        INTEGER, ALLOCATABLE, PUBLIC :: TAGXCNT( : ) ! table size for each table (number of entries of each type)

!.........  List of species that are tagged
        INTEGER, PUBLIC :: NTAGSALL   ! count of all species that are getting tags.
        CHARACTER(IOVLEN3), ALLOCATABLE, PUBLIC :: TAGSPECIES( : )  ! list of species that get tags

!.........  Sorted groups of cross-references
!.........  For speciation, TAGT* values are the tagging labels

!.........  FIPS code=0, SCC=all (level 4)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT03( :,: )
        CHARACTER(SCCLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT03( : )

!.........  FIPS code=state default, SCC=0  
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT04( :,: )
        CHARACTER(STALEN3), ALLOCATABLE, PUBLIC :: TAGCHRT04( : )

!.........  FIPS code=state, SCC=all (level 4)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT06( :,: )
        CHARACTER(STSLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT06( : )

!.........  FIPS code=all, SCC=0
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT07( :,: )
        CHARACTER(FIPLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT07( : )

!.........  FIPS code=all, SCC=all (level 4)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT09( :,: )
        CHARACTER(FPSLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT09( : )

!.........  PLANT=non-blank, SCC=0
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT10( :,: )
        CHARACTER(FPLLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT10( : )

!.........  Plant=non-blank, SCC=all
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT11( :,: )
        CHARACTER(SS0LEN3), ALLOCATABLE, PUBLIC :: TAGCHRT11( : )

!.........  Additional groups for special SIC handling. Note that these are
!              put at the end since they were added later and added for
!              Cntlmat only at first.  Did not want to mess up the numbering
!              scheme for TXCNT
!.........  FIPS code = 0, SIC = 2-digit (Type 26)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT26( :,: )
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT26( : )

!.........  FIPS code = 0, SIC = all (Type 27)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT27( :,: )
        CHARACTER(SICLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT27( : )

!.........  FIPS code = state, SIC = 2-digit (Type 28)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT28( :,: )
        CHARACTER(STILEN3), ALLOCATABLE, PUBLIC :: TAGCHRT28( : )

!.........  FIPS code = state, SIC = all (Type 29)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT29( :,: )
        CHARACTER(STILEN3), ALLOCATABLE, PUBLIC :: TAGCHRT29( : )

!.........  FIPS code = all, SIC = 2-digit (Type 30)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT30( :,: )
        CHARACTER(FPILEN3), ALLOCATABLE, PUBLIC :: TAGCHRT30( : )

!.........  FIPS code = all, SIC = all (Type 31)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT31( :,: )
        CHARACTER(FPILEN3), ALLOCATABLE, PUBLIC :: TAGCHRT31( : )

!.........  Additional groups for MACT processing

!.........  FIPS code = 0, SCC = 0, MACT = all (Type 32)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT32( :,: )
        CHARACTER(MACLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT32( : )
        
!.........  FIPS code = 0, SCC = all, MACT = all (Type 33)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT33( :,: )
        CHARACTER(MSCLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT33( : )
        
!.........  FIPS code = state, SCC = 0, MACT = all (Type 34)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT34( :,: )
        CHARACTER(MSTLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT34( : )
        
!.........  FIPS code = state, SCC = all, MACT = all (Type 35)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT35( :,: )
        CHARACTER(MSSLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT35( : )
        
!.........  FIPS code = all, SCC = 0, MACT = all (Type 36)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT36( :,: )
        CHARACTER(MFPLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT36( : )
        
!.........  FIPS code = all, SCC = all, MACT = all (Type 37)
        CHARACTER(TAGLEN3), ALLOCATABLE, PUBLIC :: TAGT37( :,: )
        CHARACTER(MFSLEN3), ALLOCATABLE, PUBLIC :: TAGCHRT37( : )

        END MODULE MODTAG
