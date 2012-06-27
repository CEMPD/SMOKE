
        SUBROUTINE BLDMRGIDX

C***********************************************************************
C  subroutine BLDMRGIDX body starts at line
C
C  DESCRIPTION:
C      The purpose of this subroutine is to allocate and populate indicator
C      arrays that say which pollutants and species are present for each
C      type of input file (inventory, speciation matrix, multiplicative
C      control matrix, reactivity matrix, etc.) for each source category
C      (area, biogenic, mobile, point).  The indicator arrays store the
C      position in the list of i/o api file variables that the given pollutant
C      or species matches.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
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

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: NIPPA,
     &                      NSMATV, NMSPC,
     &                      TSVDESC,
     &                      SIINDEX, SPINDEX,
     &                      EMNAM, EANAM

C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: EANAMREP

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  EXTERNAL FUNCTIONS
        CHARACTER(2)    CRLF   
        INTEGER         INDEX1
        EXTERNAL        CRLF, INDEX1

C...........   Call allocated arrays
C...........   Group index counter for each source-category-specific list of 
C              pollutants and activities.
        INTEGER  KM( NIPPA )    !  mobile

C...........   Allocatable arrays
        INTEGER, ALLOCATABLE :: EANAMIDX( : )  ! index from EANAM to TSVDESC

C...........   Other local variables
        INTEGER         J, K, L1, L2, V    !  counters and indices

        INTEGER         IOS      ! i/o error status
        INTEGER         MCNT     ! mobile src var counter
        INTEGER         TGRP     ! tmp group number

        LOGICAL      :: EFLAG = .FALSE.  ! true: error found
        LOGICAL         NEXTGRP  ! true: time to increment group number cntr

        CHARACTER(300)     :: MESG ! message buffer
        CHARACTER(IODLEN3) :: CBUF ! tmp pol-to-species buffer
        CHARACTER(IOVLEN3) :: CPOL ! tmp pollutant buffer 
        CHARACTER(IOVLEN3) :: CSPC ! tmp species buffer 
        CHARACTER(IOVLEN3) :: PSPC ! tmp previous species  
        CHARACTER(IOVLEN3) :: VBUF ! tmp variable name buffer 
        CHARACTER(PLSLEN3) :: SVBUF ! tmp speciation name buffer 

        CHARACTER(16) :: PROGNAME = 'BLDMRGIDX' ! program name

C***********************************************************************
C   begin body of subroutine BLDMRGIDX

C.........  Allocate memory for building reporting flag array
        ALLOCATE( EANAMIDX( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAMIDX', PROGNAME )
        ALLOCATE( EANAMREP( NSMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAMREP', PROGNAME )
        EANAMIDX = 0
        EANAMREP = .FALSE.

C.........  Allocate memory for variable to pollutant and variable to species
C           indexes.
        ALLOCATE( SIINDEX( NSMATV, 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SIINDEX', PROGNAME )
        ALLOCATE( SPINDEX( NSMATV, 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPINDEX', PROGNAME )
        SIINDEX = 0  ! array
        SPINDEX = 0  ! array

C.........  Loop through pol-to-species names
        PSPC = ' '
        DO V = 1, NSMATV

C.............  Extract pollutant name and species name
            CBUF = TSVDESC( V )
            L1 = INDEX   ( CBUF, SPJOIN )
            L2 = LEN_TRIM( CBUF )
            CPOL = CBUF(    1:L1-1 )
            CSPC = CBUF( L1+1:L2   )

C.............  Determine location in the list of pollutants/activities
            K = INDEX1( CPOL, NIPPA, EANAM )
            J = INDEX1( CSPC, NMSPC, EMNAM )

C.............  Store pol/act and species indices
            SIINDEX( V, 1 ) = K             ! store pol/act index
            SPINDEX( V, 1 ) = J             ! store species index

C.............  Save index into TSVDESC array from EANAM
            IF( EANAMIDX( K ) .EQ. 0 ) THEN
                EANAMIDX( K ) = V
            END IF
            PSPC  = CSPC

        END DO

C.........  Build array to indicate when report values should be saved
        DO V = 1, NIPPA
            IF ( EANAMIDX( V ) > 0 ) THEN
                EANAMREP( EANAMIDX( V ) ) = .TRUE.
            END IF
        END DO
        
        DEALLOCATE( EANAMIDX )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE BLDMRGIDX
