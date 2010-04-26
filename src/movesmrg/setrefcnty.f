
        SUBROUTINE SETREFCNTY

C***********************************************************************
C  subroutine SETREFCNTY body starts at line
C
C  DESCRIPTION:
C      This subroutine sets up reference county information for Movesmrg.
C      It reads the county cross-reference file, fuel month file, and
C      reference county emission factors list.
C
C  PRECONDITIONS REQUIRED:  
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 4/10 by C. Seppanen
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
C.........  This module contains data structures and flags specific to Movesmrg
        USE MODMVSMRG, ONLY: XDEV, MDEV, FDEV, NSRCCELLS, 
     &                       NREFSRCS, REFSRCS,
     &                       RPPFLAG, MMDEV

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC

C...........   This module is the source inventory arrays
        USE MODSOURC, ONLY: IFIP

C.........  This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: NINVIFIP, INVIFIP

C.........  This module is used for reference county information
        USE MODMBSET, ONLY: NINVC, NREFC, MCREFSORT, MCREFIDX

        IMPLICIT NONE

C.........  INCLUDES:

C.........  EXTERNAL FUNCTIONS and their descriptions:
        INTEGER   INDEXINT1
        
        EXTERNAL  INDEXINT1

C.........  LOCAL VARIABLES and their descriptions:

C.........  Array to store counties inside grid
        INTEGER, ALLOCATABLE :: GRDFIP( : )

C.........  Arrays for building list of reference county sources
        INTEGER, ALLOCATABLE :: INVCNTY( : )    ! inventory counties sorted by reference county
        INTEGER, ALLOCATABLE :: REFCNTY( : )    ! list of sorted reference counties
        INTEGER, ALLOCATABLE :: SRCREFIDX( : )  ! ref. county index for each source

C.........  Other local variables
        INTEGER   I, J, K, S  ! indexes and counters
        INTEGER   MXNREFSRCS  ! max. no. sources per reference county
        INTEGER   NGRDFIP     ! no. counties inside grid
        INTEGER   PFIP        ! tmp previous iteration fip
        INTEGER   REFFIP      ! current reference county fip
        INTEGER   REFIDX      ! current reference county index
        INTEGER   IOS         ! error status

        CHARACTER(300)     MESG    ! message buffer

        CHARACTER(16) :: PROGNAME = 'SETREFCNTY' ! program name

C***********************************************************************
C   begin body of subroutine SETREFCNTY

C.........  Build list of counties within the grid
        ALLOCATE( GRDFIP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRDFIP', PROGNAME )
        GRDFIP = 0   ! array
        
        PFIP = -9
        NGRDFIP = 0
        DO S = 1, NSRC

            IF( NSRCCELLS( S ) == 0 ) CYCLE
            
            IF( IFIP( S ) .NE. PFIP ) THEN
                NGRDFIP = NGRDFIP + 1
                GRDFIP( NGRDFIP ) = IFIP( S )
                PFIP = IFIP( S )
            END IF
        
        END DO

C.........  Read county cross-reference file
        CALL RDMXREF( XDEV, NGRDFIP, GRDFIP )
        
        DEALLOCATE( GRDFIP )

C.........  Build list of inventory counties sorted by reference county
        ALLOCATE( INVCNTY( NINVC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INVCNTY', PROGNAME )
        
        DO I = 1, NINVC
            INVCNTY( I ) = MCREFSORT( I,1 )
        END DO

        ALLOCATE( REFCNTY( NREFC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REFCNTY', PROGNAME )
        
        DO I = 1, NREFC
            REFCNTY( I ) = MCREFIDX( I,1 )
        END DO

C.........  Build list of sources for each reference county
C           Start by counting number of sources for each reference county
        ALLOCATE( NREFSRCS( NREFC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NREFSRCS', PROGNAME )
        ALLOCATE( SRCREFIDX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRCREFIDX', PROGNAME )
        NREFSRCS = 0      ! array
        SRCREFIDX = 0   ! array
        
        PFIP = -9
        REFFIP = -9
        REFIDX = 0
        DO S = 1, NSRC
        
            IF( IFIP( S ) .NE. PFIP ) THEN
                J = INDEXINT1( IFIP( S ), NINVC, INVCNTY )
                IF( J <= 0 ) THEN
                    WRITE( MESG, 94010 ) 'ERROR: Inventory county ',
     &                IFIP( S ), 'not listed in county cross-reference '
     &                // 'file'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                IF( MCREFSORT( J,2 ) .NE. REFFIP ) THEN
                    REFFIP = MCREFSORT( J,2 )
                    
                    REFIDX = INDEXINT1( REFFIP, NREFC, REFCNTY )
                    IF( REFIDX <= 0 ) THEN
                        WRITE( MESG, 94010 ) 'INTERNAL ERROR: ' //
     &                    'Problem with reference county mapping for '
     &                    // 'county ', IFIP( S )
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                END IF
                
                PFIP = IFIP( S )
            END IF
            
            NREFSRCS( REFIDX ) = NREFSRCS( REFIDX ) + 1
            
            SRCREFIDX( S ) = REFIDX

        END DO

C.........  Get maximum number of sources per reference county
        MXNREFSRCS = 0
        DO I = 1, NREFC
            
            IF( NREFSRCS( I ) > MXNREFSRCS ) MXNREFSRCS = NREFSRCS( I )
            
        END DO
        
        ALLOCATE( REFSRCS( NREFC, MXNREFSRCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REFSRCS', PROGNAME )
        REFSRCS = 0   ! array
        NREFSRCS = 0  ! array
        
        REFIDX = 0
        DO S = 1, NSRC
        
            REFIDX = SRCREFIDX( S )
            K = NREFSRCS( REFIDX ) + 1
            REFSRCS( REFIDX, K ) = S
            NREFSRCS( REFIDX ) = K
        
        END DO
        
        DEALLOCATE( INVCNTY, REFCNTY, SRCREFIDX )

C.........  Read fuel month reference file
        CALL RDFMREF( MDEV )

C.........  Read emission factors file list
        CALL RDMRCLIST( FDEV )

C.........  Read Met4moves output file
        IF( RPPFLAG ) CALL RDMETMOVES( MMDEV )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE SETREFCNTY
