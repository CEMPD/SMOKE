
        SUBROUTINE WRSPDSUM( MDEV, REFIDX, SPDFLAG )

C***********************************************************************
C  subroutine body starts at line 104
C
C  DESCRIPTION:
C       Writes the SPDSUM output file
C
C  PRECONDITIONS REQUIRED:
C       MDEV has been opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***********************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC, ONLY: IFIP, IRCLAS, SPEED, VMT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC

C.........  This module is used for MOBILE6 setup information         
        USE MODMBSET, ONLY: NINVC, NREFC, MCREFIDX, MCREFSORT, MVREFSORT
        
C.........  This module is for cross reference tables
        USE MODXREF, ONLY: SPDPROFID
        
        IMPLICIT NONE

C.........  INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'M6CNST3.EXT'   !  MOBILE6 constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER        FIND1FIRST
        INTEGER        FINDR1FIRST
        INTEGER        CVTRDTYPE
        CHARACTER*2    CRLF    
        
        EXTERNAL  FINDR1FIRST, FIND1FIRST, CVTRDTYPE, CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: MDEV     ! SPDSUM file unit no.
        INTEGER, INTENT (IN) :: REFIDX   ! index of reference county
        LOGICAL, INTENT (IN) :: SPDFLAG  ! true: use speed profiles
        
C...........   Local allocatable arrays
        REAL,    ALLOCATABLE :: COUNTYSPDS ( :,: )  ! sources, speeds, and road types per county
        REAL,    ALLOCATABLE :: SPDSORT    ( :,: )  ! COUNTYSPDS sorted by road type and speed
        
        INTEGER, ALLOCATABLE :: SRCIDX( :,:)        ! starting and ending indices into FIPS 
                                                    !    for each county in ref. county
        INTEGER, ALLOCATABLE :: IDX( : )            ! index to sort COUNTYSPDS

C...........   Other local variables
        INTEGER I, J, K                   ! counters and indices

        INTEGER IOS                       ! I/O status

        INTEGER REFCOUNTY                 ! ref. county FIPS code
        INTEGER INVCOUNTY                 ! inv. county FIPS code
        INTEGER COUNTY                    ! stores ref. or inv. county depending on SPATFLAG

        INTEGER STINVIDX                  ! starting index of inv. counties in MCREFSORT
        INTEGER ENDINVIDX                 ! ending index of inv. counties

        INTEGER STSRCIDX                  ! starting index of sources in IFIP
        INTEGER ENDSRCIDX                 ! ending index of sources

        INTEGER STARTPOS                  ! starting position for storing source info 

        INTEGER NCOUNTYREF                ! no. counties in ref. county
        INTEGER NSRCCOUNTY                ! no. sources in a county
        
        LOGICAL :: SPATFLAG = .FALSE.     ! true: group speeds by ref. county
        LOGICAL :: LASAFLAG = .FALSE.     ! true: treat local roads as arterial
        
        CHARACTER*300          MESG      !  message buffer 

        CHARACTER*16 :: PROGNAME = 'WRSPDSUM'   ! program name
        
C***********************************************************************
C   begin body of subroutine WRSPDSUM
                
        REFCOUNTY = MCREFIDX( REFIDX,1 )

        LASAFLAG = .FALSE.
        SPATFLAG = .FALSE.

C.........  Check if local roads should be treated as arterial
        IF( MVREFSORT( REFIDX,4 ) /= 1 ) THEN
            LASAFLAG = .TRUE.
        END IF

C.........  Check if spatial averaging is requested
        IF( MVREFSORT( REFIDX,2 ) /= 1 ) THEN
            SPATFLAG = .TRUE.
        END IF
        
C.........  Get starting and ending indices into MCREFSORT array
        STINVIDX = MCREFIDX( REFIDX,2 )
        
        IF( REFIDX == NREFC ) THEN
            ENDINVIDX = NINVC
        ELSE
            ENDINVIDX = MCREFIDX( REFIDX + 1,2 ) - 1
        END IF

C.........  If grouping, calculate no. counties in this ref. county and allocate
C           array to store source indices
        IF( SPATFLAG ) THEN       
            NCOUNTYREF = ENDINVIDX - STINVIDX + 1
            
            ALLOCATE( SRCIDX( NCOUNTYREF,2 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRCIDX', PROGNAME )
            
            SRCIDX = 0.
        ELSE
            NCOUNTYREF = 1
        END IF

        NSRCCOUNTY = 0
        
C.........  Loop through inventory counties 
        DO J = STINVIDX, ENDINVIDX
            
            INVCOUNTY = MCREFSORT( J,1 )
            
C.............  Find starting index for this county in IFIP array
            K = FIND1FIRST( INVCOUNTY, NSRC, IFIP )
                
            STSRCIDX = K
            
C.............  Find end index for this county                                
            DO                
                K = K + 1
                IF( K > NSRC .OR. IFIP( K ) /= INVCOUNTY ) EXIT
            END DO

            ENDSRCIDX = K - 1

C.............  If grouping, store start and end indices for current county            
            IF( SPATFLAG ) THEN
                SRCIDX( J - STINVIDX + 1,1 ) = STSRCIDX
                SRCIDX( J - STINVIDX + 1,2 ) = ENDSRCIDX
            ELSE
                NSRCCOUNTY = 0
            END IF

C.............  Calculate no. sources in this county (running total if grouping)            
            NSRCCOUNTY = NSRCCOUNTY + ( ENDSRCIDX - STSRCIDX + 1 )

C.............  If grouping by ref. county, skip rest of loop except on last time through
            IF( SPATFLAG .AND. J /= ENDINVIDX ) CYCLE
            
C.............  Allocate arrays for storing sources                
            ALLOCATE( COUNTYSPDS( NSRCCOUNTY,3 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'COUNTYSPDS', PROGNAME )
            ALLOCATE( SPDSORT( NSRCCOUNTY,3 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SPDSORT', PROGNAME )
            ALLOCATE( IDX( NSRCCOUNTY ), STAT=IOS )
            CALL CHECKMEM( IOS, 'IDX', PROGNAME )

C.............  Initialize arrays
            COUNTYSPDS = 0.        
            SPDSORT    = 0.
            
            STARTPOS = 1
            
C.............  If grouping, we need to store source info for all counties, so loop
C               through all counties.  If not grouping, we only need to store info for
C               the current county, so only loop once           
            DO I = 1, NCOUNTYREF
            
C.................  Set starting and ending indices for each county            
                IF( SPATFLAG ) THEN
                    STSRCIDX  = SRCIDX( I,1 )
                    ENDSRCIDX = SRCIDX( I,2 )
                END IF
                    
C.................  Store source information in speed array
                CALL STORE_SOURCE_INFO( COUNTYSPDS, STSRCIDX, ENDSRCIDX,
     &                                  STARTPOS, NSRCCOUNTY, LASAFLAG )                                 

C.................  Update starting position so source info is stored at correct 
C                   place in speed array
                STARTPOS = STARTPOS + ( ENDSRCIDX - STSRCIDX + 1 )

            END DO

C.............  Initialize index array for sorting
            DO K = 1, NSRCCOUNTY
                IDX( K ) = K
            END DO

C.............  Sort speeds array by road type, then speed, then source
            CALL SORTR3( NSRCCOUNTY, IDX, COUNTYSPDS( :,3 ), 
     &                   COUNTYSPDS( :,2 ), COUNTYSPDS( :,1 ) )

            DO K = 1, NSRCCOUNTY
                SPDSORT( K,: ) = COUNTYSPDS( IDX( K ),: )
            END DO                    

C.............  Set which county to output           
            IF( SPATFLAG ) THEN
                COUNTY = REFCOUNTY
            ELSE
                COUNTY = INVCOUNTY
            END IF
            
C.............  Write information to SPDSUM file
            WRITE( MESG,94010 ) 'Writing to speed summary file ' //
     &             'for county', COUNTY, '...'
            CALL M3MESG( MESG )
            
            CALL WRITE_ARRAY_TO_SPDSUM( MDEV, SPDSORT, NSRCCOUNTY, 
     &                                  COUNTY )
            
            DEALLOCATE( COUNTYSPDS, SPDSORT, IDX )
            
            IF( SPATFLAG ) DEALLOCATE( SRCIDX )        
        
        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
      
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
94020   FORMAT( 3( A, 1X ), I8, 1X, A, 1X )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

            SUBROUTINE STORE_SOURCE_INFO( SPDARRAY, STIDX, ENDIDX, 
     &                                    STPOS, NUMSRC, LASAFLAG )

C.............  Subprogram arguments
            REAL,    INTENT (OUT) :: SPDARRAY( NUMSRC,3 ) ! unsorted array of spds, roads, sources
            INTEGER, INTENT (IN)  :: STIDX                ! starting index into FIPS array
            INTEGER, INTENT (IN)  :: ENDIDX               ! ending index into FIPS array
            INTEGER, INTENT (IN)  :: STPOS                ! starting position in SPDARRAY
            INTEGER, INTENT (IN)  :: NUMSRC               ! no. sources in speed array
            LOGICAL, INTENT (IN)  :: LASAFLAG             ! true: treat local roads as arterial
            
C.............  Local subprogram variables
            INTEGER :: K, N                ! counters and indices

C.............................................................................
            
            N = STPOS

C.............  Loop through sources and store information                
            DO K = STIDX, ENDIDX
                
                IF( N > NUMSRC ) EXIT
                
C.................  Store source number                    
                SPDARRAY( N,1 ) = K
                                
C.................  Convert facility type to road class and store
                SPDARRAY( N,3 ) = CVTRDTYPE( IRCLAS( K ), LASAFLAG )

C.................  Make sure we need to calculate emission factors for this source
                IF( VMT( K ) == 0 ) THEN
                    SPDARRAY( N,2 ) = IMISS3
                    N = N + 1
                    CYCLE
                END IF

C.................  If current source is not a local road, then store speed info
                IF( SPDARRAY( N,3 ) /= M6LOCAL ) THEN

C.....................  Store speed from inventory
C                       Rounding of speeds could be added here       
                    SPDARRAY( N,2 ) = SPEED( K )

C.....................  If current source uses a speed profile, store the negative ID
                    IF( SPDFLAG ) THEN
                        IF( SPDPROFID( K ) > 0 ) THEN
                            SPDARRAY( N,2 ) = -SPDPROFID( K )
                        END IF
                    END IF
                    
                END IF

                N = N + 1

            END DO
     
            END SUBROUTINE STORE_SOURCE_INFO

C.............................................................................
C.............................................................................
            
            SUBROUTINE WRITE_ARRAY_TO_SPDSUM( MDEV, SPDARRAY, NUMSRC, 
     &                                        COUNTY )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: MDEV                 ! SPDSUM file unit no.
            REAL,    INTENT (IN) :: SPDARRAY( NUMSRC,3 ) ! sorted array of spds, roads, sources 
            INTEGER, INTENT (IN) :: NUMSRC               ! no. sources in speed array
            INTEGER, INTENT (IN) :: COUNTY               ! county FIPS code

C.............  Local subprogram variables
            INTEGER :: K, L, M                ! counters and indices
            
            INTEGER NSRCSPEED                 ! no. sources for a particular speed (and road)
            INTEGER STFREEWAY                 ! starting index for freeway sources
            INTEGER STARTERIAL                ! starting index for arterial sources
            INTEGER STLOCAL                   ! starting index for local sources
            
            INTEGER CURRFREE                  ! current freeway position
            INTEGER CURRART                   ! current arterial position
            
            INTEGER ROADTYPE                  ! current road type for alternating printing
            
            REAL    CURRROAD                  ! current road type for county
            REAL    CURRSPEED                 ! current speed for county

C.............................................................................
            
C.............  Find the start of each road type
            STFREEWAY  = FINDR1FIRST( REAL( M6FREEWAY ), NUMSRC, 
     &                                SPDARRAY( :,3 ) )         
            STARTERIAL = FINDR1FIRST( REAL( M6ARTERIAL ), NUMSRC, 
     &                                SPDARRAY( :,3 ) )            
            STLOCAL    = FINDR1FIRST( REAL( M6LOCAL ), NUMSRC, 
     &                                SPDARRAY( :,3 ) )

C.............  If find comes back -1, we have no sources of that type, so set
C               starting position to next source starting position
            IF( STLOCAL <= 0 ) THEN
                STLOCAL = NUMSRC + 1
            END IF
            
            IF( STARTERIAL <= 0 ) THEN
                STARTERIAL = STLOCAL
            END IF
            
            IF( STFREEWAY <= 0 ) THEN
                STFREEWAY = STARTERIAL
            END IF
            
            CURRFREE = STFREEWAY
            CURRART  = STARTERIAL
            
            ROADTYPE = M6FREEWAY

C.............  Loop through sources for current road type and speed            
            DO            
                
C.................  Set index into speed array based on desired road type

                IF( ROADTYPE == M6FREEWAY ) THEN
                    K = CURRFREE
                    
C.....................  If we're done with freeway sources, go to arterial sources                    
                    IF( K >= STARTERIAL ) THEN
                        ROADTYPE = M6ARTERIAL
                    END IF
                END IF
                
                IF( ROADTYPE == M6ARTERIAL ) THEN
                    K = CURRART
                    
C.....................  If we're done with arterial sources, check if we still have
C                       freeway sources to output.  Otherwise, go to local sources
                    IF( K >= STLOCAL ) THEN
                        IF( CURRFREE >= STARTERIAL ) THEN
                            ROADTYPE = M6LOCAL
                        ELSE
                            ROADTYPE = M6FREEWAY
                            K = CURRFREE
                        END IF
                    END IF
                END IF
                
                IF( ROADTYPE == M6LOCAL ) THEN
                    K = STLOCAL
                END IF

C.................  Make sure we don't go out of bounds in the speed array
                IF( K > NUMSRC ) RETURN
                
C.................  Get speed of current source and make sure it's valid               
                CURRSPEED = SPDARRAY( K,2 )

                IF( CURRSPEED == IMISS3 ) THEN
                    DO
                        K = K + 1
                        IF( K > NUMSRC ) RETURN
                        CURRSPEED = SPDARRAY( K,2 )
                        IF( CURRSPEED /= IMISS3 ) EXIT
                    END DO
                END IF

C.................  Get road type of current source
                CURRROAD = SPDARRAY( K,3 )
                
                L = K + 1

C.................  Calculate the number of sources for this speed and road type                    
                DO
                    IF( L > NUMSRC .OR.
     &                  SPDARRAY( L,2 ) /= CURRSPEED .OR.
     &                  SPDARRAY( L,3 ) /= CURRROAD ) EXIT      
                    L = L + 1
                END DO

                NSRCSPEED = L - K

C.................  Write current speed to SPDSUM file                
                CALL WRITE_SPD_TO_SPDSUM( MDEV, SPDARRAY( K:L-1,: ), 
     &                                    NSRCSPEED, COUNTY )                           

C.................  Increment correct index based on current road type and set
C                   road type to opposite of current
                SELECT CASE( ROADTYPE )
                
                CASE( M6FREEWAY )
                    CURRFREE = K + NSRCSPEED
                    ROADTYPE = M6ARTERIAL
                    
                CASE( M6ARTERIAL )
                    CURRART = K + NSRCSPEED
                    ROADTYPE = M6FREEWAY
                    
C.................  Since there is only one line of local sources, if the road type
C                   is local, we're done                    
                CASE( M6LOCAL )
                    EXIT
                
                END SELECT              
            
            END DO  
                    
            END SUBROUTINE WRITE_ARRAY_TO_SPDSUM

C.............................................................................
C.............................................................................

            SUBROUTINE WRITE_SPD_TO_SPDSUM( MDEV, SPDARRAY, NUMSRC, 
     &                                        COUNTY )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: MDEV                 ! SPDSUM file unit no.
            REAL,    INTENT (IN) :: SPDARRAY( NUMSRC,3 ) ! sorted array of sources for one speed and road type 
            INTEGER, INTENT (IN) :: NUMSRC               ! no. sources in speed array
            INTEGER, INTENT (IN) :: COUNTY               ! county FIPS code

C.............  Local subprogram variables
            INTEGER :: M                      ! counters and indices
            
            INTEGER :: NSRCLINE = 7           ! no. sources per line in SPDSUM
            INTEGER NLINES                    ! no. lines in SPDSUM for a speed and road type
            
            INTEGER STSRCWR                   ! starting index of sources for line of SPDSUM
            INTEGER ENDSRCWR                  ! ending index of sources

            LOGICAL IDFLAG                    ! true: current spd is a profile ID

C.............................................................................

C.............  Check if current speed is actually a profile ID
            IF( SPDARRAY( 1,2 ) < 0 ) THEN
                IDFLAG = .TRUE.
            ELSE
                IDFLAG = .FALSE.
            END IF
            
C.............  Set starting and ending indices for first line
            STSRCWR  = 1
            ENDSRCWR = NSRCLINE
                
            IF( ENDSRCWR > NUMSRC ) THEN
                ENDSRCWR = NUMSRC
            END IF

C.............  If no. sources is more than NSRCLINE, use multiple lines
            NLINES = ( NUMSRC - 1 ) / NSRCLINE

C.............  Loop through all lines except last one (to write continuation character)
            DO M = 1, NLINES

C.................  If writing a profile ID, use different format (to avoid decimal point)
C                   and convert to positive number
                IF( IDFLAG ) THEN

#ifdef DBLBSLASH                
                    WRITE( MDEV,93020 ) COUNTY, 
     &                 INT( SPDARRAY( 1,3 ) ), -INT( SPDARRAY( 1,2 ) ), 
     &                 INT( SPDARRAY( STSRCWR:ENDSRCWR,1 ) ), '\\'
#else                
                    WRITE( MDEV,93020 ) COUNTY, 
     &                 INT( SPDARRAY( 1,3 ) ), -INT( SPDARRAY( 1,2 ) ), 
     &                 INT( SPDARRAY( STSRCWR:ENDSRCWR,1 ) ), '\'
#endif

                ELSE                

#ifdef DBLBSLASH
                    WRITE( MDEV,93010 ) COUNTY, 
     &                 INT( SPDARRAY( 1,3 ) ), SPDARRAY( 1,2 ), 
     &                 INT( SPDARRAY( STSRCWR:ENDSRCWR,1 ) ), '\\'
#else
                    WRITE( MDEV,93010 ) COUNTY, 
     &                 INT( SPDARRAY( 1,3 ) ), SPDARRAY( 1,2 ), 
     &                 INT( SPDARRAY( STSRCWR:ENDSRCWR,1 ) ), '\'
#endif

                END IF
                    
                STSRCWR  = STSRCWR  + NSRCLINE
                ENDSRCWR = ENDSRCWR + NSRCLINE
                    
            END DO

C.............  Check that ending index is not greater than actual no. sources                
            IF( ENDSRCWR > NUMSRC ) THEN
                ENDSRCWR = NUMSRC
            END IF
                
C.............  Write last line (may be only line) for this speed
            IF( STSRCWR <= ENDSRCWR ) THEN

                IF( IDFLAG ) THEN
                    WRITE( MDEV,93020 ) COUNTY,
     &                 INT( SPDARRAY( 1,3 ) ), -INT( SPDARRAY( 1,2 ) ), 
     &                 INT( SPDARRAY( STSRCWR:ENDSRCWR,1 ) )
                ELSE          	
                    WRITE( MDEV,93010 ) COUNTY,
     &                 INT( SPDARRAY( 1,3 ) ), SPDARRAY( 1,2 ), 
     &                 INT( SPDARRAY( STSRCWR:ENDSRCWR,1 ) )
                END IF
     
            END IF

C--------------  SUBPROGRAM FORMAT  STATEMENTS   --------------------------

93010       FORMAT( I6, 1X, I1, 1X, F6.2, 7( 1X, I6 ), 1X, 1A )  

93020       FORMAT( I6, 1X, I1, 1X, I6, 7( 1X, I6 ), 1X, 1A )
            
            END SUBROUTINE WRITE_SPD_TO_SPDSUM
        
        END SUBROUTINE WRSPDSUM
        
