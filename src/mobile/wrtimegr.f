
        SUBROUTINE WRTIMEGR

C***********************************************************************
C  subroutine body starts at line 82
C
C  DESCRIPTION:
C       Writes the time group output files
C
C  PRECONDITIONS REQUIRED: none
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
C***********************************************************************

C.........  MODULES for public variables
C.........  This module is used for MOBILE6 setup information 
        USE MODMBSET, ONLY: NINVC, NREFC, MCREFIDX, MCREFSORT, 
     &                      MVREFSORT, DAILY, WEEKLY, MONTHLY, EPISLEN
        
        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER, EXTERNAL :: PROMPTFFILE
        
C...........   Other local variables         
        INTEGER I, J                      ! counters and indices
        INTEGER CURRTIME                  ! current time period
        INTEGER LASAFLAG                  ! treat local as arterial setting

        INTEGER CURRDEV                   ! current file unit no.
        INTEGER DAYDEV                    ! unit number for daily time period group file 
        INTEGER WEEKDEV                   ! unit number for weekly group file
        INTEGER MONTHDEV                  ! unit number for monthly group file
        INTEGER EPISDEV                   ! unit number for episode length time period group
        
        LOGICAL :: OPENDAY   = .FALSE.    ! true: opened day group file
        LOGICAL :: OPENWEEK  = .FALSE.    ! true: opened week group file
        LOGICAL :: OPENMONTH = .FALSE.    ! true: opened month group file
        LOGICAL :: OPENEPIS  = .FALSE.    ! true: opened episode length group file

        CHARACTER(LEN=FIPLEN3) REFCOUNTY                 ! ref. county FIPS code
        CHARACTER(LEN=FIPLEN3) INVCOUNTY                 ! inv. county FIPS code
                
        CHARACTER*16 :: PROGNAME = 'WRTIMEGR'   ! program name

C***********************************************************************
C   begin body of subroutine WRTIMEGR

C.........  Loop through ref. counties in MVREFSORT
        DO I = 1, NREFC
            WRITE( REFCOUNTY, '(I6)' ) MVREFSORT( I,1 )
            CALL PADZERO( REFCOUNTY )

C.............  Get time period for current ref. county
            CURRTIME = MVREFSORT( I,3 )

C.............  Set file unit based on current time period
C               Open appropriate files if necessary            
            SELECT CASE( CURRTIME )
            
            CASE( DAILY )
            	IF( .NOT. OPENDAY ) THEN
            	    DAYDEV = PROMPTFFILE(
     &                       'Enter logical name for DAILYGROUP file',
     &                       .FALSE., .TRUE., 'DAILYGROUP', PROGNAME )
                    OPENDAY = .TRUE.
                END IF
                
                CURRDEV = DAYDEV
                
            CASE( WEEKLY )
                IF( .NOT. OPENWEEK ) THEN
                    WEEKDEV = PROMPTFFILE(
     &                        'Enter logical name for WEEKLYGROUP file',
     &                        .FALSE., .TRUE., 'WEEKLYGROUP', PROGNAME )
                    OPENWEEK = .TRUE.
                END IF
                
                CURRDEV = WEEKDEV
                
            CASE( MONTHLY )
                IF( .NOT. OPENMONTH ) THEN
                    MONTHDEV = PROMPTFFILE(
     &                       'Enter logical name for MONTHLYGROUP file',
     &                       .FALSE., .TRUE., 'MONTHLYGROUP', PROGNAME )
                    OPENMONTH = .TRUE.
                END IF
                
                CURRDEV = MONTHDEV
                
            CASE( EPISLEN )
                IF( .NOT. OPENEPIS ) THEN
                    EPISDEV = PROMPTFFILE(
     &                       'Enter logical name for EPISODEGROUP file',
     &                       .FALSE., .TRUE., 'EPISODEGROUP', PROGNAME )
                    OPENEPIS = .TRUE.
            	END IF
            	
            	CURRDEV = EPISDEV
            	
            END SELECT

C.............  Get local-to-arterial setting for this ref. county
            LASAFLAG = MVREFSORT( I,4 )
        	
C.............  Check if this ref. county is not spatially averaged            	
            IF( MVREFSORT( I,2 ) == 1 ) THEN

C.................  Loop through inv. counties using this ref. county
                J = MCREFIDX( I,2 )

                DO
                    WRITE( INVCOUNTY, '(I6)' ) MCREFSORT( J,1 )
                    CALL PADZERO( INVCOUNTY )
                    
                    WRITE( CURRDEV,93010 ) 
     &                     INVCOUNTY, REFCOUNTY, LASAFLAG
                    
                    J = J + 1
                    
                    IF( J > NINVC ) EXIT
                    IF( I /= NREFC ) THEN
                        IF( J == MCREFIDX( I + 1,2 ) ) EXIT
                    END IF

                END DO

C.............  Otherwise, write the ref. county number
            ELSE
                WRITE( CURRDEV,93010 ) 
     &                 REFCOUNTY, REFCOUNTY, LASAFLAG
            END IF
                    
        END DO ! end ref. county loop

        IF( OPENDAY )   CLOSE( DAYDEV )
        IF( OPENWEEK )  CLOSE( WEEKDEV )
        IF( OPENMONTH ) CLOSE( MONTHDEV )
        IF( OPENEPIS )  CLOSE( EPISDEV )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )  
93010   FORMAT( A6, 1X, A6, 1X, I1 )  
        
        END SUBROUTINE WRTIMEGR
        