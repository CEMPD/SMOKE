
        SUBROUTINE WRUCOUT( UONAME, UPNAME, UENAME )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine writes to the uncertainty output files.
C
C  PRECONDITIONS REQUIRED:
C      Logical file name UPNAME and UENAME opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C      Functions: I/O API functions
C
C  REVISION  HISTORY:
C      Created 9/2001 by A. Holland
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains uncertainty-specific settings
        USE MODUNCERT
        
C.........  This module contains the information about the source category
        USE MODINFO        
        
        
      IMPLICIT NONE

C.........  INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C.........  Subroutine arguments and their descriptions:
        CHARACTER(*), INTENT (IN) :: UONAME     ! I/O API output file name
        CHARACTER(*), INTENT (IN) :: UPNAME     ! parametric output file name
        CHARACTER(*), INTENT (IN) :: UENAME     ! empirical output file name
        

C...........   Other local variables
        INTEGER                 I, J, K, L        ! counters and indicies
        
        LOGICAL, SAVE  ::  FIRSTIME = .TRUE.      ! true if first time 
        LOGICAL, SAVE  ::  SECONDTIME = .TRUE.    ! true if second time

        CHARACTER*300           MESG              !  message buffer
        CHARACTER(LEN=IOVLEN3)  CBUF              !  tmp pollutant name
 
        CHARACTER*16 :: PROGNAME = 'WRUCOUT'      !  program name

C***********************************************************************
C   begin body of program WRUCOUT


        IF( FIRSTIME ) THEN
        
            MESG = 'Writing uncertainty output file...'
            CALL M3MSG2( MESG )
        
C.........  Create message to use in case there is an error
            MESG = 'Error writing uncertainty output file "' //
     &             UONAME( 1:LEN_TRIM( UONAME ) ) // '"'        
        
C.........  Write the first I/O API file, one variable at a time

            IF( .NOT. WRITE3( UONAME, 'SRCNUM', 0, 0, SRCNUM ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            DO I = 1, NUOVAR
            
                CBUF = UONAMES( I )
                L = LEN_TRIM( CBUF )
        
                IF( .NOT. WRITE3( UONAME, 'MTH_'//CBUF(1:L),
     &                            0, 0, METHOD(1,I) ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                IF( .NOT. WRITE3( UONAME, 'TYP_'//CBUF(1:L),
     &                            0, 0, EPTYP(1,I) ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                IF( .NOT. WRITE3( UONAME, 'NEP_'//CBUF(1:L),
     &                            0, 0, NUMEP(1,I) ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
                
                IF( .NOT. WRITE3( UONAME, 'UIX_'//CBUF(1:L),
     &                            0, 0, UNCIDX(1,I) ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                
            END DO
            
            FIRSTIME = .FALSE.
            
        ELSE IF( SECONDTIME ) THEN
        
            MESG = 'Writing uncertainty output parametric file...'
            CALL M3MSG2( MESG )
        
C.........  Write the second I/O API file, one variable at a time            
        
C.........  Create message to use in case there is an error
            MESG = 'Error writing parametric output file "' //
     &             UPNAME( 1:LEN_TRIM( UPNAME ) ) // '"'        


            IF( .NOT. WRITE3( UPNAME, 'PARMS', 0, 0, PARMS ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
            SECONDTIME = .FALSE.
            
            
        ELSE
        
            MESG = 'Writing uncertainty output empirical file...'
            CALL M3MSG2( MESG )
       
            
C.........  Write the third I/O API file, one variable at a time            
        
C.........  Create message to use in case there is an error
            MESG = 'Error writing empirical output file "' //
     &             UENAME( 1:LEN_TRIM( UENAME ) ) // '"'        


            IF( .NOT. WRITE3( UENAME, 'EMFVAL', 0, 0, EMFVAL ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        
            IF( .NOT. WRITE3( UENAME, 'PROBVAL', 0, 0, PROBVAL ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93100   FORMAT( I2 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94100   FORMAT( 9( A, I2 ) )

        END SUBROUTINE WRUCOUT
