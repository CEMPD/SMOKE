
        SUBROUTINE PRCLINUC( IREC, NSEGS, LINE, SEGMENT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    The PRCLINUC processes the uncertainty file by setting variables
C    that can be used by the calling routine(s) for futher processing.
C    This routine is used for several cases.  It is able to interact
C    with various other routines at various stages of processing. 
C    Settings are passed to other routines with modules.
C
C  PRECONDITIONS REQUIRED:
C    Uncertainty file is opened
C    LINE read, separated into SEGMENTS, and left-justified
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C    Subroutines:  I/O API subroutines
C
C    Functions:  I/O API functions, CRLF, INDEX1, STR2INT, STR2REAL
C
C  REVISION  HISTORY:
C     Created 8/2001 by A Holland
C
C***********************************************************************
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
C***********************************************************************

C...........   MODULES for public variables
C.........  This module contains uncertainty-specific settings
        USE MODUNCERT
        
C.........  This module contains the information about the source category
        USE MODINFO        


        IMPLICIT NONE
	
C...........   INCLUDES
    	INCLUDE  'EMCNST3.EXT'	


C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2   CRLF
        INTEGER       INDEX1
        INTEGER       STR2INT
	REAL	      STR2REAL

        EXTERNAL   CRLF, INDEX1, STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: IREC        ! line counter
        INTEGER     , INTENT (IN) :: NSEGS       ! no. line segments
        CHARACTER(*), INTENT (IN) :: LINE        ! full input record less commnt
        CHARACTER(*), INTENT (IN) :: SEGMENT( * )! parsed input record

C...........   Other local variables
        INTEGER          I, J, L, L2         ! counters and indices
        INTEGER          IOS                 ! i/o status

        INTEGER, SAVE :: IREC_MIN = 99999999 ! minimum record number
        INTEGER, SAVE :: PREC     = -9       ! previous call record number

        LOGICAL, SAVE :: FIRSTIME  = .TRUE.   ! true: first time routine called

        CHARACTER*300    MESG              !  message buffer

        CHARACTER*16 :: PROGNAME = 'PRCLINUC' ! program name

C***********************************************************************
C   begin body of subroutine PRCLINUC

C.........  Initialize for start of file
        IF( FIRSTIME .OR. IREC .EQ. IREC_MIN ) THEN

            INPACKET  = .FALSE.
            INFAPKT   = .FALSE.
            INEMPPKT  = .FALSE.
	    INPARPKT  = .FALSE.
	    
	    EMPSTART  = .FALSE.
	    FASTART   = .FALSE.

            PKTCOUNT  = 0         ! array
            PKTSTATUS = .FALSE.   ! array

            FIRSTIME = .FALSE.

        END IF


C.........  Initializations needed for every line
        INPARPKT  = .FALSE.
	EMPSTART  = .FALSE.
	FASTART   = .FALSE.

C.........  Set length of full line
        L2 = LEN_TRIM( LINE )

C.........  If a packet is not currently being read, determine what packet
C           starts on the current line
C.........  Will skip line if no packet has been found and the line does not
C           start with "/"
        IF( .NOT. INPACKET ) THEN

            J = INDEX( LINE( 2:L2 ), '/' )

C.............  Error if packet started but does not end with a slash
            IF( LINE( 1:1 ) .EQ. '/' .AND.
     &                J     .LE. 0         ) THEN

                UC_ERROR = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Packet at line', IREC,
     &                 'started but not finished.'
                CALL M3MSG2( MESG )

C.............  Packet found. Interpret it...
            ELSE IF( J .GT. 0 ) THEN
                J = J + 1

C.................  Set packet name and extract group label
                PCKTNAM = LINE( 1:J )

C.................  Find packet name in list of valid names
                PKT_IDX = INDEX1( PCKTNAM, NUCPCKT, UCPCKTS )

C.................  Warning if packet is not recognized
                IF( PKT_IDX .LE. 0 ) THEN
                    INPACKET = .TRUE.
                    WRITE( MESG,94010 ) 'WARNING: Packet "' //
     &                     PCKTNAM( 1:J ) // '" in uncertainty file'
     &                     // CRLF() // BLANK10 //' file at line', IREC,
     &                     'is not recognized. Packet will be ignored.'
                    CALL M3MESG( MESG )

C.................  Packet recognized
                ELSE

C.....................  Count the number of groups and the number of reports
                    PKTCOUNT ( PKT_IDX ) = PKTCOUNT( PKT_IDX ) + 1
                    PKTSTATUS( PKT_IDX ) = .TRUE.

                    INPACKET = .TRUE.

C.....................  Set first line number of packet
                    PKTSTART = IREC

C.....................  Packet-specific processing
                    SELECT CASE( PKT_IDX )

C.....................  A factor assignment packet
                    CASE( FA_IDX )
                       
                        INPACKET = .TRUE.            
			INFAPKT  = .TRUE.
			FASTART  = .TRUE.
			FAENTN   = 0

C.....................  An empirical packet
                    CASE( EMP_IDX )

			EMP_%ETYPE  = SEGMENT( 2 )                        
                        IF( SEGMENT( 2 ) .NE. 'S' .AND. 
     &                       SEGMENT( 2 ) .NE. 'L') THEN
                          L = LEN_TRIM( SEGMENT( 2 ) )
                          WRITE( MESG,94010 ) 'ERROR: 
     &   	            Empirical type ' //
     &                      'setting "' // SEGMENT( 2 )( 1:L ) // 
     &                      '" at line', IREC, 'is invalid. 
     &                      Must be ' //
     &                      '"S" or "L".'
                          CALL M3MSG2( MESG )
                        END IF
			    
			  EMP_%EMPNAM  = SEGMENT( 3 )
		    
			
			  EMPENTN = 0
			  INEMPPKT  = .TRUE.
			  EMPSTART  = .TRUE.
                          INPACKET  = .TRUE.    
			          

C.....................  A parametric packet
                    CASE( PAR_IDX )

			PAR_%PARNAM  = SEGMENT( 2 )
			PAR_%PTYPE   = SEGMENT( 3 )
                        IF( SEGMENT( 3 ).NE.'N'.AND.SEGMENT( 3 )
     &                      .NE.'L'.AND.SEGMENT( 3 ).NE.'G'.AND.
     &                      SEGMENT( 3 ).NE.'W'.AND.SEGMENT( 3 )
     &                      .NE.'B') THEN
			  L = LEN_TRIM( SEGMENT( 3 ) )
                          WRITE( MESG,94010 ) 'ERROR: 
     &			    Parametric type ' //
     &                      'setting "' // SEGMENT( 3 )( 1:L ) // 
     &                      '" at line', IREC, 'is invalid. 
     &                      Must be ' //
     &                      '"N", "L", "G","W" or "B".'
                          CALL M3MSG2( MESG )
                         END IF
			    
			PAR_%NUMP    = STR2INT( SEGMENT( 4 ) )
			
			DO I = 1, PAR_%NUMP
			
			  PARA( I )   = STR2REAL( SEGMENT( 4 + I ) )
			  
			END DO

                        PKTEND   = IREC
                        INPACKET = .FALSE.             ! end implied
			INPARPKT = .TRUE.


                     END SELECT

                END IF
 
            END IF

            GO TO 999    ! to end of routine

C........  End of packet...
       ELSE IF( SEGMENT( 1 ) .EQ. '/END/' ) THEN

C.............  Reset all packet-specific settings
            PKTEND     = IREC
            INPACKET   = .FALSE.
            INFAPKT    = .FALSE.
            INEMPPKT   = .FALSE.


            GO TO 999          ! to end of routine

C........  Missing end of packet and already in a packet
       ELSE IF( LINE( 1:1 ) .EQ. '/' ) THEN
            UC_ERROR = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Missing end of packet ' //
     &             'at line', IREC
            CALL M3MSG2( MESG )

            GO TO 999          ! to end of routine

       END IF

C.........  General packet processing
        IF( INFAPKT ) THEN

		FA_%CFIP   = SEGMENT( 1 ) 
		FA_%TSCC   = SEGMENT( 2 )
                FA_%CPOL   = SEGMENT( 3 )                                
		FA_%METH   = SEGMENT( 4 )
                IF( SEGMENT( 4 ).NE.'E'.AND.SEGMENT( 4 )
     &              .NE.'P') THEN
		  L = LEN_TRIM( SEGMENT( 4 ) )
                  WRITE( MESG,94010 ) 'ERROR: 
     &		    Factor assignment method ' //
     &              'setting "' // SEGMENT( 4 )( 1:L ) // 
     &              '" at line', IREC, 'is invalid. 
     &              Must be ' //
     &              '"E" or "P".'
                  CALL M3MSG2( MESG )
                 END IF
		    
		FA_%NDIST  = SEGMENT( 5 )
		FA_%APRCH  = SEGMENT( 6 )
                IF( SEGMENT( 6 ).NE.'S'.AND.SEGMENT( 6 )
     &              .NE.'I'.AND.SEGMENT( 6 ).NE.'ST'.AND.
     &              SEGMENT( 6 ).NE.'IT') THEN
		L = LEN_TRIM( SEGMENT( 6 ) )
                  WRITE( MESG,94010 ) 'ERROR: 
     &		    Factor assignment approach ' //
     &              'setting "' // SEGMENT( 6 )( 1:L ) // 
     &              '" at line', IREC, 'is invalid. 
     &              Must be ' //
     &              '"S", "I", "ST" or "IT".'
                  CALL M3MSG2( MESG )
                END IF
		     
		FA_%PLT    = SEGMENT( 7 )
		FA_%CHAR1  = SEGMENT( 8 )
		FA_%CHAR2  = SEGMENT( 9 )
		FA_%CHAR3  = SEGMENT( 10 )
		FA_%CHAR4  = SEGMENT( 11 )
		FA_%CHAR5  = SEGMENT( 12 )
		
		FAENTN = FAENTN + 1

            
        END IF      ! End factor assignment packet


        IF( INEMPPKT ) THEN
	
		EMISFAC  = STR2REAL( SEGMENT( 1 ) )
		EMPROB   = STR2REAL( SEGMENT( 2 ) )
		
		EMPENTN = EMPENTN + 1

        END IF              ! Empirical packet

C.........  Save variables from current call and exit from subroutine
999     CONTINUE

        PREC = IREC
        IREC_MIN  = MIN( IREC, IREC_MIN )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE PRCLINUC

