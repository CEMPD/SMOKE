
        SUBROUTINE RDEPROC( FDEV )

C***********************************************************************
C  subroutine body starts at line 103
C
C  DESCRIPTION:
C     This subroutine reads the emission processes file, which contains columns
C     for the activity, associated process, and associated pollutants.  If
C     there is more than one process per activity, then these are listed on
C     separate lines in the file. If there is more than one pollutant per
C     activity and process, then these are listed in additional columns.
C     During the read, the column number is set dynamically.  Only processes
C     for activities that are in the inventory are read.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 10/99 by M. Houyoux
C
C****************************************************************************/
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
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER     GETFLINE
        INTEGER     GETNLIST
        INTEGER     INDEX1

        EXTERNAL    GETFLINE, GETNLIST, INDEX1

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! cross-reference file unit no.

C...........   Local allocatable arrays
        INTEGER               , ALLOCATABLE :: INDX   ( : )  ! POLA sorting indx
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: SEGMENT( : )  ! line segments
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: POLNAM ( : )  ! pol names
        CHARACTER(LEN=IOVLEN3), ALLOCATABLE :: POLA   ( : )  ! all unsorted pols
        CHARACTER(LEN=IOULEN3), ALLOCATABLE :: RAWSUBS( :,: ) ! raw subtract pollutant list
        CHARACTER(LEN=IOULEN3), ALLOCATABlE :: SUBHCS ( :,: ) ! hydrocarbon names from packet

C...........   Local parametes
        CHARACTER, PARAMETER :: CONTCHAR = '\'  ! line continuation character

C...........   Other local variables
        INTEGER         I, J, K, L, M, V    !  counters and indices

        INTEGER         IOS     !  i/o status
        INTEGER         L1, L2  !  tmp string lengths
        INTEGER         LJ      !  string length for emis type joiner
        INTEGER         MXCOLS  !  maximum number of columns in the file
        INTEGER         MXPOL   !  maximum number of pollutants per process
        INTEGER         MXSUBPOL!  maximum number of pollutants in subtract packet
        INTEGER         NCOLS   !  no. columns in a row
        INTEGER         NLINES  !  number of lines in file
        INTEGER         NPOL    !  no. pollutants in a row
        INTEGER         NPUNSRT !  no. pols 
        INTEGER         NPCKTS  !  no. packets
        INTEGER         NSUBS   !  no. pols in subtract packet
        INTEGER         USEPCKT !  no. of packet to use
        INTEGER         HCIDX   !  index of hydrocarbon pollutant
        INTEGER         NTLINES !  no. of lines taking into account continuation lines 

        LOGICAL      :: INPACKET = .FALSE.   ! true: inside packet
        LOGICAL      :: FNDPOL   = .FALSE.   ! true: found pollutant on master list
        LOGICAL      :: NEWLINE  = .TRUE.    ! true: current line is new (not continued)
        LOGICAL      :: EFLAG    = .FALSE.   ! true: error found

        CHARACTER*300          LINE     !  line buffer
        CHARACTER*300          MESG     !  message buffer
        CHARACTER(LEN=IOVLEN3) ACT      !  tmp activity name
        CHARACTER(LEN=IOVLEN3) CPOL     !  tmp pollutant name
        CHARACTER(LEN=IOVLEN3) LPOL     !  tmp pollutant name from previous iter
        CHARACTER(LEN=IOVLEN3) PRC      !  tmp process name
        CHARACTER(LEN=IOVLEN3) SUBLINE( 4 ) ! pieces of subtract packet

        CHARACTER*16 :: PROGNAME = 'RDEPROC' ! program name

C***********************************************************************
C   begin body of subroutine RDEPROC

C.........  Set input and output hydrocarbon names to default, 
C           in case there are no packets
        INPUTHC  = ' ' 
        OUTPUTHC = ' '

C.........  Get the number of lines for the file and allocate array so that
C           the type of the line can be stored
        NLINES = GETFLINE( FDEV, 'Emission processes file' )
        NTLINES = NLINES

C.........  Read through file to determine the maximum no. of columns and
C           check packet information
        MXCOLS   = 0
        MXSUBPOL = 0
        NPCKTS   = 0
        NSUBS    = 0
        
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading emission processes file at line', I
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip any blank lines
            IF( LINE == ' ' ) THEN
                NTLINES = NTLINES - 1
                CYCLE
            END IF

            LINE = ADJUSTL( LINE )
            L1 = LEN_TRIM( LINE )
            
C.............  Check if this line is a packet
            IF( LINE( 1:1 ) == '/' ) THEN
                NTLINES = NTLINES - 1
            
C.................  Make sure packet ends correctly
                IF( INDEX( LINE( 2:L1 ), '/' ) <= 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Packet at line', I,
     &                     'started but not finished.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF
                
C.................  Check if this is an ending packet
                IF( LINE( 1:5 ) == '/END/' ) THEN
                    NSUBS = 0
                    INPACKET = .FALSE.
                    CYCLE
                END IF
                
C.................  Make sure we're not already in a packet
                IF( INPACKET ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: No end packet found' //
     &                     'before starting new packet at line', 
     &                     I, '.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF
                
                NPCKTS = NPCKTS + 1
                INPACKET = .TRUE.
                CYCLE
            END IF  

C.............  Check if inside a packet
            IF( INPACKET ) THEN
            	NTLINES = NTLINES - 1
                NSUBS = NSUBS + 1
                IF( NSUBS > MXSUBPOL ) THEN
                    MXSUBPOL = NSUBS
                END IF
                CYCLE
            END IF

C.............  If previous line was continued, add this lines columns to number
C               from previous line; otherwise, count columns in just this line
            IF( .NOT. NEWLINE ) THEN
                NCOLS = NCOLS + GETNLIST( L1, LINE )
            ELSE
                NCOLS = GETNLIST( L1, LINE )
            END IF

C.............  Check for continuation character in current line
            IF( LINE( L1:L1 ) == CONTCHAR ) THEN
                NEWLINE = .FALSE.
                NTLINES = NTLINES - 1
                NCOLS = NCOLS - 1
            ELSE
                NEWLINE = .TRUE.
            END IF
            
            IF( NCOLS > MXCOLS ) THEN
                MXCOLS = NCOLS
            END IF

        END DO

C.........  Check for errors so far
        IF( EFLAG ) THEN
            MESG = 'Problem reading emission processes file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MXPOL = MXCOLS - 2

C.........  Allocate memory for parsing line segements and storing pollutants
        ALLOCATE( SEGMENT( MXCOLS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )
        ALLOCATE( POLNAM( MXPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLNAM', PROGNAME )

C.........  Allocate memory for toxic pollutants
        IF( NPCKTS > 0 ) THEN
            ALLOCATE( RAWSUBS( MXSUBPOL,NPCKTS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'RAWSUBS', PROGNAME )
            ALLOCATE( SUBHCS( NPCKTS,2 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SUBHCS', PROGNAME )
            
            RAWSUBS = ' '  ! array
            SUBHCS = ' '   ! array
        END IF

C.........  Allocate memory for emission types and count for each per activity
        ALLOCATE( EMTNAM( NTLINES*MXPOL, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTNAM', PROGNAME )
        ALLOCATE( EMTIDX( NTLINES*MXPOL, NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTIDX', PROGNAME )
        ALLOCATE( NETYPE( NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NETYPE', PROGNAME )

        EMTNAM = ' '  ! array
        EMTIDX = 0    ! array
        NETYPE = 0    ! array

C.........  Rewind file
        REWIND( FDEV )

        INPACKET = .FALSE.
        NEWLINE  = .TRUE.

C.........  Store contents of emissions processes file in output order
        J = 0
        L = 0
        M = 0
        DO I = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

C.............  Skip blank line
            IF( LINE .EQ. ' ' ) CYCLE

            LINE = ADJUSTL( LINE )
            L1 = LEN_TRIM( LINE )
            
C.............  Check for packet
            IF( LINE( 1:1 ) == '/' ) THEN

C.................  Check if this is an ending packet
                IF( LINE( 1:5 ) == '/END/' ) THEN
                    L = 0
                    INPACKET = .FALSE.
                    CYCLE
                END IF

                M = M + 1
                
C.................  Separate line into segments
                CALL PARSLINE( LINE, 4, SUBLINE )

C.................  Store hydrocarbon information
                SUBHCS( M,1 ) = SUBLINE( 2 )
                SUBHCS( M,2 ) = SUBLINE( 4 )
                
                INPACKET = .TRUE.
                CYCLE
            END IF
           
            IF( INPACKET ) THEN
                L = L + 1
                RAWSUBS( L,M ) = LINE( 1:L1 )
                CYCLE
            END IF
                        
C.............  Separate line into segments
            NCOLS = GETNLIST( L1, LINE )
            CALL PARSLINE( LINE, NCOLS, SEGMENT )

            IF( NEWLINE ) THEN
                ACT = SEGMENT( 1 )
                PRC = SEGMENT( 2 )

                IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
                    NPOL = NCOLS - 3
                    POLNAM( 1:NPOL ) = SEGMENT( 3:NCOLS-1 )
                ELSE
                    NPOL = NCOLS - 2
                    POLNAM( 1:NPOL ) = SEGMENT( 3:NCOLS )
                END IF
            ELSE
            	IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
            	    NPOL = NCOLS - 1
            	    POLNAM( 1:NPOL ) = SEGMENT( 1:NCOLS-1 )
            	ELSE
            	    NPOL = NCOLS
            	    POLNAM( 1:NPOL ) = SEGMENT( 1:NCOLS )
            	END IF
            END IF

C.............  Make sure activity is in the inventory
            K = INDEX1( ACT, NIACT, ACTVTY )

C.............  Store emission processes and associated pollutants
            IF( K .GT. 0 ) THEN

                DO V = 1, NPOL
                    J = J + 1
                    L1 = LEN_TRIM( PRC )
                    L2 = LEN_TRIM( POLNAM( V ) )
                    EMTNAM( J,K ) = PRC( 1:L1 ) // ETJOIN // 
     &                              POLNAM( V )( 1:L2 )
                END DO

                NETYPE( K ) = NETYPE( K ) + NPOL

            END IF

C.............  Check if current line is continued
            IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
                NEWLINE = .FALSE.
            ELSE
                NEWLINE = .TRUE.
            END IF
    
        END DO      ! End loop over file

C.........  Set the maximum number of emission types
        MXETYPE = MAXVAL( NETYPE )

C.........  Create a list of pollutants associated with the emission types...

C.........  Allocate memory for unsorted pollutants list
        ALLOCATE( POLA( MXETYPE * NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLA', PROGNAME )
        ALLOCATE( INDX( MXETYPE * NIACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDX', PROGNAME )

C.........  Create unsorted pollutants list
        J = 0
        LJ = LEN_TRIM( ETJOIN )
        DO I = 1, NIACT

            DO K = 1, NETYPE( I )
                J = J + 1

                L = INDEX( EMTNAM( K,I ), ETJOIN )
                L2 = LEN_TRIM( EMTNAM( K,I ) )

                POLA( J ) = EMTNAM( K,I )( L+LJ:L2 )
                INDX( J ) = J

            END DO

        END DO
        NPUNSRT = J

C.........  Sort pollutants
        CALL SORTIC( NPUNSRT, INDX, POLA )

C.........  Determine number of actual pollutants
        LPOL = '-9'
        K = 0
        DO I = 1, NPUNSRT

            J = INDX( I )
            CPOL = POLA( J )

            IF( CPOL .NE. LPOL ) THEN
                K = K + 1
                LPOL = CPOL
            END IF

        END DO

        NEPOL = K

C.........  Allocate memory for sorted pollutants
        ALLOCATE( EMTPOL( NEPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMTPOL', PROGNAME )

C.........  Store sorted pollutants in a unique list
        LPOL = '-9'
        K = 0
        DO I = 1, NPUNSRT

            J = INDX( I )
            CPOL = POLA( J )

            IF( CPOL .NE. LPOL ) THEN
                K = K + 1
                EMTPOL( K ) = CPOL
                LPOL = CPOL
            END IF

        END DO

C.........  Create an index that references between the emission types and the
C           list of unique pollutants
        DO I = 1, NIACT

            DO K = 1, NETYPE( I )

                L = INDEX( EMTNAM( K,I ), ETJOIN )
                L2 = LEN_TRIM( EMTNAM( K,I ) )

                CPOL = EMTNAM( K,I )( L+LJ:L2 )
                J = INDEX1( CPOL, NEPOL, EMTPOL )
                EMTIDX( K,I ) = J

            END DO

        END DO

C.........  If we have subtraction packets, determine which set to use
        IF( NPCKTS > 0 ) THEN
            USEPCKT = 0

C.............  Loop through packets to find matching hydrocarbon type            
            DO I = 1, NPCKTS
                DO J = 1, NEPOL
                    IF( SUBHCS( I,1 ) == EMTPOL( J ) ) THEN
                    	HCIDX = J
                        USEPCKT = I
                        EXIT
                    END IF
                END DO
                
                IF( USEPCKT > 0 ) EXIT
            END DO
            
            IF( USEPCKT == 0 ) THEN
                MESG = 'ERROR: No match found between /SUBTRACT/ ' //
     &                 'packets and pollutants in emission ' //
     &                 'processes file.'
                CALL M3MSG2( MESG )
            END IF

C.............  Store input and output hydrocarbon names
            INPUTHC  = SUBHCS( USEPCKT,1 )
            OUTPUTHC = SUBHCS( USEPCKT,2 )
            EMTPOL( HCIDX ) = OUTPUTHC
            
            DO I = 1, NIACT
                DO K = 1, NETYPE( I )
                    IF( EMTIDX( K,I ) == HCIDX ) THEN
                        J = INDEX( EMTNAM( K,I ), TRIM( INPUTHC ) )
                        L = LEN( EMTNAM( K,I ) )
                        EMTNAM( K,I )( J:L ) = OUTPUTHC
                    END IF
                END DO
            END DO

C.............  Determine total number of pollutants for this packet
            NSUBS = 0
            DO I = 1, MXSUBPOL
                IF( RAWSUBS( I,USEPCKT ) /= ' ' ) THEN

                    FNDPOL = .FALSE.
C.....................  Check that pollutant is in master pollutant list
                    DO J = 1, NEPOL
                        IF( RAWSUBS( I,USEPCKT ) == EMTPOL( J ) ) THEN
                            FNDPOL = .TRUE.
                            EXIT
                        END IF
                    END DO
                    
                    IF( FNDPOL ) THEN
                        NSUBS = NSUBS + 1
                    ELSE
                        MESG = 'WARNING: Skipping pollutant ' //
     &                         TRIM( RAWSUBS( I,USEPCKT ) ) // ' in ' //
     &                         '/SUBTRACT/ packet since it is not ' //
     &                         'listed in the emission processes file.'
                        CALL M3MSG2( MESG )
                    END IF
                END IF
            END DO

C.............  Save total number of pollutants
            NSUBPOL = NSUBS

            IF( NSUBS > 0 ) THEN
                ALLOCATE( SUBPOLS( NSUBS ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SUBPOLS', PROGNAME )

C.................  Store names of pollutants            
                L = 0
                DO I = 1, MXSUBPOL
                    IF( RAWSUBS( I,USEPCKT ) /= ' ' ) THEN
                    	FNDPOL = .FALSE.
                    	
                    	DO J = 1, NEPOL
                    	    IF( RAWSUBS(I,USEPCKT) == EMTPOL( J ) ) THEN
                    	        FNDPOL = .TRUE.
                    	        EXIT
                    	    END IF
                    	END DO
                    		
                    	IF( FNDPOL ) THEN
                            L = L + 1
                            SUBPOLS( L ) = RAWSUBS( I,USEPCKT )
                        END IF
                    END IF
                END DO
            END IF
            	
        END IF

C.........  Rewind file
        REWIND( FDEV )

C.........  Deallocate memory for local arrays
        DEALLOCATE( SEGMENT, POLNAM, POLA, INDX )
        IF( ALLOCATED( RAWSUBS ) ) THEN
            DEALLOCATE( RAWSUBS, SUBHCS )
        END IF

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of emission processes file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDEPROC
