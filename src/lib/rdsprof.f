
        SUBROUTINE RDSPROF( FDEV, POLNAM, NPROFMAX, NPROF, NMSPC,
     &                      INPRF, SPECID, MOLEFACT, MASSFACT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the speciation profile, sorts it, and returns
C      the sorted data
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
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
C***************************************************************************

        IMPLICIT NONE
        
C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2   CRLF
        REAL          STR2REAL

        EXTERNAL      CRLF, STR2REAL

C...........   Subroutine arguments

        INTEGER     , INTENT    (IN) :: FDEV        ! file unit number
        CHARACTER(*), INTENT    (IN) :: POLNAM      ! pol name of interest
        INTEGER     , INTENT    (IN) :: NPROFMAX    ! max no. expected profiles
        INTEGER     , INTENT   (OUT) :: NPROF       ! no. profiles found
        INTEGER     , INTENT   (OUT) :: NMSPC       ! no. unique species IDs
        CHARACTER(*), INTENT   (OUT) :: INPRF   ( NPROFMAX ) ! spec prof numbers
        CHARACTER(*), INTENT(IN OUT) :: SPECID  ( NPROFMAX ) ! names of species
        REAL        , INTENT(IN OUT) :: MOLEFACT( NPROFMAX ) ! mole-based facs
        REAL        , INTENT(IN OUT) :: MASSFACT( NPROFMAX ) ! mass-based facs
                                
C...........   Local unsorted arrays

        INTEGER        INDXA ( NPROFMAX )   ! sorting index
        REAL           DIVISA( NPROFMAX )   ! unsorted divisors
        REAL           FACTRA( NPROFMAX )   ! unsorted split factors
        REAL           XMFA  ( NPROFMAX )   ! unsorted mass fraction
        CHARACTER*21   INPSPA( NPROFMAX )   ! unsorted profile no. // species ID
        CHARACTER*16   SPCIDA( NPROFMAX )   ! unsorted species IDs
        
C...........   Other local variables

        INTEGER         I, J      ! counters and indices
        INTEGER         IPS       ! counter for records in which POLID = POLNAM
        INTEGER         IOS       ! i/o error status
        INTEGER         IREC      ! record counter

        REAL            DIVISATP          ! tmp divisor
        REAL            FACTRATP          ! tmp split factor
        REAL            XMFATP            ! tmp mass fraction

        LOGICAL      :: EFLAG = .FALSE.   ! error flag

        CHARACTER*16    POLID             ! tmp pollutant ID
        CHARACTER*16    SPCIDTMP          ! tmp species ID
        CHARACTER*200   LINE              ! read buffer for a line
        CHARACTER*300   MESG              ! text for M3EXIT()

        CHARACTER*16 :: PROGNAME = 'RDSPROF' ! program name
       
C***********************************************************************
C   Begin body of subroutine RDSPROF

C...........  Read in speciation profile
        IREC  = 0
        IPS   = 0
        DO
        
            READ( FDEV, 93000, END=12, IOSTAT=IOS ) LINE
     
            IREC = IREC + 1
             
            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &              'I/O error', IOS, 'reading speciation profile '//
     &              'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            POLID = ADJUSTL( LINE( 7:22 ) )           
            IF ( POLID .NE. POLNAM ) CYCLE
            
            IPS = IPS + 1
            
            IF ( IPS .LE. NPROFMAX ) THEN
            
                INDXA ( IPS ) = IPS
                SPCIDTMP      = ADJUSTL ( LINE( 24:39 ) )
                INPSPA( IPS ) = ADJUSTR ( LINE(  1:5  ) ) // SPCIDTMP
                SPCIDA( IPS ) = SPCIDTMP
                FACTRA( IPS ) = STR2REAL( LINE( 41:53 ) )
                DIVISA( IPS ) = STR2REAL( LINE( 55:67 ) )
                XMFA  ( IPS ) = STR2REAL( LINE( 69:77 ) )

            END IF

        END DO
            
12      CONTINUE    !End of read on input file
     
        NPROF = IPS

        IF( NPROF .GT. NPROFMAX ) THEN  ! Check for memory overflow

            WRITE( MESG, 94010 )
     &        'INTERNAL ERROR: Number of profiles ' //
     &        'encountered: ', NPROF, CRLF() // BLANK5 //
     &        'Number of profiles expected: ', NPROFMAX

            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )      

        ENDIF
       
        CALL SORTIC( IPS, INDXA, INPSPA )    ! Sort on INPSPA
        
        DO I = 1, NPROF

            J = INDXA( I )        
            IF ( DIVISA( J ) .EQ. 0 ) THEN
                EFLAG = .TRUE.
                CYCLE
            END IF
        
            J = INDXA( I )
            
            INPRF   ( I ) = INPSPA( J )( 1:SPNLEN3 )
            SPECID  ( I ) = SPCIDA( J )
            MOLEFACT( I ) = TON2GM * FACTRA( J ) / DIVISA( J )
            MASSFACT( I ) = TON2GM * XMFA( J )

        ENDDO
        
        IF( EFLAG ) THEN
            MESG = 'At least one of the divisors was zero.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
 
C........... Sort the species names in order to loop through and count
C            the number of times the name changes. This will give us the
C            number of unique species names, NMSPC 
        
        CALL SORTIC( IPS, INDXA, SPECID )       ! Sort on SPECID
        
        SPCIDTMP = EMCMISS3  ! Initialize temporary 
        NMSPC   = 0
        DO I = 1, NPROF
            
            J = INDXA ( I )
            IF ( SPECID( J ) .NE. SPCIDTMP ) THEN
                NMSPC = NMSPC + 1
                SPCIDTMP = SPECID( J )
            END IF
            
        ENDDO
        
C......... Rewind file

       REWIND( FDEV )
       
       RETURN
       
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010  FORMAT( 10( A, :, I8, :, 1X ) )
        
       END SUBROUTINE RDSPROF                                                                            
