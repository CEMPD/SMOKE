
        SUBROUTINE RDB4FAC_CSV( BEISV, NLINES, NSEF, FDEV, NVEG, VGLIST,
     &                      BIOTYPES,LINDX, LFAC, WNTF, LWGT, FACS  ) 

C***********************************************************************
C  subroutine body starts at line XX 
C
C  DESCRIPTION:
C       Reads in the BEIS3 emissions factors from the BFAC file formatted as a csv file
C       for BELD4 only!
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C       07/13 protoype by G. Pouliot
C       07/14 cleaned up the code and add cycle commands for warning cases
C
C***********************************************************************
C
C***********************************************************************

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions:
 
        INTEGER         STR2INT
        REAL            STR2REAL

        EXTERNAL        STR2INT, STR2REAL

C...........   ARGUMENTS and their descriptions: actually-occurring ASC table


        CHARACTER(11),INTENT (IN):: BEISV   !  version of BEIS
        INTEGER, INTENT (IN)  :: NSEF    !  no. biogenic emission factors
        INTEGER, INTENT (IN)  :: FDEV    !  unit number input file 
        INTEGER, INTENT (IN)  :: NVEG  !  no. veg types
        INTEGER, INTENT (IN)  :: NLINES  !  no. lines in input file
	CHARACTER(16), INTENT(IN) :: VGLIST(NVEG)
	CHARACTER(5), INTENT(IN) :: BIOTYPES(NSEF)	

        INTEGER, INTENT (OUT)        :: LINDX( NVEG )      ! leaf area index
        REAL, INTENT (OUT)           :: LFAC( NVEG )       ! leaf biomass
        REAL, INTENT (OUT)           :: WNTF( NVEG )       ! winter factor
        REAL, INTENT (OUT)           :: LWGT( NVEG )       ! specific leaf wgt
        REAL, INTENT (OUT)           :: FACS( NVEG, NSEF ) ! emis facs

        LOGICAL,SAVE :: CHKHDR= .FALSE.   ! check the header BEISFAC4BELD5 
        LOGICAL      :: EFLAG = .FALSE.  !  error flag
        INTEGER      :: MXSEG            ! # of potential line segments

        INTEGER       I, J               !  counters
        INTEGER       ISTAT              !  iostat error
	INTEGER     :: INDEX, SINDEX

        CHARACTER(16) :: VTYP,VGID,UNIT
	REAL :: VALU
        CHARACTER(50), ALLOCATABLE :: SEGMENT( : )   ! Segments of parsed lines
        CHARACTER(300)  MESG             !  message buffer
        CHARACTER(300)  LINE             !  buffer for variables

        CHARACTER(16) :: PROGNAME = 'RDB4FAC' ! program name

C***********************************************************************
C   begin body of subroutine RDB4FAC_CSV

C.........  Set number of potential line segments
        MXSEG = 4
 
        ALLOCATE( SEGMENT( MXSEG ), STAT=ISTAT )
        CALL CHECKMEM( ISTAT, 'SEGMENT', PROGNAME )

C.......... Read in emissions factors for each veg id
        DO I = 1, NLINES

          READ( FDEV, 93030, IOSTAT=ISTAT )  LINE

          IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Error', ISTAT,
     &              'reading EMISSION FACTOR file at line', I
               CALL M3MESG( MESG )
          END IF

C............   Check the header
          IF( BEISV == 'NORMBEIS360' ) THEN   ! BEIS v3.6 does not need a header from B360FAC file
              CHKHDR = .TRUE.
          ELSE
              IF( LINE( 2:14 ) == 'BEISFAC4BELD5' ) THEN
                  CHKHDR = .TRUE.
                  CYCLE
              END IF
          END IF

          IF( I > 1 .AND. .NOT. CHKHDR ) THEN
              MESG = 'ERROR: #BEISFAC4BELD5 header is missing ' //
     &               'from BEISFAC input file'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
          END IF

C.............  Separate the line of data into each part
          CALL PARSLINE( LINE, MXSEG, SEGMENT )

          VGID = SEGMENT( 1 )
	  VTYP = SEGMENT( 2 )
	  UNIT = SEGMENT( 3 )
	  VALU = STR2REAL( SEGMENT ( 4 )) 
	  
	  INDEX = 0
	  DO J = 1, NVEG
	      IF (TRIM(VGID) .eq. TRIM(VGLIST(J)) ) THEN
	         INDEX = J

	      ENDIF
	  ENDDO
	  IF (INDEX .eq. 0) THEN
                MESG = 'WARNING: VEGETATION NAME: '//TRIM(VGID)//
     &          'found in emission factor file but not in land use file'		
               CALL M3MESG( MESG )
	       CYCLE	     
	  ENDIF

          SELECT CASE (TRIM(VTYP))

             CASE ('LAI')
	        LINDX(INDEX) = VALU
             CASE ('LFBIO')
	        LFAC(INDEX) = VALU
				
             CASE ('WFAC')
	        WNTF(INDEX) = VALU
             CASE ('SLW')
	        LWGT(INDEX) = VALU
             CASE DEFAULT
	        SINDEX = 0
	        DO J = 1, NSEF

	           IF (TRIM(VTYP) .eq. TRIM(BIOTYPES(J)) ) THEN
	              SINDEX = J
	           ENDIF
	        ENDDO

		
	        IF (SINDEX .eq. 0) THEN
                   EFLAG = .TRUE.
                   MESG = 'Error: SPECIES NAME: '//TRIM(VTYP)//
     &              'not found'
                   CALL M3MESG( MESG )	  
		   CYCLE   
	        ENDIF

  	         	    
		FACS(INDEX,SINDEX) = VALU

				
	  END SELECT  

        ENDDO

        IF( EFLAG ) THEN
            MESG = 'Problem reading biogenic emissions factors file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        DEALLOCATE(SEGMENT)

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( 8X, A16, A )
93020   FORMAT( A16, A )
93030   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )
94020   FORMAT( A )

        END SUBROUTINE RDB4FAC_CSV
