
        SUBROUTINE WREPINVEN( ADEV, CDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine writes out the REPINVEN file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 11/2002 by A. Holland
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel gsions (SMOKE) Modeling
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
C...........   This module is the inventory arrays
        USE MODSOURC
        
C.........  This module contains the information about the source category
        USE MODINFO
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS
        
C.........  This module contains the arrays for the area-to-point x-form
        USE MODAR2PT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        
C...........   EXTERNAL FUNCTIONS and their descriptions:

	INTEGER		FIND1
	INTEGER		INDEX1
        
        EXTERNAL	FIND1, INDEX1

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT (IN) :: ADEV  ! file unit no. for REPINVEN file
        INTEGER     , INTENT (IN) :: CDEV

C...........   Local variables

	CHARACTER*1	KEEP
        CHARACTER(LEN=DDSLEN3)	DESC
        CHARACTER(LEN=SCCLEN3)  SCC
        CHARACTER(LEN=IOVLEN3)  DNAME
        CHARACTER(LEN=SDSLEN3)	SCCDC
        CHARACTER(LEN=SCCLEN3)  PSCC
        
        INTEGER		I, J, K, L, S, IOS
        INTEGER		STATE
        INTEGER		NFIPS
        INTEGER		POLL
        INTEGER		PFIP
        
        INTEGER, ALLOCATABLE :: ASSIGNED( : )
        INTEGER, ALLOCATABLE :: UNASSIGN( : )
        
        REAL		VALCHECK
        REAL		DIFF
        REAL		OEMIS
        REAL		SEMIS

C...........   Other local variables

	CHARACTER*300	MESG
        CHARACTER(LEN=SDSLEN3)	CBUF

        CHARACTER*16  :: PROGNAME = 'WREPINVEN' ! program name

C***********************************************************************
C   begin body of subroutine WREPINVEN

C.........  Write out first report to REPINVEN file

        WRITE( ADEV, 93010 ) 'CAS Code', 'Keep', 'Nrecs', 
     &         'Emissions', 'CAS Description'
     
        WRITE( ADEV, 93020 ) '[tons/year]'
        
        WRITE( ADEV, 93000 ) REPEAT( '-', 85 )
          
        DO I = 1, NUNIQCAS
          IF( UCASNPOL( I ) .EQ. UCASNKEP( I ) ) THEN
            KEEP = 'Y'
          ELSE IF( UCASNPOL( I ) .NE. UCASNKEP( I ) .AND.
     &             UCASNKEP( I ) .NE. 0 ) THEN
            KEEP = 'P'
          ELSE IF( UCASNKEP( I ) .EQ. 0 ) THEN
            KEEP = 'N'
          END IF
            
          K = INDEX1( UNIQCAS( I ), NINVTBL, ITCASA )
          IF( K .GT. 0 ) THEN
            DESC = ITCASDSCA( K )
          ELSE
            MESG = 'CAS code not found in raw list'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 1 )
          END IF
          
          WRITE( ADEV, 93030 ) UNIQCAS( I ), KEEP, RECSBYCAS( I ),
     &             EMISBYCAS( I ), DESC

	END DO
        
        WRITE( ADEV, 93000 ) REPEAT( '-', 85 )
        WRITE( ADEV, 93000 ) ' '
        WRITE( ADEV, 93000 ) ' '
        
C.........  Write out second report to REPINVEN file

	WRITE( ADEV, 93040 ) 'CAS Code', 'CAS Emissions', 'Factor',
     &         'Data Name', 'Data Emissions', 'Data Description', 
     &         'CAS Description'
     
     	WRITE( ADEV, 93050 ) '[tons/year]', '[tons/year]'
        
        WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
        
        DO I = 1, NINVTBL

          IF( SORTCAS( I ) .EQ. '        ' ) CYCLE
          
          K = INDEX1( SORTCAS( I ), NUNIQCAS, UNIQCAS )
          IF( K .LE. 0 ) THEN
            WRITE( MESG, 94010 )
     &             'Sorted CAS code, ', SORTCAS( I ), ' ,was not '//
     &             'found in list of unique CAS codes.'
	    CALL M3WARN( PROGNAME, 0, 0, MESG )
          END IF
          
          J = SCASIDX( I )
          
          IF( .NOT. ITKEEPA( J ) ) CYCLE
          
          VALCHECK = EMISBYCAS( K ) * ITFACA( J )
          DIFF = VALCHECK - EMISBYPOL( I )
          IF( ABS( DIFF ) .NE. 0.0 ) THEN
            WRITE( MESG, 94020 )
     &         'WARNING: Summed emissions for pollutant, ',
     &         ITNAMA( J ), ', differ from factored emissions '//
     &         'by ', DIFF
            CALL M3MESG( MESG )
          END IF

          WRITE( ADEV, 93060 ) SORTCAS( I ), EMISBYCAS( K ),
     &           ITFACA( J ), ITNAMA( J ), EMISBYPOL( I ),
     &           ITDSCA( J ), ITCASDSCA( J )
     
        END DO

        WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
        WRITE( ADEV, 93000 ) ' '
        WRITE( ADEV, 93000 ) ' '
        

C.........  The next four reports are only created when area
C           to point allocation is occuring
       
        IF( ALLOCATED( REPAR2PT ) ) THEN
        
          CALL RDSCCDSC( CDEV )
          
C.........  Write out third report to REPINVEN file

	  WRITE( ADEV, 93000 ) ' SCCs in area-to-point factors '//
     &           'file not found in the inventory:'
          WRITE( ADEV, 93000 ) ' '
          
          DO I = 1, NA2PSCC
          
            K = INDEX1( A2PSCC( I ), NINVSCC, INVSCC )
            IF( K .LE. 0 ) THEN
              WRITE( ADEV, 93130 ) A2PSCC( I )
            END IF
            
          END DO
          
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
          WRITE( ADEV, 93000 ) ' '
          WRITE( ADEV, 93000 ) ' '
        
C.........  Write out fourth report to REPINVEN file

	  WRITE( ADEV, 93070 ) 'SCC Code', 'Data Name',
     &         'FIPS count', 'Emissions before', 'Emissions after',
     &         'SCC Description'
     
     	  WRITE( ADEV, 93080 ) '[tons/year]', '[tons/year]'
        
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )

          DO I = 1, NCONDSRC
        
            STATE = REPAR2PT( I )%STATE
            SCC   = REPAR2PT( I )%SCC
            POLL  = REPAR2PT( I )%POLL
            DNAME = INVDNAM( REPAR2PT( I )%POLL )
            NFIPS = REPAR2PT( I )%NFIPS
            OEMIS = REPAR2PT( I )%ORIGEMIS
            SEMIS = REPAR2PT( I )%SUMEMIS
           
            K = INDEX1( SCC, NINVSCC, INVSCC )
            IF( K .LE. 0 ) THEN
              WRITE( MESG, 94010 )
     &             'SCC code, ', SCC, ' ,was not '//
     &             'found in master list of SCC codes.'
	      CALL M3WARN( PROGNAME, 0, 0, MESG )
            END IF
            
            CBUF = SCCDESC( K )
            L = LEN_TRIM( CBUF )
            SCCDC = CBUF( 1:L )
            
            DO J = 1, NCONDSRC
            
              IF( REPAR2PT( J )%STATE .EQ. STATE ) CYCLE
              
              IF( REPAR2PT( J )%SCC .EQ. SCC .AND.
     &            REPAR2PT( J )%POLL .EQ. POLL ) THEN
                NFIPS = NFIPS + REPAR2PT( J )%NFIPS
                OEMIS = OEMIS + REPAR2PT( J )%ORIGEMIS
                SEMIS = SEMIS + REPAR2PT( J )%SUMEMIS
              END IF
              
            END DO
            
            WRITE( ADEV, 93090 ) SCC, DNAME, NFIPS,
     &            OEMIS, SEMIS, SCCDC
            
          END DO
          
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
          WRITE( ADEV, 93000 ) ' '
          WRITE( ADEV, 93000 ) ' '
        
C.........  Write out fifth report to REPINVEN file

	  WRITE( ADEV, 93100 ) 'State', 'SCC Code', 'Data Name',
     &         'FIPS count', 'Emissions before', 'Emissions after',
     &         'SCC Description'
     
          WRITE( ADEV, 93110 ) '[tons/year]', '[tons/year]'
        
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )

          DO I = 1, NCONDSRC
        
            STATE = REPAR2PT( I )%STATE
            SCC   = REPAR2PT( I )%SCC
            DNAME = INVDNAM( REPAR2PT( I )%POLL )
            NFIPS = REPAR2PT( I )%NFIPS
            OEMIS = REPAR2PT( I )%ORIGEMIS
            SEMIS = REPAR2PT( I )%SUMEMIS
           
            K = INDEX1( SCC, NINVSCC, INVSCC )
            IF( K .LE. 0 ) THEN
              WRITE( MESG, 94010 )
     &             'SCC code, ', SCC, ' ,was not '//
     &             'found in master list of SCC codes.'
	      CALL M3WARN( PROGNAME, 0, 0, MESG )
            END IF
            
            CBUF = SCCDESC( K )
            L = LEN_TRIM( CBUF )
            SCCDC = CBUF( 1:L )
           
            WRITE( ADEV, 93120 ) STATE, SCC, DNAME, NFIPS,
     &            OEMIS, SEMIS, SCCDC
     
          END DO
          
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
          WRITE( ADEV, 93000 ) ' '
          WRITE( ADEV, 93000 ) ' ' 
        
        
C.........  Write out sixth report to REPINVEN file

	  WRITE( ADEV, 93140 ) 'FIPS count'
	  WRITE( ADEV, 93150 ) 'SCC Code', 'Assigned',
     &         'Unassigned', 'SCC Description'
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
          
          ALLOCATE( ASSIGNED( NINVSCC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'ASSIGNED', PROGNAME )
          ALLOCATE( UNASSIGN( NINVSCC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'UNASSIGN', PROGNAME )
          ASSIGNED = 0
          UNASSIGN = 0
          
          PFIP = 99999
          PSCC = '9999999999'
          
          DO S = 1, NSRC
          
            IF( IFIP( S ) .EQ. PFIP .AND. CSCC( S ) .EQ. PSCC ) CYCLE
            
            J = INDEX1( CSCC( S ), NA2PSCC, A2PSCC )
            IF( J .GT. 0 ) THEN
            
              K = INDEX1( CSCC( S ), NINVSCC, INVSCC )
              IF( XLOCA( S ) .GT. AMISS3 .AND. 
     &             YLOCA( S ) .GT. AMISS3 ) THEN
                ASSIGNED( K ) = ASSIGNED( K ) + 1
              ELSE
                UNASSIGN( K ) = UNASSIGN( K ) + 1
              END IF
              
            END IF
            
            PFIP = IFIP( S )
            PSCC = CSCC( S )
            
          END DO
          
          DO I = 1, NINVSCC
          
            WRITE( ADEV, 93160 ) INVSCC( I ), ASSIGNED( I ),
     &             UNASSIGN( I ), SCCDESC( I )
     
          END DO
          
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
     
     
        END IF
          
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( 1X, A8, 4X, A4, 4X, A5, 6X, A9, 9X, A15 )

93020   FORMAT( 32X, A11 )

93030   FORMAT( 1X, A8, 4X, A1, 2X, I10, 4X, F16.10, 4X, A40 )

93040	FORMAT( 1X, A8, 6X, A13, 4X, A6, 2X, A9, 13X, A14, 4X,
     &          A16, 28X, A15 )
     
93050	FORMAT( 15X, A11, 36X, A11 )

93060	FORMAT( 1X, A8, 4X, F16.10, 4X, F3.1, 4X, A16, 4X, F16.10,
     &          4X, A40, 4X, A40 )
     
93070	FORMAT( 1X, A8, 5X, A9, 6X, A10, 4X, A16, 4X, A15,
     &          4X, A16 )
     
93080	FORMAT( 43X, A11, 9X, A11 )

93090	FORMAT( 1X, A10, 4X, A16,4X, I3, 4X, F16.10, 4X,
     &          F16.10, 4X, A )
     
93100	FORMAT( 1X, A5, 2X, A8, 5X, A9, 6X, A10, 4X, A16, 4X, A15,
     &          4X, A16 )
     
93110	FORMAT( 50X, A11, 9X, A11 )

93120	FORMAT( 1X, I2.2, 4X, A10, 4X, A16,4X, I3, 4X, F16.10, 4X,
     &          F16.10, 4X, A )
     
93130	FORMAT( 1X, A10 )

93140	FORMAT( 20X, A10 )

93150	FORMAT( 1X, A8, 4X, A8, 2X, A10, 4X, A15 )

93160	FORMAT( 1X, A10, I10, 2X, I10, 4X, A )
     
94010	FORMAT( 10( A, :, A8, :, 1X ) )

94020	FORMAT( A, :, A16, :, A, :, F16.10 )



        END SUBROUTINE WREPINVEN
