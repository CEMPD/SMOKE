
        SUBROUTINE WREPINVEN( ADEV, CDEV )

C***********************************************************************
C  subroutine body starts at line 108
C
C  DESCRIPTION:
C      This subroutine writes out the REPINVEN file, which contains
C      reports on the toxics emissions inventory and the area-to-point
C      sources.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      I/O API functions
C      RDSCCDESC
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C.........  MODULES for public variables
C...........   This module is the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CIFIP, CSCC, XLOCA, YLOCA
        
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NIPOL
        
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: NUNIQCAS, NINVTBL, ITCASA, UNIQCAS,
     &                      NINVSCC, INVSCC, UCASNPOL, UCASNKEP,
     &                      ITCASDSCA, RECSBYCAS, EMISBYCAS,
     &                      SORTCAS, SCASIDX, ITKEEPA, ITFACA,
     &                      EMISBYPOL, ITNAMA, INVDNAM, ITDSCA,
     &                      SCCDESC
        
C.........  This module contains the arrays for the area-to-point x-form
        USE MODAR2PT, ONLY: REPAR2PT, NA2PSCC, A2PSCC, NCONDSRC

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        
C...........   EXTERNAL FUNCTIONS and their descriptions:

C       INTEGER         INDEX1
        
C        EXTERNAL        INDEX1

C...........   SUBROUTINE ARGUMENTS

        INTEGER     , INTENT (IN) :: ADEV  ! file unit no. for REPINVEN file
        INTEGER     , INTENT (IN) :: CDEV  ! file unit no. for SCC desc. file

C...........   Local variables

        CHARACTER       KEEP               ! determines if CAS number is kept
        CHARACTER(36)   TSCCPOLL           ! temp. SCC/pollutant combination
        CHARACTER(DDSLEN3)  DESC       ! temp. SCC description
        CHARACTER(SCCLEN3)  SCC        ! temp. SCC code
        CHARACTER(IOVLEN3)  DNAME      ! temp. data name
        CHARACTER(SDSLEN3)  SCCDC      ! temp. SCC description
        CHARACTER(SCCLEN3)  PSCC       ! previous SCC code
        CHARACTER(SCCLEN3)  TSCC       ! temp. SCC code
        CHARACTER(IOVLEN3)  TDNAME     ! temp. data name
        CHARACTER(FIPLEN3)  PFIP       ! previous FIPS code

        CHARACTER(36), ALLOCATABLE :: SCCPOLL( : ) ! SCC/pollutant combination
        
        INTEGER         I, J, K, L, S, IOS ! counters and indicies
        INTEGER         STATE              ! temp. state code
        INTEGER         NFIPS              ! temp. number of FIPS codes
        INTEGER         NSCCPOLL           ! total number of SCC // pollutant combinations
        INTEGER         POLL               ! temp. pollutant
        
        INTEGER, ALLOCATABLE :: ASSIGNED( : ) ! number of FIPS codes assigned
        INTEGER, ALLOCATABLE :: UNASSIGN( : ) ! number of FIPS codes unassigned
        
        REAL            VALCHECK           ! temp. value check
        REAL            DIFF               ! temp. difference
        REAL            PDFF               ! temp. percent difference
        REAL            OEMIS              ! temp. original emissions
        REAL            SEMIS              ! temp. summed emissions

C...........   Other local variables

        CHARACTER(300)  MESG               ! temp. message buffer
        CHARACTER(SDSLEN3)  CBUF       ! temp. buffer

        CHARACTER(16) :: PROGNAME = 'WREPINVEN' ! program name

C***********************************************************************
C   begin body of subroutine WREPINVEN

C.........  Write out first report to REPINVEN file
C           This report summarizies by CAS code the emissions,
C           the number of inventory records, whether all, some, or none
C           of the pollutants are kept and the CAS description.

C............  Write out header
        WRITE( ADEV, 93010 ) 'CAS Code', 'Keep', 'Nrecs', 
     &         'Emissions', 'CAS Description'
     
        WRITE( ADEV, 93020 ) '[tons/year]'
        
        WRITE( ADEV, 93000 ) REPEAT( '-', 85 )
        
C............  Determine whether all, some, or none of the pollutants
C              are kept            
        DO I = 1, NUNIQCAS
          IF( UCASNPOL( I ) .EQ. UCASNKEP( I ) ) THEN
            KEEP = 'Y'
          ELSE IF( UCASNPOL( I ) .NE. UCASNKEP( I ) .AND.
     &             UCASNKEP( I ) .NE. 0 ) THEN
            KEEP = 'P'
          ELSE IF( UCASNKEP( I ) .EQ. 0 ) THEN
            KEEP = 'N'
          END IF
          
C............  Find CAS code in the original list to get the description            
          K = INDEX1( UNIQCAS( I ), NINVTBL, ITCASA )
          IF( K .GT. 0 ) THEN
            DESC = ITCASDSCA( K )
          ELSE
            MESG = 'CAS code not found in raw list'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 1 )
          END IF
          
C............  Write out report data fields          
          WRITE( ADEV, 93030 ) UNIQCAS( I ), KEEP, RECSBYCAS( I ),
     &             EMISBYCAS( I ), DESC

        END DO
        
        WRITE( ADEV, 93000 ) REPEAT( '-', 85 )
        WRITE( ADEV, 93000 ) ' '
        WRITE( ADEV, 93000 ) ' '
        
C.........  Write out second report to REPINVEN file
C           The report contains emissions before and after disaggregation.
C           Only pollutants that have KEEP=Y will be reported.

C............  Write out header
        WRITE( ADEV, 93040 ) 'CAS Code', 'CAS Emissions', 'Factor',
     &         'Data Name', 'Data Emissions', 'Data Description', 
     &         'CAS Description'
     
        WRITE( ADEV, 93050 ) '[tons/year]', '[tons/year]'
        
        WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
        
        DO I = 1, NINVTBL

C............  If CAS code is blank, skip
          IF( SORTCAS( I ) .EQ. '        ' ) CYCLE
          
C............  Find sorted CAS code in unique list
          K = INDEX1( SORTCAS( I ), NUNIQCAS, UNIQCAS )
          IF( K .LE. 0 ) THEN
            WRITE( MESG, 94010 )
     &             'Sorted CAS code, ', SORTCAS( I ), ' ,was not '//
     &             'found in list of unique CAS codes.'
            CALL M3WARN( PROGNAME, 0, 0, MESG )
          END IF
          
          J = SCASIDX( I )
          
C............  If pollutant is not kept then skip
          IF( .NOT. ITKEEPA( J ) ) CYCLE
          
C............  Check value of factored emissions
          VALCHECK = EMISBYCAS( K ) * ITFACA( J )
          DIFF = VALCHECK - EMISBYPOL( I )
          PDFF = ( ABS( DIFF ) / EMISBYPOL( I ) ) * 100
          IF( PDFF >  0.1  ) THEN
            WRITE( MESG, 94020 )
     &         'WARNING: Summed emissions of ', EMISBYPOL( I ),
     &         ' for pollutant, ', ITNAMA( J ),
     &         ', differ from component total '//
     &         'by ', DIFF
            CALL M3MESG( MESG )
          END IF

C............  Write out report data fields
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
        
C.........  Read in SCC description file
          CALL RDSCCDSC( CDEV )
          
C.........  Write out third report to REPINVEN file
C           This report lists the SCCs that have area-to-point source
C           factor file assignments but are not in the inventory.

C............  Write out header
          WRITE( ADEV, 93000 ) ' SCCs in area-to-point factors '//
     &           'file not found in the inventory:'
          WRITE( ADEV, 93000 ) ' '
          
          DO I = 1, NA2PSCC
          
C............  Find area-to-point SCC in inventory SCC list.  If
C              not found, write SCC to REPINVEN file
            K = INDEX1( A2PSCC( I ), NINVSCC, INVSCC )
            IF( K .LE. 0 ) THEN
              WRITE( ADEV, 93130 ) A2PSCC( I )
            END IF
            
          END DO
          
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
          WRITE( ADEV, 93000 ) ' '
          WRITE( ADEV, 93000 ) ' '
        
C.........  Write out fourth report to REPINVEN file
C           Reports the pollutant name, emissions total before and 
C           after factors are applied, the total number of FIPS codes
C           affected by SCC and the SCC description for area emissions
C           assigned to point sources.

C............  Write out header
          WRITE( ADEV, 93070 ) 'SCC Code', 'Data Name',
     &         'FIPS count', 'Emissions before', 'Emissions after',
     &         'SCC Description'
     
          WRITE( ADEV, 93080 ) '[tons/year]', '[tons/year]'
        
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )

          ALLOCATE( SCCPOLL( NCONDSRC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'SCCPOLL', PROGNAME )
          SCCPOLL = ''
          NSCCPOLL = 0

          DO I = 1, NCONDSRC

            STATE = REPAR2PT( I )%STATE
            SCC   = REPAR2PT( I )%SCC
            POLL  = REPAR2PT( I )%POLL
            DNAME = INVDNAM( REPAR2PT( I )%POLL )
            NFIPS = REPAR2PT( I )%NFIPS
            OEMIS = REPAR2PT( I )%ORIGEMIS
            SEMIS = REPAR2PT( I )%SUMEMIS

            TSCCPOLL = SCC//DNAME
            K = INDEX1( TSCCPOLL, NSCCPOLL, SCCPOLL )
            IF( K .LE. 0 ) THEN
              SCCPOLL( I ) = TSCCPOLL
              NSCCPOLL = I
            ELSE
              CYCLE
            END IF
           
C............  Find SCC in master list of SCC codes
            K = INDEX1( SCC, NINVSCC, INVSCC )
            IF( K .LE. 0 ) THEN
              WRITE( MESG, 94010 )
     &             'SCC code, ', SCC, ' ,was not '//
     &             'found in master list of SCC codes.'
              CALL M3WARN( PROGNAME, 0, 0, MESG )
              CBUF = 'Unknown SCC'
            ELSE
              CBUF = SCCDESC( K )
            END IF
            
            L = LEN_TRIM( CBUF )
            SCCDC = CBUF( 1:L )
            
            DO J = 1, NCONDSRC
            
              IF( REPAR2PT( J )%STATE .EQ. STATE ) CYCLE
        
C............  Sum emissions and count up FIPS codes by SCC      
              TSCC = REPAR2PT( J )%SCC
              TDNAME = INVDNAM( REPAR2PT( J )%POLL )
              TSCCPOLL = TSCC//TDNAME
              K = INDEX1( TSCCPOLL, NA2PSCC * NIPOL, SCCPOLL )
              IF( K .NE. 0 ) THEN
                NFIPS = NFIPS + REPAR2PT( J )%NFIPS
                OEMIS = OEMIS + REPAR2PT( J )%ORIGEMIS
                SEMIS = SEMIS + REPAR2PT( J )%SUMEMIS
              END IF
              
            END DO
            
C............  Write out report data fields
            WRITE( ADEV, 93090 ) SCC, DNAME, NFIPS,
     &            OEMIS, SEMIS, SCCDC
            
          END DO
          
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
          WRITE( ADEV, 93000 ) ' '
          WRITE( ADEV, 93000 ) ' '
        
C.........  Write out fifth report to REPINVEN file
C           This report is similar to the above one separated by state.

C............  Write out header
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

C............  Find SCC in master list of SCC codes        
            K = INDEX1( SCC, NINVSCC, INVSCC )
            IF( K .LE. 0 ) THEN
              WRITE( MESG, 94010 )
     &             'SCC code, ', SCC, ' ,was not '//
     &             'found in master list of SCC codes.'
              CALL M3WARN( PROGNAME, 0, 0, MESG )
              CYCLE
            END IF
            
            CBUF = SCCDESC( K )
            L = LEN_TRIM( CBUF )
            SCCDC = CBUF( 1:L )
           
C............  Write out report data fields
            WRITE( ADEV, 93120 ) STATE, SCC, DNAME, NFIPS,
     &            OEMIS, SEMIS, SCCDC
     
          END DO
          
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
          WRITE( ADEV, 93000 ) ' '
          WRITE( ADEV, 93000 ) ' ' 
        
        
C.........  Write out sixth report to REPINVEN file
C           Reports for each SCC, the number of FIPS codes getting assigned
C           or unassigned to point locations.

C............  Write out header
          WRITE( ADEV, 93140 ) 'FIPS count'
          WRITE( ADEV, 93150 ) 'SCC Code', 'Assigned',
     &         'Unassigned', 'SCC Description'
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
          
C............  Allocate and initialize arrays
          ALLOCATE( ASSIGNED( NINVSCC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'ASSIGNED', PROGNAME )
          ALLOCATE( UNASSIGN( NINVSCC ), STAT=IOS )
          CALL CHECKMEM( IOS, 'UNASSIGN', PROGNAME )
          ASSIGNED = 0
          UNASSIGN = 0
          
C............  Initialize previous FIPS and SCC codes
          PFIP = ' '
          PSCC = REPEAT( '9', SCCLEN3 )
          
C............  Loop through sources
          DO S = 1, NSRC
        
C............  If FIPS and SCC codes are equal to previous, cycle  
            IF( CIFIP( S ) .EQ. PFIP .AND. CSCC( S ) .EQ. PSCC ) CYCLE

C............  If SCC is in the area-to-point SCC list then determine
C              if it is assigned or unassigned     
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
     
C............  Reset previous FIPS and SCC codes to current ones       
            PFIP = CIFIP( S )
            PSCC = CSCC( S )
            
          END DO
          
C............  Write out report data fields
          DO I = 1, NINVSCC
          
            WRITE( ADEV, 93160 ) INVSCC( I ), ASSIGNED( I ),
     &             UNASSIGN( I ), SCCDESC( I )
     
          END DO
          
          WRITE( ADEV, 93000 ) REPEAT( '-', 150 )
     
C.........  Deallocate local memory
    	  IF( ALLOCATED( ASSIGNED ) ) DEALLOCATE( ASSIGNED )  
    	  IF( ALLOCATED( UNASSIGN ) ) DEALLOCATE( UNASSIGN )  
     
        END IF
          
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( 1X, A8, 12X, A4, 4X, A5, 6X, A9, 9X, A15 )

93020   FORMAT( 32X, A11 )

93030   FORMAT( 1X, A16, 4X, A1, 2X, I10, 4X, E16.10, 4X, A40 )

93040   FORMAT( 1X, A8, 14X, A13, 4X, A6, 2X, A9, 13X, A14, 4X,
     &          A16, 28X, A15 )
     
93050   FORMAT( 15X, A11, 36X, A11 )

93060   FORMAT( 1X, A16, 4X, E16.10, 3X, F6.4, 2X, A16, 4X, E16.10,
     &          4X, A40, 4X, A40 )
     
93070   FORMAT( 1X, A8, 15X, A9, 6X, A10, 4X, A16, 4X, A15,
     &          4X, A16 )
     
93080   FORMAT( 53X, A11, 9X, A11 )

93090   FORMAT( 1X, A20, 4X, A16,4X, I6, 4X, E16.10, 4X,
     &          E16.10, 4X, A )
     
93100   FORMAT( 1X, A5, 2X, A8, 15X, A9, 6X, A10, 4X, A16, 4X, A15,
     &          4X, A16 )
     
93110   FORMAT( 60X, A11, 9X, A11 )

93120   FORMAT( 1X, I2.2, 4X, A20, 4X, A16,4X, I6, 4X, E16.10, 4X,
     &          E16.10, 4X, A )
     
93130   FORMAT( 1X, A20 )

93140   FORMAT( 30X, A10 )

93150   FORMAT( 1X, A8, 14X, A8, 2X, A10, 4X, A15 )

93160   FORMAT( 1X, A20, I10, 2X, I10, 4X, A )
     
94010   FORMAT( 10( A, :, A16, :, 1X ) )

94020   FORMAT( A, :, E10.3, :, A, :, A16, :, A, :, E10.3 )



        END SUBROUTINE WREPINVEN
