
        SUBROUTINE RDSTKEMSPT( FDEV, CFLAG, WFLAG, UTMZONE,
     &                         NCNTY, FIPTOCSRC )

C***********************************************************************
C  subroutine body starts at line 232
C
C  DESCRIPTION:
C      This subroutine reads the EMS-95 point stack file.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Copied from rdemspt.f by C. Seppanen (2/03)
C
C****************************************************************************
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

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC
        
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: CSOURC, STKHT, STKDM, STKTK, STKVE, 
     &                      XLOCA, YLOCA

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'CONST3.EXT'    !  physical constants
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        INTEGER         FINDCFIRST
        INTEGER         FIND1
        INTEGER         STR2INT
        REAL            STR2REAL
        
        EXTERNAL        CRLF, FINDCFIRST, FIND1, STR2INT, STR2REAL

C...........   SUBROUTINE ARGUMENTS
        INTEGER,          INTENT (IN) :: FDEV   ! unit no. of inv file
        LOGICAL,          INTENT (IN) :: CFLAG  ! true: recalc vel w/ flow and diam
        LOGICAL,          INTENT (IN) :: WFLAG  ! true: convert lat-lons to western hemisphere
        INTEGER,          INTENT (IN) :: UTMZONE( NSRC )   ! UTM zone by source
        INTEGER,          INTENT (IN) :: NCNTY  ! no. counties in inventory
        INTEGER,          INTENT (IN) :: FIPTOCSRC( NCNTY+1,2 )  ! index into CSOURC

C...........   Other local variables
        INTEGER              I, K1                  ! counters and indices
        
        INTEGER              ACTCNT                 ! actual count of processed records
        INTEGER              IOS                    ! I/O status
        INTEGER              IREC                   ! no. records in file
        INTEGER              IFIP                   ! integer FIPS code
        INTEGER              STIDX, ENDIDX          ! start and end idx into CSOURC
        INTEGER              KEYLEN                 ! length of source key
        INTEGER              ZONE                   ! UTM zone
        
        REAL                 FSAV                   ! temporary flow value
        REAL                 STKD                   ! stack diameter
        REAL                 STKH                   ! stack height
        REAL                 STKT                   ! stack temperature
        REAL                 STKV                   ! stack velocity
        REAL                 STKF                   ! stack flow
        REAL                 XVAL                   ! UTM X coordinate
        REAL                 XX                     ! lat-lon X coordinate
        REAL                 YVAL                   ! UTM Y coordinate
        REAL                 YY                     ! lat-lon Y coordinate
        
        CHARACTER(LEN=FIPLEN3) CFIP      ! fip code
        CHARACTER(LEN=PLTLEN3) FCID      ! facility ID
        CHARACTER(LEN=CHRLEN3) SKID      ! stack ID
        
        CHARACTER(LEN=ALLLEN3) SRCKEY               ! source key
        CHARACTER(LEN=300)     LINE                 ! line from file
        CHARACTER(LEN=300)     MESG                 ! message buffer

        CHARACTER(LEN=16) :: PROGNAME = 'RDSTKEMSPT' ! program name

C***********************************************************************
C   begin body of subroutine RDSTKEMSPT

        ACTCNT = 0
        IREC = 0

C.........  Loop through lines in file
        DO
    
            READ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1
            
C.............  Check for I/O errors
            IF( IOS > 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,
     &              'reading stack file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
    
C.............  Check for end of file
            IF( IOS < 0 ) EXIT
    
C.............  Read FIPS code from file
            CFIP( 1:1 ) = '0'
            CFIP( 2:3 ) = LINE( 1:2 )
            CFIP( 4:6 ) = LINE( 3:5 )
                        
C.............  Replace blanks with zeros        
            DO I = 1,FIPLEN3
                IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
            END DO
            
C.............  Find county in FIPTOCSRC array
            IFIP = STR2INT( CFIP( 1:6 ) )
            K1 = FIND1( IFIP, NCNTY+1, FIPTOCSRC( :,1 ) )

C.............  Make sure county is in the inventory
            IF( K1 <= 0 ) CYCLE
            
C.............  Set starting and ending indices            
            STIDX = FIPTOCSRC( K1,2 )
            ENDIDX = FIPTOCSRC( K1+1,2 ) - 1
            
C.............  Build source key
            FCID = ADJUSTL( LINE(  6:20 ) )
            SKID = ADJUSTL( LINE( 21:32 ) )
            
            CALL BLDCSRC( CFIP, FCID, SKID, CHRBLNK3, CHRBLNK3,
     &                    CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                    SRCKEY )

            KEYLEN = LEN_TRIM( SRCKEY )
            
C.............  Find source key in CSOURC array
            K1 = FINDCFIRST( SRCKEY, ENDIDX-STIDX+1, 
     &                       CSOURC( STIDX:ENDIDX )( 1:KEYLEN ) )
    
C.............  If key not found, go to next line
            IF( K1 <= 0 ) CYCLE

C.............  Shift index to account for only searching part of CSOURC array
            K1 = K1 + STIDX - 1
                
C.............  Read and check stack parameters
            STKD = STR2REAL( LINE( 33:40 ) )
            STKH = STR2REAL( LINE( 41:47 ) )
            STKT = STR2REAL( LINE( 48:54 ) )
            STKV = STR2REAL( LINE( 55:61 ) )
            STKF = STR2REAL( LINE( 62:71 ) )
            
            IF( STKD > 0.0 ) STKD = STKD * FT2M   ! diameter - ft to m
            IF( STKH > 0.0 ) STKH = STKH * FT2M   ! height - ft to m
            IF( STKT > 0.0 ) STKT = (STKT-32.)*FTOC + CTOK  ! temperature - F to K
            IF( STKV > 0.0 ) STKV = STKV * FT2M   ! velocity - ft/s to m/s
            IF( STKF > 0.0 ) STKF = STKF * FLWE2M ! flow - ft^3/min to m^3/s
    
C.............  Calculate velocity from flow and diameter if needed
            IF( CFLAG        .OR.
     &        ( STKV <= 0.0  .AND.
     &          STKD >  0.0  .AND.
     &          STKF >  0.0 )      ) THEN
                
                IF( STKD > 0.0 ) THEN
                    STKV = STKF / ( 0.25 * PI * STKD * STKD )
                ELSE
                    WRITE( MESG,94010 ) 'Zero or negative stack ' //
     &                     'diameter at line', IREC, 'in stack ' //
     &                     'file.' // CRLF() // BLANK10 // 'Using ' //
     &                     'velocity from file rather than ' //
     &                     'recalculating from stack flow'
                    CALL M3MESG( MESG )
                END IF
    
C.............  Compare flow to velocity and diameter. Set to exact flow input
C               value only if it is consistent with velocity and diameter
            ELSE IF( STKF > 0 ) THEN
                FSAV = STKV * 0.25 * PI * STKD * STKD
                
                IF( ( STKF - FSAV ) / STKF > 0.001 ) THEN
                    FSAV = STKF
                END IF
                
            END IF
    
C.............  Read and check coordinates
            XVAL = STR2REAL( LINE( 72:80 ) )
            YVAL = STR2REAL( LINE( 81:89 ) )
    
C.............  If invalid values, pull previously stored values from facility file
            IF( XVAL <= 0.0 .OR. YVAL <= 0.0 ) THEN
                XVAL = XLOCA( K1 )
                YVAL = YLOCA( K1 )
    
C.................  Still no values, skip to next line
                IF( XVAL <= 0.0 .OR. YVAL <= 0.0 ) CYCLE
            END IF
    
C.............  Convert coordinates from UTM to lat-lon
            ZONE = UTMZONE( K1 )
            
            IF( ZONE > 0 ) THEN
                CALL UTM2LL( XVAL, YVAL, ZONE, XX, YY )
            
            ELSE
C.................  Make sure coordinates are within lat-lon range
                IF( ABS( XVAL ) > 180. .OR. ABS( YVAL ) > 180. ) THEN
                    WRITE( MESG,94010 ) 'Invalid (X,Y) coordinates ' //
     &                     'at line', IREC, 'in facility file'
                    CALL M3MESG( MESG )
                    CYCLE
                ELSE
                    XX = XVAL
                    YY = YVAL
                    IF( WFLAG .AND. XX > 0 ) XX = -XX  ! convert to western hemisphere
                END IF
            END IF
            
C.............  Increment count of actual records            
            ACTCNT = ACTCNT + 1
            
C.............  Loop through matching sources and set values
            DO
            
C.................  Check for end of array or end of matching sources
                IF( K1 > NSRC ) EXIT
                IF( CSOURC( K1 )( 1:KEYLEN ) /= SRCKEY ) EXIT
                
                STKHT( K1 ) = STKH
                STKDM( K1 ) = STKD
                STKTK( K1 ) = STKT
                STKVE( K1 ) = STKV
                XLOCA( K1 ) = XX
                YLOCA( K1 ) = YY
                
                K1 = K1 + 1
            END DO
    
        END DO

C.........  Write count of total records processed    
        WRITE( MESG,94010 )
     &         'STACK FILE processed' // CRLF() // BLANK10 //
     &         'Actual STACK record-count ', ACTCNT
        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94125   FORMAT( I5 )

94300   FORMAT( A, I2.2, A, I2.2, A )

        END SUBROUTINE RDSTKEMSPT
