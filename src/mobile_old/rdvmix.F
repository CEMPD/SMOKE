
C Version "@(#)$Id$ $Source$ $Date$ 


        SUBROUTINE RDVMIX( VDEV, JYEAR, IFIP, IRCLAS, ILINK, VMT,
     &                     MXMAT )

C***********************************************************************
C  program body starts at line 177
C
C  DESCRIPTION:
C       Read MVTVM vehicle mix factors file and map vehicle
C       mixes into sources, to produce VMT mix transform matrix 
C       MXMAT
C
C  PRECONDITIONS REQUIRED:  
C       VDEV already opened and positioned at beginning of file.
C	IFIP, IRCLAS, ILINK from MOBL
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       FIND1, FIND2, FIND3 GETEFILE, GETNUM, GETYN, DAYMON, 
C       TIME2SEC, TRIMLEN, WKDAY
C
C  REVISION  HISTORY:
C       Prototype  1/96 by CJC.
C***********************************************************************

      IMPLICIT NONE


C...........   INCLUDES:
        
        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters


C...........   ARGUMENTS:

        INTEGER         VDEV                   ! unit number for VMT mix file
        INTEGER         JYEAR                  ! simulation year
        INTEGER         IFIP  ( NMSRC )        ! FIP codes for sources
        INTEGER         IRCLAS( NMSRC )        ! EPA road classification code
        INTEGER         ILINK ( NMSRC )        ! link number or 0 for non-link
        REAL            VMT   ( NMSRC )        !  VMT totals by source
        REAL            MXMAT( NMSRC, NVTYPE )  !  VMT mix transform matrix

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1, FIND2, FIND3
        INTEGER         INDEX1
        INTEGER         TRIMLEN
        
        EXTERNAL    FIND1, FIND2, FIND3, INDEX1, TRIMLEN

                                            
C...........   LOCAL VARIABLES and their descriptions:

C...........   Vehicle Mix Factor arrays:
C.......   Before sorting:

        INTEGER         NVMIX                   !  # of mixes from MVMTM file
        INTEGER         INDMA( MXAMIX )         !  index array
        INTEGER         VFIPA( MXAMIX )         !  FIP         for this mix
        INTEGER         VRDTA( MXAMIX )         !  Road type   for this mix
        INTEGER         VLNKA( MXAMIX )         !  link number for this mix
        REAL            VFACA( NVTYPE, MXAMIX ) !  VMT mix profile
        REAL            AFAC ( NVTYPE )         !  VMT mix prof scratch array
                          
C.......   After sorting and filtering into different classes:

        REAL            MFAC0( NVTYPE ) !  Ultimate fallback veh mixes

        INTEGER         NVMXR                   !  roadclass-only veh mixes
        INTEGER         VRDTR( NRCLAS )
        REAL            MFACR( NVTYPE, NRCLAS )

        INTEGER         NVMXS                   !  state&roadclass veh mixes
        INTEGER         VFIPS( MXSTRC )
        INTEGER         VRDTS( MXSTRC )
        REAL            MFACS( NVTYPE, MXSTRC )

        INTEGER         NVMXF                   !  FIP&roadclass veh mixes
        INTEGER         VFIPF( MXFRC )
        INTEGER         VRDTF( MXFRC )
        REAL            MFACF( NVTYPE, MXFRC )

        INTEGER         NVMXL                   !  FIP&roadclass&link veh mixes
        INTEGER         VFIPL( MXFLK )
        INTEGER         VRDTL( MXFLK )
        INTEGER         VLNKL( MXFLK )
        REAL            MFACL( NVTYPE, MXFLK )


C...........   Other local variables

        LOGICAL         EFLAG   !  error-flag
        LOGICAL         DFLAG   !  true if defaults exist
        DATA            DFLAG / .FALSE. /

        INTEGER         IOS
        INTEGER         IREC
        INTEGER         K, F, S, V
        INTEGER         YEAR, FIP, RDT, LNK
        REAL            SCR
        CHARACTER*256   MESG    !  buffer for M3EXIT()


C***********************************************************************
C   begin body of program TMPMOBIL
        
C..........   Read VMT mix profiles file:

        IREC  = 0
        S     = 0
        EFLAG = .FALSE.
        
22      CONTINUE        !  head of loop reading VMT MIX file

            IREC  = IREC + 1
            READ( VDEV, *, IOSTAT=IOS, END = 25 )
     &          YEAR,
     &          FIP,
     &          RDT, 
     &          LNK, 
     &        ( AFAC( V ), V = 1, NVTYPE )


            IF ( YEAR .EQ. JYEAR ) THEN

                S = S + 1

                IF ( S .LE. MXAMIX ) THEN

                    INDMA( S ) = S
                    VFIPA( S ) = FIP
                    VRDTA( S ) = RDT
                    VLNKA( S ) = LNK

                    SCR = 0.0
                    DO  23  V = 1, NVTYPE
                        SCR = SCR + AFAC( V )
23                  CONTINUE
                    IF ( SCR .EQ. 0.0 ) THEN
                        WRITE( MESG,94010 ) 
     &                  'Vehicle mix weights sum to 0 at line', IREC
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    ELSE
                        SCR = 1.0 / SCR
                        DO  24  V = 1, NVTYPE
                            VFACA( V , S ) = SCR * AFAC( V )
24                      CONTINUE
                    END IF

                END IF

            ELSE IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 )
     &              'I/O error', IOS, 
     &              'reading VMT MIX file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                        
            END IF
            
            GO TO 22 !  to head of VMT MIX file loop

25      CONTINUE        !  end of loop reading VMT MIX file

        NVMIX = S
        IF ( S .GT. MXAMIX ) THEN
            WRITE( MESG,94010 )
     &          'Vehicle mix table overflow.  Max', MXAMIX,
     &          'actual', S
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        ELSE IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading VMT MIX file', 2 )
        ELSE IF ( S .EQ. 0 ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'No relevant mixes in VMT MIX file', 2 )
        END IF


C...........   Sort and process the vehicle mix factor arrays:


        NVMXR = 0
        NVMXS = 0
        NVMXF = 0
        NVMXL = 0
        
        DO  39  S = 1, NVMIX
        
            K   = INDMA( S )
            FIP = VFIPA( K )
            RDT = VRDTA( K )
            LNK = VLNKA( K )
            
            IF ( FIP .EQ. 0 ) THEN
                
                IF ( RDT .EQ. 0 ) THEN

                    DFLAG = .TRUE.
                    
                    DO  31  V = 1, NVTYPE
                        MFAC0( V ) = VFACA( V,K )
31                  CONTINUE
                
                ELSE
                    
                    NVMXR = NVMXR + 1

                    IF ( NVMXR .LE. NRCLAS ) THEN

                        VRDTR( NVMXR ) = RDT
                        
                        DO  32  V = 1, NVTYPE
                            MFACR( V,NVMXR ) = VFACA( V,K )
32                      CONTINUE

                    ELSE
                        
                        WRITE( MESG,94010 ) 
     &                  'Vehicle Mix RDTY overflow at subscript', NVMXR
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.

                    END IF
                
                END IF
            
            ELSE IF ( MOD( FIP, 1000 ) .EQ. 0 ) THEN

                NVMXS = NVMXS + 1

                IF ( NVMXS .LE. MXSTRC ) THEN

                    VFIPS( NVMXS ) = FIP / 1000
                    VRDTS( NVMXS ) = RDT
                    DO  33  V = 1, NVTYPE
                        MFACS( V,NVMXS ) = VFACA( V,K )
33                  CONTINUE

                ELSE

                    WRITE( MESG,94010 )
     &              'Vehicle Mix STATE overflow at subscript', NVMXS
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.

                END IF

            ELSE IF ( LNK .EQ. 0 ) THEN

                NVMXF = NVMXF + 1

                IF ( NVMXF .LE. MXFRC ) THEN

                    VFIPF( NVMXF ) = FIP
                    VRDTF( NVMXF ) = RDT
                    DO  34  V = 1, NVTYPE
                        MFACF( V, NVMXF ) = VFACA( V,K )
34                  CONTINUE

                ELSE

                    WRITE( MESG,94010 )
     &              'Vehicle Mix FIP overflow at subscript', NVMXF
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.

                END IF

            ELSE

                NVMXL = NVMXL + 1

                IF ( NVMXL .LE. MXFLK ) THEN

                    VFIPL( NVMXL ) = FIP
                    VRDTL( NVMXL ) = RDT
                    VLNKL( NVMXL ) = LNK
                    DO  35  V = 1, NVTYPE
                        MFACL( V, NVMXL ) = VFACA( V,K )
35                  CONTINUE

                ELSE

                    WRITE( MESG,94010 )
     &              'Vehicle Mix LNK overflow at subscript', NVMXL
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.

                END IF

            END IF
            
39      CONTINUE
        
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error processing VEHICLE MIX file', 2 )
        END IF
        

C...........   Construct mapping of source number to vehicle mixes.

        WRITE( *,92000 ) ' ', 'Building vehicle mix source indexes.'
        EFLAG = .FALSE.

        DO  299  S = 1, NMSRC

            FIP = IFIP  ( S )
            RDT = IRCLAS( S )
            LNK = ILINK ( S )
        
            F = FIND3( FIP, RDT, LNK, NVMXL, VFIPL, VRDTL, VLNKL )

            IF ( F .GT. 0 ) THEN        !  fip-rdt-lnk match
                DO  281  V = 1, NVTYPE
                    MXMAT( S,V ) = MFACL( V,F ) * VMT( S )
281             CONTINUE
                GO TO 299               !  to end of mix file processing
            END IF

            F = FIND2( FIP, RDT, NVMXF, VFIPF, VRDTF )

            IF ( F .GT. 0 ) THEN        !  fip-rdt match
                DO  282  V = 1, NVTYPE
                    MXMAT( S,V ) = MFACF( V,F ) * VMT( S )
282             CONTINUE
                GO TO 299               !  to end of mix file processing
            END IF

            F = FIND2( FIP/1000, RDT, NVMXS, VFIPS, VRDTS )

            IF ( F .GT. 0 ) THEN        !  state-rdty match
                DO  283  V = 1, NVTYPE
                    MXMAT( S,V ) = MFACS( V,F ) * VMT( S )
283             CONTINUE
                GO TO 299               !  to end of mix file processing
            END IF

            F = FIND1( RDT, NVMXR, VRDTR )

            IF ( F .GT. 0 ) THEN        !  rdt match

                DO  284  V = 1, NVTYPE
                    MXMAT( S,V ) = MFACR( V,F ) * VMT( S )
284             CONTINUE

            ELSEIF( DFLAG ) THEN        !  no match:  use fallback fractions

                DO  285  V = 1, NVTYPE
                    MXMAT( S,V ) = MFAC0( V ) * VMT( S )
285             CONTINUE

            ELSE                        !  Error

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &                 'No VMT MIX match found for FIP:', FIP, 
     &                 ' Road class:', RDT
                CALL M3MESG( MESG )

            END IF

                        
299     CONTINUE        !  end loop on sources S, constructing TNDX
        
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error applying VMT mix from MVMTM', 2 )
        END IF
       
        RETURN


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

        END

