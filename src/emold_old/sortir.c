
/***********************************************************************
   void qsortr1() for    C   starts at line 116
   void qsortr2() for    C   starts at line 219
   void qsortr3() for    C   starts at line 324
   void qsortr4() for    C   starts at line 431
   void SORTR1() for Fortran starts at line 542
   void SORTR2() for Fortran starts at line 570
   void SORTR3() for Fortran starts at line 598
   void SORTR4() for Fortran starts at line 628

    PURPOSE:
        qsortr<k>() and SORTR<k>() sort index-tables ind[] on the basis of 
        <k> parallel arrays holding <k>-tuple keys.  The SORTR<k>() are 
        designed to be called from Fortran, so their subscripting from ind[] 
        is offset by 1 for indexing into key tables  tbl<k>[]

    ALGORITHM:
        quick sort (q.v. Sedgwick, _Algorithms_, 2nd. ed., chapter 9)
       
   PRECONDITIONS:
        ind[ N ] initialized with 1, ..., N  
        (Fortran-style subscripts to the key tables tbl<k>[])
    
    REVISION HISTORY:
        Prototypes 9/95 by CJC
************************************************************************/

/****************************************************************************
*
* Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
*                System
* File: @(#)$Id$
*
* COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
* All Rights Reserved
*
* See file COPYRIGHT for conditions of use.
*
* Environmental Programs Group
* MCNC--North Carolina Supercomputing Center
* P.O. Box 12889
* Research Triangle Park, NC  27709-2889
*
* env_progs@mcnc.org
*
* Pathname: $Source$
* Last updated: $Date$ 
*
****************************************************************************/
 
                         /**  DEAL WITH  FELDMANISMS  OF MANY UN*X F77'S   **/
 
#if __sgi || __sun || __osf__ || __mips__
 
#define  SORTR1   sortr1_
#define  SORTR2   sortr2_
#define  SORTR3   sortr3_
#define  SORTR4   sortr4_

#elif __hpux || _AIX

#define  SORTR1   sortr1
#define  SORTR2   sortr2
#define  SORTR3   sortr3
#define  SORTR4   sortr4

#endif                                   /** #IF SGI OR SUN OR OSF OR MIPS **/


                    /** MACROS FOR COMPARING INDICES INTO KEY-TUPLE TABLES **/
        /** CMPk( X,Y ) iff index X > index Y in terms of k-tuple tables   **/
        /** LCMPk( I,K1, ...,Kk ) iff k-tuple[ ind[ I ] ] > (Y1, ..., Yk)  **/
        /** i.e.,  1  iff *out*of order,  -1 iff *in*order,  0 if *equal*  **/

#define CMP1( X,Y ) \
 ( tbl1[ X ] > tbl1[ Y ] ? 1 : ( tbl1[ X ] < tbl1[ Y ] ? -1 : 0 ) )

#define CMP2( X,Y ) \
 ( tbl1[ X ] > tbl1[ Y ] ? 1 : ( tbl1[ X ] < tbl1[ Y ] ? -1 : \
 ( tbl2[ X ] > tbl2[ Y ] ? 1 : ( tbl2[ X ] < tbl2[ Y ] ? -1 : 0 ) ) ) )

#define CMP3( X,Y ) \
 ( tbl1[ X ] > tbl1[ Y ] ? 1 : ( tbl1[ X ] < tbl1[ Y ] ? -1 : \
 ( tbl2[ X ] > tbl2[ Y ] ? 1 : ( tbl2[ X ] < tbl2[ Y ] ? -1 : \
 ( tbl3[ X ] > tbl3[ Y ] ? 1 : ( tbl3[ X ] < tbl3[ Y ] ? -1 : 0 ) ) ) ) ) )

#define CMP4( X,Y ) \
 ( tbl1[ X ] > tbl1[ Y ] ? 1 : ( tbl1[ X ] < tbl1[ Y ] ? -1 : \
 ( tbl2[ X ] > tbl2[ Y ] ? 1 : ( tbl2[ X ] < tbl2[ Y ] ? -1 : \
 ( tbl3[ X ] > tbl3[ Y ] ? 1 : ( tbl3[ X ] < tbl3[ Y ] ? -1 : \
 ( tbl4[ X ] > tbl4[ Y ] ? 1 : ( tbl4[ X ] < tbl4[ Y ] ? -1 : 0 ) ) ) ) ) ) ) )


#define LCMP1( I,K1 ) \
 ( tbl1[ ind[ I ] ] > K1 ? 1 : ( tbl1[ ind[ I ] ] < K1 ? -1 : 0 ) )

#define LCMP2( I,K1,K2 ) \
 ( tbl1[ ind[ I ] ] > K1 ? 1 : ( tbl1[ ind[ I ] ] < K1 ? -1 : \
 ( tbl2[ ind[ I ] ] > K2 ? 1 : ( tbl2[ ind[ I ] ] < K2 ? -1 : 0 ) ) ) )

#define LCMP3( I,K1,K2,K3 ) \
 ( tbl1[ ind[ I ] ] > K1 ? 1 : ( tbl1[ ind[ I ] ] < K1 ? -1 : \
 ( tbl2[ ind[ I ] ] > K2 ? 1 : ( tbl2[ ind[ I ] ] < K2 ? -1 : \
 ( tbl3[ ind[ I ] ] > K3 ? 1 : ( tbl3[ ind[ I ] ] < K3 ? -1 : 0 ) ) ) ) ) )

#define LCMP4( I,K1,K2,K3,K4 ) \
 ( tbl1[ ind[ I ] ] > K1 ? 1 : ( tbl1[ ind[ I ] ] < K1 ? -1 : \
 ( tbl2[ ind[ I ] ] > K2 ? 1 : ( tbl2[ ind[ I ] ] < K2 ? -1 : \
 ( tbl3[ ind[ I ] ] > K3 ? 1 : ( tbl3[ ind[ I ] ] < K3 ? -1 : \
 ( tbl4[ ind[ I ] ] > K4 ? 1 : ( tbl4[ ind[ I ] ] < K4 ? -1 : 0 ) ) ) ) ) ) ) )


/********************* BODIES OF THE PRIVATE SORT-ROUTINES *******************/

void qsortr1( int         n,          /** number of elements             **/
              int         ind[],      /** index-array                    **/
              const float tbl1[] )    /** first  key-component in tuple  **/
    {
    int   i, j ;
    float k1 ;
    int   a, b, c ;
    int   l, r , p, q ;
    int   t, u, v, w ;
    
    if ( n > 2 )
        {
                /** DO SORT-3 TO GET MEDIAN-OF-3  -- Q.V.KNUTH VOL3 P. 182 **/
        p = n / 2 ;
        a = ind[ 0 ] ;
        b = ind[ p ] ;
        c = ind[ n-1 ] ;
        u = CMP1( a,b ) ;
        v = CMP1( b,c ) ;
        w = CMP1( a,c ) ;

        if ( u > 0 )                                       /** A,B REVERSED **/
            {
            if ( v > 0 )                                   /** ABC ~~~> CBA **/
                {
                ind[ 0   ] = c ;
                ind[ n-1 ] = a ;
                }
            else{
                if ( w > 0 )                               /** ABC ~~~> BCA **/
                    {
                    ind[ 0   ] = b ;
                    ind[ p   ] = c ;
                    ind[ n-1 ] = a ;
                    }
                else{                                      /** ABC ~~~> BAC **/
                    ind[ 0   ] = b ;
                    ind[ p   ] = a ;
                    }
                }
            }
        else if ( v > 0 )                           /** A,B OK; BC REVERSED **/
            {
            if ( w > 0 )                                   /** ABC ~~~> CAB **/
                {
                ind[ 0   ] = c ;
                ind[ p   ] = a ;
                ind[ n-1 ] = b ;
                }
            else{                                          /** ABC ~~~> ACB **/
                ind[ p   ] = c ;
                ind[ n-1 ] = b ;
                }
            }
        
                                        /** IF N > 3, PARTITION AND RECURSE **/
        if ( n > 3 ) 
            {
            b  = ind[  p ] ;
            k1 = tbl1[ b ] ;
            
            t        = ind[ p ] ;
            ind[ p ] = ind[ 0 ] ;
            ind[ 0 ] = t ;

            for( l = 1, r = n-1 ; ; l++, r-- )
                {
                for( ; LCMP1( l,k1 ) < 0 ; l++ ) ;           /** EMPTY BODY **/

                for( ; LCMP1( r,k1 ) > 0 ; r-- ) ;           /** EMPTY BODY **/

                if ( l < r )
                    {
                    t        = ind[ l ] ;
                    ind[ l ] = ind[ r ] ;
                    ind[ r ] = t ;
                    }
                else break ;
                }

            t        = ind[ r ] ;
            ind[ r ] = ind[ 0 ] ;
            ind[ 0 ] = t ;

            qsortr1( r,   ind,   tbl1 ) ;
            qsortr1( n-l, ind+l, tbl1 ) ;
            }                                     /** END IF-CLAUSE:  N > 3 **/
        }                                         /** END IF-CLAUSE:  N > 2 **/

    else if ( n == 2 )
        {
        a = ind[ 0 ] ;
        b = ind[ 1 ] ;
        if ( CMP1( a,b ) > 0 )
            {
            ind[ 0 ] = b ;
            ind[ 1 ] = a ;
            }
        }                                     /** END ELSE-IF CLAUSE:  N==2 **/

    }   /** .............................................END VOID qsortr1() **/


void qsortr2( int         n,           /** number of elements             **/
              int         ind[],       /** index-array                    **/
              const float tbl1[],      /** first  key-component in tuple  **/
              const float tbl2[] )     /** second key-component in tuple  **/
    {
    int   i, j ;
    float k1, k2 ;
    int   a, b, c ;
    int   l, r , p, q ;
    int   t, u, v, w ;
    
    if ( n > 2 )
        {
                /** DO SORT-3 TO GET MEDIAN-OF-3  -- Q.V.KNUTH VOL3 P. 182 **/
        p = n / 2 ;
        a = ind[ 0 ] ;
        b = ind[ p ] ;
        c = ind[ n-1 ] ;
        u = CMP2( a,b ) ;
        v = CMP2( b,c ) ;
        w = CMP2( a,c ) ;

        if ( u > 0 )                                       /** A,B REVERSED **/
            {
            if ( v > 0 )                                   /** ABC ~~~> CBA **/
                {
                ind[ 0   ] = c ;
                ind[ n-1 ] = a ;
                }
            else{
                if ( w > 0 )                               /** ABC ~~~> BCA **/
                    {
                    ind[ 0   ] = b ;
                    ind[ p   ] = c ;
                    ind[ n-1 ] = a ;
                    }
                else{                                      /** ABC ~~~> BAC **/
                    ind[ 0   ] = b ;
                    ind[ p   ] = a ;
                    }
                }
            }
        else if ( v > 0 )                           /** A,B OK; BC REVERSED **/
            {
            if ( w > 0 )                                   /** ABC ~~~> CAB **/
                {
                ind[ 0   ] = c ;
                ind[ p   ] = a ;
                ind[ n-1 ] = b ;
                }
            else{                                          /** ABC ~~~> ACB **/
                ind[ p   ] = c ;
                ind[ n-1 ] = b ;
                }
            }
        
                                        /** IF N > 3, PARTITION AND RECURSE **/
        if ( n > 3 ) 
            {
            b  = ind[  p ] ;
            k1 = tbl1[ b ] ;
            k2 = tbl2[ b ] ;
            
            t        = ind[ p ] ;
            ind[ p ] = ind[ 0 ] ;
            ind[ 0 ] = t ;

            for( l = 1, r = n-1 ; ; l++, r-- )
                {
                for( ; LCMP2( l,k1,k2 ) < 0 ; l++ ) ;        /** EMPTY BODY **/

                for( ; LCMP2( r,k1,k2 ) > 0 ; r-- ) ;        /** EMPTY BODY **/

                if ( l < r )
                    {
                    t        = ind[ l ] ;
                    ind[ l ] = ind[ r ] ;
                    ind[ r ] = t ;
                    }
                else break ;
                }

            t        = ind[ r ] ;
            ind[ r ] = ind[ 0 ] ;
            ind[ 0 ] = t ;

            qsortr2( r,   ind,   tbl1, tbl2 ) ;
            qsortr2( n-l, ind+l, tbl1, tbl2 ) ;
            }                                     /** END IF-CLAUSE:  N > 3 **/
        }                                         /** END IF-CLAUSE:  N > 2 **/

    else if ( n == 2 )
        {
        a = ind[ 0 ] ;
        b = ind[ 1 ] ;
        if ( CMP2( a,b ) > 0 )
            {
            ind[ 0 ] = b ;
            ind[ 1 ] = a ;
            }
        }                                     /** END ELSE-IF CLAUSE:  N==2 **/

    }   /** .............................................END VOID qsortr2() **/


void qsortr3( int         n,           /** number of elements             **/
              int         ind[],       /** index-array                    **/
              const float tbl1[],      /** first  key-component in tuple  **/
              const float tbl2[],      /** second key-component in tuple  **/
              const float tbl3[] )     /** third  key-component in tuple  **/
    {
    int   i, j ;
    float k1, k2, k3 ;
    int   a, b, c ;
    int   l, r , p, q ;
    int   t, u, v, w ;
    
    if ( n > 2 )
        {
                /** DO SORT-3 TO GET MEDIAN-OF-3  -- Q.V.KNUTH VOL3 P. 182 **/
        p = n / 2 ;
        a = ind[ 0 ] ;
        b = ind[ p ] ;
        c = ind[ n-1 ] ;
        u = CMP3( a,b ) ;
        v = CMP3( b,c ) ;
        w = CMP3( a,c ) ;

        if ( u > 0 )                                       /** A,B REVERSED **/
            {
            if ( v > 0 )                                   /** ABC ~~~> CBA **/
                {
                ind[ 0   ] = c ;
                ind[ n-1 ] = a ;
                }
            else{
                if ( w > 0 )                               /** ABC ~~~> BCA **/
                    {
                    ind[ 0   ] = b ;
                    ind[ p   ] = c ;
                    ind[ n-1 ] = a ;
                    }
                else{                                      /** ABC ~~~> BAC **/
                    ind[ 0   ] = b ;
                    ind[ p   ] = a ;
                    }
                }
            }
        else if ( v > 0 )                           /** A,B OK; BC REVERSED **/
            {
            if ( w > 0 )                                   /** ABC ~~~> CAB **/
                {
                ind[ 0   ] = c ;
                ind[ p   ] = a ;
                ind[ n-1 ] = b ;
                }
            else{                                          /** ABC ~~~> ACB **/
                ind[ p   ] = c ;
                ind[ n-1 ] = b ;
                }
            }
        
                                        /** IF N > 3, PARTITION AND RECURSE **/
        if ( n > 3 ) 
            {
            b  = ind[  p ] ;
            k1 = tbl1[ b ] ;
            k2 = tbl2[ b ] ;
            k3 = tbl3[ b ] ;
            
            t        = ind[ p ] ;
            ind[ p ] = ind[ 0 ] ;
            ind[ 0 ] = t ;

            for( l = 1, r = n-1 ; ; l++, r-- )
                {
                for( ; LCMP3( l,k1,k2,k3 ) < 0 ; l++ ) ;     /** EMPTY BODY **/

                for( ; LCMP3( r,k1,k2,k3 ) > 0 ; r-- ) ;     /** EMPTY BODY **/

                if ( l < r )
                    {
                    t        = ind[ l ] ;
                    ind[ l ] = ind[ r ] ;
                    ind[ r ] = t ;
                    }
                else break ;
                }

            t        = ind[ r ] ;
            ind[ r ] = ind[ 0 ] ;
            ind[ 0 ] = t ;

            qsortr3( r,   ind,   tbl1, tbl2, tbl3 ) ;
            qsortr3( n-l, ind+l, tbl1, tbl2, tbl3 ) ;
            }                                     /** END IF-CLAUSE:  N > 3 **/
        }                                         /** END IF-CLAUSE:  N > 2 **/

    else if ( n == 2 )
        {
        a = ind[ 0 ] ;
        b = ind[ 1 ] ;
        if ( CMP3( a,b ) > 0 )
            {
            ind[ 0 ] = b ;
            ind[ 1 ] = a ;
            }
        }                                     /** END ELSE-IF CLAUSE:  N==2 **/

    }  /** ..............................................END VOID qsortr3() **/


void qsortr4( int         n,           /** number of elements             **/
              int         ind[],       /** index-array                    **/
              const float tbl1[],      /** first  key-component in tuple  **/
              const float tbl2[],      /** second key-component in tuple  **/
              const float tbl3[],      /** third  key-component in tuple  **/
              const float tbl4[] )     /** fourth key-component in tuple  **/
    {
    int   i, j ;
    float k1, k2, k3, k4 ;
    int   a, b, c ;
    int   l, r , p, q ;
    int   t, u, v, w ;
    
    if ( n > 2 )
        {
                /** DO SORT-3 TO GET MEDIAN-OF-3  -- Q.V.KNUTH VOL3 P. 182 **/
        p = n / 2 ;
        a = ind[ 0 ] ;
        b = ind[ p ] ;
        c = ind[ n-1 ] ;
        u = CMP4( a,b ) ;
        v = CMP4( b,c ) ;
        w = CMP4( a,c ) ;

        if ( u > 0 )                                       /** A,B REVERSED **/
            {
            if ( v > 0 )                                   /** ABC ~~~> CBA **/
                {
                ind[ 0   ] = c ;
                ind[ n-1 ] = a ;
                }
            else{
                if ( w > 0 )                               /** ABC ~~~> BCA **/
                    {
                    ind[ 0   ] = b ;
                    ind[ p   ] = c ;
                    ind[ n-1 ] = a ;
                    }
                else{                                      /** ABC ~~~> BAC **/
                    ind[ 0   ] = b ;
                    ind[ p   ] = a ;
                    }
                }
            }
        else if ( v > 0 )                           /** A,B OK; BC REVERSED **/
            {
            if ( w > 0 )                                   /** ABC ~~~> CAB **/
                {
                ind[ 0   ] = c ;
                ind[ p   ] = a ;
                ind[ n-1 ] = b ;
                }
            else{                                          /** ABC ~~~> ACB **/
                ind[ p   ] = c ;
                ind[ n-1 ] = b ;
                }
            }
        
                                        /** IF N > 3, PARTITION AND RECURSE **/
        if ( n > 3 ) 
            {
            b  = ind[  p ] ;
            k1 = tbl1[ b ] ;
            k2 = tbl2[ b ] ;
            k3 = tbl3[ b ] ;
            k4 = tbl4[ b ] ;
            
            t        = ind[ p ] ;
            ind[ p ] = ind[ 0 ] ;
            ind[ 0 ] = t ;

            for( l = 1, r = n-1 ; ; l++, r-- )
                {
                for( ; LCMP4( l,k1,k2,k3,k4 ) < 0 ; l++ ) ;  /** EMPTY BODY **/

                for( ; LCMP4( r,k1,k2,k3,k4 ) > 0 ; r-- ) ;  /** EMPTY BODY **/

                if ( l < r )
                    {
                    t        = ind[ l ] ;
                    ind[ l ] = ind[ r ] ;
                    ind[ r ] = t ;
                    }
                else break ;
                }

            t        = ind[ r ] ;
            ind[ r ] = ind[ 0 ] ;
            ind[ 0 ] = t ;

            qsortr4( r,   ind,   tbl1, tbl2, tbl3, tbl4 ) ;
            qsortr4( n-l, ind+l, tbl1, tbl2, tbl3, tbl4 ) ;
            }                                     /** END IF-CLAUSE:  N > 3 **/
        }                                         /** END IF-CLAUSE:  N > 2 **/

    else if ( n == 2 )
        {
        a = ind[ 0 ] ;
        b = ind[ 1 ] ;
        if ( CMP4( a,b ) > 0 )
            {
            ind[ 0 ] = b ;
            ind[ 1 ] = a ;
            }
        }                                     /** END ELSE-IF CLAUSE:  N==2 **/

    }  /**...............................................END VOID qsortr4() **/


/**********************   BODIES OF THE F77 SORT-ROUTINES   *******************/

void SORTR1( const int * nelts,          /** number of elements              **/
             int         ind[],          /** index-array                     **/
             const float tbl1[] )        /** table: key-component in 1-tuple **/
    {
    int n, i ;

    n = *nelts ;

                              /** CONVERT FROM FORTRAN SUBSCRIPTS TO C SUBS **/
    for ( i = 0 ; i < n ; i++ )
        {
        ind[ i ]-- ;
        } ;

                                             /** CALL C QSORT ROUTINE **/
    qsortr1( n, ind, tbl1 ) ;
        
                         /** CONVERT FROM C SUBSCRIPTS BACK TO FORTRAN SUBS **/
    for ( i = 0 ; i < n ; i++ )
        {
        ind[ i ]++ ;
        }

    return ;

    } /***********************************************  END FUNCTION SORTR) **/


void SORTR2( const int * nelts,    /** number of elements                     **/
             int         ind[],    /** index-array                            **/
             const float tbl1[],   /** table:  first  key-component in tuple  **/
             const float tbl2[] )  /** table:  second key-component in tuple  **/
    {    
    int n, i, u ;

    n = *nelts ;

                              /** CONVERT FROM FORTRAN SUBSCRIPTS TO C SUBS **/
    for ( i = 0 ; i < n ; i++ )
        {
        ind[ i ]-- ;
        } ;

                                           /** CALL C-BINDING QSORT ROUTINE **/
    qsortr2( n, ind, tbl1, tbl2 ) ;
        
                         /** CONVERT FROM C SUBSCRIPTS BACK TO FORTRAN SUBS **/
    for ( i = 0 ; i < n ; i++ )
        {
        ind[ i ]++ ;
        }

    return ;

    } /*********************************************  END FUNCTION SORTR2() **/

void SORTR3( const int * nelts,   /** number of elements                     **/
             int         ind[],   /** index-array                            **/
             const float tbl1[],  /** table:  first  key-component in tuple  **/
             const float tbl2[],  /** table:  second key-component in tuple  **/
             const float tbl3[] ) /** table:  third  key-component in tuple  **/
    {
    int n, i ;

    n = *nelts ;

                              /** CONVERT FROM FORTRAN SUBSCRIPTS TO C SUBS **/
    for ( i = 0 ; i < n ; i++ )
        {
        ind[ i ]-- ;
        } ;

                                           /** CALL C-BINDING QSORT ROUTINE **/
    qsortr3( n, ind, tbl1, tbl2, tbl3 ) ;
        
                         /** CONVERT FROM C SUBSCRIPTS BACK TO FORTRAN SUBS **/
    for ( i = 0 ; i < n ; i++ )
        {
        ind[ i ]++ ;
        }

    return ;

    } /*********************************************  END FUNCTION SORTR3() **/


void SORTR4( const int * nelts,    /** number of elements                     **/
             int         ind[],    /** index-array                            **/
             const float tbl1[],   /** table:  first  key-component in tuple  **/
             const float tbl2[],   /** table:  second key-component in tuple  **/
             const float tbl3[],   /** table:  third  key-component in tuple  **/
             const float tbl4[] )  /** table:  fourth key-component in tuple  **/
    {
    int n, i ;

    n = *nelts ;

                              /** CONVERT FROM FORTRAN SUBSCRIPTS TO C SUBS **/
    for ( i = 0 ; i < n ; i++ )
        {
        ind[ i ]-- ;
        } ;

                                           /** CALL C-BINDING QSORT ROUTINE **/
    qsortr4( n, ind, tbl1, tbl2, tbl3, tbl4 ) ;
        
                         /** CONVERT FROM C SUBSCRIPTS BACK TO FORTRAN SUBS **/
    for ( i = 0 ; i < n ; i++ )
        {
        ind[ i ]++ ;
        }

    return ;

    } /*********************************************  END FUNCTION SORTR4() **/

