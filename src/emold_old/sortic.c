
/***********************************************************************
   void qsortic() for    C   starts at line   70
   void SORTIC () for Fortran starts at line 173

    PURPOSE:
        qsortic() and SORTIC() sort index-tables ind[] on the basis of 
        Fortran character-string arrays holding keys.  SORTIC() is
        designed to be called from Fortran, so its subscripting from ind[] 
        is offset by 1 for indexing into key tables  tblc[]

    ALGORITHM:
        quick sort (q.v. Sedgwick, _Algorithms_, 2nd. ed., chapter 9)
       
   PRECONDITIONS:
        ind[ N ] initialized with 1, 2, ..., N      for SORTIC(),
                                  0, k, ..., (N-1)k for qsortic()
    
    REVISION HISTORY:
        Prototypes 8/96 by CJC
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

#include  <string.h>

                         /**  DEAL WITH  FELDMANISMS  OF MANY UN*X F77'S   **/
 
#if __sgi || __sun || __osf__ || __mips__
 
#define  SORTIC   sortic_

#elif __hpux || _AIX

#define  SORTIC   sortic

#endif                                   /** #IF SGI OR SUN OR OSF OR MIPS **/


                    /** MACROS FOR COMPARING INDICES INTO KEY-TUPLE TABLES **/
        /** CMPk( X,Y ) iff index X > index Y in terms of k-tuple tables   **/
        /** i.e.,  1  iff *out*of order,  -1 iff *in*order,  0 if *equal*  **/

#define CMPC( X,Y )   strncmp( tblc + (X) , tblc + (Y) , (size_t) k )


/********************* BODIES OF THE PRIVATE SORT-ROUTINES *******************/

void qsortic( int        n,          /** number of elements              **/
              int        ind[],      /** index-array                     **/
              const char tblc[],     /** first  key-component in tuple   **/
              const int  k )         /** key-length as a Fortran string  **/
    {
    int  i, j ;
    int  k1 ;
    int  a, b, c ;
    int  l, r , p, q ;
    int  t, u, v, w ;
    
    if ( n > 2 )
        {
                /** DO SORT-3 TO GET MEDIAN-OF-3  -- Q.V.KNUTH VOL3 P. 182 **/
        p = n / 2 ;
        a = ind[ 0 ] ;
        b = ind[ p ] ;
        c = ind[ n-1 ] ;
        u = CMPC( a,b ) ;
        v = CMPC( b,c ) ;

        if ( u > 0 )                                       /** A,B REVERSED **/
            {
            if ( v > 0 )                                   /** ABC ~~~> CBA **/
                {
                ind[ 0   ] = c ;
                ind[ n-1 ] = a ;
                }
            else{
                w = CMPC( a,c ) ;
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
            w = CMPC( a,c ) ;
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
            b        = ind[ p ] ;
            ind[ p ] = ind[ 0 ] ;
            ind[ 0 ] = b ;

            for( l = 1, r = n-1 ; ; l++, r-- )
                {
                for( ; CMPC( ind[l], b ) < 0 ; l++ ) ;       /** EMPTY BODY **/

                for( ; CMPC( ind[r], b ) > 0 ; r-- ) ;       /** EMPTY BODY **/

                if ( l < r )
                    {
                    t        = ind[ l ] ;
                    ind[ l ] = ind[ r ] ;
                    ind[ r ] = t ;
                    }
                else break ;
                }

            ind[ 0 ] = ind[ r ] ;
            ind[ r ] = b ;

            qsortic( r,   ind,   tblc, k ) ;
            qsortic( n-l, ind+l, tblc, k ) ;
            }                                     /** END IF-CLAUSE:  N > 3 **/
        }                                         /** END IF-CLAUSE:  N > 2 **/

    else if ( n == 2 )
        {
        a = ind[ 0 ] ;
        b = ind[ 1 ] ;
        if ( CMPC( a,b ) > 0 )
            {
            ind[ 0 ] = b ;
            ind[ 1 ] = a ;
            }
        }                                     /** END ELSE-IF CLAUSE:  N==2 **/

    }   /** .............................................END VOID qsortic() **/


/**********************   BODY OF THE F77 SORT-ROUTINE   *******************/

void SORTIC( const int *nelts,          /** number of elements              **/
             int        ind[],          /** index-array                     **/
             const char tblc[],         /** table: key-component in 1-tuple **/
             const int  k )             /** string length in table          **/
    {
    int n, i, t ;

    n = *nelts ;

                              /** CONVERT FROM FORTRAN SUBSCRIPTS TO C SUBS **/
    for ( i = 0 ; i < n ; i++ )           /** ALSO ADJUST FOR STRING LENGTH **/
        {
        ind[ i ] = k * ( ind[ i ] - 1 ) ;
        } ;

                                             /** CALL C QSORT ROUTINE **/
    qsortic( n, ind, tblc, k ) ;
        
                         /** CONVERT FROM C SUBSCRIPTS BACK TO FORTRAN SUBS **/
    for ( i = 0 ; i < n ; i++ )        /** ALSO DE-ADJUST FOR STRING LENGTH **/
        {
        ind[ i ] = 1 + ( ind[ i ] / k ) ;
        }

    return ;

    } /*********************************************  END FUNCTION SORTIC() **/


