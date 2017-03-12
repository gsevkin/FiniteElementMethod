//==============================================================================
//
//  $Id: Gauss_Points.C,v 1.9 2006/07/14 10:51:30 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
template <typename TYPE, int nr_>
Gauss_Init<TYPE,nr_>::Gauss_Init ( )
   {
     gaussianPoints = new TYPE*[numberGaussianPoints+1] ;
     for ( int i = 0 ;i < numberGaussianPoints ; ++ i )
        {
          gaussianPoints[i] = new TYPE [4] ;
#ifdef INITIALIZE_ARRAYS
          for (int k=0; k<4; ++k)
             {
               gaussianPoints[i][k] = 0.;
             }
#endif
        }
     stress = new TYPE [3] ;
#ifdef INITIALIZE_ARRAYS
     for (int k=0; k<3; ++k)
        {
          stress[k] = 0.;
        }
#endif
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int nr_>
Gauss_Init<TYPE,nr_>::~Gauss_Init ( )
   {
      delete [] stress ;
      for ( int gaussIter = 0 ; gaussIter < numberGaussianPoints ; ++gaussIter)
         {
           delete [] gaussianPoints [gaussIter] ;
         }
      delete [] gaussianPoints ;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int nr_>
inline TYPE**
Gauss_Init<TYPE,nr_>::getGaussianPoints() const
   {
      return gaussianPoints ;
   }
//-----------------------------------------------------------------------------
template <typename TYPE>
Gaussian_Points<TYPE, Gauss3 , cuboid > :: Gaussian_Points()
   {
     this->stress [0] = 0.5 ;
     this->stress [1] = 0.5 ;
     this->stress [2] = 0.5 ;
     TYPE x [3] =
        {
          0.5 + sqrt ( 0.15 ) , 
          0.5 ,
          0.5 - sqrt ( 0.15 )
         } ;
     TYPE w [3] =
        {
          5./18. ,
          4./9. ,
          5./18.
        };
     int c = 0;
     for ( int i = 0; i < 3 ; ++ i )
        {
          for ( int j = 0 ; j < 3 ; ++ j )
              {
                for ( int k = 0 ; k < 3 ; ++ k )
                   {
                     this->gaussianPoints [c] [0] = x [i] ;
                     this->gaussianPoints [c] [1] = x [j] ;
                     this->gaussianPoints [c] [2] = x [k] ;
                     this->gaussianPoints [c++] [3] = w [i] * w [j] * w [k] ;
                   }
              }
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss2 , cuboid >::Gaussian_Points()
   {
     this->stress [0] = 0.5 ; 
     this->stress [1] = 0.5 ;  
     this->stress [2] = 0.5 ;
     TYPE x [2] = 
       {
         0.5 * ( 1. + 1.0 / sqrt ( 3. ) ) ,
         0.5 * ( 1. - 1.0 / sqrt ( 3. ) ) 
       } ;
     TYPE w [2] = 
       {
         1. / 2. ,
         1. / 2.
       } ; 
     int c = 0 ;
     for ( int i = 0 ; i < 2 ; ++ i )
        { 
          for ( int j = 0 ; j < 2 ; ++ j )
             {
               for ( int k = 0 ; k < 2 ; ++ k )
                  {
                    this->gaussianPoints [c] [0] = x [i] ;
                    this->gaussianPoints [c] [1] = x [j] ;
                    this->gaussianPoints [c] [2] = x [k] ;
                    this->gaussianPoints [c++] [3] = w [i] * w [j] * w[k] ;
                  }
             }
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss1 , cuboid >::Gaussian_Points() 
   {
     this->stress [0] = 0.5 ; 
     this->stress [1] = 0.5 ;  
     this->stress [2] = 0.5 ;
     this->gaussianPoints [0] [0] = 0.5 ;
     this->gaussianPoints [0] [1] = 0.5 ;
     this->gaussianPoints [0] [2] = 0.5 ;
     this->gaussianPoints [0] [3] = 1. ;
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss2 , tetrahedron >::Gaussian_Points()
   {
     this->stress [0] = 0.25 ; 
     this->stress [1] = 0.25 ;  
     this->stress [2] = 0.25 ;
     TYPE x [2] = 
        {
          ( 5. - sqrt ( 5. ) ) / 20. ,
          ( 5. + 3. * sqrt ( 5. ) ) / 20.
        } ;
     TYPE w [2] = 
        {
          1. / 24. ,
          1. / 24. 
        } ;
     int c = 0 ;
     for ( int i = 0 ; i < 2 ; ++ i )
        { 
          for ( int j = 0 ; j < 2 - i ; ++ j )
             { 
               for ( int k = 0 ; k < 2 - i - j ; ++ k )
                  {
                    this->gaussianPoints [c] [0] = x [i] ;
                    this->gaussianPoints [c] [1] = x [j] ;
                    this->gaussianPoints [c] [2] = x [k] ;
                    this->gaussianPoints [c++] [3] = w [0];
                  }  
             }
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss1 , tetrahedron >::Gaussian_Points()
   {
     this->stress [0] = 0.25 ; 
     this->stress [1] = 0.25 ;  
     this->stress [2] = 0.25 ;
     this->gaussianPoints [0] [0] = 0.25 ; 
     this->gaussianPoints [0] [1] = 0.25 ;
     this->gaussianPoints [0] [2] = 0.25 ;
     this->gaussianPoints [0] [3] = 1. / 6. ;
  }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss2 , pyramid >::Gaussian_Points()
   {
     this->stress [0] = 0.4 ; 
     this->stress [1] = 0.4 ;  
     this->stress [2] = 0.2 ;
     TYPE x [2] = 
        {
          ( 5. - sqrt ( 5. ) ) / 20. ,
          ( 5. + 3. * sqrt ( 5. ) ) / 20.
        } ;
     TYPE w [2] = 
        {
          1. / 24. ,
          1. / 24. 
        } ;
     int c = 0 ;
     for ( int i = 0 ; i < 2 ; ++ i )
        { 
          for ( int j = 0 ; j < 2 - i ; ++ j )
             { 
               for ( int k = 0 ; k < 2 - i - j ; ++ k )
                  {
                    this->gaussianPoints [c] [0] = x [i] + x[j];
                    this->gaussianPoints [c] [1] = x [j] ;
                    this->gaussianPoints [c] [2] = x [k] ;
                    this->gaussianPoints [c++] [3] = w [0];
                  }  
             }
        }
     for ( int i = 0 ; i < 2 ; ++ i )
        { 
          for ( int j = 0 ; j < 2 - i ; ++ j )
             { 
               for ( int k = 0 ; k < 2 - i - j ; ++ k )
                  {
                    this->gaussianPoints [c] [0] = x [i] ;
                    this->gaussianPoints [c] [1] = x [j] + x[i];
                    this->gaussianPoints [c] [2] = x [k] ;
                    this->gaussianPoints [c++] [3] = w [0];
                  }  
             }
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss1 , pyramid >::Gaussian_Points()
   {
     this->stress [0] = 0.4 ; 
     this->stress [1] = 0.4 ;  
     this->stress [2] = 0.2 ;
     this->gaussianPoints [0] [0] = 0.5 ; 
     this->gaussianPoints [0] [1] = 0.25 ;
     this->gaussianPoints [0] [2] = 0.25 ;
     this->gaussianPoints [0] [3] = 1. / 6. ;
     this->gaussianPoints [1] [0] = 0.25 ; 
     this->gaussianPoints [1] [1] = 0.5 ;
     this->gaussianPoints [1] [2] = 0.25 ;
     this->gaussianPoints [1] [3] = 1. / 6. ;
  }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss3 , prism >::Gaussian_Points()
   {
     this->stress [0] = 1. / 3. ; 
     this->stress [1] = 1. / 3. ;  
     this->stress [2] = 1. / 2. ;
     TYPE w [6] = 
        {
          ( 155. - sqrt ( 15. ) ) / 2400. , 
          ( 155. + sqrt ( 15. ) ) / 2400. , 
          9. / 80.,
          5. / 18. , 
          4. / 9. , 
          5. / 18. 
        } ;
     TYPE x [8] = 
        {
          ( 6. - sqrt ( 15. ) ) / 21. , 
          ( 9. + 2. * sqrt ( 15. ) ) / 21. , 
          ( 6. + sqrt ( 15. ) ) / 21.,
          ( 9. - 2. * sqrt ( 15. ) ) / 21. , 
          1. / 3. ,
          0.5 + sqrt ( 0.15 ) , 
          0.5 , 
          0.5 - sqrt ( 0.15 ) 
        } ;
     for (int i= 0; i < 3; ++i)
        {
          this->gaussianPoints [0] [0] = x [0] ; 
          this->gaussianPoints [0] [1] = x [0] ; 
          this->gaussianPoints [0] [2] = x [5+i] ; 
          this->gaussianPoints [0] [3] = w [0] * w[4+i];
          this->gaussianPoints [1] [0] = x [1] ; 
          this->gaussianPoints [1] [1] = x [0] ; 
          this->gaussianPoints [1] [2] = x [5+i] ; 
          this->gaussianPoints [1] [3] = w [0] * w[4+i];
          this->gaussianPoints [2] [0] = x [0] ; 
          this->gaussianPoints [2] [1] = x [1] ; 
          this->gaussianPoints [2] [2] = x [5+i] ; 
          this->gaussianPoints [2] [3] = w [0] * w[4+i];
          this->gaussianPoints [3] [0] = x [2] ;  
          this->gaussianPoints [3] [1] = x [2] ; 
          this->gaussianPoints [3] [2] = x [5+i] ; 
          this->gaussianPoints [3] [3] = w [1] * w[4+i];
          this->gaussianPoints [4] [0] = x [3] ; 
          this->gaussianPoints [4] [1] = x [2] ; 
          this->gaussianPoints [4] [2] = x [5+i] ; 
          this->gaussianPoints [4] [3] = w [1] * w[4+i];
          this->gaussianPoints [5] [0] = x [2] ; 
          this->gaussianPoints [5] [1] = x [3] ; 
          this->gaussianPoints [5] [2] = x [5+i] ; 
          this->gaussianPoints [5] [3] = w [1] * w[4+i];
          this->gaussianPoints [6] [0] = x [4] ; 
          this->gaussianPoints [6] [1] = x [4] ; 
          this->gaussianPoints [6] [2] = x [5+i] ; 
          this->gaussianPoints [6] [3] = w [2] * w[4+i];  
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss2 , prism >::Gaussian_Points()
   {
     this->stress [0] = 1. / 3. ; 
     this->stress [1] = 1. / 3. ;  
     this->stress [2] = 1. / 2. ;
     TYPE x [4] = 
        { 
          1. / 6. , 
          4. / 6. ,
          0.5 * ( 1 + 1./sqrt(3.) ),
          0.5 * ( 1 - 1./sqrt(3.) )
        } ; 
     TYPE w = 1./12. ; 
     this->gaussianPoints [0] [0] = x [0] ;  
     this->gaussianPoints [0] [1] = x [0] ;  
     this->gaussianPoints [0] [2] = x [3] ;   
     this->gaussianPoints [0] [3] = w ;
     this->gaussianPoints [1] [0] = x [0] ;  
     this->gaussianPoints [1] [1] = x [1] ;  
     this->gaussianPoints [1] [2] = x [3] ;   
     this->gaussianPoints [1] [3] = w ;
     this->gaussianPoints [2] [0] = x [1] ;  
     this->gaussianPoints [2] [1] = x [0] ;  
     this->gaussianPoints [2] [2] = x [3] ;   
     this->gaussianPoints [2] [3] = w ;    
     this->gaussianPoints [3] [0] = x [0] ;  
     this->gaussianPoints [3] [1] = x [0] ;  
     this->gaussianPoints [3] [2] = x [4] ;   
     this->gaussianPoints [3] [3] = w ;
     this->gaussianPoints [4] [0] = x [0] ;  
     this->gaussianPoints [4] [1] = x [1] ;  
     this->gaussianPoints [4] [2] = x [4] ;   
     this->gaussianPoints [4] [3] = w ;
     this->gaussianPoints [5] [0] = x [1] ;  
     this->gaussianPoints [5] [1] = x [0] ;  
     this->gaussianPoints [5] [2] = x [4] ;   
     this->gaussianPoints [5] [3] = w ;    
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss1 , prism >::Gaussian_Points()
   {
     this->stress [0] = 1. / 3. ; 
     this->stress [1] = 1. / 3. ;  
     this->stress [2] = 1. / 2. ;
     this->gaussianPoints [0] [0] = 1. / 3. ;
     this->gaussianPoints [0] [1] = 1. / 3. ; 
     this->gaussianPoints [0] [2] = 1. / 2. ;
     this->gaussianPoints [0] [3] = 0.5 ; 
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss3 , quadrangle >::Gaussian_Points()
   {
     this->stress [0] = 0.5 ; 
     this->stress [1] = 0.5 ;  
     this->stress [2] = 0. ;
     TYPE x [3] = 
        {
          0.5 + sqrt ( 0.15 ) , 
          0.5 , 
          0.5 - sqrt ( 0.15 )
        } ;
     TYPE w [3] = 
        {
          5. / 18. , 
          4. / 9. , 
          5. / 18.
        } ;
     int c = 0 ;
     for ( int i = 0 ; i < 3 ; ++ i )
        { 
          for ( int j = 0 ; j < 3 ; ++ j )
             {
               this->gaussianPoints [c] [0] = x [i] ;
               this->gaussianPoints [c] [1] = x [j] ;
               this->gaussianPoints [c] [2] = 0. ;
               this->gaussianPoints [c++] [3] = w [i] * w [j] ;
             }  
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss2 , quadrangle >::Gaussian_Points()
   {
     this->stress [0] = 0.5 ; 
     this->stress [1] = 0.5 ;  
     this->stress [2] = 0. ;
     TYPE x [2] = 
        { 
          0.5 * ( 1. + 1.0 / sqrt ( 3. ) ) , 
          0.5 * ( 1. - 1.0 / sqrt ( 3. ) ) 
        } ;
     TYPE w [2] = 
        {
          1. / 2. , 
          1. / 2. 
        } ;
     int c = 0 ;
     for ( int i = 0 ; i < 2 ; ++ i )
        { 
          for ( int j = 0 ;j < 2 ; ++ j )
             {
               this->gaussianPoints [c] [0] = x [i] ;
               this->gaussianPoints [c] [1] = x [j] ;
               this->gaussianPoints [c] [2] = 0. ; 
               this->gaussianPoints [c++] [3] = w [i] * w [j] ;
             } 
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss1 , quadrangle >::Gaussian_Points()
   {
     this->stress [0] = 0.5 ; 
     this->stress [1] = 0.5 ;  
     this->stress [2] = 0. ;
     this->gaussianPoints [0] [0] = 0.5 ; 
     this->gaussianPoints [0] [1] = 0.5 ;
     this->gaussianPoints [0] [2] = 0. ; 
     this->gaussianPoints [0] [3] = 2. ;
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss3 , triangle >::Gaussian_Points()
   {
     this->stress [0] = 1. / 3. ; 
     this->stress [1] = 1. / 3. ;  
     this->stress [2] = 0. ;
     TYPE w [3] = 
        {
          ( 155. - sqrt ( 15. ) ) / 2400. , 
          ( 155. + sqrt ( 15. ) ) / 2400. , 
          9. / 80.
        } ;
     TYPE x [5] = 
        {
          ( 6. - sqrt ( 15. ) ) / 21. , 
          ( 9. + 2. * sqrt ( 15. ) ) / 21. , 
          ( 6. + sqrt ( 15. ) ) / 21.,
          ( 9. - 2. * sqrt ( 15. ) ) / 21. , 
          1. / 3. 
        } ;
     this->gaussianPoints [0] [0] = x [0] ; 
     this->gaussianPoints [0] [1] = x [0] ; 
     this->gaussianPoints [0] [2] = 0.; 
     this->gaussianPoints [0] [3] = w [0] ;
     this->gaussianPoints [1] [0] = x [1] ; 
     this->gaussianPoints [1] [1] = x [0] ; 
     this->gaussianPoints [1] [2] = 0. ; 
     this->gaussianPoints [1] [3] = w [0] ;
     this->gaussianPoints [2] [0] = x [0] ; 
     this->gaussianPoints [2] [1] = x [1] ; 
     this->gaussianPoints [2] [2] = 0. ; 
     this->gaussianPoints [2] [3] = w [0] ;
     this->gaussianPoints [3] [0] = x [2] ;  
     this->gaussianPoints [3] [1] = x [2] ; 
     this->gaussianPoints [3] [2] = 0. ; 
     this->gaussianPoints [3] [3] = w [1] ;
     this->gaussianPoints [4] [0] = x [3] ; 
     this->gaussianPoints [4] [1] = x [2] ; 
     this->gaussianPoints [4] [2] = 0. ; 
     this->gaussianPoints [4] [3] = w [1] ;
     this->gaussianPoints [5] [0] = x [2] ; 
     this->gaussianPoints [5] [1] = x[3] ; 
     this->gaussianPoints [5] [2] = 0. ; 
     this->gaussianPoints [5] [3] = w [1] ;
     this->gaussianPoints [6] [0] = x [4] ; 
     this->gaussianPoints [6] [1] = x [4] ; 
     this->gaussianPoints [6] [2] = 0. ; 
     this->gaussianPoints [6] [3] = w [2] ;  
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss2 , triangle >::Gaussian_Points()
   {
     this->stress [0] = 1. / 3. ; 
     this->stress [1] = 1. / 3. ;  
     this->stress [2] = 0. ;
     TYPE x [2] = 
        { 
          1. / 6. , 
          4. / 6. 
        } ; 
     TYPE w = 1./6. ; 
     this->gaussianPoints [0] [0] = x [0] ;  
     this->gaussianPoints [0] [1] = x [0] ;  
     this->gaussianPoints [0] [2] = 0. ;   
     this->gaussianPoints [0] [3] = w ;
     this->gaussianPoints [1] [0] = x [0] ;  
     this->gaussianPoints [1] [1] = x [1] ;  
     this->gaussianPoints [1] [2] = 0. ;   
     this->gaussianPoints [1] [3] = w ;
     this->gaussianPoints [2] [0] = x [1] ;  
     this->gaussianPoints [2] [1] = x [0] ;  
     this->gaussianPoints [2] [2] = 0. ;   
     this->gaussianPoints [2] [3] = w ;    
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss1 , triangle >::Gaussian_Points()
   {
     this->stress [0] = 1. / 3. ; 
     this->stress [1] = 1. / 3. ;  
     this->stress [2] = 0. ;
     this->gaussianPoints [0] [0] = 1. / 3. ;
     this->gaussianPoints [0] [1] = 1. / 3. ; 
     this->gaussianPoints [0] [2] = 0. ;
     this->gaussianPoints [0] [3] = 0.5 ; 
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss3 , interval >::Gaussian_Points()
   {
     this->stress [0] = 0.5 ; 
     this->stress [1] = 0. ;  
     this->stress [2] = 0. ;
     TYPE x [3] = 
        {
          0.5 + sqrt ( 0.15 ) , 
          0.5 , 
          0.5 - sqrt ( 0.15 ) 
        } ;
     TYPE w [3] = 
        {
          5. / 18. , 
          4. / 9. , 
          5. / 18. 
        } ;
     for ( int i = 0 ; i < 3 ; ++ i )
        {
          this->gaussianPoints [i] [0] = x [i] ;
          this->gaussianPoints [i] [1] = 0. ;
          this->gaussianPoints [i] [2] = 0. ;
          this->gaussianPoints [i] [3] = w [i] ;
        }   
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
Gaussian_Points<TYPE, Gauss2 , interval >::Gaussian_Points()
   {
     this->stress [0] = 0.5 ; 
     this->stress [1] = 0. ;  
     this->stress [2] = 0. ;
     TYPE x [2] = 
        {
          0.5 * ( 1. + 1.0 / sqrt ( 3. ) ) , 
          0.5 * ( 1. - 1.0 / sqrt ( 3. ) ) 
        } ;
     TYPE w [2] = 
        {
          0.5, 
          0.5
        } ;
     for ( int i = 0 ; i < 2 ; ++ i )
        {
          this->gaussianPoints [i] [0] = x [i] ; 
          this->gaussianPoints [i] [1] = 0. ; 
          this->gaussianPoints [i] [2] = 0. ; 
          this->gaussianPoints [i] [3] = w [i] ;
        }  
   }
//-----------------------------------------------------------------------------
template <typename TYPE> 
 Gaussian_Points<TYPE, Gauss1 , interval >::Gaussian_Points()
   {
     this->stress [0] = 0.5 ; 
     this->stress [1] = 0. ;  
     this->stress [2] = 0. ;
     this->gaussianPoints [0] [0] = 0.5 ; 
     this->gaussianPoints [0] [1] = 0. ; 
     this->gaussianPoints [0] [2] = 0. ;
     this->gaussianPoints [0] [3] = 1. ;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
CornerManagement<TYPE,dimension,numberCorners>::CornerManagement()
   {
     enum { sizeOfBaryArray = numberCorners + 5,
            sizeOfPointArray = (numberCorners+1)*dimension };
     barys = new TYPE [sizeOfBaryArray] ;
     coordinatesOfVertices = new TYPE [sizeOfPointArray];
#ifdef INITIALIZE_ARRAYS
     for (int k=0; k<sizeOfBaryArray; ++k)
         {
           barys[k] = 0.;
         }
     for (int k=0; k<sizeOfPointArray; ++k)
         {
           coordinatesOfVertices[k] = 0.;
         }
#endif
   }
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
CornerManagement<TYPE,dimension,numberCorners>::~CornerManagement()
   {
     delete [] barys;
     delete [] coordinatesOfVertices;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
inline TYPE* 
CornerManagement<TYPE,dimension,numberCorners>::getCornerList() const
   {
     return coordinatesOfVertices;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
inline TYPE* 
CornerManagement<TYPE,dimension,numberCorners>::getArrayOfBarycentricCoordinates() const
   {
     return barys;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
inline TYPE
CornerManagement <TYPE,dimension,numberCorners> :: edge_length (const unsigned int cnr_1, const unsigned int cnr_2) const
   {
     TYPE result = pow(this->coordinatesOfVertices[cnr_1*dimension] - this->coordinatesOfVertices[cnr_2*dimension],2);
     for (int i=1; i< dimension; ++i)
        {
          result += pow(this->coordinatesOfVertices[cnr_1*dimension+i] - this->coordinatesOfVertices[cnr_2*dimension+i],2);
        }
     return sqrt( result );
   }
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
inline TYPE 
CornerManagement<TYPE,dimension,numberCorners>::edge_(SpaceDirection direction) const
   {
     return coordinatesOfVertices[direction+dimension] - coordinatesOfVertices[direction];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberCorners>
template <class Element>
inline void
Corner_Classes <TYPE,D3,numberCorners,cuboid> :: compute_Barycentric_Coordinates (const Element& element) 
   {
     this->barys [4] =
        PyramidVolume<1,2,6,5,0>::compute(element) *
        PyramidVolume<4,5,6,7,0>::compute(element) *
        PyramidVolume<2,6,7,3,0>::compute(element) ;
     this->barys [5] =
        PyramidVolume<0,3,7,4,1>::compute(element) *
        PyramidVolume<2,6,7,3,1>::compute(element) *
        PyramidVolume<4,5,6,7,1>::compute(element) ;
     this->barys [6] =
        PyramidVolume<0,3,7,4,2>::compute(element) *
        PyramidVolume<0,1,5,4,2>::compute(element) *
        PyramidVolume<4,5,6,7,2>::compute(element) ;
     this->barys [7] =
        PyramidVolume<1,2,6,5,3>::compute(element) *
        PyramidVolume<0,1,5,4,3>::compute(element) *
        PyramidVolume<4,5,6,7,3>::compute(element) ;
     this->barys [8] =
        PyramidVolume<1,2,6,5,4>::compute(element) *
        PyramidVolume<2,6,7,3,4>::compute(element) *
        PyramidVolume<0,1,2,3,4>::compute(element) ;
     this->barys [9] =
        PyramidVolume<0,3,7,4,5>::compute(element) *
        PyramidVolume<2,6,7,3,5>::compute(element) *
        PyramidVolume<0,1,2,3,5>::compute(element) ;
     this->barys [10] =
        PyramidVolume<0,3,7,4,6>::compute(element) *
        PyramidVolume<0,1,5,4,6>::compute(element) *
        PyramidVolume<0,1,2,3,6>::compute(element) ;
     this->barys [11] =
        PyramidVolume<1,2,6,5,7>::compute(element) *
        PyramidVolume<0,1,5,4,7>::compute(element) *
        PyramidVolume<0,1,2,3,7>::compute(element) ;

     this->barys [dirX] =  this->barys [5] + this->barys [6] + this->barys [9] + this->barys [10] ;
     this->barys [dirY] =  this->barys [6] + this->barys [7] + this->barys [10] + this->barys [11] ;
     this->barys [dirZ] =  this->barys [8] + this->barys [9] + this->barys [10] + this->barys [11] ;
     this->barys [weight] =  1. ;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberCorners>
template <class Element>
inline void
Corner_Classes <TYPE,D3,numberCorners,pyramid> :: compute_Barycentric_Coordinates (const Element& element) 
   {
    this->barys [4] =  TetrahedronVolume<1,2,4,0>::compute(element) * TetrahedronVolume<2,3,4,0>::compute(element) ;
    this->barys [5] =  TetrahedronVolume<2,3,4,1>::compute(element) * TetrahedronVolume<3,4,0,1>::compute(element) ;
    this->barys [6] =  TetrahedronVolume<3,4,0,2>::compute(element) * TetrahedronVolume<0,1,4,2>::compute(element) ;
    this->barys [7] =  TetrahedronVolume<0,1,4,3>::compute(element) * TetrahedronVolume<1,2,4,3>::compute(element) ;
    this->barys [8] =  PyramidVolume<0,1,2,3,4>::compute(element) ;

    this->barys[0] =  this->barys[5] + this->barys[6];
    this->barys[1] =  this->barys[6] + this->barys[7];
    this->barys[2] =  this->barys[8];
    this->barys[3] =  1.;
    this->barys[3] =  1.;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberCorners>
template <class Element>
inline void
Corner_Classes <TYPE,D3,numberCorners,prism> :: compute_Barycentric_Coordinates (const Element& element) 
   {
     this->barys [4] =  PyramidVolume<1,4,5,2,0>::compute(element) * TetrahedronVolume<3,4,5,0>::compute(element) ;
     this->barys [5] =  PyramidVolume<0,2,5,3,1>::compute(element) * TetrahedronVolume<3,4,5,1>::compute(element) ;
     this->barys [6] =  PyramidVolume<1,0,3,4,2>::compute(element) * TetrahedronVolume<3,4,5,2>::compute(element) ;
     this->barys [7] =  PyramidVolume<1,4,5,2,3>::compute(element) * TetrahedronVolume<0,1,2,3>::compute(element) ;
     this->barys [8] =  PyramidVolume<0,2,5,3,4>::compute(element) * TetrahedronVolume<0,1,2,4>::compute(element) ;
     this->barys [9] =  PyramidVolume<1,0,3,4,5>::compute(element) * TetrahedronVolume<0,1,2,5>::compute(element) ;

     this->barys [0] =  this->barys[5] + this->barys[8];
     this->barys [1] =  this->barys[6] + this->barys[9];
     this->barys [2] =  this->barys[7] + this->barys[8] + this->barys[9];
     this->barys [3] =  1.;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberCorners>
template <class Element>
inline void
Corner_Classes <TYPE,D3,numberCorners,tetrahedron> :: compute_Barycentric_Coordinates (const Element& element) 
   {
    this->barys [4] =  TetrahedronVolume<1,2,3,0>::compute(element) ;
    this->barys [5] =  TetrahedronVolume<0,2,3,1>::compute(element) ;
    this->barys [6] =  TetrahedronVolume<0,1,3,2>::compute(element) ;
    this->barys [7] =  TetrahedronVolume<0,1,2,3>::compute(element) ;

    this->barys[0] =  this->barys[5];
    this->barys[1] =  this->barys[6];
    this->barys[2] =  this->barys[7];
    this->barys[3] =  1.;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
template <class Element>
inline void
Corner_Classes <TYPE,dimension,numberCorners,quadrangle> :: compute_Barycentric_Coordinates (const Element& element) 
   {
      assert (dimension == D2);
      this->barys [5] = TriangleVolume<1,2,0>::compute(element) * TriangleVolume<2,3,0>::compute(element) ;
      this->barys [5] = TriangleVolume<2,3,1>::compute(element) * TriangleVolume<3,0,1>::compute(element) ;
      this->barys [6] = TriangleVolume<0,1,2>::compute(element) * TriangleVolume<3,0,2>::compute(element) ;
      this->barys [7] = TriangleVolume<0,1,3>::compute(element) * TriangleVolume<1,2,3>::compute(element) ;

      this->barys [0] =  this->barys [5] + this->barys [6] ;
      this->barys [1] =  this->barys [7] + this->barys [6] ;
      this->barys [2] =  0. ;
      this->barys [3] =  1. ;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
template <class Element>
inline void
Corner_Classes <TYPE,dimension,numberCorners,triangle> :: compute_Barycentric_Coordinates (const Element& element) 
   {
     assert (dimension == D2);
     this->barys [4] =  TriangleVolume<1,2,0>::compute(element) ;
     this->barys [5] =  TriangleVolume<0,2,1>::compute(element) ;
     this->barys [6] =  TriangleVolume<0,1,2>::compute(element) ;
     this->barys [0] =  this->barys[5];
     this->barys [1] =  this->barys[6];
     this->barys [2] =  0.;
     this->barys [3] =  1.;
   }
//==============================================================================
