//==============================================================================
//
//  $Id: Mapping.C,v 1.23 2006/07/25 13:52:27 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension_, int _numberGaussianPoints>
MappingManagement<TYPE,dimension_,_numberGaussianPoints>::MappingManagement() 
   {
     mappingValues = new TYPE**[_numberGaussianPoints+1];
     extensionValues = new TYPE*[_numberGaussianPoints+1];
     for (int i=0; i < _numberGaussianPoints+1; ++i)
        {
          mappingValues[i] = new TYPE* [dimension_];
          for (int j=0; j < dimension_; ++j)
             {
               mappingValues[i][j] = new TYPE [dimMap];
#ifdef INITIALIZE_ARRAYS
               for (int k=0; k < dimMap; ++k)
                  {
                    mappingValues[i][j][k] = 0.;
                  }
#endif
             }
          extensionValues[i] = new TYPE [dimExt];
#ifdef INITIALIZE_ARRAYS
          for (int k=0; k < dimExt; ++k)
             {
               extensionValues[i][k] = 0.;
             }
#endif
        }
     jacobianMatrixEntries = new TYPE[(_numberGaussianPoints+1)*dimVal];
#ifdef INITIALIZE_ARRAYS
     for (int k=0; k < (_numberGaussianPoints+1)*dimVal; ++k)
        {
          jacobianMatrixEntries[k] = 0.;
        }
#endif
   }
// ---------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension_, int _numberGaussianPoints>
MappingManagement<TYPE,dimension_,_numberGaussianPoints>::~MappingManagement()
   {
     for (int i=0; i < _numberGaussianPoints+1; ++i)
        {
          for (int j=0; j < dimension_; ++j)
           {
             delete [] mappingValues[i][j];
           }
          delete [] mappingValues[i];
          delete [] extensionValues[i];
        }
     delete [] mappingValues;
     delete [] jacobianMatrixEntries;
     delete [] extensionValues;
   }
// ---------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension_, int _numberGaussianPoints>
inline TYPE**
MappingManagement<TYPE,dimension_,_numberGaussianPoints>:: getMappingValues(int gaussLevel) const
    {
      return mappingValues[gaussLevel];
    }
// ---------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension_, int _numberGaussianPoints>
inline TYPE*
MappingManagement<TYPE,dimension_,_numberGaussianPoints>:: getValuesOfJacobianMatrix(int gaussLevel) const
    {
      return jacobianMatrixEntries+(gaussLevel*dimVal);
    }
// ---------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension_, int _numberGaussianPoints>
inline const TYPE&
MappingManagement<TYPE,dimension_,_numberGaussianPoints>:: getDeterminant(int gaussLevel) const
    {
      return extensionValues[gaussLevel][determinant];
    }
// ---------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension_, int _numberGaussianPoints>
inline const TYPE&
MappingManagement<TYPE,dimension_,_numberGaussianPoints>:: getNormalLength(int gaussLevel) const
    {
      return extensionValues[gaussLevel][normalLength];
    }
// ---------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension_, int _numberGaussianPoints>
inline const TYPE&
MappingManagement<TYPE,dimension_,_numberGaussianPoints>:: getTransformationFactorTimesWeight(int gaussLevel) const
    {
      return extensionValues[gaussLevel][trafoFactor];
    }
//-----------------------------------------------------------------------------
template <typename TYPE, class TRAFO, int _numberGaussianPoints>
template <class Element>
void
Mappings<TYPE,D3,interior,TRAFO,_numberGaussianPoints>::point(const Element& element, const unsigned int gaussLevel)
  { 
    TYPE* val_ = this->getValuesOfJacobianMatrix(gaussLevel);
    TYPE* map_dirX = this->getMappingValues(gaussLevel)[dirX];
    TYPE* map_dirY = this->getMappingValues(gaussLevel)[dirY];
    TYPE* map_dirZ = this->getMappingValues(gaussLevel)[dirZ];
    map_dirX[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirX>(element,gaussLevel);
    map_dirX[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirX>(element,gaussLevel);
    map_dirX[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirX>(element,gaussLevel);
    map_dirX[dz] = trafo.template evaluateGaussianPoints<0,0,1,dirX>(element,gaussLevel);
    map_dirY[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirY>(element,gaussLevel);
    map_dirY[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirY>(element,gaussLevel);
    map_dirY[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirY>(element,gaussLevel);
    map_dirY[dz] = trafo.template evaluateGaussianPoints<0,0,1,dirY>(element,gaussLevel);
    map_dirZ[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirZ>(element,gaussLevel);
    map_dirZ[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirZ>(element,gaussLevel);
    map_dirZ[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirZ>(element,gaussLevel);
    map_dirZ[dz] = trafo.template evaluateGaussianPoints<0,0,1,dirZ>(element,gaussLevel);
 // Precomputations concerning the derivative with respect to x
    val_[0] = map_dirZ[dz] * map_dirY[dy] - map_dirY[dz] * map_dirZ[dy];
    val_[1] = map_dirZ[dx] * map_dirY[dz] - map_dirY[dx] * map_dirZ[dz];
    val_[2] = map_dirZ[dy] * map_dirY[dx] - map_dirY[dy] * map_dirZ[dx];
 // Precomputations concerning the derivative with respect to y
    val_[3] = map_dirZ[dy] * map_dirX[dz] - map_dirZ[dz] * map_dirX[dy];
    val_[4] = map_dirZ[dz] * map_dirX[dx] - map_dirZ[dx] * map_dirX[dz];
    val_[5] = map_dirZ[dx] * map_dirX[dy] - map_dirZ[dy] * map_dirX[dx];
 // Precomputations concerning the derivative with respect to z
    val_[6] = map_dirY[dz] * map_dirX[dy] - map_dirY[dy] * map_dirX[dz];
    val_[7] = map_dirY[dx] * map_dirX[dz] - map_dirY[dz] * map_dirX[dx];
    val_[8] = map_dirY[dy] * map_dirX[dx] - map_dirY[dx] * map_dirX[dy];

 /* The following is the determinant, just written in a different manner. 
    But in principle just the rule of Sarrus, only factored out the terms 
    conerning the x direction. The advantage is, we achieve a formula, 
    that used parts that are already computed!
 */
    this->extensionValues[gaussLevel][determinant] =
      Colsamm_Internal_Functions::ABS<TYPE>::abs_(map_dirX[dx] * val_[0] + map_dirX[dy] * val_[1] + map_dirX[dz] * val_[2] );
     assert ( 1.e-10 < Colsamm_Internal_Functions::ABS<TYPE>::abs_(this->extensionValues[gaussLevel][determinant])) ;
    this->extensionValues[gaussLevel][trafoFactor] = this->extensionValues[gaussLevel][determinant];
  }
//-----------------------------------------------------------------------------
template <typename TYPE, class TRAFO, int _numberGaussianPoints>
template <class Element>
inline void 
Mappings<TYPE,D2,interior,TRAFO,_numberGaussianPoints>::point(const Element& element, const unsigned int gaussLevel )
    {
      TYPE* val_ = this->getValuesOfJacobianMatrix(gaussLevel);
      TYPE* map_dirX = this->getMappingValues(gaussLevel)[dirX];
      TYPE* map_dirY = this->getMappingValues(gaussLevel)[dirY];
   // Auswertungen fr die Transformation
      map_dirX[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirX>(element,gaussLevel);
      map_dirX[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirX>(element,gaussLevel);
      map_dirX[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirX>(element,gaussLevel);
      map_dirY[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirY>(element,gaussLevel);
      map_dirY[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirY>(element,gaussLevel);
      map_dirY[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirY>(element,gaussLevel);
   // Determinante, Determinante * Gewicht fr Integration, Ableitungen  
      val_[0] = map_dirY[dy];
      val_[1] = -map_dirY[dx];
      val_[2] = -map_dirX[dy];
      val_[3] = map_dirX[dx];
      this->extensionValues[gaussLevel][determinant] =
                          Colsamm_Internal_Functions::ABS<TYPE>::abs_( map_dirX[dx]*map_dirY[dy] - map_dirX[dy]*map_dirY[dx]);
      this->extensionValues[gaussLevel][trafoFactor] = this->extensionValues[gaussLevel][determinant];
    }
//-----------------------------------------------------------------------------
template <typename TYPE, class TRAFO, int _numberGaussianPoints>
template <class Element>
inline void 
Mappings<TYPE,D1,interior,TRAFO,_numberGaussianPoints>::point ( const Element& element, const unsigned int gaussLevel )
   {
     TYPE* val_ = this->getValuesOfJacobianMatrix(gaussLevel);
     TYPE* map_dirX = this->getMappingValues(gaussLevel)[dirX];
     map_dirX[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirX>(element,gaussLevel);
     map_dirX[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirX>(element,gaussLevel);
     val_[0] = 1.;
     std::cout << val_[0] << std::endl;
     this->extensionValues[gaussLevel][determinant] = Colsamm_Internal_Functions::ABS<TYPE>::abs_(map_dirX[dx]);
     this->extensionValues[gaussLevel][trafoFactor] = this->extensionValues[gaussLevel][determinant];
  }
// -----------------------------------------------------------------------------
template <typename TYPE, class TRAFO, int _numberGaussianPoints>
template <class Element>
inline void 
Mappings<TYPE,D2,boundary,TRAFO,_numberGaussianPoints>::point ( const Element& element, const unsigned int gaussLevel )
   {
     TYPE* map_dirX = this->getMappingValues(gaussLevel)[dirX];
     TYPE* map_dirY = this->getMappingValues(gaussLevel)[dirY];
     TYPE* val_ = this->getValuesOfJacobianMatrix(gaussLevel);
     map_dirX[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirX>(element,gaussLevel);
     map_dirX[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirX>(element,gaussLevel);
     map_dirX[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirX>(element,gaussLevel);
     map_dirY[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirY>(element,gaussLevel);
     map_dirY[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirY>(element,gaussLevel);
     map_dirY[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirY>(element,gaussLevel);
     val_[0] = map_dirY[dy];
     val_[1] = -map_dirY[dx];
     val_[2] = -map_dirX[dy];
     val_[3] = map_dirX[dx];
#if 0
     for (int i=0; i < 3; ++i)
        std::cout << "map_dirX["<< i<< "] :" << map_dirX[i] << " ";
     std::cout << std::endl;
     for (int i=0; i < 3; ++i)
        std::cout << "map_dirY["<< i<< "] :" << map_dirY[i] << " ";
     std::cout << std::endl;
#endif
     this->extensionValues[gaussLevel][determinant] =      
                     Colsamm_Internal_Functions::ABS<TYPE>::abs_( map_dirX[dx]*map_dirY[dy] - map_dirX[dy]*map_dirY[dx]);
     this->extensionValues[gaussLevel][normalLength] = sqrt(map_dirX[dx]*map_dirX[dx] + map_dirY[dx]*map_dirY[dx]);
     assert (1.e-16 <  fabs(this->extensionValues[gaussLevel][normalLength] ) );
 //  std::cout << "Length of Normal : " <<  this->extensionValues[gaussLevel][normalLength] << std::endl;
     this->extensionValues[gaussLevel][determinant] = this->extensionValues[gaussLevel][normalLength];
     this->extensionValues[gaussLevel][trafoFactor] = this->extensionValues[gaussLevel][normalLength];
   }
// -----------------------------------------------------------------------------
template <typename TYPE, class TRAFO, int _numberGaussianPoints>
template <class Element>
inline void 
Mappings<TYPE,D3,boundary,TRAFO,_numberGaussianPoints>::point ( const Element& element, const unsigned int gaussLevel )
   {
     TYPE* val_ = this->getValuesOfJacobianMatrix(gaussLevel);
     TYPE* map_dirX = this->getMappingValues(gaussLevel)[dirX];
     TYPE* map_dirY = this->getMappingValues(gaussLevel)[dirY];
     TYPE* map_dirZ = this->getMappingValues(gaussLevel)[dirZ];
      map_dirX[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirX>(element,gaussLevel);
      map_dirX[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirX>(element,gaussLevel);
      map_dirX[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirX>(element,gaussLevel);
      map_dirX[dz] = trafo.template evaluateGaussianPoints<0,0,1,dirX>(element,gaussLevel); 
      map_dirY[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirY>(element,gaussLevel);
      map_dirY[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirY>(element,gaussLevel);
      map_dirY[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirY>(element,gaussLevel);
      map_dirY[dz] = trafo.template evaluateGaussianPoints<0,0,1,dirY>(element,gaussLevel);
      map_dirZ[fn] = trafo.template evaluateGaussianPoints<0,0,0,dirZ>(element,gaussLevel);
      map_dirZ[dx] = trafo.template evaluateGaussianPoints<1,0,0,dirZ>(element,gaussLevel);
      map_dirZ[dy] = trafo.template evaluateGaussianPoints<0,1,0,dirZ>(element,gaussLevel);
      map_dirZ[dz] = trafo.template evaluateGaussianPoints<0,0,1,dirZ>(element,gaussLevel);
  // Computations for the derivetion with respect to x
      val_[0] = map_dirZ[dz] * map_dirY[dy] - map_dirY[dz] * map_dirZ[dy];
      val_[1] = map_dirZ[dx] * map_dirY[dz] - map_dirY[dx] * map_dirZ[dz];
      val_[2] = map_dirZ[dy] * map_dirY[dx] - map_dirY[dy] * map_dirZ[dx];
  // Computations for the derivetion with respect to y
      val_[3] = map_dirZ[dy] * map_dirX[dz] - map_dirZ[dz] * map_dirX[dy];
      val_[4] = map_dirZ[dz] * map_dirX[dx] - map_dirZ[dx] * map_dirX[dz];
      val_[5] = map_dirZ[dx] * map_dirX[dy] - map_dirZ[dy] * map_dirX[dx];
  // Computations 6or the derivetion with respect to z
      val_[6] = map_dirY[dz] * map_dirX[dy] - map_dirY[dy] * map_dirX[dz];
      val_[7] = map_dirY[dx] * map_dirX[dz] - map_dirY[dz] * map_dirX[dx];
      val_[8] = map_dirY[dy] * map_dirX[dx] - map_dirY[dx] * map_dirX[dy];
      this->extensionValues[gaussLevel][determinant] =
        Colsamm_Internal_Functions::ABS<TYPE>::abs_(map_dirX[dx] * val_[0] + map_dirX[dy] * val_[1] + map_dirX[dz] * val_[2] );
   // assert ( 1.e-10 < Colsamm_Internal_Functions::ABS<TYPE>::abs_(this->extensionValues[gaussLevel][determinant])) ;
#if 0
      this->extensionValues[gaussLevel][normalLength] = sqrt(
         pow(map_dirY[dx]*map_dirZ[dy]-map_dirZ[dx]*map_dirY[dy],2) +
         pow(map_dirZ[dx]*map_dirX[dy]-map_dirZ[dy]*map_dirX[dx],2) +
         pow(map_dirX[dx]*map_dirY[dy]-map_dirX[dy]*map_dirY[dx],2) );
#else
      this->extensionValues[gaussLevel][normalLength] = 
           sqrt( pow(val_[2],2) + pow(val_[5],2) + pow(val_[8],2) );
#endif
     this->extensionValues[gaussLevel][trafoFactor] = this->extensionValues[gaussLevel][normalLength];
     std::cout << "Normale: (" << val_[2] << ", " << val_[5] << ", " << val_[8] << ")" << std::endl;
    }
// -----------------------------------------------------------------------------
/* -----------------------------------------------------------------------------
   Comments conocerning evaluateGaussianPointsutaion for boundary derivatives:
     * since the parametrization is defined on a two dimensional element, we 
       need to pass the Gaussian points for the corresponding two dimensional
       reference element. However, by our definition, the Gaussian points are 
       always three dimensional, with components 0 of those dimensions not used, 
       we can just use the common Gaussian points.
     * the components of the normal are stored in val[l][4], val[l][7] and 
       val[l][10]; this can be verified very easily!
     * the length of the normal is multiplied to the Gaussian weight. 
     * val[l][1] is used as factor in at the end of the integral for specific 
       Gaussian point l. Thus, it has to store the length of the normal instead 
       of the determinant!
     * the determinant becomes stored in val[l][15]
  
   -----------------------------------------------------------------------------
*/
// -----------------------------------------------------------------------------
template <typename TYPE, class TRAFO, int _numberGaussianPoints>
template <class Element>
void 
Mappings<TYPE,D3,boundaryPlain,TRAFO,_numberGaussianPoints>::point(const Element& element, const unsigned int gaussLevel)
   {
     TYPE* val_ = this->getValuesOfJacobianMatrix(gaussLevel);
     TYPE* map_dirX = this->getMappingValues(gaussLevel)[dirX];
     TYPE* map_dirY = this->getMappingValues(gaussLevel)[dirY];
     TYPE* map_dirZ = this->getMappingValues(gaussLevel)[dirZ];
     map_dirX[fn] = TRAFO :: template evaluateGaussianPoints <0,0,0,dirX>(element,gaussLevel);
     map_dirX[dx] = TRAFO :: template evaluateGaussianPoints <1,0,0,dirX>(element,gaussLevel);
     map_dirX[dy] = TRAFO :: template evaluateGaussianPoints <0,1,0,dirX>(element,gaussLevel);
     map_dirX[dz] = TRAFO :: template evaluateGaussianPoints <0,0,1,dirX>(element,gaussLevel);
     map_dirY[fn] = TRAFO :: template evaluateGaussianPoints <0,0,0,dirY>(element,gaussLevel);
     map_dirY[dx] = TRAFO :: template evaluateGaussianPoints <1,0,0,dirY>(element,gaussLevel);
     map_dirY[dy] = TRAFO :: template evaluateGaussianPoints <0,1,0,dirY>(element,gaussLevel);
     map_dirY[dz] = TRAFO :: template evaluateGaussianPoints <0,0,1,dirY>(element,gaussLevel);
     map_dirZ[fn] = TRAFO :: template evaluateGaussianPoints <0,0,0,dirZ>(element,gaussLevel);
     map_dirZ[dx] = TRAFO :: template evaluateGaussianPoints <1,0,0,dirZ>(element,gaussLevel);
     map_dirZ[dy] = TRAFO :: template evaluateGaussianPoints <0,1,0,dirZ>(element,gaussLevel);
     map_dirZ[dz] = TRAFO :: template evaluateGaussianPoints <0,0,1,dirZ>(element,gaussLevel);
  // Computations for the derivetion with respect to x
     val_[0] = 0.;
     val_[1] = 0.;
     val_[2] = map_dirZ[dy] * map_dirY[dx] - map_dirY[dy] * map_dirZ[dx];
  // Computations for the derivetion with respect to y
     val_[3] = 0.;
     val_[4] = 0.;
     val_[5] = map_dirZ[dx] * map_dirX[dy] - map_dirZ[dy] * map_dirX[dx];
  // Computations 6or the derivetion with respect to z
     val_[6] = 0.;
     val_[7] = 0.;
     val_[8] = map_dirY[dy] * map_dirX[dx] - map_dirY[dx] * map_dirX[dy];
  // the length of the normal 
     this->extensionValues[gaussLevel][determinant] = 1.;
  // Gaussion weight times length of the normal (parametrization)
  /* The determinant (needed for the computations cencerning the inverse of the 3d onto mapping.
     The somehow strange looking formula is correct in this manner. We modified the rule of Sarrus
     by factoring out the compunents PA::map[l][0][*], since the other terms are already computed!
  */
      assert ( 1.e-10 < Colsamm_Internal_Functions::ABS<TYPE>::abs_(this->extensionValues[gaussLevel][determinant])) ;
#if 1
      this->extensionValues[gaussLevel][normalLength] = sqrt(
         pow(map_dirY[dx]*map_dirZ[dy]-map_dirZ[dx]*map_dirY[dy],2) +
         pow(map_dirZ[dx]*map_dirX[dy]-map_dirZ[dy]*map_dirX[dx],2) +
         pow(map_dirX[dx]*map_dirY[dy]-map_dirX[dy]*map_dirY[dx],2) );
#else
      this->extensionValues[gaussLevel][normalLength] = 
           sqrt( pow(val_[2],2) + pow(val_[5],2) + pow(val_[8],2) );
#endif
     this->extensionValues[gaussLevel][trafoFactor] = this->extensionValues[gaussLevel][normalLength];
   }
//==============================================================================
