#define BASISFUNCTION_TEMPLATES \
  template <typename TYPE, \
            DIMENSION dim_,\
            int numberGaussianPoints_plus1,\
            int numberSet1_,\
            int numberSet2_,\
            int numberOfBasisFunctionSets>

#define BASISFUNCTION_TYPE\
   Basis_Functions <\
        TYPE,\
        dim_,\
        numberGaussianPoints_plus1,\
        oneDimensional,\
        numberSet1_,\
        numberSet2_,\
        numberOfBasisFunctionSets>
//-----------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
inline
BASISFUNCTION_TYPE :: Basis_Functions () 
   : numberOfBasisFunctionSets_(numberOfBasisFunctionSets),
#ifdef SECOND_DERIVATIVES
     sizeOfValues(1+dim_+(dim_*(dim_+1))/2)
#else 
     sizeOfValues(1+dim_)
#endif
   {
     numberSet[0] = numberSet1_;
     if (1 < numberOfBasisFunctionSets)
       {
         numberSet[1] = numberSet2_;
       }
  // std::cout << " numberSet[0] " << numberSet[0] << "   numberSet[1]"   << numberSet[1] << std::endl;
     int totalSize = 0;
     for (int s=0; s < numberOfBasisFunctionSets; ++s)
        {
          totalSize += numberSet[s];
        }
     int sizeOf_basisFunctionsAtGaussianPoints = numberGaussianPoints_plus1*sizeOfValues*(totalSize),
         sizeOf_derivativeValues = dim_*numberOfBasisFunctionSets;
  // std::cout << "sizeOf_derivativeValues" << sizeOf_derivativeValues << std::endl;
     basisFunctionsAtGaussianPoints = new TYPE[sizeOf_basisFunctionsAtGaussianPoints];
     derivativeValues = new TYPE**[sizeOf_derivativeValues];
     int i = 0;
     for (int s=0; s < numberOfBasisFunctionSets; ++s)
        {
          for (; i < dim_*(s+1); ++i)
             {
                    derivativeValues[i] = new TYPE* [numberSet[s]];
                     for (int j=0; j < numberSet[s]; ++j)
                        {
                         derivativeValues[i][j] = new TYPE [numberGaussianPoints_plus1];
#ifdef INITIALIZE_ARRAYS
                         for (int k=0; k < numberGaussianPoints_plus1; ++k)
                            {
                              derivativeValues[i][j][k] = 0.;
                            }
#endif
                        }
             }
        }
   }
//-----------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
template <typename baseTYPE>
inline void 
BASISFUNCTION_TYPE :: eval(baseTYPE**_gaussianPoints) 
   {
     baseTYPE* gaussianPoints_ = _gaussianPoints[numberGaussianPoints_plus1-1];
     TYPE* basisFunctionsAtGaussianPoints_ = NULL;
     for (int numberSetsIter = 0 ; numberSetsIter < numberOfBasisFunctionSets ; ++numberSetsIter )
        {
          for (unsigned int functionSetIter = 0 ; functionSetIter < Basis_Fn[numberSetsIter].size() ; ++functionSetIter )
             {
               int pos = ((numberSetsIter* numberSet1_ + functionSetIter)* (numberGaussianPoints_plus1) 
                         + numberGaussianPoints_plus1-1)*sizeOfValues;
               basisFunctionsAtGaussianPoints_ = (basisFunctionsAtGaussianPoints+pos);
               for (int evlautaionIter = 0 ; evlautaionIter < sizeOfValues ; ++evlautaionIter )
                  {
                    basisFunctionsAtGaussianPoints_[evlautaionIter] =
                      (*Basis_Fn[numberSetsIter][functionSetIter][evlautaionIter]) (gaussianPoints_,0,0,0,TYPE(0));
                  }
             } 
        }
   }
//-----------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
template <class Expression, typename baseTYPE>
inline void 
BASISFUNCTION_TYPE :: Set_BFN (BasisFunctionSet setNumber,
                               const FunctionExpr<Expression>& expr,
                               baseTYPE **_gaussianPoints)
    {
      assert (0 <= setNumber && setNumber < numberOfBasisFunctionSets);
      assert (numberSet[setNumber] != 0);
      const int cnt = Basis_Fn[setNumber].size();
      assert ( cnt < numberSet[setNumber] );
      Base** BasisFN = new Base*[ sizeOfValues ];
   // Setting the basisfunction expr and its derivatives
      BasisFN [fn] = new Expression(expr) ;
#if 1
      BasisFN [dx] = new FunctionDDX<Expression>(expr) ;
#ifdef SECOND_DERIVATIVES
      BasisFN [dxx] = new FunctionD2DX2<Expression>(expr) ;
#endif
      if (1 < dim_)
         {
           BasisFN [dy] = new FunctionDDY<Expression>(expr) ;
#ifdef SECOND_DERIVATIVES
           BasisFN [dyy] = new FunctionD2DY2<Expression>(expr) ;
           BasisFN [dxy] = new FunctionD2DXY<Expression>(expr) ;
#endif
           if (2 < dim_)
              {
                BasisFN [dz] = new FunctionDDZ<Expression>(expr) ;
#ifdef SECOND_DERIVATIVES
                BasisFN [dzz] = new FunctionD2DZ2<Expression>(expr) ;
                BasisFN [dxz] = new FunctionD2DXZ<Expression>(expr) ;
                BasisFN [dyz] = new FunctionD2DYZ<Expression>(expr) ;
#endif
              }
         }
#endif
      Basis_Fn[setNumber].push_back(BasisFN);
      baseTYPE* levelGaussianPoints = NULL;
      TYPE* basisFunctionsAtGaussianPoints_ = NULL;
      for ( int gaussIter = 0 ; gaussIter < numberGaussianPoints_plus1-1 ; ++gaussIter )
         {
           int pos = (( setNumber* numberSet1_ + cnt )* numberGaussianPoints_plus1 + gaussIter )*sizeOfValues;
           basisFunctionsAtGaussianPoints_ = (basisFunctionsAtGaussianPoints+pos);
           levelGaussianPoints = _gaussianPoints[gaussIter];
           for ( int evaluationIter = 0 ; evaluationIter < sizeOfValues ; ++evaluationIter )
              {
                basisFunctionsAtGaussianPoints_[evaluationIter] = (*BasisFN[evaluationIter])(levelGaussianPoints,0,0,0,TYPE(0));
              }
         }
    }
//-----------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
BASISFUNCTION_TYPE :: ~Basis_Functions()
   {
#if 1
     for ( int numberSetsIter = 0 ; numberSetsIter < numberOfBasisFunctionSets ; ++numberSetsIter )
        {
          for ( unsigned int functionSetIter = 0 ; functionSetIter < Basis_Fn[numberSetsIter].size() ; ++functionSetIter )
             {
               for ( int functionTypeIter = 0 ; functionTypeIter < sizeOfValues ; ++functionTypeIter )
                  {
                    delete (Basis_Fn [numberSetsIter] [functionSetIter] [functionTypeIter]) ;
                  }
               delete [] Basis_Fn [numberSetsIter] [functionSetIter] ;
             }
           Basis_Fn[numberSetsIter].clear() ;
        }
#endif
     delete [] basisFunctionsAtGaussianPoints ;
     int i = 0;
     for (int s=0; s < numberOfBasisFunctionSets; ++s)
        {
          for (; i < dim_*(s+1); ++i)
             {
               for (int j=0; j < numberSet[s]; ++j)
                  {
                    delete [] derivativeValues[i][j];
                  }
               delete [] derivativeValues[i];
             }
        }
     delete [] derivativeValues;
   }
// ---------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
inline TYPE***
BASISFUNCTION_TYPE :: getDerivativeValues() const
    {
      return derivativeValues;
    }
// ---------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
inline unsigned int
 BASISFUNCTION_TYPE:: getNumberOfBasisFunctions(BasisFunctionSet setNumber) const
    {
       assert (0 <= setNumber);
       if ( numberOfBasisFunctionSets <= setNumber)
          {
            return 0;
          }
       else 
          {
         // This was the older and acutally false version!
         // return Basis_Fn[setNumber].size();
            return numberSet[setNumber];
          }
    }
// ---------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
template <class Element>
inline void
BASISFUNCTION_TYPE:: precomputeDerivatives(Element& element, BasisFunctionSet functionSetNumber)
    {
     const int positionInDerivativeArray = dim_*functionSetNumber, basisFunctionLoopSize = Basis_Fn[functionSetNumber].size();
     TYPE*** act_derivativeValues = derivativeValues + positionInDerivativeArray;
     for (int gaussLevel = 0 ; gaussLevel < numberGaussianPoints_plus1-1; ++gaussLevel )
        {
          TYPE* value = element.getValuesOfJacobianMatrix(gaussLevel);
          const TYPE (&ext_) = element.getDeterminant(gaussLevel);
          for (int functionIter = 0; functionIter < basisFunctionLoopSize; ++functionIter)
             {
               const TYPE* basisFunctionsAtGaussianPoints_ =
                     getBasisFunctionsAtGaussianPoints_(functionSetNumber,gaussLevel,functionIter);
               for (int k = 0; k < dim_; ++k)
                  {
                    TYPE deriv = basisFunctionsAtGaussianPoints_[dx] * value[dim_*k];
                    for(int l = dx; l < dim_; ++l)
                       {
                         deriv += basisFunctionsAtGaussianPoints_[dx+l] * value [dim_*k+l];
                       }
                    act_derivativeValues[k][functionIter][gaussLevel] = deriv/ext_;
                  }
             }
        }
    }
// ---------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
template <class Element>
inline void
BASISFUNCTION_TYPE:: precomputeDerivativesAdditional(Element& element, BasisFunctionSet number)
    {
     int gaussLevel = numberGaussianPoints_plus1-1;
     const int deriPos = dim_*number, loopSize = Basis_Fn[number].size();
     TYPE* value = element.getValuesOfJacobianMatrix(gaussLevel);
     const TYPE (&ext_) = element.getDeterminant(gaussLevel);
     const TYPE* basisFunctionsAtGaussianPoints_ = NULL;
     for ( int j = 0 ; j < loopSize ; ++ j )
        {
          basisFunctionsAtGaussianPoints_ = getBasisFunctionsAtGaussianPoints_(number,gaussLevel,j);
          for (int k = 0; k < dim_; ++k)
             {
               TYPE deriv = basisFunctionsAtGaussianPoints_[dx] * value [dim_*k];
               for(int l = dx; l < dim_; ++l)
                  {
                    deriv += basisFunctionsAtGaussianPoints_[l+dx] * value [dim_*k+l];
                  }
               derivativeValues[deriPos+k][j][gaussLevel] = deriv/ext_;
             }
        }
    }
// ---------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
inline const TYPE *
 BASISFUNCTION_TYPE:: getBasisFunctionsAtGaussianPoints_(int basisFunctionSetNumber, int gaussLevel, int actualBasisFunction) const
    {
      const int pos = ((basisFunctionSetNumber*numberSet1_+actualBasisFunction)*numberGaussianPoints_plus1+gaussLevel)*sizeOfValues;
      return (this->basisFunctionsAtGaussianPoints + pos);
    }
// ---------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
template <BasisFunctionSet basisFunctionSetNumber,unsigned int functionId,unsigned int iterLength>
inline const TYPE 
 BASISFUNCTION_TYPE::getBasisFunctionsAtGaussianPointsFn(const int (&iterator)[iterLength]) const 
    {
      assert ( functionId < iterLength - 1 );
      enum {startPos = basisFunctionSetNumber*numberSet1_};
      int realFunctionId = 0; //= (basisFunctionSetNumber==1)?(iterLength-2):((functionId==iterLength-2)?0:functionId);
      if (basisFunctionSetNumber == 1 ) 
        {
          realFunctionId = iterLength-2; 
          assert (functionId == 0);
        }
      else
        {
          realFunctionId = functionId;
        }
      const int pos = ((startPos+iterator[realFunctionId])*numberGaussianPoints_plus1+iterator[iterLength-1])*sizeOfValues;
      return (this->basisFunctionsAtGaussianPoints + pos)[fn];
    }
// ---------------------------------------------------------------------------
BASISFUNCTION_TEMPLATES
template <BasisFunctionSet basisFunctionSetNumber,unsigned int functionId, what derivativeType, unsigned int iterLength>
inline const TYPE 
 BASISFUNCTION_TYPE::getElementMappingDerivativeValues(const int (&iterator)[iterLength]) const 
   { 
#if 0
      std::cout <<"BaisfuncitonSet " << basisFunctionSetNumber << "  functionID " << functionId << " It(" << iterator[0] <<","<< iterator[1] <<","<< iterator[2] <<","<< iterator[3] <<")" << std::endl;
      std::cout << basisFunctionSetNumber <<" " << derivativeType-1 << " " <<iterator[functionId] << " " <<  iterator[iterLength-1] << std::endl;
#endif
      int realFunctionId = 0; //= (basisFunctionSetNumber==1)?(iterLength-2):((functionId==iterLength-2)?0:functionId);
      if (basisFunctionSetNumber == 1 ) 
        {
          realFunctionId = iterLength-2; 
          assert (functionId == 0);
        }
      else
        {
          realFunctionId = functionId;
        }
      
#if 0
      std::cout << "basisFunctionSetNumber " << basisFunctionSetNumber << std::endl;
      std::cout << "functionId " << functionId << std::endl;
      std::cout << "realFunctionId " << realFunctionId << std::endl;
      std::cout << "iterator[realFunctionId] " << iterator[realFunctionId] << std::endl;
#endif
      assert ( functionId < iterLength - 1 );
      assert ( basisFunctionSetNumber < numberOfBasisFunctionSets_ );
      assert ( iterator[realFunctionId] < numberSet[basisFunctionSetNumber] );
      return derivativeValues[basisFunctionSetNumber*dim_+derivativeType-1][iterator[realFunctionId]][iterator[iterLength-1]];
   }
// ---------------------------------------------------------------------------

