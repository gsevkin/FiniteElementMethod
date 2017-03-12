
// ---------------------------------------------------------------------------
#define DERIVATIVE_BASISFUNCTION_TEMPLATES \
  template <typename TYPE, \
            DIMENSION dim_,\
            int numberGaussianPoints_plus1,\
            int numberSet1_,\
            int numberSet2_,\
            int numberOfBasisFunctionSets>

#define DERIVATIVE_BASISFUNCTION_TYPE\
   Basis_Functions <\
        TYPE,\
        dim_,\
        numberGaussianPoints_plus1,\
        derivativeBasisFunctions,\
        numberSet1_,\
        numberSet2_,\
        numberOfBasisFunctionSets>
//-----------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
inline
DERIVATIVE_BASISFUNCTION_TYPE :: Basis_Functions ()
   : numberOfBasisFunctionSets_(numberOfBasisFunctionSets), sizeOfValues(1+dim_)
   {
   // Zuweisung der Groesse der ersten Menge der Basisfunktionen
     numberSet[0] = numberSet1_;
   // Zuweisung der Groesse der zweiten Menge der Basisfunktionen,
   // falls sie nicht null ist.
     if (1 < numberOfBasisFunctionSets)
       {
         numberSet[1] = numberSet2_;
       }
   // Berechnung der gesamten Anzahl an Basisfunktionen ( #Menge1 + #Menge2 )
     int totalSize = 0;
     for (int s=0; s < numberOfBasisFunctionSets; ++s)
        {
          totalSize += numberSet[s];
        }
   // die Groesse der vorberechneten Arrays f"ur Basisfunktionen und Ableitungen der Basisfunktionen
     int sizeOf_basisFunctionsAtGaussianPoints = numberGaussianPoints_plus1*sizeOfValues*(totalSize) * dim_,
         sizeOf_derivativeValues = dim_*numberOfBasisFunctionSets;
   // Allokation der Felder ...
     basisFunctionsAtGaussianPoints = new TYPE[sizeOf_basisFunctionsAtGaussianPoints];
     vector_basisFunctionsAtGaussianPoints = new TYPE[sizeOf_basisFunctionsAtGaussianPoints];
     derivativeValues = new TYPE**[sizeOf_derivativeValues];
     vector_derivativeValues = new TYPE**[sizeOf_derivativeValues];

     int i = 0;
     for (int s=0; s < numberOfBasisFunctionSets; ++s)
        {
          for (; i < dim_*(s+1); ++i) // <-- This size has to change for vector functions! 
             {
               derivativeValues[i] = new TYPE* [numberSet[s] * dim_];
               vector_derivativeValues[i] = new TYPE* [numberSet[s] * dim_];
                for (int j=0; j < numberSet[s] * dim_; ++j)
                   {
                    derivativeValues[i][j] = new TYPE [numberGaussianPoints_plus1];
                    vector_derivativeValues[i][j] = new TYPE [numberGaussianPoints_plus1];
#ifdef INITIALIZE_ARRAYS
                    for (int k=0; k < numberGaussianPoints_plus1; ++k)
                       {
                         derivativeValues[i][j][k] = 0.;
                         vector_derivativeValues[i][j][k] = 0.;
                       }
#endif
                   }
             }
        }
     value = new TYPE*[numberGaussianPoints_plus1];
   }
//-----------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <typename baseTYPE>
inline void 
DERIVATIVE_BASISFUNCTION_TYPE :: eval(baseTYPE**_gaussianPoints) 
   {
   // Die Auswertung der Basisfunktionen fuer den zusaetzlichen Punkt
     baseTYPE* gaussianPoints_ = _gaussianPoints[numberGaussianPoints_plus1-1];
     TYPE* basisFunctionsAtGaussianPoints_ = NULL;
     for (int numberSetsIter = 0 ; numberSetsIter < numberOfBasisFunctionSets ; ++numberSetsIter )
        {
          for (unsigned int functionSetIter = 0 ; functionSetIter < Basis_Fn[numberSetsIter].size() ; ++functionSetIter )
             {
               int pos = ((numberSetsIter* numberSet1_*dim_ + functionSetIter)* (numberGaussianPoints_plus1) 
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
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <class Expr, typename baseTYPE>
inline void 
DERIVATIVE_BASISFUNCTION_TYPE :: Set_BFN (BasisFunctionSet setNumber,
                               const FunctionExprVector1D<Expr>& expr,
                               const unsigned int cnr_1, const unsigned int cnr_2, 
                               baseTYPE **_gaussianPoints)
    {
      assert (dim_ == D1);
      int cnt = 0;
      assert ( cnt < numberSet[setNumber] * dim_);
      baseTYPE* levelGaussianPoints = NULL;
      TYPE* basisFunctionsAtGaussianPoints_ = NULL;
      Base** BasisFN = NULL;
   // Setting the basisfunction expr and its derivatives
      cnt = Basis_Fn[setNumber].size();
      assert ( cnt < numberSet[setNumber] * dim_);
      BasisFN = new Base*[ sizeOfValues ];
      BasisFN [fn] = new Expr(expr.getExpr()) ;
      BasisFN [dx] = new FunctionDDX<Expr>(expr.getExpr()) ;
      Basis_Fn[setNumber].push_back(BasisFN);
      for ( int gaussIter = 0 ; gaussIter < numberGaussianPoints_plus1-1 ; ++gaussIter )
         {
          int pos = (( setNumber* numberSet1_*dim_ + cnt )* numberGaussianPoints_plus1 + gaussIter )*sizeOfValues;
          basisFunctionsAtGaussianPoints_ = (basisFunctionsAtGaussianPoints+pos);
          levelGaussianPoints = _gaussianPoints[gaussIter];
          for ( int evaluationIter = 0 ; evaluationIter < sizeOfValues ; ++evaluationIter )
             {
               basisFunctionsAtGaussianPoints_[evaluationIter] = (*BasisFN[evaluationIter])(levelGaussianPoints,0,0,0,TYPE(0));
             }
         }
      cnt = Basis_Fn[setNumber].size();
      edges.push_back(std::pair<unsigned int, unsigned int>(cnr_1, cnr_2) );
    }
//-----------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <class Expr1, class Expr2, typename baseTYPE>
inline void 
DERIVATIVE_BASISFUNCTION_TYPE :: Set_BFN (BasisFunctionSet setNumber,
                               const FunctionExprVector2D<Expr1,Expr2>& expr,
                               const unsigned int cnr_1, const unsigned int cnr_2, 
                               baseTYPE **_gaussianPoints)
    {
      assert (dim_ == D2);
      int cnt = 0;
      assert ( cnt < numberSet[setNumber] * dim_);
      baseTYPE* levelGaussianPoints = NULL;
      TYPE* basisFunctionsAtGaussianPoints_ = NULL;
      Base** BasisFN = NULL;
   // Setting the basisfunction expr and its derivatives
      cnt = Basis_Fn[setNumber].size();
      assert ( cnt < numberSet[setNumber] * dim_);
      BasisFN = new Base*[ sizeOfValues ];
      BasisFN [fn] = new Expr1(expr.getExpr1()) ;
      BasisFN [dx] = new FunctionDDX<Expr1>(expr.getExpr1()) ;
      BasisFN [dy] = new FunctionDDY<Expr1>(expr.getExpr1()) ;
      Basis_Fn[setNumber].push_back(BasisFN);
      for ( int gaussIter = 0 ; gaussIter < numberGaussianPoints_plus1-1 ; ++gaussIter )
         {
          int pos = (( setNumber* numberSet1_*dim_ + cnt )* numberGaussianPoints_plus1 + gaussIter )*sizeOfValues;
          basisFunctionsAtGaussianPoints_ = (basisFunctionsAtGaussianPoints+pos);
          levelGaussianPoints = _gaussianPoints[gaussIter];
          for ( int evaluationIter = 0 ; evaluationIter < sizeOfValues ; ++evaluationIter )
             {
               basisFunctionsAtGaussianPoints_[evaluationIter] = (*BasisFN[evaluationIter])(levelGaussianPoints,0,0,0,TYPE(0));
             }
         }
      cnt = Basis_Fn[setNumber].size();
      assert ( cnt < numberSet[setNumber] * dim_);
      BasisFN = new Base*[ sizeOfValues ];
      BasisFN [fn] = new Expr2(expr.getExpr2()) ;
      BasisFN [dx] = new FunctionDDX<Expr2>(expr.getExpr2()) ;
      BasisFN [dy] = new FunctionDDY<Expr2>(expr.getExpr2()) ;
      Basis_Fn[setNumber].push_back(BasisFN);
      for ( int gaussIter = 0 ; gaussIter < numberGaussianPoints_plus1-1 ; ++gaussIter )
         {
          int pos = (( setNumber* numberSet1_*dim_ + cnt )* numberGaussianPoints_plus1 + gaussIter )*sizeOfValues;
          basisFunctionsAtGaussianPoints_ = (basisFunctionsAtGaussianPoints+pos);
          levelGaussianPoints = _gaussianPoints[gaussIter];
          for ( int evaluationIter = 0 ; evaluationIter < sizeOfValues ; ++evaluationIter )
             {
               basisFunctionsAtGaussianPoints_[evaluationIter] = (*BasisFN[evaluationIter])(levelGaussianPoints,0,0,0,TYPE(0));
             }
         }
      edges.push_back(std::pair<unsigned int, unsigned int>(cnr_1, cnr_2) );
    }
//-----------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <class Expr1, class Expr2, class Expr3, typename baseTYPE>
inline void 
DERIVATIVE_BASISFUNCTION_TYPE :: Set_BFN (BasisFunctionSet setNumber,
                               const FunctionExprVector3D<Expr1,Expr2,Expr3>& expr,
                               const unsigned int cnr_1, const unsigned int cnr_2, 
                               baseTYPE **_gaussianPoints)
    {
     
      assert (dim_ == D3);
      int cnt = 0;
      baseTYPE* levelGaussianPoints = NULL;
      TYPE* basisFunctionsAtGaussianPoints_ = NULL;
      Base** BasisFN = NULL;
   // Setting the basisfunction expr and its derivatives
   //-------------------------------------------------------------
      cnt = Basis_Fn[setNumber].size();
      assert ( cnt < numberSet[setNumber] * dim_);
      BasisFN = new Base*[ sizeOfValues ];
      BasisFN [fn] = new Expr1(expr.getExpr1()) ;
      BasisFN [dx] = new FunctionDDX<Expr1>(expr.getExpr1()) ;
      if (1 < dim_)
         {
           BasisFN [dy] = new FunctionDDY<Expr1>(expr.getExpr1()) ;
           if (2 < dim_)
              {
                BasisFN [dz] = new FunctionDDZ<Expr1>(expr.getExpr1()) ;
              }
         }
      Basis_Fn[setNumber].push_back(BasisFN);
      for ( int gaussIter = 0 ; gaussIter < numberGaussianPoints_plus1-1 ; ++gaussIter )
         {
          int pos = (( setNumber* numberSet1_*dim_ + cnt )* numberGaussianPoints_plus1 + gaussIter )*sizeOfValues;
          basisFunctionsAtGaussianPoints_ = (basisFunctionsAtGaussianPoints+pos);
          levelGaussianPoints = _gaussianPoints[gaussIter];
          for ( int evaluationIter = 0 ; evaluationIter < sizeOfValues ; ++evaluationIter )
             {
             //assert ( fabs(basisFunctionsAtGaussianPoints_[evaluationIter]) < 1.e-10 );
               basisFunctionsAtGaussianPoints_[evaluationIter] = (*BasisFN[evaluationIter])(levelGaussianPoints,0,0,0,TYPE(0));
             }
         }
   //-------------------------------------------------------------
      cnt = Basis_Fn[setNumber].size();
      assert ( cnt < numberSet[setNumber] * dim_);
      BasisFN = new Base*[ sizeOfValues ];
      BasisFN [fn] = new Expr2(expr.getExpr2()) ;
      BasisFN [dx] = new FunctionDDX<Expr2>(expr.getExpr2()) ;
      if (1 < dim_)
         {
           BasisFN [dy] = new FunctionDDY<Expr2>(expr.getExpr2()) ;
           if (2 < dim_)
              {
                BasisFN [dz] = new FunctionDDZ<Expr2>(expr.getExpr2()) ;
              }
         }
      Basis_Fn[setNumber].push_back(BasisFN);
      for ( int gaussIter = 0 ; gaussIter < numberGaussianPoints_plus1-1 ; ++gaussIter )
         {
          int pos = (( setNumber* numberSet1_*dim_ + cnt )* numberGaussianPoints_plus1 + gaussIter )*sizeOfValues;
          basisFunctionsAtGaussianPoints_ = (basisFunctionsAtGaussianPoints+pos);
          levelGaussianPoints = _gaussianPoints[gaussIter];
          for ( int evaluationIter = 0 ; evaluationIter < sizeOfValues ; ++evaluationIter )
             {
            // assert ( fabs(basisFunctionsAtGaussianPoints_[evaluationIter]) < 1.e-10 );
               basisFunctionsAtGaussianPoints_[evaluationIter] = (*BasisFN[evaluationIter])(levelGaussianPoints,0,0,0,TYPE(0));
             }
         }
   //-------------------------------------------------------------
      cnt = Basis_Fn[setNumber].size();
      assert ( cnt < numberSet[setNumber] * dim_);
      BasisFN = new Base*[ sizeOfValues ];
      BasisFN [fn] = new Expr3(expr.getExpr3()) ;
      BasisFN [dx] = new FunctionDDX<Expr3>(expr.getExpr3()) ;
      if (1 < dim_)
         {
           BasisFN [dy] = new FunctionDDY<Expr3>(expr.getExpr3()) ;
           if (2 < dim_)
              {
                BasisFN [dz] = new FunctionDDZ<Expr3>(expr.getExpr3()) ;
              }
         }
      Basis_Fn[setNumber].push_back(BasisFN);
      for ( int gaussIter = 0 ; gaussIter < numberGaussianPoints_plus1-1 ; ++gaussIter )
         {
          int pos = (( setNumber* numberSet1_*dim_ + cnt )* numberGaussianPoints_plus1 + gaussIter )*sizeOfValues;
          basisFunctionsAtGaussianPoints_ = (basisFunctionsAtGaussianPoints+pos);
          levelGaussianPoints = _gaussianPoints[gaussIter];
          for ( int evaluationIter = 0 ; evaluationIter < sizeOfValues ; ++evaluationIter )
             {
             //assert ( fabs(basisFunctionsAtGaussianPoints_[evaluationIter]) < 1.e-10 );
               basisFunctionsAtGaussianPoints_[evaluationIter] = (*BasisFN[evaluationIter])(levelGaussianPoints,0,0,0,TYPE(0));
             }
         }
      edges.push_back(std::pair<unsigned int, unsigned int>(cnr_1, cnr_2) );
    }
//-----------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
DERIVATIVE_BASISFUNCTION_TYPE :: ~Basis_Functions()
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
     delete [] vector_basisFunctionsAtGaussianPoints ;
     int i = 0;
     for (int s=0; s < numberOfBasisFunctionSets; ++s)
        {
          for (; i < dim_*(s+1); ++i)
             {
               for (int j=0; j < numberSet[s]*dim_; ++j)
                  {
                    delete [] derivativeValues[i][j];
                    delete [] vector_derivativeValues[i][j];
                  }
               delete [] derivativeValues[i];
               delete [] vector_derivativeValues[i];
             }
        }
     delete [] derivativeValues;
     delete [] vector_derivativeValues;
     delete [] value;
   }
// ---------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
inline TYPE***
DERIVATIVE_BASISFUNCTION_TYPE :: getDerivativeValues() const
    {
      return derivativeValues;
    }
// ---------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
inline unsigned int
 DERIVATIVE_BASISFUNCTION_TYPE:: getNumberOfBasisFunctions(BasisFunctionSet setNumber) const
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
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <class Element>
inline void
DERIVATIVE_BASISFUNCTION_TYPE:: precomputeDerivatives(Element& element, BasisFunctionSet functionSetNumber)
    {
     const int positionInDerivativeArray = dim_*functionSetNumber, basisFunctionLoopSize = Basis_Fn[functionSetNumber].size();
     TYPE*** act_derivativeValues = derivativeValues+positionInDerivativeArray;
     for (int gaussLevel = 0 ; gaussLevel < numberGaussianPoints_plus1-1; ++gaussLevel )
        {
          value[gaussLevel] = element.getValuesOfJacobianMatrix(gaussLevel);
          const TYPE ext_ = element.getDeterminant(gaussLevel);
          for (int functionIter = 0; functionIter < basisFunctionLoopSize; ++functionIter)
             {
               if ( this->typesDOF [functionIter] == functionDOF )
                 {
                   const TYPE* basisFunctionsAtGaussianPoints_ =
                       getBasisFunctionsAtGaussianPoints_(functionSetNumber,gaussLevel,functionIter);
                    for (int k = 0; k < dim_; ++k)
                       {
                         TYPE deriv = basisFunctionsAtGaussianPoints_[dx] * value[gaussLevel][dim_*k];
                         for(int l = dx; l < dim_; ++l)
                            {
                              deriv += basisFunctionsAtGaussianPoints_[dx+l] * value[gaussLevel][dim_*k+l];
                            }
                         act_derivativeValues[k][functionIter][gaussLevel] = deriv/ext_;
                       }
                  }
                 else
                  {
                   const TYPE* basisFunctionsAtGaussianPoints_ =
                       getBasisFunctionsAtGaussianPoints_(functionSetNumber,gaussLevel,functionIter);
                    for (int k = 0; k < dim_; ++k)
                       {
                         act_derivativeValues[k][functionIter][gaussLevel] = basisFunctionsAtGaussianPoints_[k]; 
                       }
                  }
             }
        }
     // Here we want to precompute the additional calculations concerning the vector basis functions
     for (int gaussLevel = 0 ; gaussLevel < numberGaussianPoints_plus1-1; ++gaussLevel )
        {
          const TYPE ext_ = element.getDeterminant(gaussLevel);
          for (int functionIter = 0; functionIter < basisFunctionLoopSize; ++functionIter)
             {
               if ( this->typesDOF [functionIter] == functionDOF )
                 {
                    int pos_vec = ((functionSetNumber*numberSet1_*dim_+functionIter)*
                                                  numberGaussianPoints_plus1 + gaussLevel)*sizeOfValues;
                   // (vector_basisFunctionsAtGaussianPoints+pos_vec)[fn] = (this->basisFunctionsAtGaussianPoints + pos)[fn];
                 }
                else
                 {
                    int pos_vec = ((functionSetNumber*numberSet1_*dim_+functionIter)*
                                                  numberGaussianPoints_plus1 + gaussLevel)*sizeOfValues;
                    int pos = //((functionSetNumber*numberSet1_*dim_+vec_nr*dim_)*
                                           (  numberGaussianPoints_plus1 + gaussLevel)*sizeOfValues;
                    TYPE result = 0; // (this->basisFunctionsAtGaussianPoints + pos)[fn] * value[gaussLevel][vec_pos*dim_];
                    for (int i=1; i < dim_; ++i)
                       {
                          pos = (//(functionSetNumber*numberSet1_*dim_+vec_nr*dim_+i)*
                                                        numberGaussianPoints_plus1 + gaussLevel)*sizeOfValues;
                        //  result += (this->basisFunctionsAtGaussianPoints + pos)[fn] * value[gaussLevel][i+vec_pos*dim_];
                       }
                 //   (vector_basisFunctionsAtGaussianPoints+pos_vec)[fn] =
                   //      result /ext_ ; //* element.edge_length(edges[vec_nr].first, edges[vec_nr].second);
                }
             }
        }
    }
// ---------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <class Element>
inline void
DERIVATIVE_BASISFUNCTION_TYPE:: precomputeDerivativesAdditional(Element& element, BasisFunctionSet functionSetNumber)
    {
     int gaussLevel = numberGaussianPoints_plus1-1;
     const int positionInDerivativeArray = dim_*functionSetNumber, basisFunctionLoopSize = Basis_Fn[functionSetNumber].size();
     value[gaussLevel] = element.getValuesOfJacobianMatrix(gaussLevel);
     const TYPE ext_ = element.getDeterminant(gaussLevel);
     const TYPE* basisFunctionsAtGaussianPoints_ = NULL;
     for ( int j = 0 ; j < basisFunctionLoopSize ; ++ j )
        {
          basisFunctionsAtGaussianPoints_ = getBasisFunctionsAtGaussianPoints_(functionSetNumber,gaussLevel,j);
          for (int k = 0; k < dim_; ++k)
             {
               TYPE deriv = basisFunctionsAtGaussianPoints_[dx] * value[gaussLevel][dim_*k];
               for(int l = dx; l < dim_; ++l)
                  {
                    deriv += basisFunctionsAtGaussianPoints_[l+dx] * value[gaussLevel][dim_*k+l];
                  }
               derivativeValues[positionInDerivativeArray+k][j][gaussLevel] = deriv/ext_;
             }
        }
     // Here we want to precompute the additional calculations concerning the vector basis functions
     for (int functionIter = 0; functionIter < basisFunctionLoopSize; ++functionIter)
        {
          const int vec_nr = functionIter/dim_, vec_pos = functionIter - vec_nr*dim_;
          int pos = ((functionSetNumber*numberSet1_*dim_+vec_nr*dim_)*
                                        numberGaussianPoints_plus1 + gaussLevel)*sizeOfValues;
          int pos_vec = ((functionSetNumber*numberSet1_*dim_+functionIter)*
                                        numberGaussianPoints_plus1 + gaussLevel)*sizeOfValues;
          TYPE result =  (this->basisFunctionsAtGaussianPoints + pos)[fn] * value[gaussLevel][vec_pos*dim_];
          for (int i=1; i < dim_; ++i)
             {
              pos = ((functionSetNumber*numberSet1_*dim_+vec_nr*dim_+i)*
                                                 numberGaussianPoints_plus1 + gaussLevel)*sizeOfValues;
               result += (this->basisFunctionsAtGaussianPoints + pos)[fn] * value[gaussLevel][i+vec_pos*dim_];
             }
          (vector_basisFunctionsAtGaussianPoints+pos_vec)[fn] =
                result /ext_ * element.edge_length(edges[vec_nr].first, edges[vec_nr].second);
        }
     for (int functionIter = 0; functionIter < basisFunctionLoopSize; ++functionIter)
        {
          const int vec_nr = functionIter/dim_, vec_pos = functionIter - vec_nr*dim_;
          for (int dd = 0; dd < dim_; ++dd)
             {
               TYPE result = derivativeValues[functionSetNumber*dim_+dd][vec_nr*dim_][gaussLevel] * value[gaussLevel][vec_pos*dim_];
               for (int i=1; i < dim_; ++i)
                  {
                    result += derivativeValues[functionSetNumber*dim_+dd][vec_nr*dim_+i][gaussLevel] * value[gaussLevel][i+vec_pos*dim_];
                  }
               vector_derivativeValues[functionSetNumber*dim_+dd][functionIter][gaussLevel] = 
                      result / ext_ * element.edge_length(edges[vec_nr].first,edges[vec_nr].second);
             }
        }

    }
// ---------------------------------------------------------------------------
// These two functions have to change, since we have to multiply with the 
// Jacobian matrix of the actual mapping!
DERIVATIVE_BASISFUNCTION_TEMPLATES
inline const TYPE *
 DERIVATIVE_BASISFUNCTION_TYPE:: getBasisFunctionsAtGaussianPoints_(int basisFunctionSetNumber, int gaussLevel, int actualBasisFunction) const
    {
      const int pos = ((basisFunctionSetNumber*numberSet1_*dim_+actualBasisFunction)*numberGaussianPoints_plus1 
                        + gaussLevel)*sizeOfValues;
      return (this->basisFunctionsAtGaussianPoints + pos);
    }
// ---------------------------------------------------------------------------
// These two functions have to change, since we have to multiply with the 
// Jacobian matrix of the actual mapping!
// ---------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <BasisFunctionSet basisFunctionSetNumber,unsigned int functionId,unsigned int iterLength>
inline const TYPE 
 DERIVATIVE_BASISFUNCTION_TYPE::getBasisFunctionsAtGaussianPointsFn(const int (&iterator)[iterLength]) const 
    {
   // the computation of the actual function number might be wrong for mixed basis functions!
      assert ( functionId < iterLength - 1 );
      int pos = ((basisFunctionSetNumber*numberSet1_*dim_+ iterator[functionId])*
                                             numberGaussianPoints_plus1 + iterator[iterLength-1])*sizeOfValues;
      return (this->basisFunctionsAtGaussianPoints + pos)[fn]; 
    }
// ---------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <BasisFunctionSet basisFunctionSetNumber, unsigned int functionId, what derivativeType, unsigned int iterLength>
inline const TYPE 
 DERIVATIVE_BASISFUNCTION_TYPE::getElementMappingDerivativeValues(const int (&iterator)[iterLength]) const 
   { 
      assert ( functionId < iterLength - 1 );
      return derivativeValues[basisFunctionSetNumber*dim_+derivativeType-1]
                                    [iterator[functionId]]
                                    [iterator[iterLength-1]];
   }
// ---------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <BasisFunctionSet basisFunctionSetNumber,SpaceDirection direction,unsigned int functionId,unsigned int iterLength>
inline const TYPE 
 DERIVATIVE_BASISFUNCTION_TYPE::getBasisFunctionsAtGaussianPointsFn_Vector(const int (&iterator)[iterLength]) const 
    {
      if ( (unsigned int)(direction) < (unsigned int)(dim_))
         {
   // the computation of the actual function number might be wrong for mixed basis functions!
            assert ( functionId < iterLength - 1 );
            int pos = ((basisFunctionSetNumber*numberSet1_*dim_+ iterator[functionId] *dim_ + direction)*
                                             numberGaussianPoints_plus1 + iterator[iterLength-1])*sizeOfValues;
            return (this->vector_basisFunctionsAtGaussianPoints + pos)[fn]; 
         }
      else 
         {
            return (TYPE)(0.);
         }
    }
// ---------------------------------------------------------------------------
DERIVATIVE_BASISFUNCTION_TEMPLATES
template <BasisFunctionSet basisFunctionSetNumber,SpaceDirection direction,
          unsigned int functionId, what derivativeType, unsigned int iterLength>
inline const TYPE 
 DERIVATIVE_BASISFUNCTION_TYPE::getElementMappingDerivativeValues_Vector(const int (&iterator)[iterLength]) const 
   { 
      assert ( functionId < iterLength - 1 );
      return vector_derivativeValues[basisFunctionSetNumber*dim_+derivativeType-1]
                                    [iterator[functionId]*dim_+direction]
                                    [iterator[iterLength-1]];
   }


