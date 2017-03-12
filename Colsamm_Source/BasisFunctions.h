//==============================================================================
//
//  $Id: BasisFunctions.h,v 1.6 2006/07/25 13:52:27 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dim_, int numberGaussianPoints_plus1, BasisFunctionType BFnTYPE, 
          int numberSet1_, int numberSet2_ = 0, int numberOfBasisFunctionSets = (0 < numberSet2_)? 2:1>
class Basis_Functions;  
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dim_, int numberGaussianPoints_plus1, 
          int numberSet1_, int numberSet2_, int numberOfBasisFunctionSets>
class Basis_Functions <TYPE,dim_,numberGaussianPoints_plus1,oneDimensional,
                                  numberSet1_, numberSet2_, numberOfBasisFunctionSets>
   {
     protected:
      const int numberOfBasisFunctionSets_, sizeOfValues;
      int numberSet[numberOfBasisFunctionSets];
      typedef std::vector<Base**> BasisFnVector[numberOfBasisFunctionSets];
      BasisFnVector Basis_Fn;
      TYPE* basisFunctionsAtGaussianPoints, ***derivativeValues;
     public:
      inline Basis_Functions (); 
      inline ~Basis_Functions ( ); 
      template <typename baseTYPE>
      inline void eval(baseTYPE** gaussianPoints_);
      template <class Expression, typename baseTYPE>
      inline void 
       Set_BFN (BasisFunctionSet setNumber,  const FunctionExpr<Expression>& a_, baseTYPE** gaussianPoints_); 
    /** 
     Methods for returning the amounts of basis functions
     */ 
      inline unsigned int getNumberOfBasisFunctions(BasisFunctionSet setNumber) const;
      int size_Set1() const { return numberSet[BasisSet_1];} 
      int size_Set2() const { return numberSet[BasisSet_2];} 
      int numberFunctions_Set1() const { return numberSet[BasisSet_1];} 
      int numberFunctions_Set2() const { return numberSet[BasisSet_2];} 
      inline TYPE*** getDerivativeValues() const ;
      inline const TYPE* getBasisFunctionsAtGaussianPoints_(int basisFunctionSetNumber,int gaussLevel,int actualBasisFunction)const;

      template <BasisFunctionSet basisFunctionSetNumber,unsigned int functionId,unsigned int iterLength>
      inline const TYPE getBasisFunctionsAtGaussianPointsFn(const int (&iterator)[iterLength]) const ;

      template <BasisFunctionSet basisFunctionSetNumber, unsigned int functionId, what derivativeType, unsigned int iterLength>
      inline const TYPE getElementMappingDerivativeValues(const int (&iterator)[iterLength]) const ;

      template <class Element>
      inline void precomputeDerivatives(Element& element, BasisFunctionSet number);
      template <class Element>
      inline void precomputeDerivativesAdditional(Element& element, BasisFunctionSet number);
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dim_, int numberGaussianPoints_plus1,
          int numberSet1_, int numberSet2_, int numberOfBasisFunctionSets>
class Basis_Functions <TYPE,dim_,numberGaussianPoints_plus1,vectorBasisFunctions,
                                  numberSet1_, numberSet2_, numberOfBasisFunctionSets>
   {
     protected:
      const int numberOfBasisFunctionSets_, sizeOfValues;
      int numberSet[numberOfBasisFunctionSets];
      typedef std::vector<Base**> BasisFnVector[numberOfBasisFunctionSets];
      BasisFnVector Basis_Fn;
      std::vector<std::pair<unsigned int, unsigned int > > edges;
      TYPE* basisFunctionsAtGaussianPoints, ***derivativeValues, **value;
      TYPE* vector_basisFunctionsAtGaussianPoints, ***vector_derivativeValues;
     public:
      inline Basis_Functions (); 
      inline ~Basis_Functions ( ); 
      template <typename baseTYPE>
      inline void eval(baseTYPE** gaussianPoints_);
      template <class Expr, typename baseTYPE>
      inline void 
      Set_BFN (BasisFunctionSet setNumber, const FunctionExprVector1D<Expr>& expr,
               const unsigned cnr_1, const unsigned cnr_2, 
               baseTYPE **_gaussianPoints);
      template <class Expr1, class Expr2, typename baseTYPE>
      inline void 
      Set_BFN (BasisFunctionSet setNumber, const FunctionExprVector2D<Expr1,Expr2>& expr,
               const unsigned cnr_1, const unsigned cnr_2, 
               baseTYPE **_gaussianPoints);
      template <class Expr1, class Expr2, class Expr3, typename baseTYPE>
      inline void 
      Set_BFN (BasisFunctionSet setNumber, const FunctionExprVector3D<Expr1,Expr2,Expr3>& expr, 
               const unsigned cnr_1, const unsigned cnr_2, 
               baseTYPE **_gaussianPoints);
    /** 
     Methods for returning the amounts of basis functions
     */ 
      inline unsigned int getNumberOfBasisFunctions(BasisFunctionSet setNumber) const;
      int size_Set1() const { return numberSet[BasisSet_1];} 
      int size_Set2() const { return numberSet[BasisSet_2];} 
      int numberFunctions_Set1() const { return numberSet[BasisSet_1]*dim_;} 
      int numberFunctions_Set2() const { return numberSet[BasisSet_2]*dim_;} 
      inline TYPE*** getDerivativeValues() const ;
      inline const TYPE& getElementMappingDerivativeValues(const int derivPos, const int iterPos, const int gaussLevel) const ;

      inline const TYPE* getBasisFunctionsAtGaussianPoints_(int basisFunctionSetNumber,int gaussLevel,int actualBasisFunction)const;

      template <BasisFunctionSet basisFunctionSetNumber,unsigned int functionId,unsigned int iterLength>
      inline const TYPE getBasisFunctionsAtGaussianPointsFn(const int (&iterator)[iterLength]) const ;

      template <BasisFunctionSet basisFunctionSetNumber, unsigned int functionId, what derivativeType, unsigned int iterLength>
      inline const TYPE getElementMappingDerivativeValues(const int (&iterator)[iterLength]) const ;

      template <BasisFunctionSet basisFunctionSetNumber,SpaceDirection direction,unsigned int functionId,unsigned int iterLength>
      inline const TYPE getBasisFunctionsAtGaussianPointsFn_Vector(const int (&iterator)[iterLength]) const ;

      template <BasisFunctionSet basisFunctionSetNumber,SpaceDirection direction,
                unsigned int functionId, what derivativeType, unsigned int iterLength>
      inline const TYPE getElementMappingDerivativeValues_Vector(const int (&iterator)[iterLength]) const ;
      template <class Element>
      inline void precomputeDerivatives(Element& element, BasisFunctionSet number);
      template <class Element>
      inline void precomputeDerivativesAdditional(Element& element, BasisFunctionSet number);
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dim_, int numberGaussianPoints_plus1,
          int numberSet1_, int numberSet2_, int numberOfBasisFunctionSets>
class Basis_Functions <TYPE,dim_,numberGaussianPoints_plus1,derivativeBasisFunctions,
                                  numberSet1_, numberSet2_, numberOfBasisFunctionSets>
   {
     protected:
      const int numberOfBasisFunctionSets_, sizeOfValues;
      int numberSet[numberOfBasisFunctionSets];
      typedef std::vector<Base**> BasisFnVector[numberOfBasisFunctionSets];
      BasisFnVector Basis_Fn;
      std::vector<std::pair<unsigned int, unsigned int > > edges;
      TYPE* basisFunctionsAtGaussianPoints, ***derivativeValues, **value;
      TYPE* vector_basisFunctionsAtGaussianPoints, ***vector_derivativeValues;
      TypeDOF typesDOF[numberSet1_];
     public:
      inline Basis_Functions (); 
      inline ~Basis_Functions ( ); 
      template <typename baseTYPE>
      inline void eval(baseTYPE** gaussianPoints_);
      template <class Expr, typename baseTYPE>
      inline void 
      Set_BFN (BasisFunctionSet setNumber, const FunctionExprVector1D<Expr>& expr,
               const unsigned cnr_1, const unsigned cnr_2, 
               baseTYPE **_gaussianPoints);
      template <class Expr1, class Expr2, typename baseTYPE>
      inline void 
      Set_BFN (BasisFunctionSet setNumber, const FunctionExprVector2D<Expr1,Expr2>& expr,
               const unsigned cnr_1, const unsigned cnr_2, 
               baseTYPE **_gaussianPoints);
      template <class Expr1, class Expr2, class Expr3, typename baseTYPE>
      inline void 
      Set_BFN (BasisFunctionSet setNumber, const FunctionExprVector3D<Expr1,Expr2,Expr3>& expr, 
               const unsigned cnr_1, const unsigned cnr_2, 
               baseTYPE **_gaussianPoints);
    /** 
     Methods for returning the amounts of basis functions
     */ 
      inline unsigned int getNumberOfBasisFunctions(BasisFunctionSet setNumber) const;
      int size_Set1() const { return numberSet[BasisSet_1];} 
      int size_Set2() const { return numberSet[BasisSet_2];} 
      int numberFunctions_Set1() const { return numberSet[BasisSet_1]*dim_;} 
      int numberFunctions_Set2() const { return numberSet[BasisSet_2]*dim_;} 
      inline TYPE*** getDerivativeValues() const ;
      inline const TYPE& getElementMappingDerivativeValues(const int derivPos, const int iterPos, const int gaussLevel) const ;

      inline const TYPE* getBasisFunctionsAtGaussianPoints_(int basisFunctionSetNumber,int gaussLevel,int actualBasisFunction)const;

      template <BasisFunctionSet basisFunctionSetNumber,unsigned int functionId,unsigned int iterLength>
      inline const TYPE getBasisFunctionsAtGaussianPointsFn(const int (&iterator)[iterLength]) const ;

      template <BasisFunctionSet basisFunctionSetNumber, unsigned int functionId, what derivativeType, unsigned int iterLength>
      inline const TYPE getElementMappingDerivativeValues(const int (&iterator)[iterLength]) const ;

      template <BasisFunctionSet basisFunctionSetNumber,SpaceDirection direction,unsigned int functionId,unsigned int iterLength>
      inline const TYPE getBasisFunctionsAtGaussianPointsFn_Vector(const int (&iterator)[iterLength]) const ;

      template <BasisFunctionSet basisFunctionSetNumber,SpaceDirection direction,
                unsigned int functionId, what derivativeType, unsigned int iterLength>
      inline const TYPE getElementMappingDerivativeValues_Vector(const int (&iterator)[iterLength]) const ;
      template <class Element>
      inline void precomputeDerivatives(Element& element, BasisFunctionSet number);
      template <class Element>
      inline void precomputeDerivativesAdditional(Element& element, BasisFunctionSet number);
   } ;
//==============================================================================

