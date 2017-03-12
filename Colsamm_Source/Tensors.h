namespace Colsamm_LinA
  {
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints, int dimension1>
    class Matrix;
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints, int dimension>
    class Vector
       {
          private: 
            TYPE vectorData[numberOfGaussianPoints][dimension]; 
          public:
            inline Vector();
            inline Vector<TYPE,numberOfGaussianPoints,dimension>&
               operator= (const Vector<TYPE,numberOfGaussianPoints,dimension>& v);
            inline const TYPE& operator()(const int gaussLevel, int position) const;
            inline TYPE& operator()(const int gaussLevel, int position);
            template <class Element>
            inline void setAsTrafoFactor (ExtensionValues position, const Element& element);
            inline void computeNormalLength ( const Matrix<TYPE,numberOfGaussianPoints,D2>& m);
            inline void computeNormalLength ( const Matrix<TYPE,numberOfGaussianPoints,D3>& m);
       };
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints>
    class Vector<TYPE,numberOfGaussianPoints,D1>
       {
          private: 
            TYPE vectorData[numberOfGaussianPoints]; 
          public:
            inline Vector();
            inline Vector<TYPE,numberOfGaussianPoints,D1>&
               operator= (const Vector<TYPE,numberOfGaussianPoints,D1>& v);
            template <class Expr, class Element>
            inline void operator() (const Expr& expr, const Element& element);
            inline const TYPE& operator()(const int gaussLevel, const int direction) const;
            inline TYPE& operator()(const int gaussLevel, const int direction);
       };
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints>
    class Vector<TYPE,numberOfGaussianPoints,D2>
       {
          private: 
            TYPE vectorData[numberOfGaussianPoints][2]; 
          public:
            inline Vector();
            inline Vector<TYPE,numberOfGaussianPoints,D2>& 
               operator= (const Vector<TYPE,numberOfGaussianPoints,D2>& v);
            template <class Expr, class Element>
            inline void operator() (const Expr& expr, const Element& element);
            inline const TYPE& operator()(const int gaussLevel,const int direction) const;
            inline TYPE& operator()(const int gaussLevel, const int direction);
       };
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints>
    class Vector<TYPE,numberOfGaussianPoints,D3>
       {
          private: 
            TYPE vectorData[numberOfGaussianPoints][3]; 
          public:
            inline Vector();
            inline Vector<TYPE,numberOfGaussianPoints,D3>& 
               operator= (const Vector<TYPE,numberOfGaussianPoints,D3>& v);
            template <class Expr, class Element>
            inline void operator() (const Expr& expr, const Element& element);
            inline const TYPE& operator()(const int gaussLevel,const int direction) const;
            inline TYPE& operator()(const int gaussLevel, const int direction);
            inline void computeNormalLength ( const Matrix<TYPE,numberOfGaussianPoints,D3>& m);
       };
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints>
    class Matrix<TYPE,numberOfGaussianPoints,D1>
       {
          private: 
            TYPE matrixData[numberOfGaussianPoints]; 
          public:
            inline Matrix();
            inline Matrix<TYPE,numberOfGaussianPoints,D1>& 
               operator= (const Matrix<TYPE,numberOfGaussianPoints,D1>& v);
            template <class Expr, class Element>
            inline void operator() (const Expr& expr, const Element& element);
            inline const TYPE& 
            operator()(const int gaussLevel,const int direction1 = 0,const int direction2 = 0) const;
            inline TYPE& 
            operator()(const int gaussLevel,const int direction1 = 0,const int direction2 = 0); 
            
            inline void computeInverseAndDeterminant( const Matrix<TYPE,numberOfGaussianPoints,D1>& m, 
                                                      Vector<TYPE,numberOfGaussianPoints,_extensionSize>& v);
       };
    
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints>
    class Matrix<TYPE,numberOfGaussianPoints,D2>
       {
          private: 
            TYPE matrixData[numberOfGaussianPoints][2][2]; 
          public:
            inline Matrix();
            inline Matrix<TYPE,numberOfGaussianPoints,D2>& 
               operator= (const Matrix<TYPE,numberOfGaussianPoints,D2>& v);
            template <class Expr, class Element>
            inline void operator() (const Expr& expr, const Element& element);
            inline const TYPE& 
            operator()(const int gaussLevel,const int direction1,const int direction2) const;
            inline TYPE& 
            operator()(const int gaussLevel,const int direction1,const int direction2); 
            
            inline void computeInverseAndDeterminant( const Matrix<TYPE,numberOfGaussianPoints,D2>& m, 
                                                      Vector<TYPE,numberOfGaussianPoints,_extensionSize>& v);
       };
    
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints>
    class Matrix<TYPE,numberOfGaussianPoints,D3>
       {
          private: 
            TYPE matrixData[numberOfGaussianPoints][3][3]; 
          public:
            inline Matrix();
            inline Matrix<TYPE,numberOfGaussianPoints,D3>& 
               operator= (const Matrix<TYPE,numberOfGaussianPoints,D3>& v);
            template <class Expr, class Element>
            inline void operator() (const Expr& expr, const Element& element);
            inline const TYPE& 
            operator()(const int gaussLevel,const int direction1,const int direction2) const;
            inline TYPE& 
            operator()(const int gaussLevel,const int direction1,const int direction2); 
            
            inline void computeInverseAndDeterminant( const Matrix<TYPE,numberOfGaussianPoints,D3>& m, 
                                                      Vector<TYPE,numberOfGaussianPoints,_extensionSize>& v);
       };
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
    class Function_Values
       {
          private:
            TYPE functionData[numberOfGaussianPoints][numberOfBasisFunctions]; 
          public:
            Function_Values();
 
            template <class Expr, typename baseTYPE>
            inline void operator()(const int functionIter, const FunctionExpr<Expr>& expr, baseTYPE** gaussianPoints);
            template <typename baseTYPE>
            inline void operator()(const std::vector<Base**>& expr, baseTYPE** gaussianPoints);
            inline const TYPE& 
            operator()(const int gaussLevel,const int fucntionIter) const;
            inline TYPE& 
            operator()(const int gaussLevel,const int fucntionIter); 

            void applyVectorialChainRule(const Matrix<TYPE,numberOfGaussianPoints-1,D1>& m,
                                         const Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>& f,
                                         const TYPE* edgeLength );

            void applyVectorialChainRule(const Matrix<TYPE,numberOfGaussianPoints-1,D2>& m,
                                         const Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>& f,
                                         const TYPE* edgeLength );

            void applyVectorialChainRule(const Matrix<TYPE,numberOfGaussianPoints-1,D3>& m,
                                         const Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>& f,
                                         const TYPE* edgeLength );
       };
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions, int dimension>
    class Function_First_Derivatives;
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
    class Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions, D1>
       {
          private:
            TYPE functionData[numberOfGaussianPoints][numberOfBasisFunctions]; 
          public:
            Function_First_Derivatives();
 
            template <class Expr, typename baseTYPE>
            inline void operator()(const int functionIter, const FunctionExpr<Expr>& expr, baseTYPE** gaussianPoints);
            template <typename baseTYPE>
            inline void operator()(const std::vector<Base**>& expr, baseTYPE** gaussianPoints);
            inline const TYPE& 
            operator()(const int gaussLevel,const int fucntionIter, int direction=0) const;
            inline TYPE& 
            operator()(const int gaussLevel,const int fucntionIter, int direction=0); 

            inline void applyChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D1>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>& v);

            inline void applyVectorialChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D1>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>& v,
                const TYPE* edgeLength );
       };
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
    class Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
       {
          private:
            TYPE functionData[numberOfGaussianPoints][numberOfBasisFunctions][2]; 
          public:
            Function_First_Derivatives();
 
            template <class Expr, typename baseTYPE>
            inline void operator() (const int functionIter, const FunctionExpr<Expr>& expr, baseTYPE** gaussianPoints);
            template <typename baseTYPE>
            inline void operator()(const std::vector<Base**>& expr, baseTYPE** gaussianPoints);
            inline const TYPE& 
            operator()(const int gaussLevel,const int fucntionIter, int direction) const;
            inline TYPE& 
            operator()(const int gaussLevel,const int fucntionIter, int direction); 

            inline void applyChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D2>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>& v);
            inline void applyChainRuleLast( const Matrix<TYPE,numberOfGaussianPoints-1,D2>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>& v);

            inline void applyVectorialChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D2>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>& v,
                const TYPE* edgeLength );
       };
    //-----------------------------------------------------------------------------
    template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
    class Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
       {
          private:
            TYPE functionData[numberOfGaussianPoints][numberOfBasisFunctions][3]; 
          public:
            Function_First_Derivatives();
 
            template <class Expr, typename baseTYPE>
            inline void operator() (const int functionIter, const FunctionExpr<Expr>& expr, baseTYPE** gaussianPoints);
            template <typename baseTYPE>
            inline void operator()(const std::vector<Base**>& expr, baseTYPE** gaussianPoints);
            inline const TYPE& 
            operator()(const int gaussLevel,const int fucntionIter, int direction) const;
            inline TYPE& 
            operator()(const int gaussLevel,const int fucntionIter, int direction); 

            inline void applyChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D3>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>& v);

            inline void applyChainRuleLast( const Matrix<TYPE,numberOfGaussianPoints-1,D3>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>& v);

            inline void applyVectorialChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D3>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>& v,
                const TYPE* edgeLength );
       };
  }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int dimension>
inline std::ostream& operator << (std::ostream& os, const Colsamm_LinA::Vector<TYPE,numberOfGaussianPoints,dimension>& v);
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int dimension>
inline std::ostream& operator << (std::ostream& os, const Colsamm_LinA::Matrix<TYPE,numberOfGaussianPoints,dimension>& m);
//-----------------------------------------------------------------------------

