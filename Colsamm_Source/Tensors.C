//-----------------------------------------------------------------------------
#define VEC_TEMPLATE template <typename TYPE, int numberOfGaussianPoints, int dimension>
#define VEC_CLASS_TEMPLATE Colsamm_LinA::Vector<TYPE,numberOfGaussianPoints,dimension>
#define MAT_TEMPLATE template <typename TYPE, int numberOfGaussianPoints, int dimension>
#define MAT_CLASS_TEMPLATE Colsamm_LinA::Matrix<TYPE,numberOfGaussianPoints,dimension>
#define RES_VEC_TEMPLATE template <typename TYPE, int numberOfGaussianPoints>
#define RES_VEC_CLASS_TEMPLATE(dim) Colsamm_LinA::Vector<TYPE,numberOfGaussianPoints,dim>
#define RES_MAT_TEMPLATE template <typename TYPE, int numberOfGaussianPoints>
#define RES_MAT_CLASS_TEMPLATE(dim) Colsamm_LinA::Matrix<TYPE,numberOfGaussianPoints,dim>
//-----------------------------------------------------------------------------
VEC_TEMPLATE
VEC_CLASS_TEMPLATE::Vector()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           for (int dim=0; dim < dimension; ++dim)
              {
                 vectorData[gaussLevel][dim] = (TYPE)0;
             }
        }
#endif
   }
//-----------------------------------------------------------------------------
VEC_TEMPLATE
VEC_CLASS_TEMPLATE&
VEC_CLASS_TEMPLATE::operator= (const VEC_CLASS_TEMPLATE& v)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           for (int dim=0; dim < numberOfGaussianPoints; ++dim)
              {
                 vectorData[gaussLevel][dim] = v.vectorData[gaussLevel][dim];
             }
        }
   }
//-----------------------------------------------------------------------------
VEC_TEMPLATE
const TYPE&
VEC_CLASS_TEMPLATE::operator()(const int gaussLevel,int position) const
   {
      assert( position < dimension ); 
      return vectorData[gaussLevel][position];
   }
//-----------------------------------------------------------------------------
VEC_TEMPLATE
TYPE&
VEC_CLASS_TEMPLATE::operator()(const int gaussLevel,int position)
   {
      assert( position < dimension ); 
      return vectorData[gaussLevel][position];
   }
//-----------------------------------------------------------------------------
VEC_TEMPLATE
template <class Element>
void
VEC_CLASS_TEMPLATE
   :: setAsTrafoFactor (ExtensionValues position, const Element& element)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel][trafoFactor] = 
                          vectorData[gaussLevel][position] * element.getGaussianPoints()[gaussLevel][weight];
        }
   }
//-----------------------------------------------------------------------------
VEC_TEMPLATE
void
VEC_CLASS_TEMPLATE
   :: computeNormalLength ( const Matrix<TYPE,numberOfGaussianPoints,D2>& m)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel][normalLength] = 
            sqrt(m(gaussLevel,dirX,dirX)*m(gaussLevel,dirX,dirX) + m(gaussLevel,dirY,dirX)*m(gaussLevel,dirY,dirX));
        }
   }
//-----------------------------------------------------------------------------
VEC_TEMPLATE
void
VEC_CLASS_TEMPLATE
   :: computeNormalLength ( const Matrix<TYPE,numberOfGaussianPoints,D3>& m)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel][normalLength] = 
            sqrt( 
                  pow(m(gaussLevel,dirY,dirX)*m(gaussLevel,dirZ,dirY)-m(gaussLevel,dirZ,dirX)*m(gaussLevel,dirY,dirY),2) +
                  pow(m(gaussLevel,dirZ,dirX)*m(gaussLevel,dirX,dirY)-m(gaussLevel,dirZ,dirY)*m(gaussLevel,dirX,dirX),2) +
                  pow(m(gaussLevel,dirX,dirX)*m(gaussLevel,dirY,dirY)-m(gaussLevel,dirX,dirY)*m(gaussLevel,dirY,dirX),2) );
        }
   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
RES_VEC_CLASS_TEMPLATE(D1)
   ::Vector ()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel] = (TYPE)0;
        }
#endif
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
RES_VEC_CLASS_TEMPLATE(D1)&
RES_VEC_CLASS_TEMPLATE(D1)
   ::operator= (const RES_VEC_CLASS_TEMPLATE(D1)& v)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel] = v.vectorData[gaussLevel];
        }
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
const TYPE& 
RES_VEC_CLASS_TEMPLATE(D1)
   ::operator()(const int gaussLevel, const int direction) const
   {
      assert( (int) direction < 1 ); 
      return vectorData[gaussLevel];
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
TYPE& 
RES_VEC_CLASS_TEMPLATE(D1)
   ::operator()(const int gaussLevel, const int direction)
   {
      assert( (int) direction < 1 ); 
      return vectorData[gaussLevel];
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
template <class Expr, class Element>
void
RES_VEC_CLASS_TEMPLATE(D1)  
   ::operator() (const Expr& expr, const Element& element)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel] = expr.template evaluateGaussianPoints<0,0,0,dirX>(element,gaussLevel);
        }
   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
RES_VEC_CLASS_TEMPLATE(D2)
   ::Vector()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel][0] = (TYPE)0;
          vectorData[gaussLevel][1] = (TYPE)0;
        }
#endif
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
RES_VEC_CLASS_TEMPLATE(D2)&
RES_VEC_CLASS_TEMPLATE(D2)
   ::operator= (const RES_VEC_CLASS_TEMPLATE(D2)& v)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel][0] = v.vectorData[gaussLevel][0];
          vectorData[gaussLevel][1] = v.vectorData[gaussLevel][1];
        }
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
TYPE& 
RES_VEC_CLASS_TEMPLATE(D2)
   ::operator()(const int gaussLevel, const int direction)
   {
      assert( (int) direction < 2 ); 
      return vectorData[gaussLevel][direction];
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
const TYPE& 
RES_VEC_CLASS_TEMPLATE(D2)
   ::operator()(const int gaussLevel, const int direction) const
   {
      assert( (int) direction < 2 ); 
      return vectorData[gaussLevel][direction];
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
template <class Expr, class Element>
void
RES_VEC_CLASS_TEMPLATE(D2)  
   ::operator() (const Expr& expr, const Element& element)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel][dirX] = expr.template evaluateGaussianPoints<0,0,0,dirX>(element,gaussLevel);
          vectorData[gaussLevel][dirY] = expr.template evaluateGaussianPoints<0,0,0,dirY>(element,gaussLevel);
        }
   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
RES_VEC_CLASS_TEMPLATE(D3)
   ::Vector()
   {
#ifdef INITIALIZE_ARRAYS
     for (int i=0; i < numberOfGaussianPoints; ++i)
        {
           vectorData[i][0] = (TYPE)0;
           vectorData[i][1] = (TYPE)0;
           vectorData[i][2] = (TYPE)0;
        }
#endif
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
RES_VEC_CLASS_TEMPLATE(D3)&
RES_VEC_CLASS_TEMPLATE(D3)
   ::operator= (const RES_VEC_CLASS_TEMPLATE(D3)& v)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel][0] = v.vectorData[gaussLevel][0];
          vectorData[gaussLevel][1] = v.vectorData[gaussLevel][1];
          vectorData[gaussLevel][2] = v.vectorData[gaussLevel][2];
        }
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
const TYPE& 
RES_VEC_CLASS_TEMPLATE(D3)
   ::operator()(const int gaussLevel, const int direction) const
   {
      assert( (int) direction < 3 ); 
      return vectorData[gaussLevel][direction];
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
TYPE& 
RES_VEC_CLASS_TEMPLATE(D3)
   ::operator()(const int gaussLevel, const int direction)
   {
      assert( (int) direction < 3 ); 
      return vectorData[gaussLevel][direction];
   }
//-----------------------------------------------------------------------------
RES_VEC_TEMPLATE
template <class Expr, class Element>
void
RES_VEC_CLASS_TEMPLATE(D3)  
   ::operator() (const Expr& expr, const Element& element)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          vectorData[gaussLevel][dirX] = expr.template evaluateGaussianPoints<0,0,0,dirX>(element,gaussLevel);
          vectorData[gaussLevel][dirY] = expr.template evaluateGaussianPoints<0,0,0,dirY>(element,gaussLevel);
          vectorData[gaussLevel][dirZ] = expr.template evaluateGaussianPoints<0,0,0,dirZ>(element,gaussLevel);
        }
   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
RES_MAT_CLASS_TEMPLATE(D1)
   ::Matrix()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          matrixData[gaussLevel] = (TYPE)0;
        }
#endif
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
RES_MAT_CLASS_TEMPLATE(D1)&
RES_MAT_CLASS_TEMPLATE(D1)
   ::operator= (const RES_MAT_CLASS_TEMPLATE(D1)& m)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          matrixData[gaussLevel] = m.matrixData[gaussLevel];
        }
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
const TYPE& 
RES_MAT_CLASS_TEMPLATE(D1)
   ::operator()(const int gaussLevel, const int direction1, const int direction2) const
   {
      return matrixData[gaussLevel];
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
TYPE& 
RES_MAT_CLASS_TEMPLATE(D1)
   ::operator()(const int gaussLevel, const int direction1, const int direction2)
   {
      return matrixData[gaussLevel];
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
template <class Expr, class Element>
void
RES_MAT_CLASS_TEMPLATE(D1)  
   ::operator() (const Expr& expr, const Element& element)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          matrixData[gaussLevel] = expr.template evaluateGaussianPoints<1,0,0,dirX>(element,gaussLevel);
        }
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
void
RES_MAT_CLASS_TEMPLATE(D1)  
   ::computeInverseAndDeterminant( const RES_MAT_CLASS_TEMPLATE(D1)& m, 
                                         RES_VEC_CLASS_TEMPLATE(_extensionSize)&v)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           v(gaussLevel,determinant) = fabs(m(gaussLevel));
#ifdef ASSERTS
           assert ( 1.e-10 < v(gaussLevel,determinant) );
#endif
           matrixData[gaussLevel] = 1./m(gaussLevel);
        }
   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
RES_MAT_CLASS_TEMPLATE(D2)
   ::Matrix()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          matrixData[gaussLevel][0][0] = (TYPE)0;
          matrixData[gaussLevel][0][1] = (TYPE)0;
          matrixData[gaussLevel][1][0] = (TYPE)0;
          matrixData[gaussLevel][1][1] = (TYPE)0;
        }
#endif
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
RES_MAT_CLASS_TEMPLATE(D2)&
RES_MAT_CLASS_TEMPLATE(D2)
   ::operator= (const RES_MAT_CLASS_TEMPLATE(D2)& m)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          matrixData[gaussLevel][0][0] = m.matrixData[gaussLevel][0][0];
          matrixData[gaussLevel][0][1] = m.matrixData[gaussLevel][0][1];
          matrixData[gaussLevel][1][0] = m.matrixData[gaussLevel][1][0];
          matrixData[gaussLevel][1][1] = m.matrixData[gaussLevel][1][1];
        }
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
const TYPE& 
RES_MAT_CLASS_TEMPLATE(D2)
   ::operator()(const int gaussLevel, const int direction1, const int direction2) const
   {
      assert( (int) direction1 < 2 && (int)direction2 < 2 ); 
      return matrixData[gaussLevel][direction1][direction2];
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
TYPE& 
RES_MAT_CLASS_TEMPLATE(D2)
   ::operator()(const int gaussLevel, const int direction1, const int direction2)
   {
      assert( (int) direction1 < 2 && (int)direction2 < 2 ); 
      return matrixData[gaussLevel][direction1][direction2];
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
template <class Expr, class Element>
void
RES_MAT_CLASS_TEMPLATE(D2)  
   ::operator() (const Expr& expr, const Element& element)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          matrixData[gaussLevel][dirX][dirX] = expr.template evaluateGaussianPoints<1,0,0,dirX>(element,gaussLevel);
          matrixData[gaussLevel][dirX][dirY] = expr.template evaluateGaussianPoints<0,1,0,dirX>(element,gaussLevel);
          matrixData[gaussLevel][dirY][dirX] = expr.template evaluateGaussianPoints<1,0,0,dirY>(element,gaussLevel);
          matrixData[gaussLevel][dirY][dirY] = expr.template evaluateGaussianPoints<0,1,0,dirY>(element,gaussLevel);
        }
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
void
RES_MAT_CLASS_TEMPLATE(D2)  
   ::computeInverseAndDeterminant( const RES_MAT_CLASS_TEMPLATE(D2)& m, 
                                         RES_VEC_CLASS_TEMPLATE(_extensionSize)&v)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           v(gaussLevel,determinant) = fabs(
               m(gaussLevel,dirX,dirX)*m(gaussLevel,dirY,dirY) - 
               m(gaussLevel,dirY,dirX)*m(gaussLevel,dirX,dirY) );
#ifdef ASSERTS
           assert ( 1.e-10 < v(gaussLevel,determinant) );
#endif
           matrixData[gaussLevel][dirX][dirX] = m(gaussLevel,dirY,dirY) / v(gaussLevel,determinant) ;
           matrixData[gaussLevel][dirX][dirY] = -m(gaussLevel,dirX,dirY) / v(gaussLevel,determinant) ;
           matrixData[gaussLevel][dirY][dirX] = -m(gaussLevel,dirY,dirX) / v(gaussLevel,determinant) ;
           matrixData[gaussLevel][dirY][dirY] = m(gaussLevel,dirX,dirX) / v(gaussLevel,determinant) ;
        }
   }
//-----------------------------------------------------------------------------
#if 0
RES_MAT_TEMPLATE
template <class Expr, typename baseTYPE>
void 
RES_MAT_CLASS_TEMPLATE(D2)  
   ::operator() (const FunctionExpr<Expr>& expr, baseTYPE **_gaussianPoints)
   {
      for (int gaussLevel = 0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
         {
            matrixData[gaussLevel][0] = (d_dx(expr))(_gaussianPoints[gaussLevel]);
            matrixData[gaussLevel][1] = (d_dy(expr))(_gaussianPoints[gaussLevel]);
         }
   }
#endif
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
RES_MAT_CLASS_TEMPLATE(D3)
   ::Matrix()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           matrixData[gaussLevel][0][0] = (TYPE)0;
           matrixData[gaussLevel][0][1] = (TYPE)0;
           matrixData[gaussLevel][0][2] = (TYPE)0;
           matrixData[gaussLevel][1][0] = (TYPE)0;
           matrixData[gaussLevel][1][1] = (TYPE)0;
           matrixData[gaussLevel][1][2] = (TYPE)0;
           matrixData[gaussLevel][2][0] = (TYPE)0;
           matrixData[gaussLevel][2][1] = (TYPE)0;
           matrixData[gaussLevel][2][2] = (TYPE)0;
        }
#endif
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
RES_MAT_CLASS_TEMPLATE(D3)&
RES_MAT_CLASS_TEMPLATE(D3)
   ::operator= (const RES_MAT_CLASS_TEMPLATE(D3)& m)
   {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          matrixData[gaussLevel][0][0] = m.matrixData[gaussLevel][0][0];
          matrixData[gaussLevel][0][1] = m.matrixData[gaussLevel][0][1];
          matrixData[gaussLevel][0][2] = m.matrixData[gaussLevel][0][2];
          matrixData[gaussLevel][1][0] = m.matrixData[gaussLevel][1][0];
          matrixData[gaussLevel][1][1] = m.matrixData[gaussLevel][1][1];
          matrixData[gaussLevel][1][2] = m.matrixData[gaussLevel][1][2];
          matrixData[gaussLevel][2][0] = m.matrixData[gaussLevel][2][0];
          matrixData[gaussLevel][2][1] = m.matrixData[gaussLevel][2][1];
          matrixData[gaussLevel][2][2] = m.matrixData[gaussLevel][2][2];
        }
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
const TYPE& 
RES_MAT_CLASS_TEMPLATE(D3)
   ::operator()(const int gaussLevel, const int direction1, const int direction2) const
   {
      assert( (int) direction1 < 3 && (int)direction2 < 3 ); 
      return matrixData[gaussLevel][direction1][direction2];
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
TYPE& 
RES_MAT_CLASS_TEMPLATE(D3)
   ::operator()(const int gaussLevel, const int direction1, const int direction2)
   {
      assert( (int) direction1 < 3 && (int)direction2 < 3 ); 
      return matrixData[gaussLevel][direction1][direction2];
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
template <class Expr, class Element>
void
RES_MAT_CLASS_TEMPLATE(D3)  
   ::operator() (const Expr& expr, const Element& element)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          matrixData[gaussLevel][dirX][dirX] = expr.template evaluateGaussianPoints<1,0,0,dirX>(element,gaussLevel);
          matrixData[gaussLevel][dirX][dirY] = expr.template evaluateGaussianPoints<0,1,0,dirX>(element,gaussLevel);
          matrixData[gaussLevel][dirX][dirZ] = expr.template evaluateGaussianPoints<0,0,1,dirX>(element,gaussLevel);

          matrixData[gaussLevel][dirY][dirX] = expr.template evaluateGaussianPoints<1,0,0,dirY>(element,gaussLevel);
          matrixData[gaussLevel][dirY][dirY] = expr.template evaluateGaussianPoints<0,1,0,dirY>(element,gaussLevel);
          matrixData[gaussLevel][dirY][dirZ] = expr.template evaluateGaussianPoints<0,0,1,dirY>(element,gaussLevel);

          matrixData[gaussLevel][dirZ][dirX] = expr.template evaluateGaussianPoints<1,0,0,dirZ>(element,gaussLevel);
          matrixData[gaussLevel][dirZ][dirY] = expr.template evaluateGaussianPoints<0,1,0,dirZ>(element,gaussLevel);
          matrixData[gaussLevel][dirZ][dirZ] = expr.template evaluateGaussianPoints<0,0,1,dirZ>(element,gaussLevel);
        }
   }
//-----------------------------------------------------------------------------
RES_MAT_TEMPLATE
void
RES_MAT_CLASS_TEMPLATE(D3)  
   ::computeInverseAndDeterminant( const RES_MAT_CLASS_TEMPLATE(D3)& m, 
                                         RES_VEC_CLASS_TEMPLATE(_extensionSize)&v)
   { 
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           matrixData[gaussLevel][dirX][dirX] = 
               m(gaussLevel,dirY,dirY)*m(gaussLevel,dirZ,dirZ) - m(gaussLevel,dirY,dirZ)*m(gaussLevel,dirZ,dirY);
           matrixData[gaussLevel][dirX][dirY] = 
               m(gaussLevel,dirZ,dirY)*m(gaussLevel,dirX,dirZ) - m(gaussLevel,dirX,dirY)*m(gaussLevel,dirZ,dirZ);
           matrixData[gaussLevel][dirX][dirZ] = 
               m(gaussLevel,dirY,dirZ)*m(gaussLevel,dirX,dirY) - m(gaussLevel,dirX,dirZ)*m(gaussLevel,dirY,dirY);

           matrixData[gaussLevel][dirY][dirX] = 
               m(gaussLevel,dirY,dirZ)*m(gaussLevel,dirZ,dirX) - m(gaussLevel,dirY,dirX)*m(gaussLevel,dirZ,dirZ);
           matrixData[gaussLevel][dirY][dirY] = 
               m(gaussLevel,dirX,dirX)*m(gaussLevel,dirZ,dirZ) - m(gaussLevel,dirX,dirZ)*m(gaussLevel,dirZ,dirX);
           matrixData[gaussLevel][dirY][dirZ] = 
               m(gaussLevel,dirX,dirZ)*m(gaussLevel,dirY,dirX) - m(gaussLevel,dirX,dirX)*m(gaussLevel,dirY,dirZ);

           matrixData[gaussLevel][dirZ][dirX] = 
               m(gaussLevel,dirY,dirX)*m(gaussLevel,dirZ,dirY) - m(gaussLevel,dirY,dirY)*m(gaussLevel,dirZ,dirX);
           matrixData[gaussLevel][dirZ][dirY] = 
               m(gaussLevel,dirX,dirY)*m(gaussLevel,dirZ,dirX) - m(gaussLevel,dirX,dirX)*m(gaussLevel,dirZ,dirY);
           matrixData[gaussLevel][dirZ][dirZ] = 
               m(gaussLevel,dirX,dirX)*m(gaussLevel,dirY,dirY) - m(gaussLevel,dirX,dirY)*m(gaussLevel,dirY,dirX);

          v(gaussLevel,determinant) = fabs(
              m(gaussLevel,dirX,dirX) * matrixData[gaussLevel][dirX][dirX] +
              m(gaussLevel,dirX,dirY) * matrixData[gaussLevel][dirY][dirX] +
              m(gaussLevel,dirX,dirZ) * matrixData[gaussLevel][dirZ][dirX] );

#ifdef ASSERTS
           assert ( 1.e-10 < v(gaussLevel,determinant) );
#endif

           matrixData[gaussLevel][dirX][dirX] /= v(gaussLevel,determinant);
           matrixData[gaussLevel][dirX][dirY] /= v(gaussLevel,determinant);
           matrixData[gaussLevel][dirX][dirZ] /= v(gaussLevel,determinant);
           matrixData[gaussLevel][dirY][dirX] /= v(gaussLevel,determinant);
           matrixData[gaussLevel][dirY][dirY] /= v(gaussLevel,determinant);
           matrixData[gaussLevel][dirY][dirZ] /= v(gaussLevel,determinant);
           matrixData[gaussLevel][dirZ][dirX] /= v(gaussLevel,determinant);
           matrixData[gaussLevel][dirZ][dirY] /= v(gaussLevel,determinant);
           matrixData[gaussLevel][dirZ][dirZ] /= v(gaussLevel,determinant);
        }
   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
#undef VEC_TEMPLATE
#undef VEC_CLASS_TEMPLATE
#undef MAT_TEMPLATE
#undef MAT_CLASS_TEMPLATE
#undef RES_VEC_TEMPLATE
#undef RES_VEC_CLASS_TEMPLATE
#undef RES_MAT_TEMPLATE
#undef RES_MAT_CLASS_TEMPLATE
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
Colsamm_LinA::Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>
   ::Function_Values()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
              {
                functionData[gaussLevel][functionIter] = (TYPE)0;
              }
        }
#endif
 
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
template <typename baseTYPE>
void
Colsamm_LinA::Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>
   ::operator() (const std::vector<Base**> &expr, baseTYPE** gaussianPoints)
   {
     for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
        {
          functionData[numberOfGaussianPoints-1][functionIter] = (*expr[functionIter][fn])(gaussianPoints[numberOfGaussianPoints-1]);
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
template <class Expr, typename baseTYPE>
void
Colsamm_LinA::Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>
   ::operator() (const int functionIter, const FunctionExpr<Expr>& expr, baseTYPE** gaussianPoints)
   {
     const Expr& expr_(expr);
     assert (functionIter < numberOfBasisFunctions);
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints - 1; ++gaussLevel)
        {
           functionData[gaussLevel][functionIter] = expr_(gaussianPoints[gaussLevel]);
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
const TYPE& 
Colsamm_LinA::Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>
   ::operator()(const int gaussLevel, const int functionIter) const
   {
      assert (functionIter < numberOfBasisFunctions);
      return functionData[gaussLevel][functionIter];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
TYPE& 
Colsamm_LinA::Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>
   ::operator()(const int gaussLevel,const int functionIter)
   {
      assert (functionIter < numberOfBasisFunctions);
      return functionData[gaussLevel][functionIter];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>
   ::applyVectorialChainRule(const Matrix<TYPE,numberOfGaussianPoints-1,D1>& m,
		     const Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>& f,
		     const TYPE* edgeLength)
   {
      for (int gaussLevel = 0 ; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel )
         {
           for (int functionIter = 0; functionIter < numberOfBasisFunctions; ++functionIter)
              {
                functionData[gaussLevel][functionIter] = f(gaussLevel,functionIter) * m(gaussLevel,dirX,dirX) *
                                                           edgeLength[functionIter];
              }
         }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>
   ::applyVectorialChainRule(const Matrix<TYPE,numberOfGaussianPoints-1,D2>& m,
		     const Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>& f,
		     const TYPE* edgeLength)
   {
      for (int gaussLevel = 0 ; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel )
         {
           for (int functionIter = 0; functionIter < numberOfBasisFunctions/2; ++functionIter)
              {
                 functionData[gaussLevel][2*functionIter+dirX] = 
                         ( f(gaussLevel,2*functionIter+dirX) * m(gaussLevel,dirX,dirX) + 
                           f(gaussLevel,2*functionIter+dirY) * m(gaussLevel,dirY,dirX) ) * edgeLength[functionIter];

                 functionData[gaussLevel][2*functionIter+dirY] = 
                         ( f(gaussLevel,2*functionIter+dirX) * m(gaussLevel,dirX,dirY) + 
                           f(gaussLevel,2*functionIter+dirY) * m(gaussLevel,dirY,dirY) ) * edgeLength[functionIter];
              }
         }
    }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>
   ::applyVectorialChainRule(const Matrix<TYPE,numberOfGaussianPoints-1,D3>& m,
		     const Function_Values<TYPE,numberOfGaussianPoints,numberOfBasisFunctions>& f,
		     const TYPE* edgeLength)
   {
      for (int gaussLevel = 0 ; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel )
         {
           for (int functionIter = 0; functionIter < numberOfBasisFunctions/3; ++functionIter)
              {
                 functionData[gaussLevel][3*functionIter+dirX] = 
                          (f(gaussLevel,3*functionIter+dirX) * m(gaussLevel,dirX,dirX) + 
                           f(gaussLevel,3*functionIter+dirY) * m(gaussLevel,dirX,dirY) + 
                           f(gaussLevel,3*functionIter+dirZ) * m(gaussLevel,dirX,dirZ) ) * edgeLength[functionIter];

                 functionData[gaussLevel][3*functionIter+dirY] = 
                          (f(gaussLevel,3*functionIter+dirX) * m(gaussLevel,dirY,dirX) + 
                           f(gaussLevel,3*functionIter+dirY) * m(gaussLevel,dirY,dirY) + 
                           f(gaussLevel,3*functionIter+dirZ) * m(gaussLevel,dirY,dirZ) ) * edgeLength[functionIter];

                 functionData[gaussLevel][3*functionIter+dirZ] = 
                          (f(gaussLevel,3*functionIter+dirX) * m(gaussLevel,dirZ,dirX) + 
                           f(gaussLevel,3*functionIter+dirY) * m(gaussLevel,dirZ,dirY) + 
                           f(gaussLevel,3*functionIter+dirZ) * m(gaussLevel,dirZ,dirZ) ) * edgeLength[functionIter];
              }
         }
   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>
   ::Function_First_Derivatives()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
              {
                functionData[gaussLevel][functionIter] = (TYPE)0;
              }
        }
#endif
 
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
template <class Expr, typename baseTYPE>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>
   ::operator() (const int functionIter, const FunctionExpr<Expr>& expr, baseTYPE** gaussianPoints)
   {
     const Expr& expr_(expr);
     assert (functionIter < numberOfBasisFunctions);
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints - 1; ++gaussLevel)
        {
           functionData[gaussLevel][functionIter] = (d_dx(expr_))(gaussianPoints[gaussLevel]);
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
template <typename baseTYPE>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>
   ::operator() (const std::vector<Base**> &expr, baseTYPE** gaussianPoints)
   {
     for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
        {
          functionData[numberOfGaussianPoints-1][functionIter] = ((*expr[functionIter])[dx])(gaussianPoints[numberOfGaussianPoints-1]);
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
const TYPE& 
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>
   ::operator()(const int gaussLevel, const int functionIter, int direction) const
   {
      assert (functionIter < numberOfBasisFunctions);
      return functionData[gaussLevel][functionIter];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
TYPE& 
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>
   ::operator()(const int gaussLevel,const int functionIter, int direction)
   {
      assert (functionIter < numberOfBasisFunctions);
      return functionData[gaussLevel][functionIter];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>
   ::applyChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D1>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>& v)
  {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel)
        {
           for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
              {
                functionData[gaussLevel][functionIter] = 
                       m(gaussLevel,dirX) * v(gaussLevel,functionIter,dirX) ;
              }
        }

  }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>
   :: applyVectorialChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D1>& m,
                               const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D1>& v,
                               const TYPE* edgeLength )
   {
     for (int gaussLevel = 0 ; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel )
        {
          for (int functionIter = 0; functionIter < numberOfBasisFunctions; ++functionIter)
             {
                functionData[gaussLevel][functionIter] = v(gaussLevel,functionIter,dirX) * m(gaussLevel,dirX,dirX) *
                                                           edgeLength[functionIter];
             }
        }

   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
   ::Function_First_Derivatives()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
              {
                functionData[gaussLevel][functionIter][0] = (TYPE)0;
                functionData[gaussLevel][functionIter][1] = (TYPE)0;
              }
        }
#endif
 
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
template <class Expr, typename baseTYPE>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
   ::operator() (const int functionIter, const FunctionExpr<Expr>& expr, baseTYPE** gaussianPoints)
   {
     const Expr& expr_(expr);
     assert (functionIter < numberOfBasisFunctions);
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints - 1; ++gaussLevel)
        {
           functionData[gaussLevel][functionIter][dirX] = (d_dx(expr_))(gaussianPoints[gaussLevel]);
           functionData[gaussLevel][functionIter][dirY] = (d_dy(expr_))(gaussianPoints[gaussLevel]);
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
template <typename baseTYPE>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
   ::operator() (const std::vector<Base**> &expr, baseTYPE** gaussianPoints)
   {
     for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
        {
          functionData[numberOfGaussianPoints-1][functionIter][dirX] = ((*expr[functionIter][dx]))(gaussianPoints[numberOfGaussianPoints-1]);
          functionData[numberOfGaussianPoints-1][functionIter][dirY] = ((*expr[functionIter][dy]))(gaussianPoints[numberOfGaussianPoints-1]);
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
const TYPE& 
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
   ::operator()(const int gaussLevel, const int functionIter, int direction) const
   {
      assert (functionIter < numberOfBasisFunctions);
      assert (direction < 2);
      return functionData[gaussLevel][functionIter][direction];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
TYPE& 
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
   ::operator()(const int gaussLevel,const int functionIter, int direction)
   {
      assert (functionIter < numberOfBasisFunctions);
      assert (direction < 2);
      return functionData[gaussLevel][functionIter][direction];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
   ::applyChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D2>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>& v)
  {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel)
        {
           for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
              {
                functionData[gaussLevel][functionIter][dirX] = 
                       m(gaussLevel,dirX,dirX) * v(gaussLevel,functionIter,dirX) + 
                       m(gaussLevel,dirY,dirX) * v(gaussLevel,functionIter,dirY) ; 
                functionData[gaussLevel][functionIter][dirY] = 
                       m(gaussLevel,dirX,dirY) * v(gaussLevel,functionIter,dirX) + 
                       m(gaussLevel,dirY,dirY) * v(gaussLevel,functionIter,dirY) ; 
              }
        }

  }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
   ::applyChainRuleLast( const Matrix<TYPE,numberOfGaussianPoints-1,D2>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>& v)
  {
      for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
         {
           functionData[numberOfGaussianPoints-1][functionIter][dirX] = 
                  m(numberOfGaussianPoints-1,dirX,dirX) * v(numberOfGaussianPoints-1,functionIter,dirX) + 
                  m(numberOfGaussianPoints-1,dirY,dirX) * v(numberOfGaussianPoints-1,functionIter,dirY) ; 
           functionData[numberOfGaussianPoints-1][functionIter][dirY] = 
                  m(numberOfGaussianPoints-1,dirX,dirY) * v(numberOfGaussianPoints-1,functionIter,dirX) + 
                  m(numberOfGaussianPoints-1,dirY,dirY) * v(numberOfGaussianPoints-1,functionIter,dirY) ; 
         }
  }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>
   :: applyVectorialChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D2>& m,
                               const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D2>& v,
                               const TYPE* edgeLength )
   {
     for (int gaussLevel = 0 ; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel )
        {
          for (int functionIter = 0; functionIter < numberOfBasisFunctions/2; ++functionIter)
             {
                functionData[gaussLevel][2*functionIter+dirX][dirX] = 
                        ( v(gaussLevel,2*functionIter+dirX,dirX) * m(gaussLevel,dirX,dirX) +
                          v(gaussLevel,2*functionIter+dirY,dirX) * m(gaussLevel,dirY,dirX) ) * edgeLength[functionIter];

                functionData[gaussLevel][2*functionIter+dirX][dirY] = 
                        ( v(gaussLevel,2*functionIter+dirX,dirX) * m(gaussLevel,dirX,dirY) +
                          v(gaussLevel,2*functionIter+dirY,dirX) * m(gaussLevel,dirY,dirY) ) * edgeLength[functionIter];

                functionData[gaussLevel][2*functionIter+dirY][dirX] = 
                        ( v(gaussLevel,2*functionIter+dirX,dirY) * m(gaussLevel,dirX,dirX) +
                          v(gaussLevel,2*functionIter+dirY,dirY) * m(gaussLevel,dirY,dirX) ) * edgeLength[functionIter];

                functionData[gaussLevel][2*functionIter+dirY][dirY] = 
                        ( v(gaussLevel,2*functionIter+dirX,dirY) * m(gaussLevel,dirX,dirY) +
                          v(gaussLevel,2*functionIter+dirY,dirY) * m(gaussLevel,dirY,dirY) ) * edgeLength[functionIter];

             }
        }

   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
   ::Function_First_Derivatives()
   {
#ifdef INITIALIZE_ARRAYS
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
           for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
              {
                functionData[gaussLevel][functionIter][0] = (TYPE)0;
                functionData[gaussLevel][functionIter][1] = (TYPE)0;
                functionData[gaussLevel][functionIter][2] = (TYPE)0;
              }
        }
#endif
 
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
template <class Expr, typename baseTYPE>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
   ::operator() (const int functionIter, const FunctionExpr<Expr>& expr, baseTYPE** gaussianPoints)
   {
     const Expr& expr_(expr);
     assert (functionIter < numberOfBasisFunctions);
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints - 1; ++gaussLevel)
        {
           functionData[gaussLevel][functionIter][dirX] = (d_dx(expr_))(gaussianPoints[gaussLevel]);
           functionData[gaussLevel][functionIter][dirY] = (d_dy(expr_))(gaussianPoints[gaussLevel]);
           functionData[gaussLevel][functionIter][dirZ] = (d_dz(expr_))(gaussianPoints[gaussLevel]);
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
template <typename baseTYPE>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
   ::operator() (const std::vector<Base**> &expr, baseTYPE** gaussianPoints)
   {
     for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
        {
          functionData[numberOfGaussianPoints-1][functionIter][dirX] = ((*expr[functionIter][dx]))(gaussianPoints[numberOfGaussianPoints-1]);
          functionData[numberOfGaussianPoints-1][functionIter][dirY] = ((*expr[functionIter][dy]))(gaussianPoints[numberOfGaussianPoints-1]);
          functionData[numberOfGaussianPoints-1][functionIter][dirZ] = ((*expr[functionIter][dz]))(gaussianPoints[numberOfGaussianPoints-1]);
        }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
const TYPE& 
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
   ::operator()(const int gaussLevel, const int functionIter, int direction) const
   {
      assert (functionIter < numberOfBasisFunctions);
      assert (direction < 3);
      return functionData[gaussLevel][functionIter][direction];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
TYPE& 
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
   ::operator()(const int gaussLevel,const int functionIter, int direction)
   {
      assert (functionIter < numberOfBasisFunctions);
      assert (direction < 3);
      return functionData[gaussLevel][functionIter][direction];
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
   ::applyChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D3>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>& v)
  {
     for (int gaussLevel=0; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel)
        {
           for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
              {
                functionData[gaussLevel][functionIter][dirX] = 
                       m(gaussLevel,dirX,dirX) * v(gaussLevel,functionIter,dirX) + 
                       m(gaussLevel,dirY,dirX) * v(gaussLevel,functionIter,dirY) + 
                       m(gaussLevel,dirZ,dirX) * v(gaussLevel,functionIter,dirZ) ; 
                functionData[gaussLevel][functionIter][dirY] = 
                       m(gaussLevel,dirX,dirY) * v(gaussLevel,functionIter,dirX) + 
                       m(gaussLevel,dirY,dirY) * v(gaussLevel,functionIter,dirY) + 
                       m(gaussLevel,dirZ,dirY) * v(gaussLevel,functionIter,dirZ) ; 
                functionData[gaussLevel][functionIter][dirZ] = 
                       m(gaussLevel,dirX,dirZ) * v(gaussLevel,functionIter,dirX) + 
                       m(gaussLevel,dirY,dirZ) * v(gaussLevel,functionIter,dirY) + 
                       m(gaussLevel,dirZ,dirZ) * v(gaussLevel,functionIter,dirZ) ; 
              }
        }

  }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
   ::applyChainRuleLast( const Matrix<TYPE,numberOfGaussianPoints-1,D3>& m,
                const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>& v)
  {
      for (int functionIter=0; functionIter < numberOfBasisFunctions; ++functionIter)
         {
           functionData[numberOfGaussianPoints-1][functionIter][dirX] = 
                  m(numberOfGaussianPoints-1,dirX,dirX) * v(numberOfGaussianPoints-1,functionIter,dirX) + 
                  m(numberOfGaussianPoints-1,dirY,dirX) * v(numberOfGaussianPoints-1,functionIter,dirY) + 
                  m(numberOfGaussianPoints-1,dirZ,dirX) * v(numberOfGaussianPoints-1,functionIter,dirZ) ; 
           functionData[numberOfGaussianPoints-1][functionIter][dirY] = 
                  m(numberOfGaussianPoints-1,dirX,dirY) * v(numberOfGaussianPoints-1,functionIter,dirX) + 
                  m(numberOfGaussianPoints-1,dirY,dirY) * v(numberOfGaussianPoints-1,functionIter,dirY) + 
                  m(numberOfGaussianPoints-1,dirZ,dirY) * v(numberOfGaussianPoints-1,functionIter,dirZ) ; 
           functionData[numberOfGaussianPoints-1][functionIter][dirZ] = 
                  m(numberOfGaussianPoints-1,dirX,dirZ) * v(numberOfGaussianPoints-1,functionIter,dirX) + 
                  m(numberOfGaussianPoints-1,dirY,dirZ) * v(numberOfGaussianPoints-1,functionIter,dirY) + 
                  m(numberOfGaussianPoints-1,dirZ,dirZ) * v(numberOfGaussianPoints-1,functionIter,dirZ) ; 
         }
  }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int numberOfBasisFunctions>
void
Colsamm_LinA::Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>
   :: applyVectorialChainRule( const Matrix<TYPE,numberOfGaussianPoints-1,D3>& m,
                               const Function_First_Derivatives<TYPE,numberOfGaussianPoints,numberOfBasisFunctions,D3>& v,
                               const TYPE* edgeLength )
   {
     std::cout << m << std::endl;
     for (int gaussLevel = 0 ; gaussLevel < numberOfGaussianPoints-1; ++gaussLevel )
        {
          for (int functionIter = 0; functionIter < numberOfBasisFunctions/3; ++functionIter)
             {
                functionData[gaussLevel][3*functionIter+dirX][dirX] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirX) * m(gaussLevel,dirX,dirX) +
                          v(gaussLevel,3*functionIter+dirY,dirX) * m(gaussLevel,dirY,dirX) +
                          v(gaussLevel,3*functionIter+dirZ,dirX) * m(gaussLevel,dirZ,dirX) ) * edgeLength[functionIter];

                functionData[gaussLevel][3*functionIter+dirX][dirY] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirX) * m(gaussLevel,dirX,dirY) +
                          v(gaussLevel,3*functionIter+dirY,dirX) * m(gaussLevel,dirY,dirY) +
                          v(gaussLevel,3*functionIter+dirZ,dirX) * m(gaussLevel,dirZ,dirY) ) * edgeLength[functionIter];

                functionData[gaussLevel][3*functionIter+dirX][dirZ] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirX) * m(gaussLevel,dirX,dirZ) +
                          v(gaussLevel,3*functionIter+dirY,dirX) * m(gaussLevel,dirY,dirZ) +
                          v(gaussLevel,3*functionIter+dirZ,dirX) * m(gaussLevel,dirZ,dirZ) ) * edgeLength[functionIter];

                functionData[gaussLevel][3*functionIter+dirY][dirX] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirY) * m(gaussLevel,dirX,dirX) +
                          v(gaussLevel,3*functionIter+dirY,dirY) * m(gaussLevel,dirY,dirX) +
                          v(gaussLevel,3*functionIter+dirZ,dirY) * m(gaussLevel,dirZ,dirX) ) * edgeLength[functionIter];

                functionData[gaussLevel][3*functionIter+dirY][dirY] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirY) * m(gaussLevel,dirX,dirY) +
                          v(gaussLevel,3*functionIter+dirY,dirY) * m(gaussLevel,dirY,dirY) +
                          v(gaussLevel,3*functionIter+dirZ,dirY) * m(gaussLevel,dirZ,dirY) ) * edgeLength[functionIter];

                functionData[gaussLevel][3*functionIter+dirY][dirZ] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirY) * m(gaussLevel,dirX,dirZ) +
                          v(gaussLevel,3*functionIter+dirY,dirY) * m(gaussLevel,dirY,dirZ) +
                          v(gaussLevel,3*functionIter+dirZ,dirY) * m(gaussLevel,dirZ,dirZ) ) * edgeLength[functionIter];

                functionData[gaussLevel][3*functionIter+dirZ][dirX] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirZ) * m(gaussLevel,dirX,dirX) +
                          v(gaussLevel,3*functionIter+dirY,dirZ) * m(gaussLevel,dirY,dirX) +
                          v(gaussLevel,3*functionIter+dirZ,dirZ) * m(gaussLevel,dirZ,dirX) ) * edgeLength[functionIter];

                functionData[gaussLevel][3*functionIter+dirZ][dirY] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirZ) * m(gaussLevel,dirX,dirY) +
                          v(gaussLevel,3*functionIter+dirY,dirZ) * m(gaussLevel,dirY,dirY) +
                          v(gaussLevel,3*functionIter+dirZ,dirZ) * m(gaussLevel,dirZ,dirY) ) * edgeLength[functionIter];

                functionData[gaussLevel][3*functionIter+dirZ][dirZ] = 
                        ( v(gaussLevel,3*functionIter+dirX,dirZ) * m(gaussLevel,dirX,dirZ) +
                          v(gaussLevel,3*functionIter+dirY,dirZ) * m(gaussLevel,dirY,dirZ) +
                          v(gaussLevel,3*functionIter+dirZ,dirZ) * m(gaussLevel,dirZ,dirZ) ) * edgeLength[functionIter];
                for (int i=0; i< 3; ++i)
                   {
                      for (int j=0; j< 3; ++j)
                         {
                            std::cout << functionData[gaussLevel][3*functionIter+i][j] << " " ;
                         }
                       std::cout << std::endl;
                   }
             }
        }

   }
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int dimension>
inline std::ostream& operator << (std::ostream& os, const Colsamm_LinA::Vector<TYPE,numberOfGaussianPoints,dimension>& v)
   {
     for (int gaussLevel= 0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          os << "Gauss-Level " << gaussLevel << ": "<< std::endl;
          os << "    (";
          for ( int dim=0; dim < dimension-1; ++dim)
             {
                os << v(gaussLevel,dim) << ", " ;
             }
           os << v(gaussLevel,dimension-1) << " ) " << std::endl;
        } 
      return os;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints, int dimension>
inline std::ostream& operator << (std::ostream& os, const Colsamm_LinA::Matrix<TYPE,numberOfGaussianPoints,dimension>& m)
   {
     for (int gaussLevel= 0; gaussLevel < numberOfGaussianPoints; ++gaussLevel)
        {
          os << "Gauss-Level " << gaussLevel << ": "<< std::endl;
          for ( int dim1=0; dim1 < dimension; ++dim1)
             {
                os << "    | ";
                for ( int dim2=0; dim2 < dimension-1; ++dim2)
                   {
                      os << m(gaussLevel,dim1,dim2) << ", " ;
                   }
                os << m(gaussLevel,dim1,dimension-1) << " | " << std::endl;
             }
        } 
      return os;
   }
//-----------------------------------------------------------------------------

