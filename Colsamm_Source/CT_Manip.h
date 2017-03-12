//==============================================================================
//
//  $Id: CT_Manip.h,v 1.21 2006/07/19 09:05:18 jochen Exp $
//
//==============================================================================
#define Define_Integrand(A)  (A); \
		        __typeof__(A)
template <class A>
inline void 
 print(const A& a); 
//-----------------------------------------------------------------------------
// Data structures for rearranging the expressions. This is only applied for 
// the derivatives so far. But there might be a huge potential for optimizing, 
// e.g. the multiplication of a constant with a vector, where we should not
// immediately multiply every compontent with the constant, but only apply the 
// product to the result!
//-----------------------------------------------------------------------------
template <class A> 
struct MANIP;
//-----------------------------------------------------------------------------
// Constructor functions to build the corresponding expression objects!
//-----------------------------------------------------------------------------
typedef BasisFunction<BasisSet_1,0> v_;
typedef BasisFunction<BasisSet_1,1> w_;
typedef BasisFunction<BasisSet_2,0> p_;
//-----------------------------------------------------------------------------
  template <int k,class T>
  inline VertexVector<k,T> 
   vec(const T(&d));
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC1_<k,float> 
   func(const typename FunctionInterfaces<float>::functionInterfaceOne& f_);
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC2_<k,float> 
   func(const typename FunctionInterfaces<float>::functionInterfaceTwo& f_);
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC3_<k,float> 
   func(const typename FunctionInterfaces<float>::functionInterfaceThree& f_);
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC1_<k,double> 
   func(const typename FunctionInterfaces<double>::functionInterfaceOne& f_);
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC2_<k,double> 
   func(const typename FunctionInterfaces<double>::functionInterfaceTwo& f_);
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC3_<k,double> 
   func(const typename FunctionInterfaces<double>::functionInterfaceThree& f_);
//-----------------------------------------------------------------------------
  template <int k, class T>
  inline FUNC1_<k,T> 
   func(const typename FunctionInterfaces<T>::functionInterfaceOne& f_);
//-----------------------------------------------------------------------------
  template <int k, class T>
  inline FUNC2_<k,T> 
   func(const typename FunctionInterfaces<T>::functionInterfaceTwo& f_);
//-----------------------------------------------------------------------------
 template <int k, class T>
 inline FUNC3_<k,T> 
  func(const typename FunctionInterfaces<T>::functionInterfaceThree& f_);
//-----------------------------------------------------------------------------
 template <class A> 
 inline Conjugate<A> 
  conj(const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template <class A, class B> inline Add_<A,B>
 operator + (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b);
//-----------------------------------------------------------------------------
 template <int A, int B> inline Double_<-(A+B)>
 operator + (const Double_<A>& a, const Double_<B>& b); 
//-----------------------------------------------------------------------------
 template <class A, class B> inline Sub_<A,B>
 operator - (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b);
//-----------------------------------------------------------------------------
 template <class A> inline Min_<A>
 operator - (const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template <class A, class B> inline Mult_<A,B>
 operator * (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b);
//-----------------------------------------------------------------------------
 template <class A, class B> inline Div_<A,B>
 operator / (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b);
//-----------------------------------------------------------------------------
 template <class A> 
 inline SQRT_<A>
  Sqrt (const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template <int p, class A, typename eTYPE> 
 inline POW_<p,A,eTYPE>
  Pow (const BasisFunctionExpr<A>& a, eTYPE exponent);
//-----------------------------------------------------------------------------
 template <class A> 
 inline EXP_<A>
  Exp (const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template<class A> 
 inline typename MANIP<Deri_<A,dx> >::I 
  d_dx(const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template<class A> 
 inline typename MANIP<Deri_<A,dy> >::I 
  d_dy(const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template<class A> 
 inline typename MANIP<Deri_<A,dz> >::I 
  d_dz(const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
/* A vector is a type list of its components, namely a vector of two  
 * components or a component and vector (which is a component itself)!
 * The evaluation of the operator '&' is from left to right (hopefully),
 * and thereby, we build a template type of the structure
 *   VECTOR<VECTOR<VECTOR<....>,.>,.> 
 * This means, we can only access the VECTOR from the end to the top! 
 * A vector (so far) is always a line or a column, not differentiated 
 * between them! 
*/
template <class A, class B>
struct VECTOR : public BasisFunctionExpr<VECTOR<A,B> >
   {
     enum {hasV = _MAX_(A::hasV,B::hasV),
           hasW = _MAX_(A::hasW,B::hasW), 
           hasP = _MAX_(A::hasP,B::hasP)} ;
     enum {length = A::length+1, size = 1};
     inline static void 
      print(std::ostream &os);
};
//-----------------------------------------------------------------------------
template <class A, class B>
VECTOR<A,B> 
operator &(const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b);
//-----------------------------------------------------------------------------
/* A matrix is a concatenation of two vectors or an vector and a matrix. 
 * The result of the operator '|' is from left to right. Thus, the template 
 * type we achieve is something of the type 
 * MATRIX< MATRIX < MATRIX < ..., VECTOR <..> >, VECTOR <...> >, VECTOR <...> >
*/

template <class A, class B>
struct MATRIX : public BasisFunctionExpr<MATRIX<A,B> >
   {
     enum {hasV = _MAX_(A::hasV,B::hasV),
           hasW = _MAX_(A::hasW,B::hasW), 
           hasP = _MAX_(A::hasP,B::hasP)} ;
     enum {length = B::length ,size = A::size + 1 };
     inline static void 
      print(std::ostream& os);
   };
//-----------------------------------------------------------------------------
template <class A, class B, class C, class D>
MATRIX<VECTOR<A,B>,VECTOR<C,D> > 
 operator | (const VECTOR<A,B>& a, const VECTOR<C,D>& b);
//-----------------------------------------------------------------------------
template <class A, class B, class C, class D>
MATRIX<MATRIX<A,B>,VECTOR<C,D> > 
 operator |(const MATRIX<A,B>& a, const VECTOR<C,D>& b);
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a vector with a one dimensional
 * value. The result is the vector, muliplied with the value in each component.
*/
template <class A, class B>
struct BasisFunctionExpr_VECTOR_MULT; 
//-----------------------------------------------------------------------------
template <class A, class V1, class V2>
typename BasisFunctionExpr_VECTOR_MULT< A, VECTOR<V1,V2> > :: Erg
 operator * (const BasisFunctionExpr<A>& a, const VECTOR<V1,V2>& b);
//-----------------------------------------------------------------------------
template <class A, class V1, class V2>
typename BasisFunctionExpr_VECTOR_MULT< A, VECTOR<V1,V2> > :: Erg
 operator * (const VECTOR<V1,V2>& a, const BasisFunctionExpr<A>& b);
//-----------------------------------------------------------------------------
/* This struct describes the addition of two vectors. The result is the vector, 
 * that holds the sum of the vector components. 
*/
template <class A, class B>
struct CT_VECTOR_Add;  
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename CT_VECTOR_Add<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg 
 operator+(const VECTOR<A1,A2>& a, const VECTOR<B1,B2>& b);
//-----------------------------------------------------------------------------
template <class A, class B>
struct CT_VECTOR_Sub;  
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename CT_VECTOR_Sub<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg 
 operator-(const VECTOR<A1,A2>& a, const VECTOR<B1,B2>& b);
//-----------------------------------------------------------------------------
template <class A>
struct CT_VECTOR_Min; 
//-----------------------------------------------------------------------------
template <class A, class B>
typename CT_VECTOR_Min<VECTOR<A,B> >::Erg 
 operator-(const VECTOR<A,B>& a);
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of two vectors. The result is  
 * the sum that holds the products of the vectors components. 
*/
template <class A, class B>
struct VECTOR_VECTOR_MULT; 
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename VECTOR_VECTOR_MULT<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg 
 operator*(const VECTOR<A1,A2>& a, const VECTOR<B1,B2>& b);
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a vector with a matrix!  
 * The result is a vector of the line sums . 
*/
template <class A, class B>
struct VECTOR_MATRIX_MULT;
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename VECTOR_MATRIX_MULT<VECTOR<A1,A2>,MATRIX<B1,B2> >::Erg 
 operator*(const VECTOR<A1,A2>& a, const MATRIX<B1,B2>& b);
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a matrix with a vector!  
 * The result is a vector of the line sums . 
*/
template <class A, class B>
struct MATRIX_VECTOR_MULT;
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename MATRIX_VECTOR_MULT<MATRIX<A1,A2>,VECTOR<B1,B2> >::Erg 
 operator*(const MATRIX<A1,A2>& a, const VECTOR<B1,B2>& b); 
// ----------------------------------------------------------------------------
 template <class A> 
 inline VECTOR< VECTOR<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>,
                       typename MANIP<Deri_<A,dz> >::I >
  grad(const BasisFunctionExpr<A>& a); 
// ----------------------------------------------------------------------------
 template <class A> 
 inline  VECTOR<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>
  grad2(const BasisFunctionExpr<A>& a); 
// ----------------------------------------------------------------------------
 template <class A, class B, class C> 
 inline Add_<Add_<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>,
                  typename MANIP<Deri_<A,dz> >::I >
  div(const VECTOR< VECTOR <A,B>, C >& a); 
// ----------------------------------------------------------------------------
 inline VECTOR<VECTOR< NormalComponent<0>,NormalComponent<1> >, NormalComponent<2> >
  N(); 
// ----------------------------------------------------------------------------
 inline VECTOR<VECTOR< UnitNormalComponent<0>,UnitNormalComponent<1> >, UnitNormalComponent<2> >
  N_Unit(); 
// ----------------------------------------------------------------------------

/* Concluding comments
 * At this point there is still more to introduce. A compile time data structure
 * which holds tensor data. Therefore we would need an additional operator to be 
 * overloaded, whicht has less priority than operator '&'. 
*/

//==============================================================================

