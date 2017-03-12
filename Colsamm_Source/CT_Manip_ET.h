//==============================================================================
//
//  $Id: CT_Manip_ET.h,v 1.12 2006/07/19 09:05:18 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
#define Define_Integrand(A, var)  \
		        __typeof__(A) var (A)
//-----------------------------------------------------------------------------
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
// Constructor functions to build the corresponding expression objects!
//-----------------------------------------------------------------------------
typedef BasisFunction<BasisSet_1,0> v_;
typedef BasisFunction<BasisSet_1,1> w_;
typedef BasisFunction<BasisSet_2,0> p_;
//-----------------------------------------------------------------------------
  template <class T>
  inline VertexVector<T> 
   vec(const T(&d));
//-----------------------------------------------------------------------------
  template <class T>
  inline FUNC1_<T> 
   func(const typename FunctionInterfaces<T>::functionInterfaceOne& f_);
//-----------------------------------------------------------------------------
  template <class T>
  inline FUNC2_<T> 
   func(const typename FunctionInterfaces<T>::functionInterfaceTwo& f_);
//-----------------------------------------------------------------------------
 template <class T>
 inline FUNC3_<T> 
  func(const typename FunctionInterfaces<T>::functionInterfaceThree& f_);
//-----------------------------------------------------------------------------
 template <class A> 
 inline Conjugate<A> 
  conj(const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template <class A, class B> inline Add_<A,B>
 operator + (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b);
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
 inline Sinus_<A>
  Sin (const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template <class A> 
 inline Cosinus_<A>
  Cos (const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template <class A> 
 inline SQRT_<A>
  Sqrt (const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template <class A, typename eTYPE> 
 inline POW_<A,eTYPE>
  Pow (const BasisFunctionExpr<A>& a, eTYPE exponent);
//-----------------------------------------------------------------------------
 template <class A> 
 inline EXP_<A>
  Exp (const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template<class A> 
 inline Deri_<A,dx> 
  d_dx(const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template<class A> 
 inline Deri_<A,dy>
  d_dy(const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
 template<class A> 
 inline Deri_<A,dz> 
  d_dz(const BasisFunctionExpr<A>& a);
//-----------------------------------------------------------------------------
//==============================================================================

