//==============================================================================
//
//  $Id: VertexComputations.h,v 1.1 2006/06/09 10:12:46 jochen Exp $
//
//==============================================================================
//----------------------------------------------------------------------------
template <class A>
const A& performFunc(const A& a) {return a;}
//----------------------------------------------------------------------------
template < class A > 
struct PolynomialExprForVertexExpr {};
//----------------------------------------------------------------------------
struct ConstantForVertexExpression ;
//----------------------------------------------------------------------------
template < int p = 1 >
struct PolynomialXForVertexExpr; 
//----------------------------------------------------------------------------
typedef PolynomialXForVertexExpr<1> _U;
typedef PolynomialXForVertexExpr<2> _U_2;
typedef PolynomialXForVertexExpr<3> _U_3;
//----------------------------------------------------------------------------
template < int p = 1 >
struct PolynomialYForVertexExpr; 
//----------------------------------------------------------------------------
typedef PolynomialYForVertexExpr<1> _V;
typedef PolynomialYForVertexExpr<2> _V_2;
typedef PolynomialYForVertexExpr<3> _V_3;
//----------------------------------------------------------------------------
template < int p = 1 >
struct PolynomialZForVertexExpr; 
//----------------------------------------------------------------------------
typedef PolynomialZForVertexExpr<1> _W;
typedef PolynomialZForVertexExpr<2> _W_2;
typedef PolynomialZForVertexExpr<3> _W_3;
//----------------------------------------------------------------------------
template < class A , class B >
 struct PolynomialExprForVertexExprAddition; 
//----------------------------------------------------------------------------
template < class A, class B >
inline PolynomialExprForVertexExprAddition <A,B>
operator + ( const PolynomialExprForVertexExpr <A>& a , 
             const PolynomialExprForVertexExpr <B>& b ); 
//----------------------------------------------------------------------------
template < class A , class B > 
struct PolynomialExprForVertexExprSubtraction; 
//-------------------------------------------------------------------------------------
template < class A , class B > 
inline PolynomialExprForVertexExprSubtraction <A,B>
operator - ( const PolynomialExprForVertexExpr <A>& a , 
             const PolynomialExprForVertexExpr <B>& b ); 
//-------------------------------------------------------------------------------------
template < class A >
struct PolynomialExprForVertexExprUnaryMinus; 
//-------------------------------------------------------------------------------------
template < class A >
inline PolynomialExprForVertexExprUnaryMinus <A> 
operator - ( const PolynomialExprForVertexExpr <A>& a ); 
//-------------------------------------------------------------------------------------
template <class A, class B, int ds, int dx, int dy, int dz, class Element>
struct Multiplication;
//-------------------------------------------------------------------------------------
template<class A, class B>
struct PolynomialExprForVertexExpr_Mult; 
//-------------------------------------------------------------------------------------
template<class A, class B>
inline PolynomialExprForVertexExpr_Mult <A,B>
operator * ( const PolynomialExprForVertexExpr <A>& a , 
             const PolynomialExprForVertexExpr <B>& b ); 
//-------------------------------------------------------------------------------------
template<class A, class D>
class PolynomialExprForVertexExpr_CMult; 
//-------------------------------------------------------------------------------------
template < class B > 
inline PolynomialExprForVertexExpr_CMult <B,float>
 operator * (float a, const PolynomialExprForVertexExpr <B> & b ); 
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CMult <A,float>
 operator * (const PolynomialExprForVertexExpr <A>& a , float b); 
//-------------------------------------------------------------------------------------
template < class B > 
inline PolynomialExprForVertexExpr_CMult <B,double>
 operator * (double a, const PolynomialExprForVertexExpr <B> & b ); 
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CMult <A,double>
 operator * (const PolynomialExprForVertexExpr <A>& a , double b); 
//-------------------------------------------------------------------------------------
template < class B > 
inline PolynomialExprForVertexExpr_CMult < B, std::complex<double> >
 operator * (const std::complex<double>& a, const PolynomialExprForVertexExpr <B>& b ); 
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CMult < A, std::complex<double> >
 operator * (const PolynomialExprForVertexExpr <A>& a, const std::complex<double>& b ); 
//-------------------------------------------------------------------------------------
template<class A, class B>
struct PolynomialExprForVertexExpr_Div; 
//-------------------------------------------------------------------------------------
template<class A, class B>
inline PolynomialExprForVertexExpr_Div <A,B>
 operator / ( const PolynomialExprForVertexExpr <A>& a , 
              const PolynomialExprForVertexExpr <B>& b ); 
//-------------------------------------------------------------------------------------
template<class A, class D>
class PolynomialExprForVertexExpr_CDiv; 
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CDiv <A,float>
 operator / ( const PolynomialExprForVertexExpr<A> & a , float b ); 
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CDiv <A,double>
 operator / ( const PolynomialExprForVertexExpr<A> & a , double b ); 
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CDiv < A, std::complex<double> >
 operator / (const PolynomialExprForVertexExpr <A>& a, const std::complex<double>& b ); 
//-----------------------------------------------------------------------------
template < class A >
struct EXP_F ;
//-----------------------------------------------------------------------------
template < class A > 
inline EXP_F <A> 
 Exp ( const PolynomialExprForVertexExpr <A> & a ); 
//-----------------------------------------------------------------------------
template < class A > 
struct VertexExpr {};
//-----------------------------------------------------------------------------
template < int vertexNumber >
struct Vertex_ ;
//-----------------------------------------------------------------------------
typedef Vertex_<0> P_0;
typedef Vertex_<1> P_1;
typedef Vertex_<2> P_2;
typedef Vertex_<3> P_3;
typedef Vertex_<4> P_4;
typedef Vertex_<5> P_5;
typedef Vertex_<6> P_6;
typedef Vertex_<7> P_7;
typedef Vertex_<8> P_8;
typedef Vertex_<9> P_9;
typedef Vertex_<10> P_10;
typedef Vertex_<11> P_11;
typedef Vertex_<12> P_12;
typedef Vertex_<13> P_13;
typedef Vertex_<14> P_14;
typedef Vertex_<15> P_15;
typedef Vertex_<16> P_16;
typedef Vertex_<17> P_17;
typedef Vertex_<18> P_18;
typedef Vertex_<19> P_19;
typedef Vertex_<20> P_20;
typedef Vertex_<21> P_21;
//-----------------------------------------------------------------------------
template <class A, class B>
struct VertexAdd ;
//-----------------------------------------------------------------------------
template <class A, class B>
struct VertexSub ;
//-----------------------------------------------------------------------------
template <typename Type, class B>
struct ConstVertexMult ;
//-----------------------------------------------------------------------------
/* This cross product (initialized by the operator *) exists only for elements
   with dimension 3. However, this is not caught yet. 
*/
template <class A, class B>
struct VertexCrossProduct ;
//-----------------------------------------------------------------------------
template <class A, class B>
struct VertexPolynomialProduct; 
//-----------------------------------------------------------------------------
template  < class A >
struct VertexMin ;
//-----------------------------------------------------------------------------
#if 1
template <class B>
inline ConstVertexMult <double,B>
 operator * (double a , const VertexExpr <B>& b){
     return  ConstVertexMult <double,B>(a);
  }
#endif
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexPolynomialProduct <A,B>
 operator * (const VertexExpr <A>& a , const PolynomialExprForVertexExpr <B>& b); 
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexPolynomialProduct <A,B>
 operator * ( const PolynomialExprForVertexExpr <B>& b , const VertexExpr <A>& a); 
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexAdd<A,B>
 operator + (const VertexExpr <A>& a , const VertexExpr <B>& b); 
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexSub <A,B>
 operator - (const VertexExpr <A>& a, const VertexExpr <B>& b); 
//-----------------------------------------------------------------------------
template < class A > 
inline VertexMin <A> 
 operator - (const VertexExpr <A>& a); 
//-----------------------------------------------------------------------------
template <class Element, class A, class B>
inline typename Element::bTYPE
 vec (const Element& element, const VertexExpr <A>& a, const VertexExpr <B>& b); 
//-----------------------------------------------------------------------------
template <class Element, class A>
inline typename Element::bTYPE
 length ( const Element& element, const VertexExpr <A>& a ); 
//-----------------------------------------------------------------------------
template <class Element, class A, class B>
inline typename Element::bTYPE
 det (const Element& element, const VertexExpr <A>& a, const VertexExpr <B>& b ); 
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexCrossProduct <A,B>
 operator * (const VertexExpr <A>& a , const VertexExpr <B>& b ); 
//-----------------------------------------------------------------------------
template < class A >
struct D_DX_Vertex;
//-----------------------------------------------------------------------------
template < class A >
struct D_DY_Vertex;
//-----------------------------------------------------------------------------
template < class A >
struct D_DZ_Vertex;
//-----------------------------------------------------------------------------
template < class A > 
inline D_DX_Vertex<A> 
 d_dx ( const VertexExpr<A>& a ); 
//-----------------------------------------------------------------------------
template < class A > 
inline D_DY_Vertex<A> 
 d_dy ( const VertexExpr<A>& a ); 
//-----------------------------------------------------------------------------
template < class A > 
inline D_DZ_Vertex<A> 
 d_dz ( const VertexExpr<A>& a ); 
//==============================================================================
