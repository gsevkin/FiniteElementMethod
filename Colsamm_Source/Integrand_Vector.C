//==============================================================================
//
//  $Id: Integrand_Vector.C,v 1.1 2006/07/19 09:34:25 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
/* A vector is a type list of its components, namely a vector of two  
 * components or a component and vector (which is a component itself)!
 * The evaluation of the operator '&' is from left to right (hopefully),
 * and thereby, we build a template type of the structure
 *   Integrand_Vector<Integrand_Vector<Integrand_Vector<....>,.>,.> 
 * This means, we can only access the Integrand_Vector from the end to the top! 
 * A vector (so far) is always a line or a column, not differentiated 
 * between them! 
*/
//-----------------------------------------------------------------------------
template <class A, class B>
void
Integrand_Vector<A,B>::print(std::ostream &os)const
   { 
#if 1
      a.print(os); 
      os << "\033[44m" << " & " << "\033[m" /*<< std::endl*/;  
      b.print(os);
#else
      os << "Integrand_Vector < " ;  
      a.print(os); 
      os << " , "; 
      b.print(os); 
      os << " > ";
#endif
   }
#if 1
//-----------------------------------------------------------------------------
template <class A, class B>
template<class Element, unsigned int iterLength>
typename Element::TYPE
Integrand_Vector<A,B>::evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
   {
      return b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template <class A1, class A2, class B>
template<class Element, unsigned int iterLength>
typename Element::TYPE
Integrand_Vector<Integrand_Vector<A1,A2>,B>::
  evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
   {
      return b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template <class A, class B>
inline Integrand_Vector<A,B> 
operator &(const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
   {
     return Integrand_Vector<A,B>(a,b);
   }
//-----------------------------------------------------------------------------
Integrand_Vector<Constant<int>,Constant<int> > 
inline operator &(const Constant<int>& a, const Constant<int>& b)
   {
     return Integrand_Vector<Constant<int>,Constant<int> >(a,b);
   }
//-----------------------------------------------------------------------------
template <class B>
Integrand_Vector<Constant<int>,B> 
inline operator &(const Constant<int>& a, const BasisFunctionExpr<B>& b)
   {
     return Integrand_Vector<Constant<int>,B>(a,b);
   }
#endif
//-----------------------------------------------------------------------------
#if 1
template <class A, class B, int l>
Integrand_Vector<A,B> 
inline operator &(const BasisFunctionVectorExpr<A,l>& a, const BasisFunctionExpr<B>& b)
   {
     return Integrand_Vector<A,B>(a,b);
   }
#endif
// ----------------------------------------------------------------------------
#if 1
 inline Integrand_Vector< NormalComponent<0>,NormalComponent<1> >
  N_2D() 
     { 
       return NormalComponent<0>() & NormalComponent<1>();
     }
// ----------------------------------------------------------------------------
 inline Integrand_Vector< UnitNormalComponent<0>,UnitNormalComponent<1> >
  N_Unit_2D() 
     { 
       return UnitNormalComponent<0> () & UnitNormalComponent<1>();
     }
#endif
// ----------------------------------------------------------------------------
#if 1
 inline Integrand_Vector<Integrand_Vector< NormalComponent<0>,NormalComponent<1> >, NormalComponent<2> >
  N() 
     { 
       return NormalComponent<0>() & NormalComponent<1>() & NormalComponent<2>();
     }
// ----------------------------------------------------------------------------
 inline Integrand_Vector<Integrand_Vector< UnitNormalComponent<0>,UnitNormalComponent<1> >, UnitNormalComponent<2> >
  N_Unit() 
     { 
       return UnitNormalComponent<0> () & UnitNormalComponent<1>() & UnitNormalComponent<2>();
     }
#endif
// ----------------------------------------------------------------------------
template<class A, class B, int l, class Element, unsigned int iterLength>
inline typename Element::TYPE 
 evaluate(const BasisFunctionVectorExpr<A,l>& a, const BasisFunctionVectorExpr<B,l>& b,
          const Element& element, const int (&iterator)[iterLength]) 
   {
     return static_cast<const A&>(a).evalLastComponent(element,iterator) * static_cast<const B&>(b).evalLastComponent(element,iterator) + 
            evaluate(static_cast<const A&>(a).getHeadOfVector(),static_cast<const B&>(b).getHeadOfVector(),element,iterator); 
   }
//-----------------------------------------------------------------------------
template<class A, class B, class Element, unsigned int iterLength>
inline typename Element::TYPE 
 evaluate(const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b,
          const Element& element, const int (&iterator)[iterLength]) 
   {
     return static_cast<const A&>(a).eval(element,iterator) * static_cast<const B&>(b).eval(element,iterator);
   }
//-----------------------------------------------------------------------------

/* Concluding comments
 * At this point there is still more to introduce. A compile time data structure
 * which holds tensor data. Therefore we would need an additional operator to be 
 * overloaded, whicht has less priority than operator '&'. 
*/

//==============================================================================

