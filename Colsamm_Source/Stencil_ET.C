//==============================================================================
//
//  $Id: Stencil_ET.C,v 1.10 2006/07/25 13:52:27 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
template <class A> 
template<class Element, unsigned int iterLength> 
const typename Element::TYPE
BasisFunctionExpr<A>::eval(const Element& element,const int (&iterator)[iterLength]) const
  {
    return static_cast<const A&>(*this).eval(element,iterator);
  }
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId>
void
BasisFunction<number,functionId>::print(std::ostream &os)  const
   {
     if (number == BasisSet_1 && functionId == 0) 
        os << " v_";
     if (number == BasisSet_1 && functionId == 1) 
        os << " w_";
     if (number == BasisSet_2 && functionId == 0) 
        os << " p_";
   }
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId>
template<class Element, unsigned int iterLength> 
const typename Element::TYPE
BasisFunction<number,functionId>::eval(const Element& element,const int (&iterator)[iterLength]) const
   {
     return element.template getBasisFunctionsAtGaussianPointsFn<number,functionId>(iterator);
   }
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId>
template<what typ, class Element, unsigned int iterLength> 
const typename Element::TYPE
BasisFunction<number,functionId>::evaluateDerivative(const Element& element,const int (&iterator)[iterLength]) const
   {
     if ((unsigned int)typ <= (unsigned int)Element::dimensionOfElement)
        {
          assert ( functionId < iterLength -1 );
          return element.template getElementMappingDerivativeValues<number,functionId,typ>(iterator);
        }
     else 
        { 
          return 0;
        }
   }
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId, SpaceDirection direction>
void
Vector_BasisFunction<number,functionId,direction>::print(std::ostream &os)  const
   {
     if (number == BasisSet_1 && functionId == 0) 
        os << " Vector_v_" << direction << "_";
     if (number == BasisSet_1 && functionId == 1) 
        os << " Vector_w_" << direction << "_";
     if (number == BasisSet_2 && functionId == 0) 
        os << " Vector_p_" << direction << "_";
   }
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId, SpaceDirection direction>
template<class Element, unsigned int iterLength> 
const typename Element::TYPE
Vector_BasisFunction<number,functionId,direction>
  ::eval(const Element& element,const int (&iterator)[iterLength]) const
   {
     return element.template getBasisFunctionsAtGaussianPointsFn_Vector<number,direction,functionId>(iterator);
   }
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId, SpaceDirection direction>
template<what type, class Element, unsigned int iterLength> 
const typename Element::TYPE
Vector_BasisFunction<number,functionId,direction>
  ::evaluateDerivative(const Element& element,const int (&iterator)[iterLength]) const
   {
     if ((unsigned int)type <= (unsigned int)Element::dimensionOfElement)
        {
          return element.template getElementMappingDerivativeValues_Vector<number,direction,functionId,type>(iterator);
        }
     else 
        { 
          return 0;
        }
   }
//-----------------------------------------------------------------------------
template <class T>
void
FUNC1_<T>::print(std::ostream &os)  const
   {
      os << "FUNC1_";
   }
//-----------------------------------------------------------------------------
template <class T> 
FUNC1_<T>::FUNC1_(const typename FunctionInterfaces<T>::functionInterfaceOne& f) 
   {
     func_ = f;
   }
//-----------------------------------------------------------------------------
template <class T> 
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
FUNC1_<T>::eval(const Element& element, const int (&iterator)[iterLength]) const
   { 
     typename Element::TYPE ** map = element.getMappingValues(iterator[iterLength-1]);
     return func_(map[dirX][fn]);
   }
//-----------------------------------------------------------------------------
template <class T>
void
FUNC2_<T>::print(std::ostream &os)  const
   { 
     os << "FUNC2_";
   }
//-----------------------------------------------------------------------------
template <class T> 
FUNC2_<T>::FUNC2_(const typename FunctionInterfaces<T>::functionInterfaceTwo& f) 
   {
     func_ = f;
   }
//-----------------------------------------------------------------------------
template <class T> 
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
FUNC2_<T>::eval(const Element& element, const int (&iterator)[iterLength]) const
   {     
     typename Element::TYPE ** map = element.getMappingValues(iterator[iterLength-1]);
     return func_(map[dirX][fn],map[dirY][fn]);
   }
//-----------------------------------------------------------------------------
template <class T> 
void
FUNC3_<T>::print(std::ostream &os)  const
   { 
     os << "FUNC3_";
   }
//-----------------------------------------------------------------------------
template <class T> 
FUNC3_<T>::FUNC3_(const typename FunctionInterfaces<T>::functionInterfaceThree& f) 
   {
      func_ = f;
    }
//-----------------------------------------------------------------------------
template <class T> 
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
FUNC3_<T>::eval(const Element& element, const int (&iterator)[iterLength]) const
   { 
     typename Element::TYPE ** map = element.getMappingValues(iterator[iterLength-1]);
     return func_( map[dirX][fn], 
                   map[dirY][fn],
                   map[dirZ][fn]);
   }
//-----------------------------------------------------------------------------
template <int ex, SpaceDirection direction>
void
BasisMonom<ex,direction>::print(std::ostream &os)  const
    { 
      os <<  " " << static_cast<char>('x' + direction) << "^(" << ex << ")";
    }
//-----------------------------------------------------------------------------
template <int ex, SpaceDirection direction>
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
BasisMonom<ex,direction>::eval(const Element& element,const int (&iterator)[iterLength]) const
    { 
      return Colsamm_Internal_Functions::POW<ex>::pow(element.getMappingValues(iterator[iterLength-1])[direction][fn]);
    }
//-----------------------------------------------------------------------------
template <int p>
BasisMonom<p,dirX> 
Interfaces::x_() 
   {
     return BasisMonom<p,dirX>();
   }
//-----------------------------------------------------------------------------
template <int p>
BasisMonom<p,dirY> 
Interfaces::y_() 
   {
     return BasisMonom<p,dirY>();
   }
//-----------------------------------------------------------------------------
template <int p>
BasisMonom<p,dirZ> 
Interfaces::z_() 
   {
     return BasisMonom<p,dirZ>();
   }
//-----------------------------------------------------------------------------
template <int position>
void
NormalComponent<position>::print(std::ostream &os) const
   { 
     os << "N_" << position ; 
   }
//-----------------------------------------------------------------------------
template <int position>
template<class Element, unsigned int iterLength>
typename Element::TYPE 
NormalComponent<position>::eval(const Element& element,const int (&iterator)[iterLength]) const
   { 
     enum {pos = position * Element::dimensionOfElement + Element::dimensionOfElement-1 };
     return -element.getValuesOfJacobianMatrix(iterator[iterLength-1])[pos]; 
   }
//-----------------------------------------------------------------------------
template <int position>
void 
UnitNormalComponent<position>::print(std::ostream &os) const
    { 
      os << "N_" << position ; 
    }
//-----------------------------------------------------------------------------
template <int position>
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
UnitNormalComponent<position>::eval(const Element& element,const int (&iterator)[iterLength]) const
   { 
     if ( position < Element::dimensionOfElement )
        {
          enum {pos = position * Element::dimensionOfElement + Element::dimensionOfElement-1};
          assert ( 1.e-10 < fabs(element.getNormalLength(iterator[iterLength-1])) );
          return -element.getValuesOfJacobianMatrix(iterator[iterLength-1])[pos]/
                  element.getNormalLength(iterator[iterLength-1]);
        }
     else 
        {
          return 0.;
        }
   }
//-----------------------------------------------------------------------------
template <class T>
void
VertexVector<T>::print(std::ostream &os) const
   { 
     os << "BFI_";
   }
//-----------------------------------------------------------------------------
template <class T>
VertexVector<T>::VertexVector(const T& d) : data(d)
   {}
//-----------------------------------------------------------------------------
template <class T>
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
VertexVector<T>::eval(const Element& element,const int (&iterator)[iterLength]) const
    { 
      return data[iterator[0]];
    }
//-----------------------------------------------------------------------------
template <class A> 
void
Conjugate<A>::print(std::ostream &os) const
   {
     os << "conj("; a.print(os); os << " )";
   }
//-----------------------------------------------------------------------------
template <class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Conjugate<A>::eval(const Element& element,const int (&iterator)[iterLength])const
   {
     return conj(a.eval(element,iterator));
   }
//-----------------------------------------------------------------------------
template<class A> 
void
Min_<A>::print(std::ostream &os) const
    {
      os << " -("; a.print(os); os << ")";
    }
//-----------------------------------------------------------------------------
template<class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE
Min_<A>::eval(const Element& element,const int (&iterator)[iterLength]) const
    {
      return -a.eval(element,iterator);
    }
//-----------------------------------------------------------------------------
template<class A, class B> 
void
Add_<A,B>::print(std::ostream & os) const
   { 
     a.print(os); 
     os << "\033[46m" << " + " << "\033[m" ; 
     b.print(os);
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
template<class Element, unsigned int iterLength>
typename Element::TYPE
Add_<A,B>::eval(const Element& element,const int (&iterator)[iterLength])const
   {
     return a.eval(element,iterator) + b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class TYPE, class B> 
void
CAdd_<TYPE,B>::print(std::ostream & os) const
   { 
     os << "\033[46m" << a << " + " ; 
     b.print(os); 
     os << "\033[m" ; 
   }
//-----------------------------------------------------------------------------
template<class TYPE, class B> 
template<class Element, unsigned int iterLength>
typename Element::TYPE
CAdd_<TYPE,B>::eval(const Element& element,const int (&iterator)[iterLength])const
   {
     return a + b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class A, class B>
void
Sub_<A,B>::print(std::ostream& os) const
   {
     a.print(os);
     os << "\033[45m" << " - " << "\033[m" ;
     b.print(os);
   }
//-----------------------------------------------------------------------------
template<class A, class B>
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Sub_<A,B>::eval(const Element& element, const int (&iterator)[iterLength]) const
   {
     return a.eval(element,iterator) - b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class TYPE, class B>
void
CSub_<TYPE,B>::print(std::ostream& os) const
   {
     os << "\033[45m" << a << " - " << "\033[m" ;
     b.print(os);
   }
//-----------------------------------------------------------------------------
template<class TYPE, class B>
template<class Element, unsigned int iterLength>
typename Element::TYPE 
CSub_<TYPE,B>::eval(const Element& element, const int (&iterator)[iterLength]) const
   {
     return a - b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
void
Mult_<A,B>::print(std::ostream &os) const
   {
     os << " ("; 
     a.print(os); 
     os << " )"<< "\033[43m" << " * " << "\033[m" << "("; 
     b.print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Mult_<A,B>::eval(const Element& element, const int (&iterator)[iterLength]) const
   {
     return a.eval(element,iterator) * b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class TYPE, class B> 
void
CMult_<TYPE,B>::print(std::ostream &os) const
   {
     os << " (" << a << " )"<< "\033[43m" << " * " << "\033[m" << "("; 
     b.print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class TYPE, class B> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
CMult_<TYPE,B>::eval(const Element& element, const int (&iterator)[iterLength]) const
   {
     return a*b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
void
Div_<A,B>::print(std::ostream &os) const
   {
     os << " ("; 
     a.print(os); 
     os << " )"<< "\033[43m" << " / " << "\033[m" << "("; 
     b.print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Div_<A,B>::eval(const Element& element,const int (&iterator)[iterLength]) const
   {
     return a.eval(element,iterator)/b.eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class A> 
void
SQRT_<A>::print(std::ostream &os) const
   {
     os << " sqrt("; 
     a.print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
SQRT_<A>::eval(const Element& element,const int (&iterator)[iterLength]) const
   {
     return sqrt ( a.eval(element,iterator) ) ;
   }
//-----------------------------------------------------------------------------
template<class A> 
void
Sinus_<A>::print(std::ostream &os) const
   {
     os << " sin("; 
     a.print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Sinus_<A>::eval(const Element& element,const int (&iterator)[iterLength]) const
   {
     return std::sin( a.eval(element,iterator) ) ;
   }
//-----------------------------------------------------------------------------
template<class A> 
void
Cosinus_<A>::print(std::ostream &os) const
   {
     os << " cos("; 
     a.print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Cosinus_<A>::eval(const Element& element,const int (&iterator)[iterLength]) const
   {
     return std::cos( a.eval(element,iterator) ) ;
   }
//-----------------------------------------------------------------------------
template<class A, typename eTYPE> 
POW_<A,eTYPE>::POW_ (eTYPE e)
   {
     exponent = e;
   }
//-----------------------------------------------------------------------------
template<class A, typename eTYPE> 
void 
POW_<A,eTYPE>::print(std::ostream &os) const
   {
     a.print(os); 
     os << " ^" << exponent << " ";
   }
//-----------------------------------------------------------------------------
template<class A, typename eTYPE> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
POW_<A,eTYPE>::eval(const Element& element, const int (&iterator)[iterLength]) const
   {
     return std::pow( a.eval(element,iterator), exponent ) ;
   }
//-----------------------------------------------------------------------------
template<class A> 
void 
EXP_<A>::print(std::ostream &os) const
   {
     os << "exp(";
     a.print(os); 
     os << ") ";
   }
//-----------------------------------------------------------------------------
template<class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
EXP_<A>::eval(const Element& element, const int (&iterator)[iterLength]) const
   {
     return std::exp(a.eval(element,iterator)) ;
   }
//-----------------------------------------------------------------------------
template <class A, what typ> 
void
Deri_<A,typ>::print(std::ostream &os) const
   {
     if (typ == dx)
         os << " d_x(";
     if (typ == dy)
         os << " d_y(";
     if (typ == dz)
         os << " d_z(";
     a.print(os); os << " )";
   }
//-----------------------------------------------------------------------------
template <class A, what typ> 
template<class Element, unsigned int iterLength>
const typename Element::TYPE
Deri_<A,typ>::eval(const Element& element, const int (&iterator)[iterLength]) const
   {
     return a.template evaluateDerivative<typ>(element,iterator);
   }
//==============================================================================

