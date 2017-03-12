//==============================================================================
//
//  $Id: Stencil.C,v 1.16 2006/07/25 08:11:31 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
template <class A> 
template<class Element, unsigned int iterLength> 
const typename Element::TYPE
BasisFunctionExpr<A>::eval(const Element& element,const int (&iterator)[iterLength]) 
  {
    return A::eval(element,iterator);
  }
//-----------------------------------------------------------------------------
template <int uniqueFloatId>
Float_<uniqueFloatId>::Float_(float constantValue_) 
   {
     constantValue = constantValue_;
   }
//-----------------------------------------------------------------------------
template <int uniqueFloatId>
void
Float_<uniqueFloatId>::print(std::ostream & os)
   {
     os << " Float_<" << uniqueFloatId << ">(" << constantValue <<")";
   }
 //-----------------------------------------------------------------------------
template <int uniqueFloatId>
Float_<uniqueFloatId> 
Float_<uniqueFloatId>::operator = (const float& constantValue_) 
    { 
      constantValue = constantValue_; 
      return (*this);
    }
//-----------------------------------------------------------------------------
template <int uniqueFloatId>
Float_<uniqueFloatId> 
Float_<uniqueFloatId>::operator ()(float constantValue_) 
    { 
      constantValue = constantValue_;  
      return (*this);
    }
//-----------------------------------------------------------------------------
template <int uniqueFloatId>
template<class Element, unsigned int iterLength> 
const typename Element::TYPE
Float_<uniqueFloatId>::eval(const Element& element,const int (&iterator)[iterLength]) 
  {
    return constantValue;
  }
//-----------------------------------------------------------------------------
template <int uniqueFloatId> 
float 
 Float_<uniqueFloatId>::constantValue;
//-----------------------------------------------------------------------------
template <int uniqueDoubleId>
Double_<uniqueDoubleId>::Double_(double constantValue_) 
   {
     constantValue = constantValue_;
   }
//-----------------------------------------------------------------------------
template <int uniqueDoubleId>
void
Double_<uniqueDoubleId>::print(std::ostream & os)
   {
     os << " Double_<" << uniqueDoubleId << ">(" << constantValue <<")";
   }
 //-----------------------------------------------------------------------------
template <int uniqueDoubleId>
Double_<uniqueDoubleId> 
Double_<uniqueDoubleId>::operator = (const double& constantValue_) 
    { 
      constantValue = constantValue_; 
      return (*this);
    }
//-----------------------------------------------------------------------------
template <int uniqueDoubleId>
Double_<uniqueDoubleId> 
Double_<uniqueDoubleId>::operator ()(double constantValue_) 
    { 
      constantValue = constantValue_;  
      return (*this);
    }
//-----------------------------------------------------------------------------
template <int uniqueDoubleId>
template<class Element, unsigned int iterLength> 
const typename Element::TYPE
Double_<uniqueDoubleId>::eval(const Element& element,const int (&iterator)[iterLength]) 
  {
    return constantValue;
  }
//-----------------------------------------------------------------------------
template <int uniqueDoubleId> 
double 
 Double_<uniqueDoubleId>::constantValue;
//-----------------------------------------------------------------------------
template <int uniqueComplexId>
void
Complex_<uniqueComplexId>::print(std::ostream & os)
   {
#if 1
     os << " Complex_<" << uniqueComplexId << ">(" << constantValue <<")";
#else
     os <<  " " << constantValue;
#endif
   }
//-----------------------------------------------------------------------------
template <int uniqueComplexId>
Complex_<uniqueComplexId>::Complex_(const std::complex<double>& constantValue_) 
   {
     constantValue = constantValue_;
   } 
//-----------------------------------------------------------------------------
template <int uniqueComplexId>
Complex_<uniqueComplexId> 
Complex_<uniqueComplexId>::operator = (const std::complex<double>& constantValue_) 
   { 
     constantValue = constantValue_; 
     return (*this);
   }
//-----------------------------------------------------------------------------
template <int uniqueComplexId>
Complex_<uniqueComplexId> 
Complex_<uniqueComplexId>::operator ()(const std::complex<double>& constantValue_) 
   { 
     constantValue = constantValue_;  
     return (*this);
   }
//-----------------------------------------------------------------------------
template <int uniqueComplexId>
template<class Element, unsigned int iterLength> 
const typename Element::TYPE&
Complex_<uniqueComplexId>::eval(const Element& element,const int (&iterator)[iterLength]) 
   {
     return constantValue;
   }
//-----------------------------------------------------------------------------
template <int uniqueComplexId> 
std::complex<double> 
 Complex_<uniqueComplexId>::constantValue;
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, int functionId>
void
BasisFunction<number,functionId>::print(std::ostream &os) 
   {
     if (number == BasisSet_1 && functionId == 0) 
        os << " v_";
     if (number == BasisSet_1 && functionId == 1) 
        os << " w_";
     if (number == BasisSet_2 && functionId == 0) 
        os << " p_";
   }
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, int functionId>
template<class Element, unsigned int iterLength> 
const typename Element::TYPE
BasisFunction<number,functionId>::eval(const Element& element,const int (&iterator)[iterLength]) 
   {
     return element.getBasisFunctionsAtGaussianPoints_(number,iterator[iterLength-1],iterator[functionId])[fn];
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T>
void
FUNC1_<uniqueFunctionId,T>::print(std::ostream &os) 
   {
      os << "FUNC1_<" << uniqueFunctionId << ">";
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
FUNC1_<uniqueFunctionId,T>::FUNC1_(const typename FunctionInterfaces<T>::functionInterfaceOne& f) 
   {
     func_ = f;
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
FUNC1_<uniqueFunctionId,T>::eval(const Element& element, const int (&iterator)[iterLength]) 
   { 
     typename Element::TYPE ** map = element.getMappingValues(iterator[iterLength-1]);
     return func_(map[dirX][fn]);
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
typename FunctionInterfaces<T>::functionInterfaceOne 
 FUNC1_<uniqueFunctionId,T>::func_;
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T>
void
FUNC2_<uniqueFunctionId,T>::print(std::ostream &os) 
   { 
     os << "FUNC2_<" << uniqueFunctionId << ">";
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
FUNC2_<uniqueFunctionId,T>::FUNC2_(const typename FunctionInterfaces<T>::functionInterfaceTwo& f) 
   {
     func_ = f;
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
FUNC2_<uniqueFunctionId,T>::eval(const Element& element, const int (&iterator)[iterLength]) 
   {     
     typename Element::TYPE ** map = element.getMappingValues(iterator[iterLength-1]);
     return func_(map[dirX][fn],map[dirY][fn]);
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
typename FunctionInterfaces<T>::functionInterfaceTwo
 FUNC2_<uniqueFunctionId,T>::func_;
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
void
FUNC3_<uniqueFunctionId,T>::print(std::ostream &os) 
   { 
     os << "FUNC3_<" << uniqueFunctionId << ">";
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
FUNC3_<uniqueFunctionId,T>::FUNC3_(const typename FunctionInterfaces<T>::functionInterfaceThree& f) 
   {
      func_ = f;
    }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
FUNC3_<uniqueFunctionId,T>::eval(const Element& element, const int (&iterator)[iterLength]) 
   { 
     typename Element::TYPE ** map = element.getMappingValues(iterator[iterLength-1]);
     return func_( map[dirX][fn], 
                   map[dirY][fn],
                   map[dirZ][fn]);
   }
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T> 
typename FunctionInterfaces<T>::functionInterfaceThree 
 FUNC3_<uniqueFunctionId,T>::func_;
//-----------------------------------------------------------------------------
template <int ex, SpaceDirection direction>
void
BasisMonom<ex,direction>::print(std::ostream &os) 
    { 
      os <<  " " << static_cast<char>('x' + direction) << "^(" << ex << ")";
    }
//-----------------------------------------------------------------------------
template <int ex, SpaceDirection direction>
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
BasisMonom<ex,direction>::eval(const Element& element,const int (&iterator)[iterLength]) 
   { 
     typename Element::TYPE ** map = element.getMappingValues(iterator[iterLength-1]);
     return Colsamm_Internal_Functions::POW<ex>::pow(map[direction][fn]);
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
NormalComponent<position>::print(std::ostream &os) 
   { 
     os << "N_" << position ; 
   }
//-----------------------------------------------------------------------------
template <int position>
template<class Element, unsigned int iterLength>
typename Element::TYPE 
NormalComponent<position>::eval(const Element& element,const int (&iterator)[iterLength]) 
   { 
     enum {pos = position * Element::dimensionOfElement + Element::dimensionOfElement-1 };
     return - element.getValuesOfJacobianMatrix(iterator[iterLength-1])[pos];
   }
//-----------------------------------------------------------------------------
template <int position>
void 
UnitNormalComponent<position>::print(std::ostream &os) 
    { 
      os << "N_" << position ; 
    }
//-----------------------------------------------------------------------------
template <int position>
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
UnitNormalComponent<position>::eval(const Element& element,const int (&iterator)[iterLength]) 
   { 
      enum {pos = position * Element::dimensionOfElement + Element::dimensionOfElement-1};
      return -element.getValuesOfJacobianMatrix(iterator[iterLength-1])[pos]/
              element.getNormalLength(iterator[iterLength-1]);
   }
//-----------------------------------------------------------------------------
template <int uniqueVertexVectorId, class T>
void
VertexVector<uniqueVertexVectorId,T>::print(std::ostream &os) 
   { 
     os << "BFI_<" << uniqueVertexVectorId << ">";
   }
//-----------------------------------------------------------------------------
template <int uniqueVertexVectorId, class T>
VertexVector<uniqueVertexVectorId,T>::VertexVector(const T& d) 
   {
     data=d;
   }
//-----------------------------------------------------------------------------
template <int uniqueVertexVectorId, class T>
template<class Element, unsigned int iterLength> 
typename Element::TYPE 
VertexVector<uniqueVertexVectorId,T>::eval(const Element& element,const int (&iterator)[iterLength]) 
    { 
      return data[iterator[0]];
    }
//-----------------------------------------------------------------------------
template <int uniqueVertexVectorId, class T> 
T VertexVector<uniqueVertexVectorId,T>::data;
//-----------------------------------------------------------------------------
template <class A> 
void
Conjugate<A>::print(std::ostream &os) 
   {
     os << "conj("; A::print(os); os << " )";
   }
//-----------------------------------------------------------------------------
template <class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Conjugate<A>::eval(const Element& element,const int (&iterator)[iterLength])
   {
     return conj(A::eval(element,iterator));
   }
//-----------------------------------------------------------------------------
template<class A> 
void
Min_<A>::print(std::ostream &os) 
    {
      os << " -("; A::print(os); os << ")";
    }
//-----------------------------------------------------------------------------
template<class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE
Min_<A>::eval(const Element& element,const int (&iterator)[iterLength]) 
    {
      return -A::eval(element,iterator);
    }
//-----------------------------------------------------------------------------
template<class A, class B> 
void
Add_<A,B>::print(std::ostream & os) 
   { 
     A::print(os); 
     os << "\033[46m" << " + " << "\033[m" ; 
     B::print(os);
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
template<class Element, unsigned int iterLength>
typename Element::TYPE
Add_<A,B>::eval(const Element& element,const int (&iterator)[iterLength])
   {
     return A::eval(element,iterator) + B::eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class A, class B>
void
Sub_<A,B>::print(std::ostream& os) 
   {
     A::print(os);
     os << "\033[45m" << " - " << "\033[m" ;
     B::print(os);
   }
//-----------------------------------------------------------------------------
template<class A, class B>
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Sub_<A,B>::eval(const Element& element, const int (&iterator)[iterLength]) 
   {
     return A::eval(element,iterator) - B::eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
void
Mult_<A,B>::print(std::ostream &os) 
   {
     os << " ("; 
     A::print(os); 
     os << " )"<< "\033[43m" << " * " << "\033[m" << "("; 
     B::print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Mult_<A,B>::eval(const Element& element, const int (&iterator)[iterLength]) 
   {
     return A::eval(element,iterator)*B::eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
void
Div_<A,B>::print(std::ostream &os) 
   {
     os << " ("; 
     A::print(os); 
     os << " )"<< "\033[43m" << " / " << "\033[m" << "("; 
     B::print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class A, class B> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
Div_<A,B>::eval(const Element& element,const int (&iterator)[iterLength]) 
   {
     return A::eval(element,iterator)/B::eval(element,iterator);
   }
//-----------------------------------------------------------------------------
template<class A> 
void
SQRT_<A>::print(std::ostream &os) 
   {
     os << " sqrt("; 
     A::print(os); 
     os << " )";
   }
//-----------------------------------------------------------------------------
template<class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
SQRT_<A>::eval(const Element& element,const int (&iterator)[iterLength]) 
   {
     return sqrt ( A::eval(element,iterator) ) ;
   }
//-----------------------------------------------------------------------------
template<int number, class A, typename eTYPE> 
POW_<number,A,eTYPE>::POW_ (eTYPE e)
   {
     exponent = e;
   }
//-----------------------------------------------------------------------------
template<int number, class A, typename eTYPE> 
void 
POW_<number,A,eTYPE>::print(std::ostream &os) 
   {
     A::print(os); 
     os << " ^" << exponent << " ";
   }
//-----------------------------------------------------------------------------
template<int number, class A, typename eTYPE> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
POW_<number,A,eTYPE>::eval(const Element& element, const int (&iterator)[iterLength]) 
   {
     return std::pow( A::eval(element,iterator), exponent ) ;
   }
//-----------------------------------------------------------------------------
template<class A> 
void 
EXP_<A>::print(std::ostream &os) 
   {
     os << "exp(";
     A::print(os); 
     os << ") ";
   }
//-----------------------------------------------------------------------------
template<class A> 
template<class Element, unsigned int iterLength>
typename Element::TYPE 
EXP_<A>::eval(const Element& element, const int (&iterator)[iterLength]) 
   {
     return std::exp(A::eval(element,iterator)) ;
   }
//-----------------------------------------------------------------------------
template<int p, class A, typename eTYPE> 
eTYPE POW_<p,A,eTYPE>::exponent = 0.;
//-----------------------------------------------------------------------------
template <class A, what typ> 
void
Deri_<A,typ>::print(std::ostream &os) 
   {
     if (typ == dx)
         os << " d_x(";
     if (typ == dy)
         os << " d_y(";
     if (typ == dz)
         os << " d_z(";
     A::print(os); os << " )";
   }
//-----------------------------------------------------------------------------
template <class A, what typ> 
template<class Element, unsigned int iterLength>
const typename Element::TYPE
Deri_<A,typ>::eval(const Element& element, const int (&iterator)[iterLength]) 
   {
     return A::template evaluateDerivative<typ>(element,iterator);
   }
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, int functionId>
template<what typ, class Element, unsigned int iterLength> 
const typename Element::TYPE
BasisFunction<number,functionId>::evaluateDerivative(const Element& element,const int (&iterator)[iterLength]) 
   {
      if ( (unsigned int) typ <= (unsigned int) Element::dimensionOfElement)
        {
          return element.template getElementMappingDerivativeValues<number,functionId,typ>(iterator);
        }
     else 
        { 
          return 0;
        }
   }
//==============================================================================

