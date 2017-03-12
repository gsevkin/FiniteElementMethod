//==============================================================================
//
//  $Id: CT_Manip_ET.C,v 1.10 2006/07/19 09:05:18 jochen Exp $
//
//==============================================================================
template <class A>
void 
 print(const A& a) 
    {
      std::cerr << std::endl; 
      a.print(std::cerr); 
      std::cerr << std::endl;
    }
//-----------------------------------------------------------------------------
// Constructor functions to build the corresponding expression objects!
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template <typename TYPE>
inline Constant<TYPE> C_(const TYPE& t){
      return Constant<TYPE>(t);
  }
//-----------------------------------------------------------------------------
inline Constant<double> C_(const double& t){
      return Constant<double>(t);
  }
//-----------------------------------------------------------------------------
  template <class T>
  inline VertexVector<T> 
   vec(const T(&d))
     {
       return VertexVector<T>(d);
     }
//-----------------------------------------------------------------------------
  template <class T>
  inline FUNC1_<T> 
   func(const typename FunctionInterfaces<T>::functionInterfaceOne& f_)
     {
       return FUNC1_<T>(f_);
     }
//-----------------------------------------------------------------------------
  template <class T>
  inline FUNC2_<T> 
   func(const typename FunctionInterfaces<T>::functionInterfaceTwo& f_)
     { 
       return FUNC2_<T>(f_);
     }
//-----------------------------------------------------------------------------
 template <class T>
 inline FUNC3_<T> 
  func(const typename FunctionInterfaces<T>::functionInterfaceThree& f_)
    {
      return FUNC3_<T>(f_);
    }
//-----------------------------------------------------------------------------
  inline FUNC1_<double> 
   func(const FunctionInterfaces<double>::functionInterfaceOne& f_)
     {
       return FUNC1_<double>(f_);
     }
//-----------------------------------------------------------------------------
  inline FUNC2_<double> 
   func(const FunctionInterfaces<double>::functionInterfaceTwo& f_)
     { 
       return FUNC2_<double>(f_);
     }
//-----------------------------------------------------------------------------
  inline FUNC3_<double> 
  func(const FunctionInterfaces<double>::functionInterfaceThree& f_)
    {
      return FUNC3_<double>(f_);
    }
//-----------------------------------------------------------------------------
  inline FUNC1_<std::complex<double> > 
   func(const FunctionInterfaces<std::complex<double> >::functionInterfaceOne& f_)
     {
       return FUNC1_<std::complex<double> >(f_);
     }
//-----------------------------------------------------------------------------
  inline FUNC2_<std::complex<double> > 
   func(const FunctionInterfaces<std::complex<double> >::functionInterfaceTwo& f_)
     { 
       return FUNC2_<std::complex<double> >(f_);
     }
//-----------------------------------------------------------------------------
  inline FUNC3_<std::complex<double> > 
  func(const FunctionInterfaces<std::complex<double> >::functionInterfaceThree& f_)
    {
      return FUNC3_<std::complex<double> >(f_);
    }
//-----------------------------------------------------------------------------
 template <class A> 
 inline Conjugate<A> 
  conj(const BasisFunctionExpr<A>& a) 
     {
       return Conjugate<A>(a);
     }
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
 template <class B> 
 inline CAdd_<int,B>
 operator + (const int&a, const BasisFunctionExpr<B>& b)
    {
      return CAdd_<int,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
 inline CAdd_<float,B>
 operator + (const float&a, const BasisFunctionExpr<B>& b)
    {
      return CAdd_<float,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
 inline CAdd_<double,B>
 operator + (const double&a, const BasisFunctionExpr<B>& b)
    {
      return CAdd_<double,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <typename TYPE, class B> 
 inline CAdd_<typename std::complex<TYPE>,B>
 operator + (const typename std::complex<TYPE>&a, const BasisFunctionExpr<B>& b)
    {
      return CAdd_<typename std::complex<TYPE>,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
 inline CAdd_<int,B>
 operator + (const BasisFunctionExpr<B>& b, const int&a)
    {
      return CAdd_<int,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
 inline CAdd_<float,B>
 operator + (const BasisFunctionExpr<B>& b, const float&a)
    {
      return CAdd_<float,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
 inline CAdd_<double,B>
 operator + (const BasisFunctionExpr<B>& b, const double&a)
    {
      return CAdd_<double,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class TYPE, class B> 
 inline CAdd_<typename std::complex<TYPE>,B>
 operator + (const BasisFunctionExpr<B>& b, const typename std::complex<TYPE>&a)
    {
      return CAdd_<typename std::complex<TYPE>,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class A, class B> inline Add_<A,B>
 operator + (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
    {
      return Add_<A,B>(a,b);
    }
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
 template <class B> inline CSub_<int,B>
 operator - (const int& a, const BasisFunctionExpr<B>& b)
    {
      return CSub_<int,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> inline CSub_<float,B>
 operator - (const float& a, const BasisFunctionExpr<B>& b)
    {
      return CSub_<float,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> inline CSub_<double,B>
 operator - (const double& a, const BasisFunctionExpr<B>& b)
    {
      return CSub_<double,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class TYPE, class B> inline CSub_<typename std::complex<TYPE>,B>
 operator - (const typename std::complex<TYPE>& a, const BasisFunctionExpr<B>& b)
    {
      return CSub_<typename std::complex<TYPE>,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> inline CAdd_<int,B>
 operator - (const BasisFunctionExpr<B>& b, const int& a)
    {
      return CAdd_<int,B>(-a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> inline CAdd_<float,B>
 operator - (const BasisFunctionExpr<B>& b, const float& a)
    {
      return CAdd_<float,B>(-a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> inline CAdd_<double,B>
 operator - (const BasisFunctionExpr<B>& b, const double& a)
    {
      return CAdd_<double,B>(-a,b);
    }
//-----------------------------------------------------------------------------
 template <class TYPE, class B> inline CAdd_<typename std::complex<TYPE>,B>
 operator - (const BasisFunctionExpr<B>& b, const typename std::complex<TYPE>& a)
    {
      return CAdd_<typename std::complex<TYPE>,B>(-a,b);
    }
//-----------------------------------------------------------------------------
 template <class A, class B> inline Sub_<A,B>
 operator - (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
    {
      return Sub_<A,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class A> inline Min_<A>
 operator - (const BasisFunctionExpr<A>& a)
    {
      return Min_<A>(a);
    }
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
 template <class B> 
  inline CMult_<int,B>
 operator * (const int& a, const BasisFunctionExpr<B>& b)
    { 
      return CMult_<int,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
  inline CMult_<float,B>
 operator * (const float& a, const BasisFunctionExpr<B>& b)
    { 
      return CMult_<float,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
  inline CMult_<double,B>
 operator * (const double& a, const BasisFunctionExpr<B>& b)
    { 
      return CMult_<double,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class TYPE, class B> 
  inline CMult_<typename std::complex<TYPE>,B>
 operator * (const typename std::complex<TYPE>& a, const BasisFunctionExpr<B>& b)
    { 
      return CMult_<typename std::complex<TYPE>,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
  inline CMult_<int,B>
 operator * (const BasisFunctionExpr<B>& b, const int& a)
    { 
      return CMult_<int,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
  inline CMult_<float,B>
 operator * (const BasisFunctionExpr<B>& b, const float& a)
    { 
      return CMult_<float,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class B> 
  inline CMult_<double,B>
 operator * (const BasisFunctionExpr<B>& b, const double& a)
    { 
      return CMult_<double,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <typename TYPE, class B> 
  inline CMult_<typename std::complex<TYPE>,B>
 operator * (const BasisFunctionExpr<B>& b, const typename std::complex<TYPE>& a)
    { 
      return CMult_<typename std::complex<TYPE>,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class A, class B>
  inline Mult_<A,B>
 operator * (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
    { 
      return Mult_<A,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class A, class B> inline Div_<A,B>
 operator / (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
    { 
      return Div_<A,B>(a,b);
    }
//-----------------------------------------------------------------------------
 template <class A> 
 inline Sinus_<A>
  Sin (const BasisFunctionExpr<A>& a)
     { 
       return Sinus_<A>(a);
     }
//-----------------------------------------------------------------------------
 template <class A> 
 inline Cosinus_<A>
  Cos (const BasisFunctionExpr<A>& a)
     { 
       return Cosinus_<A>(a);
     }
//-----------------------------------------------------------------------------
 template <class A> 
 inline SQRT_<A>
  Sqrt (const BasisFunctionExpr<A>& a)
     { 
       return SQRT_<A>(a);
     }
//-----------------------------------------------------------------------------
 template <class A, typename eTYPE> 
 inline POW_<A,eTYPE>
  Pow (const BasisFunctionExpr<A>& a, eTYPE exponent)
     { 
       return POW_<A,eTYPE>(a,exponent);
     }
//-----------------------------------------------------------------------------
 template <class A> 
 inline EXP_<A>
  Exp(const BasisFunctionExpr<A>& a)
     { 
       return EXP_<A>(a);
     }
//-----------------------------------------------------------------------------
 template<class A> 
 inline Deri_<A,dx> 
  d_dx(const BasisFunctionExpr<A>& a)
      {
        return Deri_<A,dx>(a);
      }
//-----------------------------------------------------------------------------
 template<class A> 
 inline Deri_<A,dy>
  d_dy(const BasisFunctionExpr<A>& a)
      {
        return Deri_<A,dy>(a);
      }
//-----------------------------------------------------------------------------
 template<class A> 
 inline Deri_<A,dz> 
  d_dz(const BasisFunctionExpr<A>& a)
      {
        return Deri_<A,dz>(a);
      }
//==============================================================================

