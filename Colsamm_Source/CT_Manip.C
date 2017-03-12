//==============================================================================
//
//  $Id: CT_Manip.C,v 1.4 2006/06/09 09:50:53 jochen Exp $
//
//==============================================================================
#define Define_Integrand(A)  (A); \
		        __typeof__(A)
template <class A>
void 
 print(const A& a) 
    {
      std::cerr << std::endl; 
      A::print(std::cerr); 
      std::cerr << std::endl;
    }
//-----------------------------------------------------------------------------
// Data structures for rearranging the expressions. This is only applied for 
// the derivatives so far. But there might be a huge potential for optimizing, 
// e.g. the multiplication of a constant with a vector, where we should not
// immediately multiply every compontent with the constant, but only apply the 
// product to the result!
//-----------------------------------------------------------------------------
template <class A> 
struct MANIP
   { 
     typedef A I; 
   };
//-----------------------------------------------------------------------------
template <class A, class B, what i>
struct MANIP<Deri_<Add_<A,B>,i > >
   {
     typedef Add_<typename MANIP<Deri_<A,i> >::I,typename MANIP<Deri_<B,i> >::I> I;
   };
//-----------------------------------------------------------------------------
template <class A, class B, what i>
struct MANIP<Deri_<Sub_<A,B>,i > >
   {
     typedef Sub_<typename MANIP<Deri_<A,i> >::I,typename MANIP<Deri_<B,i> >::I> I;
   };
//-----------------------------------------------------------------------------
template <class A, class B, what i>
struct MANIP<Deri_<Mult_<A,B>,i> >
   {
     typedef Add_<Mult_<typename MANIP<Deri_<A,i> >::I,typename MANIP<B>::I>,
                  Mult_<typename MANIP<B>::I,typename MANIP<Deri_<B,i> >::I> > I;
   };
//-----------------------------------------------------------------------------
template <class A, what i, int j>
struct MANIP<Deri_<Mult_<A,Double_<j> >,i> >
   {
     typedef Mult_<Double_<j>,typename MANIP<Deri_<A,i> >::I> I;
   };
//-----------------------------------------------------------------------------
template <class A, what i, int j>
struct MANIP<Deri_<Mult_<Double_<j>,A>,i> >
   {
     typedef Mult_<Double_<j>,typename MANIP<Deri_<A,i> >::I> I;
   };
//-----------------------------------------------------------------------------
template <class A, what i, int j>
struct MANIP<Deri_<Mult_<A,Complex_<j> >,i> > 
   {
     typedef Mult_<Complex_<j>,typename MANIP<Deri_<A,i> >::I> I;
   };
//-----------------------------------------------------------------------------
template <class A, what i, int j>
struct MANIP<Deri_<Mult_<Complex_<j>,A>,i> >
   {
     typedef Mult_<Complex_<j>,typename MANIP<Deri_<A,i> >::I> I;
   };
//-----------------------------------------------------------------------------
template <class A, what i>
struct MANIP<Deri_<Conjugate<A>,i> >
   { 
     typedef Conjugate<typename MANIP<Deri_<A,i> >::I> I;  
   };
//-----------------------------------------------------------------------------
// Constructor functions to build the corresponding expression objects!
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
  template <int k,class T>
  inline VertexVector<k,T> 
   vec(const T(&d))
     {
       return VertexVector<k,T>(d);
     }
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC1_<k,float> 
   func(const typename FunctionInterfaces<float>::functionInterfaceOne& f_)
     {
       return FUNC1_<k,float>(f_);
     }
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC2_<k,float> 
   func(const typename FunctionInterfaces<float>::functionInterfaceTwo& f_)
     {
       return FUNC2_<k,float>(f_);
     }
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC3_<k,float> 
   func(const typename FunctionInterfaces<float>::functionInterfaceThree& f_)
     {
       return FUNC3_<k,float>(f_);
     }
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC1_<k,double>
   func(const typename FunctionInterfaces<double>::functionInterfaceOne& f_)
     {
       return FUNC1_<k,double>(f_);
     }
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC2_<k,double> 
   func(const typename FunctionInterfaces<double>::functionInterfaceTwo& f_)
     {
       return FUNC2_<k,double>(f_);
     }
//-----------------------------------------------------------------------------
  template <int k>
  inline FUNC3_<k,double> 
   func(const typename FunctionInterfaces<double>::functionInterfaceThree& f_)
     {
       return FUNC3_<k,double>(f_);
     }
//-----------------------------------------------------------------------------
  template <int k, class T>
  inline FUNC1_<k,T> 
   func(const typename FunctionInterfaces<T>::functionInterfaceOne& f_)
     {
       return FUNC1_<k,T>(f_);
     }
//-----------------------------------------------------------------------------
  template <int k, class T>
  inline FUNC2_<k,T> 
   func(const typename FunctionInterfaces<T>::functionInterfaceTwo& f_)
     { 
       return FUNC2_<k,T>(f_);
     }
//-----------------------------------------------------------------------------
 template <int k, class T>
 inline FUNC3_<k,T> 
  func(const typename FunctionInterfaces<T>::functionInterfaceThree& f_)
    {
      return FUNC3_<k,T>(f_);
    }
//-----------------------------------------------------------------------------
 template <class A> 
 inline Conjugate<A> 
  conj(const BasisFunctionExpr<A>& a) 
     {
       return Conjugate<A>();
     }
//-----------------------------------------------------------------------------
 template <class A, class B> inline Add_<A,B>
 operator + (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
    {
      return Add_<A,B>();
    }
//-----------------------------------------------------------------------------
 template <int A, int B> inline Double_<-(A+B)>
 operator + (const Double_<A>& a, const Double_<B>& b) 
    {
      return Double_<-(A+B)>(a.b_+b.b_);
    }
//-----------------------------------------------------------------------------
 template <class A, class B> inline Sub_<A,B>
 operator - (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
    {
      return Sub_<A,B>();
    }
//-----------------------------------------------------------------------------
 template <class A> inline Min_<A>
 operator - (const BasisFunctionExpr<A>& a)
    {
      return Min_<A>();
    }
//-----------------------------------------------------------------------------
 template <class A, class B> inline Mult_<A,B>
 operator * (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
    { 
      return Mult_<A,B>();
    }
//-----------------------------------------------------------------------------
 template <class A, class B> inline Div_<A,B>
 operator / (const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
    { 
      return Div_<A,B>();
    }
//-----------------------------------------------------------------------------
 template <class A> 
 inline SQRT_<A>
  Sqrt (const BasisFunctionExpr<A>& a)
     { 
       return SQRT_<A>();
     }
//-----------------------------------------------------------------------------
 template <int p, class A, typename eTYPE> 
 inline POW_<p,A,eTYPE>
  Pow (const BasisFunctionExpr<A>& a, eTYPE exponent)
     { 
       return POW_<p,A,eTYPE>(exponent);
     }
//-----------------------------------------------------------------------------
 template <class A> 
 inline EXP_<A>
  Exp(const BasisFunctionExpr<A>& a)
     { 
       return EXP_<A>();
     }
//-----------------------------------------------------------------------------
 template<class A> 
 inline typename MANIP<Deri_<A,dx> >::I 
  d_dx(const BasisFunctionExpr<A>& a)
      {
        return typename MANIP<Deri_<A,dx> >::I();
      }
//-----------------------------------------------------------------------------
 template<class A> 
 inline typename MANIP<Deri_<A,dy> >::I 
  d_dy(const BasisFunctionExpr<A>& a)
      {
        return typename MANIP<Deri_<A,dy> >::I();
      }
//-----------------------------------------------------------------------------
 template<class A> 
 inline typename MANIP<Deri_<A,dz> >::I 
  d_dz(const BasisFunctionExpr<A>& a)
      {
        return typename MANIP<Deri_<A,dz> >::I();
      }
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
void
VECTOR<A,B>::print(std::ostream &os)
   { 
#if 1
      A::print(os); 
      os << "\033[44m" << " & " << "\033[m" /*<< std::endl*/;  
      B::print(os);
#else
      os << "VECTOR < " ;  
      A::print(os); 
      os << " , "; 
      B::print(os); 
      os << " > ";
#endif
   }
//-----------------------------------------------------------------------------
template <class A, class B>
VECTOR<A,B> 
operator &(const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b)
   {
     return VECTOR<A,B>();
   }
//-----------------------------------------------------------------------------
/* A matrix is a concatenation of two vectors or an vector and a matrix. 
 * The result of the operator '|' is from left to right. Thus, the template 
 * type we achieve is something of the type 
 * MATRIX< MATRIX < MATRIX < ..., VECTOR <..> >, VECTOR <...> >, VECTOR <...> >
*/
template <class A, class B>
void
MATRIX<A,B>::print(std::ostream& os)
   { 
#if 1
      A::print(os); 
      os << "\033[42m" << " | " << "\033[m" << std::endl; 
      B::print(os);
#else
      os << "MATRIX < "; 
      A::print(os); 
      os << " ,\n          " ; 
      B::print(os); os << " >" ;
#endif
   }
//-----------------------------------------------------------------------------
template <class A, class B, class C, class D>
MATRIX<VECTOR<A,B>,VECTOR<C,D> > 
 operator | (const VECTOR<A,B>& a, const VECTOR<C,D>& b)
    {
      return MATRIX<VECTOR<A,B>,VECTOR<C,D> >();
    }
//-----------------------------------------------------------------------------
template <class A, class B, class C, class D>
MATRIX<MATRIX<A,B>,VECTOR<C,D> > 
 operator |(const MATRIX<A,B>& a, const VECTOR<C,D>& b)
    {
      return MATRIX<MATRIX<A,B>,VECTOR<C,D> >();
    }
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a vector with a one dimensional
 * value. The result is the vector, muliplied with the value in each component.
*/
template <class A, class B>
struct BasisFunctionExpr_VECTOR_MULT 
   : public BasisFunctionExpr<BasisFunctionExpr_VECTOR_MULT<A,B> > 
     {
       typedef Mult_<A,B> Erg;
     };
//-----------------------------------------------------------------------------
template <class A, class B1, class B2>
struct BasisFunctionExpr_VECTOR_MULT<A,VECTOR<B1,B2> >
   {
     typedef VECTOR<typename BasisFunctionExpr_VECTOR_MULT<A,B1>::Erg,Mult_<A,B2> > Erg;
   };
//-----------------------------------------------------------------------------
template <class A, class V1, class V2>
typename BasisFunctionExpr_VECTOR_MULT< A, VECTOR<V1,V2> > :: Erg
 operator * (const BasisFunctionExpr<A>& a, const VECTOR<V1,V2>& b)
    {
      return typename BasisFunctionExpr_VECTOR_MULT<A,VECTOR<V1,V2> > :: Erg();
    }
//-----------------------------------------------------------------------------
template <class A, class V1, class V2>
typename BasisFunctionExpr_VECTOR_MULT< A, VECTOR<V1,V2> > :: Erg
 operator * (const VECTOR<V1,V2>& a, const BasisFunctionExpr<A>& b)
    {
       return typename BasisFunctionExpr_VECTOR_MULT< A, VECTOR<V1,V2> > :: Erg();
    }
//-----------------------------------------------------------------------------
/* This struct describes the addition of two vectors. The result is the vector, 
 * that holds the sum of the vector components. 
*/
template <class A, class B>
struct CT_VECTOR_Add  
   {
     typedef Add_<A,B> Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
struct CT_VECTOR_Add<VECTOR<A1,A2>,VECTOR<B1,B2> > 
   {
     typedef VECTOR<typename CT_VECTOR_Add<A1,B1>::Erg, Add_<A2,B2> > Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename CT_VECTOR_Add<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg 
 operator+(const VECTOR<A1,A2>& a, const VECTOR<B1,B2>& b)
    {
      return typename CT_VECTOR_Add<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg();
    }

//-----------------------------------------------------------------------------
template <class A, class B>
struct CT_VECTOR_Sub  
   {
     typedef Sub_<A,B> Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
struct CT_VECTOR_Sub<VECTOR<A1,A2>,VECTOR<B1,B2> > 
   {
     typedef VECTOR<typename CT_VECTOR_Sub<A1,B1>::Erg, Sub_<A2,B2> > Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename CT_VECTOR_Sub<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg 
 operator-(const VECTOR<A1,A2>& a, const VECTOR<B1,B2>& b)
    {
      return typename CT_VECTOR_Sub<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg();
    }

//-----------------------------------------------------------------------------
template <class A>
struct CT_VECTOR_Min  
   {
     typedef Min_<A> Erg;
   };
//-----------------------------------------------------------------------------
template <class A, class B>
struct CT_VECTOR_Min<VECTOR<A,B> > 
   {
     typedef VECTOR<typename CT_VECTOR_Min<A>::Erg, Min_<B> > Erg;
   };
//-----------------------------------------------------------------------------
template <class A, class B>
typename CT_VECTOR_Min<VECTOR<A,B> >::Erg 
 operator-(const VECTOR<A,B>& a)
    {
      return typename CT_VECTOR_Min<VECTOR<A,B> >::Erg();
    }

//-----------------------------------------------------------------------------
/* This struct describes the multiplication of two vectors. The result is  
 * the sum that holds the products of the vectors components. 
*/
template <class A, class B>
struct VECTOR_VECTOR_MULT 
   {
     typedef Mult_<A,B> Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
struct VECTOR_VECTOR_MULT<VECTOR<A1,A2>,VECTOR<B1,B2> >
   {
     typedef Add_<typename VECTOR_VECTOR_MULT<A1,B1>::Erg,Mult_<A2,B2> > Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename VECTOR_VECTOR_MULT<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg 
 operator*(const VECTOR<A1,A2>& a, const VECTOR<B1,B2>& b)
    {
    //if ((int)VECTOR<A1,A2>::length != (int)VECTOR<B1,B2>::length)
      if (VECTOR<A1,A2>::length - VECTOR<B1,B2>::length != 0)
         { 
           std::cerr << "multiplication of vectors with wrong dimensions!" << std::endl << " Vector1: " << std::endl;
           VECTOR<A1,A2>::print(std::cerr);
           std::cerr << std::endl << " Vector2: " << std::endl;
           VECTOR<B1,B2>::print(std::cerr);
           std::cerr << std::endl;  exit(0);
         }
      return typename VECTOR_VECTOR_MULT<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg();
    }

//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a vector with a matrix!  
 * The result is a vector of the line sums . 
*/
template <class A, class B>
struct VECTOR_MATRIX_MULT;
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
struct VECTOR_MATRIX_MULT<VECTOR<A1,A2>, MATRIX<B1,B2> >
   {
     typedef typename CT_VECTOR_Add<
             typename VECTOR_MATRIX_MULT<A1,B1>::Erg, 
             typename BasisFunctionExpr_VECTOR_MULT<A2,B2>::Erg >::Erg Erg;
   };
//-----------------------------------------------------------------------------
template <class A, class B1, class B2>
struct VECTOR_MATRIX_MULT<A, VECTOR<B1,B2> > 
   {
     typedef typename BasisFunctionExpr_VECTOR_MULT<A,VECTOR<B1,B2> >::Erg Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename VECTOR_MATRIX_MULT<VECTOR<A1,A2>,MATRIX<B1,B2> >::Erg 
 operator*(const VECTOR<A1,A2>& a, const MATRIX<B1,B2>& b)
    {
   // if ((int)VECTOR<A1,A2>::length != (int)MATRIX<B1,B2>::size)
      if (VECTOR<A1,A2>::length - MATRIX<B1,B2>::size != 0)
         {
           std::cerr << "multiplication of matrix and vector with wrong dimensions!\n Matrix: " << std::endl;
           MATRIX<B1,B2>::print(std::cerr);
           std::cerr << "\n Vector: " << std::endl;
           VECTOR<A1,A2>::print(std::cerr);
           std::cerr << std::endl;  exit(0);
         }
      return typename VECTOR_MATRIX_MULT<VECTOR<A1,A2>,MATRIX<B1,B2> >::Erg();
    }
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a matrix with a vector!  
 * The result is a vector of the line sums . 
*/
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
struct MATRIX_VECTOR_MULT<MATRIX<A1,A2>, VECTOR<B1,B2> >
   {
     typedef VECTOR<typename MATRIX_VECTOR_MULT<A1,VECTOR<B1,B2> >::Erg,
     typename MATRIX_VECTOR_MULT<A2,VECTOR<B1,B2> >::Erg> Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
struct MATRIX_VECTOR_MULT<VECTOR<A1,A2>, VECTOR<B1,B2> >
   {
     typedef typename VECTOR_VECTOR_MULT<VECTOR<A1,A2>,VECTOR<B1,B2> >::Erg Erg;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename MATRIX_VECTOR_MULT<MATRIX<A1,A2>,VECTOR<B1,B2> >::Erg 
 operator*(const MATRIX<A1,A2>& a, const VECTOR<B1,B2>& b) 
    {
   // if ((int)MATRIX<A1,A2>::length != (int)VECTOR<B1,B2>::length)
      if (MATRIX<A1,A2>::length - VECTOR<B1,B2>::length != 0)
         { 
           std::cerr << "multiplication of matrix and vector with wrong dimensions!" << std::endl << " Matrix: " << std::endl;
           MATRIX<A1,A2>::print(std::cerr);
           std::cerr << std::endl << " Vector: " << std::endl;
           VECTOR<B1,B2>::print(std::cerr);
           std::cerr << std::endl;  exit(0);
         }
      return typename MATRIX_VECTOR_MULT<MATRIX<A1,A2>,VECTOR<B1,B2> >::Erg();
    }

// ----------------------------------------------------------------------------
 template <class A> 
 inline VECTOR< VECTOR<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>,
                       typename MANIP<Deri_<A,dz> >::I >
  grad(const BasisFunctionExpr<A>& a) 
     { 
       return VECTOR< VECTOR<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>,
                       typename MANIP<Deri_<A,dz> >::I >(); 
     }
// ----------------------------------------------------------------------------
 template <class A> 
 inline  VECTOR<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>
  grad2(const BasisFunctionExpr<A>& a) 
     { 
       return  VECTOR<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>(); 
     }
// ----------------------------------------------------------------------------
 template <class A, class B, class C> 
 inline Add_<Add_<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>,
                  typename MANIP<Deri_<A,dz> >::I >
  div(const VECTOR< VECTOR <A,B>, C >& a) 
     { 
       return Add_<Add_<typename MANIP<Deri_<A,dx> >::I,typename MANIP<Deri_<A,dy> >::I>,
                  typename MANIP<Deri_<A,dz> >::I > ();
     }
// ----------------------------------------------------------------------------
 inline VECTOR<VECTOR< NormalComponent<0>,NormalComponent<1> >, NormalComponent<2> >
  N() 
     { 
       return VECTOR<VECTOR< NormalComponent<0>,NormalComponent<1> >, NormalComponent<2> >();
     }
// ----------------------------------------------------------------------------
 inline VECTOR<VECTOR< UnitNormalComponent<0>,UnitNormalComponent<1> >, UnitNormalComponent<2> >
  N_Unit() 
     { 
       return VECTOR<VECTOR< UnitNormalComponent<0>,UnitNormalComponent<1> >, UnitNormalComponent<2> >();
     }
// ----------------------------------------------------------------------------

/* Concluding comments
 * At this point there is still more to introduce. A compile time data structure
 * which holds tensor data. Therefore we would need an additional operator to be 
 * overloaded, whicht has less priority than operator '&'. 
*/

//==============================================================================

