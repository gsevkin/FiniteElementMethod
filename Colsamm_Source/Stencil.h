//==============================================================================
//
//  $Id: Stencil.h,v 1.29 2006/07/19 09:05:19 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
template <class TYPE> 
struct FunctionInterfaces
   {
     typedef TYPE (*functionInterfaceOne)(TYPE x);
     typedef TYPE (*functionInterfaceTwo)(TYPE x, TYPE y);
     typedef TYPE (*functionInterfaceThree)(TYPE x, TYPE y, TYPE z);
   };
//-----------------------------------------------------------------------------
template <class A> 
struct BasisFunctionExpr
   {
     enum{length = 1, size = 1};
     template<class Element, unsigned int iterLength> 
     static inline const typename Element::TYPE
      eval(const Element& element,const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <int uniqueFloatId>
class Float_  : public BasisFunctionExpr<Float_<uniqueFloatId> >
   {
    private:
     static float constantValue;
    public:
     enum {hasV = 0, hasW =0, hasP = 0} ;
     inline static void 
      print(std::ostream & os); 
     inline Float_(){ } 
     inline Float_(float constantValue_); 
     inline Float_<uniqueFloatId> 
      operator = (const float& constantValue_); 
     inline Float_<uniqueFloatId> 
      operator ()(float constantValue_); 
     template<class Element, unsigned int iterLength> 
     static inline const typename Element::TYPE
      eval(const Element& element,const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <int uniqueDoubleId>
class Double_  : public BasisFunctionExpr<Double_<uniqueDoubleId> >
   {
    private:
     static double constantValue;
    public:
     enum {hasV = 0, hasW =0, hasP = 0} ;
     inline static void 
      print(std::ostream & os); 
     inline Double_(){ } 
     inline Double_(double constantValue_); 
     inline Double_<uniqueDoubleId> 
      operator = (const double& constantValue_); 
     inline Double_<uniqueDoubleId> 
      operator ()(double constantValue_); 
     template<class Element, unsigned int iterLength> 
     static inline const typename Element::TYPE
      eval(const Element& element,const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <int uniqueComplexId>
struct Complex_  : public BasisFunctionExpr<Complex_<uniqueComplexId> >
   {
     private:
      static std::complex<double> constantValue;
     public:
      enum {hasV = 0, hasW =0, hasP = 0} ;
      inline static void 
       print(std::ostream & os);
      inline Complex_() {} 
      inline Complex_(const std::complex<double>& constantValue_); 
      inline Complex_<uniqueComplexId> 
       operator = (const std::complex<double>& constantValue_); 
      inline Complex_<uniqueComplexId>
       operator ()(const std::complex<double>& constantValue_); 
      template<class Element, unsigned int iterLength> 
      static inline const typename Element::TYPE&
      eval(const Element& element,const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, int functionId>
struct BasisFunction 
  : public BasisFunctionExpr<BasisFunction<number,functionId> >
    {
      enum { hasV = ((functionId == 0 && number == BasisSet_1)?1:0), 
             hasW = ((functionId == 1 && number == BasisSet_1)?1:0), 
             hasP = ((functionId == 0 && number == BasisSet_2)?1:0)};
      inline static void print(std::ostream &os); 
      template<class Element, unsigned int iterLength> 
      static inline const typename Element::TYPE
       eval(const Element& element,const int (&iterator)[iterLength]); 
      template<what typ, class Element, unsigned int iterLength> 
      static inline const typename Element::TYPE
       evaluateDerivative(const Element& element,const int (&iterator)[iterLength]); 
};
typedef BasisFunction<BasisSet_1,0> Ansatz_Function;
typedef BasisFunction<BasisSet_1,1> Testing_Function;
typedef BasisFunction<BasisSet_2,0> Mixed_Function;
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T>
struct FUNC1_ : public BasisFunctionExpr<FUNC1_<uniqueFunctionId,T> >
   {
     enum {hasV = 0, hasW =0, hasP = 0};
     static typename FunctionInterfaces<T>::functionInterfaceOne func_;
     inline static void 
      print(std::ostream &os); 
     inline FUNC1_(const typename FunctionInterfaces<T>::functionInterfaceOne& f); 
     template<class Element, unsigned int iterLength> 
     static inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T>
struct FUNC2_ : public BasisFunctionExpr<FUNC2_<uniqueFunctionId,T> >
   {
     enum {hasV = 0, hasW =0, hasP = 0};
     static typename FunctionInterfaces<T>::functionInterfaceTwo func_;
     inline static void 
      print(std::ostream &os); 
     inline FUNC2_(const typename FunctionInterfaces<T>::functionInterfaceTwo& f); 
     template<class Element, unsigned int iterLength> 
     static inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <int uniqueFunctionId, class T>
struct FUNC3_ : public BasisFunctionExpr<FUNC3_<uniqueFunctionId,T> >
   {
     enum {hasV = 0, hasW =0, hasP = 0};
     static typename FunctionInterfaces<T>::functionInterfaceThree func_;
     inline static void 
      print(std::ostream &os); 
     inline FUNC3_(const typename FunctionInterfaces<T>::functionInterfaceThree& f); 
     template<class Element, unsigned int iterLength> 
     static inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <int ex, SpaceDirection direction>
struct BasisMonom
  : public BasisFunctionExpr<BasisMonom<ex,direction> >
    {
      enum {hasV = 0, hasW = 0, hasP = 0};
      inline static void 
       print(std::ostream &os); 
      template<class Element, unsigned int iterLength> 
      static inline typename Element::TYPE 
       eval(const Element& element,const int (&iterator)[iterLength]); 
    };
//-----------------------------------------------------------------------------
namespace Interfaces 
   {
     template <int p>
     inline BasisMonom<p,dirX> x_(); 
     template <int p>
     inline BasisMonom<p,dirY> y_(); 
     template <int p>
     inline BasisMonom<p,dirZ> z_(); 
   }
//-----------------------------------------------------------------------------
template <int position>
struct NormalComponent 
  : public BasisFunctionExpr< NormalComponent<position> >
    {
      enum {hasV = 0, hasW =0, hasP = 0};
      inline static void 
       print(std::ostream &os); 
      template<class Element, unsigned int iterLength>
      static inline typename Element::TYPE 
       eval(const Element& element,const int (&iterator)[iterLength]); 
    };
//-----------------------------------------------------------------------------
template <int position>
struct UnitNormalComponent 
  : public BasisFunctionExpr< UnitNormalComponent<position> >
    {
      enum {hasV = 0, hasW =0, hasP = 0};
      inline static void 
       print(std::ostream &os); 
      template<class Element, unsigned int iterLength> 
      static inline typename Element::TYPE 
       eval(const Element& element,const int (&iterator)[iterLength]); 
    };
//-----------------------------------------------------------------------------
template <int uniqueVertexVectorId, class T>
struct VertexVector : public BasisFunctionExpr<VertexVector<uniqueVertexVectorId,T> >
   {
     static T data;
     enum {hasV = 1, hasW =0, hasP = 0};
     inline static void print(std::ostream &os); 
     inline VertexVector(const T& d); 
     template<class Element, unsigned int iterLength> 
     static inline typename Element::TYPE 
      eval(const Element& element,const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <class A> 
struct Conjugate : public BasisFunctionExpr<Conjugate<A> >
   {
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
    inline static void 
     print(std::ostream &os); 
    template<class Element, unsigned int iterLength>
    static inline typename Element::TYPE 
     eval(const Element& element,const int (&iterator)[iterLength]);
   };
//-----------------------------------------------------------------------------
template<class A> 
struct Min_ : public BasisFunctionExpr<Min_<A> >
   {
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline static void 
      print(std::ostream &os); 
     template<class Element, unsigned int iterLength>
     static inline typename Element::TYPE
      eval(const Element& element,const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template<class A, class B> 
struct Add_ : public BasisFunctionExpr<Add_<A,B> >
   {
     enum {hasV = _MAX_(A::hasV,B::hasV), 
           hasW = _MAX_(A::hasW,B::hasW),
           hasP = _MAX_(A::hasP,B::hasP)} ;
    inline static void 
     print(std::ostream & os); 
    template<class Element, unsigned int iterLength>
    static inline typename Element::TYPE
     eval(const Element& element,const int (&iterator)[iterLength]);
};
//-----------------------------------------------------------------------------
template<class A, class B>
struct Sub_ : public BasisFunctionExpr<Sub_<A,B> >
   {
     enum {hasV = _MAX_(A::hasV,B::hasV), 
           hasW = _MAX_(A::hasW,B::hasW),
           hasP = _MAX_(A::hasP,B::hasP)} ;
     inline static void
      print(std::ostream& os); 
     template<class Element, unsigned int iterLength>
     static inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]); 
};
//-----------------------------------------------------------------------------
template<class A, class B> 
struct Mult_: public BasisFunctionExpr<Mult_<A,B> >
   {
     enum {hasV = _MAX_(A::hasV,B::hasV), 
           hasW = _MAX_(A::hasW,B::hasW),
           hasP = _MAX_(A::hasP,B::hasP)} ;
     inline static void 
      print(std::ostream &os); 
     template<class Element, unsigned int iterLength>
     static inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template<class A, class B> 
struct Div_: public BasisFunctionExpr<Div_<A,B> >
   {
     enum {hasV = _MAX_(A::hasV,B::hasV), 
           hasW = _MAX_(A::hasW,B::hasW),
           hasP = _MAX_(A::hasP,B::hasP)} ;
     inline static void 
      print(std::ostream &os); 
     template<class Element, unsigned int iterLength>
     static inline typename Element::TYPE 
      eval(const Element& element,const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template<class A> 
struct SQRT_: public BasisFunctionExpr<SQRT_<A> >
   {
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline static void 
      print(std::ostream &os); 
     template<class Element, unsigned int iterLength>
     static inline typename Element::TYPE 
      eval(const Element& element,const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template<int number, class A, typename eTYPE> 
struct POW_: public BasisFunctionExpr<POW_<number,A,eTYPE> >
   {
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     static eTYPE exponent;
     inline POW_ (eTYPE e);
     inline static void 
      print(std::ostream &os); 
     template<class Element, unsigned int iterLength>
     static inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template<class A> 
struct EXP_: public BasisFunctionExpr<EXP_<A> >
   {
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline static void 
      print(std::ostream &os); 
     template<class Element, unsigned int iterLength>
     static inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]); 
   };
//-----------------------------------------------------------------------------
template <class A, what typ> 
struct Deri_ : public BasisFunctionExpr<Deri_<A,typ> >
   {
     inline static void 
      print(std::ostream &os); 
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     template<class Element, unsigned int iterLength>
     static inline const typename Element::TYPE
      eval(const Element& element, const int (&iterator)[iterLength]); 
   };
//==============================================================================

