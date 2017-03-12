//==============================================================================
//
//  $Id: Stencil_ET.h,v 1.8 2006/07/25 13:52:27 jochen Exp $
//
//==============================================================================
#define REFERENZEN 0
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
     operator const A&()const { return *static_cast<const A*>(this); }
     enum{length = 1, size = 1};
     template<class Element, unsigned int iterLength> 
     inline const typename Element::TYPE
      eval(const Element& element,const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template <typename TYPE>
class Constant : public BasisFunctionExpr<Constant<TYPE> >
   {
      private: 
        TYPE constant;
      public:
      enum { hasV = 0, 
             hasW = 0, 
             hasP = 0};
      Constant (const TYPE c_) : constant(c_) {}
      inline void print(std::ostream &os) const { os << "C_(" << constant << ")";} 
       template<class Element, unsigned int iterLength> 
       inline const typename Element::TYPE
        eval(const Element& element,const int (&iterator)[iterLength]) const
          {
            return constant;
          } 
   };
//-----------------------------------------------------------------------------
template <>
class Constant<int> : public BasisFunctionExpr<Constant<int> >
   {
      private: 
        int constant;
      public:
      enum { hasV = 0, 
             hasW = 0, 
             hasP = 0};
       Constant (const int c_) : constant(c_) {}
      inline void print(std::ostream &os) const { os << "C_("<<constant<<")";} 
       template<class Element, unsigned int iterLength> 
       inline const typename Element::TYPE
        eval(const Element& element,const int (&iterator)[iterLength]) const
          {
            return constant;
          } 
   };
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId>
struct BasisFunction 
  : public BasisFunctionExpr<BasisFunction<number,functionId> >
    {
      enum { hasV = ((functionId == 0 && number == BasisSet_1)?1:0), 
             hasW = ((functionId == 1 && number == BasisSet_1)?1:0), 
             hasP = ((functionId == 0 && number == BasisSet_2)?1:0)};
      inline void print(std::ostream &os) const; 
      template<class Element, unsigned int iterLength> 
      inline const typename Element::TYPE
       eval(const Element& element,const int (&iterator)[iterLength]) const; 
      template<what type, class Element, unsigned int iterLength> 
      inline const typename Element::TYPE
       evaluateDerivative(const Element& element,const int (&iterator)[iterLength]) const; 
};
typedef BasisFunction<BasisSet_1,0> Ansatz_Function;
typedef BasisFunction<BasisSet_1,1> Testing_Function;
typedef BasisFunction<BasisSet_2,0> Mixed_Function;
//-----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId, SpaceDirection direction>
struct Vector_BasisFunction 
  : public BasisFunctionExpr<Vector_BasisFunction<number,functionId,direction> >
    {
      enum { hasV = ((functionId == 0 && number == BasisSet_1)?1:0), 
             hasW = ((functionId == 1 && number == BasisSet_1)?1:0), 
             hasP = ((functionId == 0 && number == BasisSet_2)?1:0)};
      inline void print(std::ostream &os) const; 
      template<class Element, unsigned int iterLength> 
      inline const typename Element::TYPE
       eval(const Element& element,const int (&iterator)[iterLength]) const; 
      template<what type, class Element, unsigned int iterLength> 
      inline const typename Element::TYPE
       evaluateDerivative(const Element& element,const int (&iterator)[iterLength]) const; 
};
//-----------------------------------------------------------------------------
template <SpaceDirection direction>
struct Edge_Component 
  : public BasisFunctionExpr<Edge_Component<direction> >
    {
      enum { hasV = 0, hasW = 0,  hasP = 0};
      inline void print(std::ostream &os) const { os << "EDGE <" << direction << "> ";} 
      template<class Element, unsigned int iterLength> 
      inline const typename Element::TYPE
       eval(const Element& element,const int (&iterator)[iterLength]) const
         {
           return element.edge_(direction);
         } 
};
//-----------------------------------------------------------------------------
template <class T>
class FUNC1_ : public BasisFunctionExpr<FUNC1_<T> >
   {
     private:
     typename FunctionInterfaces<T>::functionInterfaceOne func_;
     public:
     enum {hasV = 0, hasW =0, hasP = 0};
     inline void 
      print(std::ostream &os) const; 
     inline FUNC1_(const typename FunctionInterfaces<T>::functionInterfaceOne& f); 
     template<class Element, unsigned int iterLength> 
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template <class T>
class FUNC2_ : public BasisFunctionExpr<FUNC2_<T> >
   {
     private:
     typename FunctionInterfaces<T>::functionInterfaceTwo func_;
     public:
     enum {hasV = 0, hasW =0, hasP = 0};
     inline void 
      print(std::ostream &os) const; 
     inline FUNC2_(const typename FunctionInterfaces<T>::functionInterfaceTwo& f); 
     template<class Element, unsigned int iterLength> 
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template <class T>
class FUNC3_ : public BasisFunctionExpr<FUNC3_<T> >
   {
     private:
     typename FunctionInterfaces<T>::functionInterfaceThree func_;
     public:
     enum {hasV = 0, hasW =0, hasP = 0};
     inline void 
      print(std::ostream &os) const; 
     inline FUNC3_(const typename FunctionInterfaces<T>::functionInterfaceThree& f); 
     template<class Element, unsigned int iterLength> 
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template <int ex, SpaceDirection direction>
struct BasisMonom
  : public BasisFunctionExpr<BasisMonom<ex,direction> >
    {
      enum {hasV = 0, hasW = 0, hasP = 0};
      inline void 
       print(std::ostream &os) const; 
      template<class Element, unsigned int iterLength> 
      inline typename Element::TYPE 
       eval(const Element& element,const int (&iterator)[iterLength]) const; 
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
      inline void 
       print(std::ostream &os) const; 
      template<class Element, unsigned int iterLength>
      inline typename Element::TYPE 
       eval(const Element& element,const int (&iterator)[iterLength]) const; 
    };
//-----------------------------------------------------------------------------
template <int position>
struct UnitNormalComponent 
  : public BasisFunctionExpr< UnitNormalComponent<position> >
    {
      enum {hasV = 0, hasW =0, hasP = 0};
      inline void 
       print(std::ostream &os) const; 
      template<class Element, unsigned int iterLength> 
      inline typename Element::TYPE 
       eval(const Element& element,const int (&iterator)[iterLength]) const; 
    };
//-----------------------------------------------------------------------------
template <class T>
class VertexVector : public BasisFunctionExpr<VertexVector<T> >
   {
     private:
       const T& data;
     public:
     enum {hasV = 1, hasW =0, hasP = 0};
     inline void print(std::ostream &os) const; 
     inline VertexVector(const T& d); 
     template<class Element, unsigned int iterLength> 
     inline typename Element::TYPE 
      eval(const Element& element,const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template <class A> 
class Conjugate : public BasisFunctionExpr<Conjugate<A> >
   {
     private: 
#if REFERENZEN
       const A& a;
#else
       A a;
#endif
     public:
       inline Conjugate(const A& a_) : a(a_) {};
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
    inline void 
     print(std::ostream &os) const; 
    template<class Element, unsigned int iterLength>
    inline typename Element::TYPE 
     eval(const Element& element,const int (&iterator)[iterLength]) const;
   };
//-----------------------------------------------------------------------------
template<class A> 
class Min_ : public BasisFunctionExpr<Min_<A> >
   {
     private: 
#if REFERENZEN
       const A& a;
#else
       A a;
#endif
     public:
       inline Min_(const A& a_) : a(a_) {};
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE
      eval(const Element& element,const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template<class A, class B> 
class Add_ : public BasisFunctionExpr<Add_<A,B> >
   {
     private: 
#if REFERENZEN
       const A& a;
       const B& b;
#else
       A a;
       B b;
#endif
     public:
       inline Add_(const A& a_, const B& b_) : a(a_), b(b_) {};
     enum {hasV = _MAX_(A::hasV,B::hasV), 
           hasW = _MAX_(A::hasW,B::hasW),
           hasP = _MAX_(A::hasP,B::hasP)} ;
    inline void 
     print(std::ostream & os) const; 
    template<class Element, unsigned int iterLength>
    inline typename Element::TYPE
     eval(const Element& element,const int (&iterator)[iterLength]) const;
};
//-----------------------------------------------------------------------------
template<class TYPE, typename B> 
class CAdd_ : public BasisFunctionExpr<CAdd_<TYPE,B> >
   {
     private: 
       const TYPE a;
#if REFERENZEN
       const B& b;
#else
       B b;
#endif
     public:
       inline CAdd_(const TYPE& a_, const B& b_) : a(a_), b(b_) {};
     enum {hasV = B::hasV, 
           hasW = B::hasW,
           hasP = B::hasP} ;
    inline void 
     print(std::ostream & os) const; 
    template<class Element, unsigned int iterLength>
    inline typename Element::TYPE
     eval(const Element& element,const int (&iterator)[iterLength]) const;
};
//-----------------------------------------------------------------------------
template<class A, class B>
class Sub_ : public BasisFunctionExpr<Sub_<A,B> >
   {
     private: 
#if REFERENZEN
       const A& a;
       const B& b;
#else
       A a;
       B b;
#endif
     public:
     inline Sub_(const A& a_, const B& b_) : a(a_), b(b_) {};
     enum {hasV = _MAX_(A::hasV,B::hasV), 
           hasW = _MAX_(A::hasW,B::hasW),
           hasP = _MAX_(A::hasP,B::hasP)} ;
     inline void
      print(std::ostream& os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
};
//-----------------------------------------------------------------------------
template<class TYPE, class B>
class CSub_ : public BasisFunctionExpr<CSub_<TYPE,B> >
   {
     private: 
       const TYPE a;
#if REFERENZEN
       const B& b;
#else
       B b;
#endif
     public:
     inline CSub_(const TYPE& a_, const B& b_) : a(a_), b(b_) {};
     enum {hasV = B::hasV, 
           hasW = B::hasW,
           hasP = B::hasP} ;
     inline void
      print(std::ostream& os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
};
//-----------------------------------------------------------------------------
template<class A, class B> 
class Mult_: public BasisFunctionExpr<Mult_<A,B> >
   {
     private: 
#if REFERENZEN
       const A& a;
       const B& b;
#else
       A a;
       B b;
#endif
     public:
     inline Mult_(const A& a_, const B& b_) : a(a_), b(b_) {};
     enum {hasV = _MAX_(A::hasV,B::hasV), 
           hasW = _MAX_(A::hasW,B::hasW),
           hasP = _MAX_(A::hasP,B::hasP)} ;
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template<class TYPE, class B> 
class CMult_: public BasisFunctionExpr<CMult_<TYPE,B> >
   {
     private: 
       const TYPE a;
#if REFERENZEN
       const B& b;
#else
       B b;
#endif
     public:
     inline CMult_(const TYPE& a_, const B& b_) : a(a_), b(b_) {};
     enum {hasV = B::hasV, 
           hasW = B::hasW,
           hasP = B::hasP} ;
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template<class A, class B> 
struct Div_: public BasisFunctionExpr<Div_<A,B> >
   {
     private: 
#if REFERENZEN
       const A& a;
       const B& b;
#else
       A a;
       B b;
#endif
     public:
     inline Div_(const A& a_, const B& b_) : a(a_), b(b_) {};
     enum {hasV = _MAX_(A::hasV,B::hasV), 
           hasW = _MAX_(A::hasW,B::hasW),
           hasP = _MAX_(A::hasP,B::hasP)} ;
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element,const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template<class A> 
class Sinus_: public BasisFunctionExpr<Sinus_<A> >
   {
     private: 
#if REFERENZEN
       const A& a;
#else
       A a;
#endif
     public:
       inline Sinus_(const A& a_) : a(a_) {};
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element,const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template<class A> 
class Cosinus_: public BasisFunctionExpr<Cosinus_<A> >
   {
     private: 
#if REFERENZEN
       const A& a;
#else
       A a;
#endif
     public:
       inline Cosinus_(const A& a_) : a(a_) {};
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element,const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template<class A> 
class SQRT_: public BasisFunctionExpr<SQRT_<A> >
   {
     private: 
#if REFERENZEN
       const A& a;
#else
       A a;
#endif
     public:
       inline SQRT_(const A& a_) : a(a_) {};
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element,const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template<class A, typename eTYPE> 
class POW_: public BasisFunctionExpr<POW_<A,eTYPE> >
   {
     private: 
#if REFERENZEN
       const A& a;
#else
       A a;
#endif
       eTYPE exponent;
     public:
       inline POW_(const A& a_, eTYPE expo_) : a(a_), exponent(expo_) {};
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline POW_ (eTYPE e);
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template<class A> 
class EXP_: public BasisFunctionExpr<EXP_<A> >
   {
     private: 
#if REFERENZEN
       const A& a;
#else
       A a;
#endif
     public:
       inline EXP_(const A& a_) : a(a_) {};
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     inline void 
      print(std::ostream &os) const; 
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
   };
//-----------------------------------------------------------------------------
template <class A, what typ> 
class Deri_ : public BasisFunctionExpr<Deri_<A,typ> >
   {
     private: 
#if REFERENZEN
       const A& a;
#else
       A a;
#endif
     public:
       inline Deri_(const A& a_) : a(a_) {};
     inline void 
      print(std::ostream &os) const; 
     enum {hasV = A::hasV,
           hasW = A::hasW,
           hasP = A::hasP};
     template<class Element, unsigned int iterLength>
     inline const typename Element::TYPE
      eval(const Element& element, const int (&iterator)[iterLength]) const; 
   };
//==============================================================================

