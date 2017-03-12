//==============================================================================
//
//  $Id: Integrand_Vector.h,v 1.2 2006/07/19 11:30:48 jochen Exp $
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
template <class A, int vLength>
struct BasisFunctionVectorExpr
   {
      operator const A&() const { return *static_cast<const A*>(this); }
      enum {length = vLength, size = 1};
   };
//-----------------------------------------------------------------------------
template <class A, class B >
class Integrand_Vector 
   : public BasisFunctionVectorExpr<Integrand_Vector<A,B>,A::length+B::length>
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
      inline Integrand_Vector(const A& a_, const B& b_) : a(a_), b(b_) {}
       enum {hasV = _MAX_(A::hasV,B::hasV),
             hasW = _MAX_(A::hasW,B::hasW), 
             hasP = _MAX_(A::hasP,B::hasP)} ;
       inline void print(std::ostream &os) const;
       inline const A& getHeadOfVector() const { return a;}
       const A& getA() const {return a;}
       const B& getB() const {return b;}
       template<class Element, unsigned int iterLength>
        inline typename Element::TYPE 
         evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const;
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B>
class Integrand_Vector<Integrand_Vector<A1,A2>,B> 
   : public BasisFunctionVectorExpr<Integrand_Vector<Integrand_Vector<A1,A2>,B>,Integrand_Vector<A1,A2>::length+B::length>
   {
     private:
#if REFERENZEN
      const Integrand_Vector<A1,A2>& a;
      const B& b;
#else
       Integrand_Vector<A1,A2> a;
       B b;
#endif
     public:
     inline Integrand_Vector(const Integrand_Vector<A1,A2>& a_, const B& b_) : a(a_), b(b_) {}
     enum {hasV = _MAX_((Integrand_Vector<A1,A2>::hasV),B::hasV),
           hasW = _MAX_((Integrand_Vector<A1,A2>::hasW),B::hasW), 
           hasP = _MAX_((Integrand_Vector<A1,A2>::hasP),B::hasP)} ;
     inline void print(std::ostream &os) const;
     inline const Integrand_Vector<A1,A2>& getHeadOfVector() const { return a;}
     const A1& getA() const {return a.getA();}
     const A2& getB() const {return a.getB();}
     const B&  getC() const {return b;}
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const;
};
//-----------------------------------------------------------------------------
template <class A, class B>
inline Integrand_Vector<A,B> 
operator &(const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b);
//-----------------------------------------------------------------------------
template <class A, class B, int lA>
inline Integrand_Vector<A,B> 
operator &(const BasisFunctionVectorExpr<A,lA>& a, const BasisFunctionExpr<B>& b);
//-----------------------------------------------------------------------------
template <typename TYPE, class B>
class ConstVectorMult_;
//-----------------------------------------------------------------------------
template<class A, class B, int l, class Element, unsigned int iterLength>
inline typename Element::TYPE 
 evaluate(const BasisFunctionVectorExpr<A,l>& a, const BasisFunctionVectorExpr<B,l>& b,
          const Element& element, const int (&iterator)[iterLength]); 
//-----------------------------------------------------------------------------
template<class A, class B, class Element, unsigned int iterLength>
inline typename Element::TYPE 
 evaluate(const BasisFunctionExpr<A>& a, const BasisFunctionExpr<B>& b,
          const Element& element, const int (&iterator)[iterLength]); 
//-----------------------------------------------------------------------------
template <class A, class B>
class VectorVectorMult_ : public BasisFunctionExpr<VectorVectorMult_<A,B> >
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
     inline VectorVectorMult_(const A& a_, const B& b_) : a(a_), b(b_) {}
     enum {hasV = _MAX_(A::hasV,B::hasV),
           hasW = _MAX_(A::hasW,B::hasW), 
           hasP = _MAX_(A::hasP,B::hasP)} ;
     inline void 
      print(std::ostream &os) const
        {
            os << " ("; a.print(os); os << ") * (" ; b.print(os) ; os << ") ";
        }
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      eval(const Element& element, const int (&iterator)[iterLength]) const
        {
          return evaluate(a,b,element,iterator);
        } 
   };
//-----------------------------------------------------------------------------
template <typename TYPE, class B1, class B2>
class ConstVectorMult_<TYPE,Integrand_Vector<B1,B2> > 
   : public BasisFunctionVectorExpr<ConstVectorMult_<TYPE,Integrand_Vector<B1,B2> >,Integrand_Vector<B1,B2>::length >
   {
     private: 
      TYPE a;
#if REFERENZEN
      const Integrand_Vector<B1,B2> & b;
#else
      Integrand_Vector<B1,B2> b;
#endif
     public:
     inline ConstVectorMult_(const TYPE& a_, const Integrand_Vector<B1,B2>& b_) : a(a_), b(b_) {}
     enum {hasV = Integrand_Vector<B1,B2>::hasV,
           hasW = Integrand_Vector<B1,B2>::hasW, 
           hasP = Integrand_Vector<B1,B2>::hasP} ;
     inline void 
      print(std::ostream &os) const
        {
            os << " " << a << " * (" ; b.print(os) ; os << ") ";
        }
     inline const CMult_<TYPE,B1> getHeadOfVector() const { return a * b.getHeadOfVector();}
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
        {
          return a*b.evalLastComponent(element,iterator);
        } 
   };
//-----------------------------------------------------------------------------
template <typename TYPE, class B1, class B2, class C>
class ConstVectorMult_<TYPE,Integrand_Vector<Integrand_Vector<B1,B2>,C> > 
    : public BasisFunctionVectorExpr<ConstVectorMult_<TYPE,Integrand_Vector<Integrand_Vector<B1,B2>,C> >,Integrand_Vector<Integrand_Vector<B1,B2>,C>::length >
   {
     private: 
      TYPE a;
#if REFERENZEN
      const Integrand_Vector<Integrand_Vector<B1,B2>,C>& b;
#else
      Integrand_Vector<Integrand_Vector<B1,B2>,C> b;
#endif
     public:
     inline ConstVectorMult_(const TYPE& a_, const Integrand_Vector<Integrand_Vector<B1,B2>,C>& b_) : a(a_), b(b_) {}
     enum {hasV = Integrand_Vector<Integrand_Vector<B1,B2>,C>::hasV,
           hasW = Integrand_Vector<Integrand_Vector<B1,B2>,C>::hasW, 
           hasP = Integrand_Vector<Integrand_Vector<B1,B2>,C>::hasP} ;
     inline void 
      print(std::ostream &os) const
        {
            os << " " << a << " * (" ; b.print(os) ; os << ") ";
        }
     inline const ConstVectorMult_<TYPE,Integrand_Vector<B1,B2> > getHeadOfVector() const { return a*b.getHeadOfVector();}
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
        {
          return a*b.evalLastComponent(element,iterator);
        } 
   };
//-----------------------------------------------------------------------------
template <class A, class B, int l>
inline VectorVectorMult_<A,B>
 operator*(const BasisFunctionVectorExpr<A,l>& a, const BasisFunctionVectorExpr<B,l>& b)
   {
     return VectorVectorMult_<A,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class B, int l>
inline ConstVectorMult_<int,B>
 operator*(int a, const BasisFunctionVectorExpr<B,l>& b)
   {
     return ConstVectorMult_<int,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class B, int l>
inline ConstVectorMult_<float,B>
 operator*(float a, const BasisFunctionVectorExpr<B,l>& b)
   {
     return  ConstVectorMult_<float,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class B, int l>
inline ConstVectorMult_<double,B>
 operator*(double a, const BasisFunctionVectorExpr<B,l>& b)
   {
     return ConstVectorMult_<double,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <typename TYPE, class B, int l>
inline ConstVectorMult_<typename std::complex<TYPE>,B>
 operator*(typename std::complex<TYPE> a, const BasisFunctionVectorExpr<B,l>& b)
   {
     return  ConstVectorMult_<typename std::complex<TYPE>,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class B, int l>
inline ConstVectorMult_<int,B>
 operator*(const BasisFunctionVectorExpr<B,l>& b, int a)
   {
     return ConstVectorMult_<int,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class B, int l>
inline ConstVectorMult_<float,B>
 operator*(const BasisFunctionVectorExpr<B,l>& b, float a)
   {
     return  ConstVectorMult_<float,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class B, int l>
inline ConstVectorMult_<double,B>
 operator*(const BasisFunctionVectorExpr<B,l>& b, double a)
   {
     return ConstVectorMult_<double,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <typename TYPE, class B, int l>
inline ConstVectorMult_<typename std::complex<TYPE>,B>
 operator*(const BasisFunctionVectorExpr<B,l>& b, typename std::complex<TYPE> a)
   {
     return  ConstVectorMult_<typename std::complex<TYPE>,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class A, class B>
class BasisVectorMult_;
//-----------------------------------------------------------------------------
template <class A, class B1, class B2>
class BasisVectorMult_<A, Integrand_Vector<B1,B2> > 
   : public BasisFunctionVectorExpr<BasisVectorMult_<A, Integrand_Vector<B1,B2> >,Integrand_Vector<B1,B2>::length>
   {
     private: 
#if REFERENZEN
      const A& a;
      const Integrand_Vector<B1,B2> & b;
#else
      A a;
      Integrand_Vector<B1,B2> b;
#endif
     public:
     inline BasisVectorMult_(const A& a_, const Integrand_Vector<B1,B2>& b_) : a(a_), b(b_) {}
     enum {hasV = _MAX_(A::hasV,(Integrand_Vector<B1,B2>::hasV)),
           hasW = _MAX_(A::hasW,(Integrand_Vector<B1,B2>::hasW)), 
           hasP = _MAX_(A::hasP,(Integrand_Vector<B1,B2>::hasP))} ;
     inline void 
      print(std::ostream &os) const
        {
            os << " ("; a.print(os); os << ") * (" ; b.print(os) ; os << ") ";
        }
     inline const Mult_<A,B1> getHeadOfVector() const { return a * b.getHeadOfVector();}
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
        {
          return a.eval(element,iterator)*b.evalLastComponent(element,iterator);
        } 
   };
//-----------------------------------------------------------------------------
template <class A, class B1, class B2, class C>
class BasisVectorMult_<A,Integrand_Vector<Integrand_Vector<B1,B2>,C> > 
    : public BasisFunctionVectorExpr<BasisVectorMult_<A,Integrand_Vector<Integrand_Vector<B1,B2>,C> >,
                                     Integrand_Vector<Integrand_Vector<B1,B2>,C>::length >
   {
     private: 
#if REFERENZEN
      const A& a;
      const Integrand_Vector<Integrand_Vector<B1,B2>,C>& b;
#else
      A a;
      Integrand_Vector<Integrand_Vector<B1,B2>,C> b;
#endif
     public:
     inline BasisVectorMult_(const A& a_, const Integrand_Vector<Integrand_Vector<B1,B2>,C>& b_) : a(a_), b(b_) {}
     enum {hasV = _MAX_(A::hasV,(Integrand_Vector<Integrand_Vector<B1,B2>,C>::hasV)),
           hasW = _MAX_(A::hasW,(Integrand_Vector<Integrand_Vector<B1,B2>,C>::hasW)), 
           hasP = _MAX_(A::hasP,(Integrand_Vector<Integrand_Vector<B1,B2>,C>::hasP))} ;
     inline void 
      print(std::ostream &os) const
        {
            os << " ("; a.print(os); os << ") * (" ; b.print(os) ; os << ") ";
        }
     inline const BasisVectorMult_<A,Integrand_Vector<B1,B2> > getHeadOfVector() const { return a*b.getHeadOfVector();}
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
        {
          return a.eval(element,iterator)*b.evalLastComponent(element,iterator);
        } 
   };
//-----------------------------------------------------------------------------
template <class A, class B, int l>
inline BasisVectorMult_<A,B>
 operator*(const BasisFunctionExpr<A>& a, const BasisFunctionVectorExpr<B,l>& b)
   {
     return BasisVectorMult_<A,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class A, class B, int l>
inline BasisVectorMult_<A,B>
 operator*(const BasisFunctionVectorExpr<B,l>& b, const BasisFunctionExpr<A>& a)
   {
     return BasisVectorMult_<A,B>(a,b);
   } 
//-----------------------------------------------------------------------------
template <class A>
class GRADIENT2D : public BasisFunctionVectorExpr<GRADIENT2D<A>,2>
   {
     private:
#if REFERENZEN
      const A& a;
#else
      A a;
#endif
     public:
        GRADIENT2D(const A& a_) : a(a_) {}
        enum {hasV = A::hasV,
              hasW = A::hasW, 
              hasP = A::hasP} ;
        inline Deri_<A,dx> getHeadOfVector() const { return d_dx(a);}
        template<class Element, unsigned int iterLength>
        inline typename Element::TYPE 
         evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
           { 
              return (d_dy(a)).eval(element,iterator);
           }
        
   };
//-----------------------------------------------------------------------------
#if 0
template <class A, class B>
inline VectorVectorMult_<A,GRADIENT2D<B> >
 operator*(const BasisFunctionVectorExpr<A,2>& a, const GRADIENT2D<B>& b)
   {
     return VectorVectorMult_<A,GRADIENT2D<B> >(a,b);
   } 
#endif
//-----------------------------------------------------------------------------
template <class A>
class GRADIENT3D : public BasisFunctionVectorExpr<GRADIENT3D<A>,3>
   {
     private:
#if REFERENZEN
      const A& a;
#else
      A a;
#endif
     public:
        GRADIENT3D(const A& a_) : a(a_){}
        enum {hasV = A::hasV,
              hasW = A::hasW, 
              hasP = A::hasP} ;
        inline GRADIENT2D<A> getHeadOfVector() const { return GRADIENT2D<A>(a);}
        template<class Element, unsigned int iterLength>
        inline typename Element::TYPE 
         evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
           { 
              return (d_dz(a)).eval(element,iterator);
           }
        
   };
//-----------------------------------------------------------------------------
template <typename TYPE, class B>
class ConstVectorMult_<TYPE,GRADIENT3D<B> > 
   : public BasisFunctionVectorExpr<ConstVectorMult_<TYPE,GRADIENT3D<B> >,3>
   {
     private: 
      TYPE a;
#if REFERENZEN
      const GRADIENT3D<B> & b;
#else
      GRADIENT3D<B> b;
#endif
     public:
     inline ConstVectorMult_(const TYPE& a_, const GRADIENT3D<B>& b_) : a(a_), b(b_) {}
     enum {hasV = GRADIENT3D<B>::hasV,
           hasW = GRADIENT3D<B>::hasW, 
           hasP = GRADIENT3D<B>::hasP} ;
     inline void 
      print(std::ostream &os) const
        {
            os << " " << a << " * (" ; b.print(os) ; os << ") ";
        }
     inline const ConstVectorMult_<TYPE,GRADIENT2D<B> > getHeadOfVector() const { return a * b.getHeadOfVector();}
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
        {
          return a*b.evalLastComponent(element,iterator);
        } 
   };
//-----------------------------------------------------------------------------
template <typename TYPE, class B>
class ConstVectorMult_<TYPE,GRADIENT2D<B> > 
   : public BasisFunctionVectorExpr<ConstVectorMult_<TYPE,GRADIENT2D<B> >,2>
   {
     private: 
      TYPE a;
#if REFERENZEN
      const GRADIENT2D<B> & b;
#else
      GRADIENT2D<B> b;
#endif
     public:
     inline ConstVectorMult_(const TYPE& a_, const GRADIENT2D<B>& b_) : a(a_), b(b_) {}
     enum {hasV = GRADIENT2D<B>::hasV,
           hasW = GRADIENT2D<B>::hasW, 
           hasP = GRADIENT2D<B>::hasP} ;
     inline void 
      print(std::ostream &os) const
        {
            os << " " << a << " * (" ; b.print(os) ; os << ") ";
        }
     inline const CMult_<TYPE,Deri_<B,dx> > getHeadOfVector() const { return a * b.getHeadOfVector();}
     template<class Element, unsigned int iterLength>
     inline typename Element::TYPE 
      evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
        {
          return a*b.evalLastComponent(element,iterator);
        } 
   };
//-----------------------------------------------------------------------------
 template <class A, int length> 
 ConstVectorMult_<double,A> operator-(const BasisFunctionVectorExpr<A,length>& a){
     return ConstVectorMult_<double,A>(-1., a); 
 }
//-----------------------------------------------------------------------------
 template <class A> 
 inline GRADIENT3D<A>
  grad(const BasisFunctionExpr<A>& a)
     {
        return GRADIENT3D<A>(a);
     } 
//-----------------------------------------------------------------------------
 template <class A> 
 inline GRADIENT2D<A>
  grad2(const BasisFunctionExpr<A>& a)
     {
        return GRADIENT2D<A>(a);
     } 
//-----------------------------------------------------------------------------
 template <class A> 
 inline Add_<Add_<Deri_<A,dx> ,Deri_<A,dy> >, Deri_<A,dz>  >
  div(const BasisFunctionExpr<A>& a) {
       return d_dx(a) +d_dy(a)+d_dz(a);
    }
//-----------------------------------------------------------------------------
 template <class A, class B, class C> 
 inline Add_<Add_<Deri_<A,dx> ,Deri_<A,dy> >, Deri_<A,dz>  >
  div(const Integrand_Vector< Integrand_Vector <A,B>, C >& a); 
// ----------------------------------------------------------------------------
inline Integrand_Vector< NormalComponent<0>,NormalComponent<1> >
N_2D(); 
// ----------------------------------------------------------------------------
inline Integrand_Vector< UnitNormalComponent<0>,UnitNormalComponent<1> >
N_Unit_2D(); 
// ----------------------------------------------------------------------------
inline Integrand_Vector<Integrand_Vector< NormalComponent<0>,NormalComponent<1> >, NormalComponent<2> >
N(); 
// ----------------------------------------------------------------------------
inline Integrand_Vector<Integrand_Vector< UnitNormalComponent<0>,UnitNormalComponent<1> >, UnitNormalComponent<2> >
N_Unit(); 
// ----------------------------------------------------------------------------
template <BasisFunctionSet number, unsigned int functionId>
struct BF_Integrand_Vector
   {
     typedef Integrand_Vector<Vector_BasisFunction<number,functionId,dirX>,
                    Vector_BasisFunction<number,functionId,dirY> >  
                                                                      VEC_FUNC_2D;
     typedef Integrand_Vector<Integrand_Vector<Vector_BasisFunction<number,functionId,dirX>,
                           Vector_BasisFunction<number,functionId,dirY> >,
                    Vector_BasisFunction<number,functionId,dirZ> >         VEC_FUNC_3D;

   };
// ----------------------------------------------------------------------------
 inline BF_Integrand_Vector<BasisSet_1,0>::VEC_FUNC_2D
  v_vec2D ()
    {
      return ( Vector_BasisFunction<BasisSet_1,0,dirX>() & Vector_BasisFunction<BasisSet_1,0,dirY>() );
    }

// ----------------------------------------------------------------------------
 inline BF_Integrand_Vector<BasisSet_1,1>::VEC_FUNC_2D
  w_vec2D ()
    {
      return ( Vector_BasisFunction<BasisSet_1,1,dirX>() & Vector_BasisFunction<BasisSet_1,1,dirY>() );
    }

// ----------------------------------------------------------------------------
 inline BF_Integrand_Vector<BasisSet_1,0>::VEC_FUNC_3D
  v_vec3D ()
    {
      return Vector_BasisFunction<BasisSet_1,0,dirX>() & 
             Vector_BasisFunction<BasisSet_1,0,dirY>() &
             Vector_BasisFunction<BasisSet_1,0,dirZ>();

    }
// ----------------------------------------------------------------------------
 inline BF_Integrand_Vector<BasisSet_1,1>::VEC_FUNC_3D
  w_vec3D ()
    {
      return Vector_BasisFunction<BasisSet_1,1,dirX>() & 
             Vector_BasisFunction<BasisSet_1,1,dirY>() &
             Vector_BasisFunction<BasisSet_1,1,dirZ>();

    }
// ----------------------------------------------------------------------------
template <class A, class B>
class CURL2D : public BasisFunctionExpr<CURL2D<A,B> >
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
        enum {hasV = _MAX_(A::hasV,B::hasV),
              hasW = _MAX_(A::hasW,B::hasW),
              hasP = _MAX_(A::hasP,B::hasP)} ;
       CURL2D(const Integrand_Vector<A,B>& a_) : a(a_.getA()), b(a_.getB()) {}
       CURL2D(const A& a_, const B& b_) : a(a_), b(b_) {}
       template<class Element, unsigned int iterLength>
       inline typename Element::TYPE 
         eval(const Element& element, const int (&iterator)[iterLength]) const
           {
             return d_dx(b).eval(element,iterator) - d_dy(a).eval(element,iterator);
           }
   };
// ----------------------------------------------------------------------------
 template <class A, class B>
 inline CURL2D<A,B>
  curl(const Integrand_Vector<A,B>& a)
    {
      return CURL2D<A,B>(a);
    }
// ----------------------------------------------------------------------------
template <class B, class C>
class CURL3D_Help_1 : public BasisFunctionExpr<CURL3D_Help_1<B,C> >
   {
     private:
#if REFERENZEN
       const B& b;
       const C& c;
#else
       B b;
       C c;
#endif
     public:
        enum {hasV = _MAX_(C::hasV,B::hasV),
              hasW = _MAX_(C::hasW,B::hasW),
              hasP = _MAX_(C::hasP,B::hasP)} ;
        enum {length = 1, size = 1};
       CURL3D_Help_1(const B& b_, const C& c_) : b(b_), c(c_) {}
       template<class Element, unsigned int iterLength>
       inline typename Element::TYPE 
        eval(const Element& element, const int (&iterator)[iterLength]) const
           { 
             return d_dy(c).eval(element,iterator) - d_dz(b).eval(element,iterator);
           }
   };
// ----------------------------------------------------------------------------
template <class A, class B, class C>
class CURL3D_Help : public BasisFunctionVectorExpr<CURL3D_Help<A,B,C>,3>
   {
     private:
#if REFERENZEN
       const A& a;
       const B& b;
       const C& c;
#else
       A a;
       B b;
       C c;
#endif
     public:
        enum {hasV = _MAX_(_MAX_(A::hasV,B::hasV),C::hasV),
              hasW = _MAX_(_MAX_(A::hasW,B::hasW),C::hasW),
              hasP = _MAX_(_MAX_(A::hasP,B::hasP),C::hasP)} ;
       CURL3D_Help(const A& a_, const B& b_, const C& c_) : a(a_), b(b_), c(c_) {}
       inline CURL3D_Help_1<B,C> getHeadOfVector() const { return CURL3D_Help_1<B,C>(b,c);}
       template<class Element, unsigned int iterLength>
       inline typename Element::TYPE 
        evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
           { 
             return d_dx(b).eval(element,iterator) - d_dy(a).eval(element,iterator);
           }
   };
// ----------------------------------------------------------------------------
template <class A, class B, class C>
class CURL3D : public BasisFunctionVectorExpr<CURL3D<A,B,C>,3>
   {
     private:
#if REFERENZEN
       const A& a;
       const B& b;
       const C& c;
#else
       A a;
       B b;
       C c;
#endif
     public:
       enum {hasV = _MAX_(_MAX_(A::hasV,B::hasV),C::hasV),
              hasW = _MAX_(_MAX_(A::hasW,B::hasW),C::hasW),
              hasP = _MAX_(_MAX_(A::hasP,B::hasP),C::hasP)} ;
       CURL3D(const Integrand_Vector<Integrand_Vector<A,B>,C>& a_) : a(a_.getA()), b(a_.getB()) , c(a_.getC()) {}
       inline CURL3D_Help<A,B,C> getHeadOfVector() const { return CURL3D_Help<A,B,C>(a,b,c);}
       template<class Element, unsigned int iterLength>
       inline typename Element::TYPE 
        evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
           { 
             return d_dz(a).eval(element,iterator) - d_dx(c).eval(element,iterator);
           }
   };
// ----------------------------------------------------------------------------
 template <class A, class B, class C>
 inline CURL3D<A,B,C>
  curl(const Integrand_Vector<Integrand_Vector<A,B>,C>& a)
    {
      return CURL3D<A,B,C>(a);
    }
// ----------------------------------------------------------------------------
struct EDGE2D : public BasisFunctionVectorExpr<EDGE2D,2>
   {
       enum {hasV = 0, hasW = 0, hasP = 0} ;
       inline Edge_Component<dirX> getHeadOfVector() const { return Edge_Component<dirX>();}
       template<class Element, unsigned int iterLength>
       inline typename Element::TYPE 
        evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
           { 
             return element.edge_(dirY);
           }
   };
#if 1
// ----------------------------------------------------------------------------
template <class A, class B>
class CrossProduct2D : public BasisFunctionExpr<CrossProduct2D<A,B> >
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
        enum {hasV = _MAX_(A::hasV,B::hasV),
              hasW = _MAX_(A::hasW,B::hasW),
              hasP = _MAX_(A::hasP,B::hasP)} ;
       CrossProduct2D(const A& a_, const B& b_) : a(a_), b(b_) {}
       template<class Element, unsigned int iterLength>
       inline typename Element::TYPE 
         eval(const Element& element, const int (&iterator)[iterLength]) const
           {
             return a.getA().eval(element,iterator) * b.getB().eval(element,iterator) - 
                    a.getB().eval(element,iterator) * b.getA().eval(element,iterator) ;
           }
   };
// ----------------------------------------------------------------------------
template <class A, class B>
class CrossProduct2D_2;
// ----------------------------------------------------------------------------
template <class A1, class A2, class B>
class CrossProduct2D_2<Integrand_Vector<A1,A2>,B> : public BasisFunctionVectorExpr<CrossProduct2D_2<Integrand_Vector<A1,A2>,B>,2>
   {
     private:
#if REFERENZEN
       const Integrand_Vector<A1,A2>& a;
       const B& b;
#else
       Integrand_Vector<A1,A2> a;
       B b;
#endif
     public:
        enum {hasV = _MAX_((Integrand_Vector<A1,A2>::hasV),B::hasV),
              hasW = _MAX_((Integrand_Vector<A1,A2>::hasW),B::hasW),
              hasP = _MAX_((Integrand_Vector<A1,A2>::hasP),B::hasP)} ;
       CrossProduct2D_2(const Integrand_Vector<A1,A2>& a_, const B& b_) : a(a_), b(b_) {}
       inline Mult_<A2,B> getHeadOfVector() const { return a.getB()*b; }
       template<class Element, unsigned int iterLength>
       inline typename Element::TYPE 
        evalLastComponent(const Element& element, const int (&iterator)[iterLength]) const
          { 
             return -a.getA().eval(element,iterator) * b.eval(element,iterator);
          }
   };
// ----------------------------------------------------------------------------
template <class A, class B>
inline CrossProduct2D<A,B>
operator ^ (const BasisFunctionVectorExpr<A,2>& a, const BasisFunctionVectorExpr<B,2>& b) 
   {
     return CrossProduct2D<A,B>(a,b);
   }
// ----------------------------------------------------------------------------
template <class A, class B>
inline BasisVectorMult_<A,B>
operator ^ (const BasisFunctionExpr<A>& a, const BasisFunctionVectorExpr<B,2>& b) 
   {
     return BasisVectorMult_<A,B>(a,b);
   }
// ----------------------------------------------------------------------------
template <class A, class B>
inline CrossProduct2D_2<A,B>
operator ^ (const BasisFunctionVectorExpr<A,2>& a, const BasisFunctionExpr<B>& b) 
   {
     return CrossProduct2D_2<A,B>(a,b);
   }
#endif
// ----------------------------------------------------------------------------
/* Concluding comments
 * At this point there is still more to introduce. A compile time data structure
 * which holds tensor data. Therefore we would need an additional operator to be 
 * overloaded, whicht has less priority than operator '&'. 
*/
//-----------------------------------------------------------------------------
/* A matrix is a concatenation of two vectors or an vector and a matrix. 
 * The result of the operator '|' is from left to right. Thus, the template 
 * type we achieve is something of the type 
 * MATRIX< MATRIX < MATRIX < ..., Integrand_Vector <..> >, Integrand_Vector <...> >, Integrand_Vector <...> >
*/
template <class A, int mLength, int mSize>
struct BasisFunctionMatrixExpr
   {
      enum {length = mLength, size = mSize};
   };
//-----------------------------------------------------------------------------
template <class A, class B>
class Integrand_Matrix : public BasisFunctionMatrixExpr<Integrand_Matrix<A,B>,A::length,A::size+B::size>
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
      Integrand_Matrix(const A& a_, const B& b_) : a(a_), b(b_) {}
      enum {hasV = _MAX_(A::hasV,B::hasV),
            hasW = _MAX_(A::hasW,B::hasW), 
            hasP = _MAX_(A::hasP,B::hasP)} ;
      inline void print(std::ostream& os);
   };
//-----------------------------------------------------------------------------
template <class A, class B, int l>
Integrand_Matrix<A,B> 
 operator | (const BasisFunctionVectorExpr<A,l>& a, const BasisFunctionVectorExpr<B,l>& b)
   {
      return Integrand_Matrix<A,B>(a,b);
   }
#if 0
template <class A, class B>
class MATRIX ;
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
class MATRIX<Integrand_Vector<A1,A2>,Integrand_Vector<B1,B2> > : public BasisFunctionVectorExpr<MATRIX<Integrand_Vector<A1,A2>, Integrand_Vector<B1,B2> > >
   {
     private:
      const Integrand_Vector<A1,A2>& a;
      const Integrand_Vector<B1,B2>& b;
     public:
      MATRIX(const Integrand_Vector<A1,A2>& a_, const Integrand_Vector<B1,B2>& b_) : a(a_), b(b_) {}
     enum {hasV = _MAX_(Integrand_Vector<A1,A2>::hasV,Integrand_Vector<B1,B2>::hasV),
           hasW = _MAX_(Integrand_Vector<A1,A2>::hasW,Integrand_Vector<B1,B2>::hasW), 
           hasP = _MAX_(Integrand_Vector<A1,A2>::hasP,Integrand_Vector<B1,B2>::hasP)} ;
     inline static void 
      print(std::ostream& os);
   };
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
class MATRIX<MATRIX<A1,A2>,Integrand_Vector<B1,B2> > : public BasisFunctionVectorExpr<MATRIX<Integrand_Vector<A1,A2>, Integrand_Vector<B1,B2> > >
   {
     private:
      const MATRIX<A1,A2>& a;
      const Integrand_Vector<B1,B2>& b;
     public:
      MATRIX(const MATRIX<A1,A2>& a_, const Integrand_Vector<B1,B2>& b_) : a(a_), b(b_) {}
      enum {hasV = _MAX_(MATRIX<A1,A2>::hasV,Integrand_Vector<B1,B2>::hasV),
           hasW = _MAX_(MATRIX<A1,A2>::hasW,Integrand_Vector<B1,B2>::hasW), 
           hasP = _MAX_(MATRIX<A1,A2>::hasP,Integrand_Vector<B1,B2>::hasP)} ;
      inline static void 
       print(std::ostream& os);
   };
//-----------------------------------------------------------------------------
template <class A, class B>
MATRIX<A,B> 
 operator | (const BasisFunctionVectorExpr<A>& a, const BasisFunctionVectorExpr<B>& b)
   {
      return MATRIX<A,B>(static_cast<const A&>(a),static_cast<const B&>(b));
   }
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a vector with a one dimensional
 * value. The result is the vector, muliplied with the value in each component.
*/
template <class A, class B>
struct BasisFunctionExpr_Integrand_Vector_MULT : public BasisFunctionExpr<BasisFunctionExpr_Integrand_Vector_MULT<A,B> >
   {
      private: 
        const A& a;
        const B& b;
      public:
       BasisFunctionExpr_Integrand_Vector_MULT(const A& a_, const B& b_) : a(a_), b(b_) 
          {
             assert (A::length == B::length);
             assert (A::length == 2);
          }
       enum {hasV = _MAX_(A::hasV,B::hasV),
             hasW = _MAX_(A::hasW,B::hasW), 
             hasP = _MAX_(A::hasP,B::hasP)} ;
     inline static void 
      print(std::ostream& os)
        {
          os << " (";   
          a.print(os);
          os << ") * (";   
          b.print(os);
          os << ") ";   
        }
   }; 
//-----------------------------------------------------------------------------
template <class A, class B, class C>
struct BasisFunctionExpr_Integrand_Vector_MULT<Integrand_Vector<A,B>,C> : public BasisFunctionExpr<BasisFunctionExpr_Integrand_Vector_MULT<Integrand_Vector<A,B>,C> >
   {
      private: 
        const Integrand_Vector<A,B>& a;
        const B& b;
      public:
       BasisFunctionExpr_Integrand_Vector_MULT(const Integrand_Vector<A,B>& a_, const B& b_) : a(a_), b(b_) 
          {
             assert (Integrand_Vector<A,B>::length == B::length);
             assert (Integrand_Vector<A,B>::length == B::length);
          }
       enum {hasV = _MAX_(Integrand_Vector<A,B>::hasV,B::hasV),
             hasW = _MAX_(Integrand_Vector<A,B>::hasW,B::hasW), 
             hasP = _MAX_(Integrand_Vector<A,B>::hasP,B::hasP)} ;
     inline static void 
      print(std::ostream& os)
        {
          os << " (";   
          a.print(os);
          os << ") * (";   
          b.print(os);
          os << ") ";   
        }
   }; 
//-----------------------------------------------------------------------------
template <class A, class V1, class V2>
typename BasisFunctionExpr_Integrand_Vector_MULT< A, Integrand_Vector<V1,V2> >
 operator * (const BasisFunctionExpr<A>& a, const Integrand_Vector<V1,V2>& b);
//-----------------------------------------------------------------------------
template <class A, class V1, class V2>
typename BasisFunctionExpr_Integrand_Vector_MULT< A, Integrand_Vector<V1,V2> > 
 operator * (const Integrand_Vector<V1,V2>& a, const BasisFunctionExpr<A>& b);
//-----------------------------------------------------------------------------
/* This struct describes the addition of two vectors. The result is the vector, 
 * that holds the sum of the vector components. 
*/
template <class A, class B>
struct CT_Integrand_Vector_Add;  
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename CT_Integrand_Vector_Add<Integrand_Vector<A1,A2>,Integrand_Vector<B1,B2> >
 operator+(const Integrand_Vector<A1,A2>& a, const Integrand_Vector<B1,B2>& b);
//-----------------------------------------------------------------------------
template <class A, class B>
struct CT_Integrand_Vector_Sub;  
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename CT_Integrand_Vector_Sub<Integrand_Vector<A1,A2>,Integrand_Vector<B1,B2> >
 operator-(const Integrand_Vector<A1,A2>& a, const Integrand_Vector<B1,B2>& b);
//-----------------------------------------------------------------------------
template <class A>
struct CT_Integrand_Vector_Min; 
//-----------------------------------------------------------------------------
template <class A, class B>
typename CT_Integrand_Vector_Min<Integrand_Vector<A,B> >
 operator-(const Integrand_Vector<A,B>& a);
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of two vectors. The result is  
 * the sum that holds the products of the vectors components. 
*/
template <class A, class B>
struct Integrand_Vector_Integrand_Vector_MULT; 
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename Integrand_Vector_Integrand_Vector_MULT<Integrand_Vector<A1,A2>,Integrand_Vector<B1,B2> >
 operator*(const Integrand_Vector<A1,A2>& a, const Integrand_Vector<B1,B2>& b);
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a vector with a matrix!  
 * The result is a vector of the line sums . 
*/
template <class A, class B>
struct Integrand_Vector_MATRIX_MULT;
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename Integrand_Vector_MATRIX_MULT<Integrand_Vector<A1,A2>,MATRIX<B1,B2> >
 operator*(const Integrand_Vector<A1,A2>& a, const MATRIX<B1,B2>& b);
//-----------------------------------------------------------------------------
/* This struct describes the multiplication of a matrix with a vector!  
 * The result is a vector of the line sums . 
*/
template <class A, class B>
struct MATRIX_Integrand_Vector_MULT;
//-----------------------------------------------------------------------------
template <class A1, class A2, class B1, class B2>
typename MATRIX_Integrand_Vector_MULT<MATRIX<A1,A2>,Integrand_Vector<B1,B2> >
 operator*(const MATRIX<A1,A2>& a, const Integrand_Vector<B1,B2>& b); 
#endif

//==============================================================================

