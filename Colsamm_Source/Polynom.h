//==============================================================================
//
//  $Id: Polynom.h,v 1.25 2006/07/25 13:52:27 jochen Exp $
//
//==============================================================================
template <typename Type1, typename Type2>
struct Type_Manage;
//----------------------------------------------------------------------------
/* Dieser Teil wird ben"otigt um die Typumwandlungen (auch die falschen), 
   welche durch die virtuellen Funktionen enstehen, abzudecken. Um sicherzustellen,
   dass keine flaschen Umwandlungen stattfinden, wird hier ein assert eingebaut.
*/
template <>
struct Type_Manage<int,int>{
  static int manage (const int v) {return v;}
};
template <>
struct Type_Manage<int,float>{
  static float manage (const int v) {return v;}
};
template <>
struct Type_Manage<float,int>{
  static int manage (const float v) {assert(false) ; return 0;}
};
template <>
struct Type_Manage<int,double>{
  static double manage (const int v) {return v;}
};
template <>
struct Type_Manage<double,int>{
  static int manage (const double v) {assert(false) ; return 0;}
};

template <typename Type>
struct Type_Manage<int,std::complex<Type> >{
  static std::complex<Type> manage (const int v) {return Type(v);}
};

template <typename Type>
struct Type_Manage<std::complex<Type> ,int>{
  static int manage (const std::complex<Type>& v) {assert(false) ; return 0;}
};

template <>
struct Type_Manage<float,float>{
  static float manage (const float v) {return v;}
};
template <>
struct Type_Manage<float,double>{
  static double manage (const float v) {return v;}
};
template <>
struct Type_Manage<double,float>{
  static float manage (const double v) { assert(false) ; return 0.;}
};
template <typename Type>
struct Type_Manage<float,std::complex<Type> >{
  static std::complex<Type> manage (const float v) {return std::complex<Type>(v);}
};
template <typename Type>
struct Type_Manage<std::complex<Type>,float>{
  static float manage (const std::complex<Type> v) { assert(false) ; return 0.;}
};
template <>
struct Type_Manage<double,double>{
  static double manage (const double v) {return v;}
};
template <typename Type>
struct Type_Manage<std::complex<Type>,double>{
  static double manage (const std::complex<Type> v) { assert(false) ; return 0.;}
};
template <typename Type>
struct Type_Manage<double,std::complex<Type> >{
  static std::complex<Type> manage (const double v) {return Type(v);}
};
template <>
struct Type_Manage<std::complex<double>,std::complex<double> >{
  static const std::complex<double>& manage (const std::complex<double>& v) {return v;}
};
//----------------------------------------------------------------------------
struct Base {
  virtual float operator() (float* (&d), int dx = 0, int dy = 0, int dz = 0, float a = 0.) const  =  0;
  virtual double operator() (float* (&d), int dx = 0, int dy = 0, int dz = 0, double a = 0.) const  =  0;
  virtual double operator() (double* (&d), int dx = 0, int dy = 0, int dz = 0, double a = 0.) const  =  0;
  virtual std::complex<double> operator()(float* (&d), int dx = 0, int dy = 0, int dz = 0, std::complex<double> a = 0.) const = 0;
  virtual std::complex<double> operator()(double* (&d), int dx = 0, int dy = 0, int dz = 0, std::complex<double> a = 0.) const = 0;
  virtual ~Base(){};
};
//----------------------------------------------------------------------------
template <class Expression>
struct FunctionExpr : public Base
   {
     operator Expression() const { return *static_cast<const Expression*>(this); }
     virtual float operator() (float* (&d), int dx = 0, int dy = 0, int dz = 0, float a = 0.) const ;
     virtual double operator() (float* (&d), int dx = 0, int dy = 0, int dz = 0, double a = 0.) const ;
     virtual double operator() (double* (&d), int dx = 0, int dy = 0, int dz = 0, double a = 0.) const ;
     virtual std::complex<double> operator()(float* (&d), int dx = 0, int dy = 0, int dz = 0, std::complex<double> a = 0.) const;
     virtual std::complex<double> operator()(double* (&d), int dx = 0, int dy = 0, int dz = 0, std::complex<double> a = 0.) const;
   };
//----------------------------------------------------------------------------
class Null_ : public FunctionExpr<Null_>
   {
      public:  
     template <typename RTYPE, typename DTYPE>
     RTYPE 
      eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return (RTYPE)(0);
        } 
   };
//----------------------------------------------------------------------------
template <typename TYPE>
class Const_ : public FunctionExpr<Const_<TYPE> >
   {
      private:
        const TYPE const_;
      public:  
        Const_(const TYPE& t) : const_(t) {};
     template <typename RTYPE, typename DTYPE>
     RTYPE 
      eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          if (ddx+ddy+ddz == 0)
             {
               return (RTYPE)(const_);
             }
          else
             {
               return (RTYPE)(0);
             } 
        } 
   };
//----------------------------------------------------------------------------
template <typename TYPE>
inline Const_<TYPE> c_(const TYPE& t)
   {
      return Const_<TYPE>(t);
   }

//----------------------------------------------------------------------------
template <typename eTYPE>
class MonomX : public FunctionExpr <MonomX<eTYPE> >
  {
    private:
     eTYPE exponent;
    public:
     MonomX (eTYPE exponent_ = 1) : exponent(exponent_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
      eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const; 
  };
//----------------------------------------------------------------------------
template <typename eTYPE>
class MonomY : public FunctionExpr <MonomY<eTYPE> >
  {
    private:
     eTYPE exponent;
    public:
     MonomY (eTYPE exponent_ = 1) : exponent(exponent_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
      eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const; 
  };
//----------------------------------------------------------------------------
template <typename eTYPE>
class MonomZ : public FunctionExpr <MonomZ<eTYPE> >
  {
    private:
     eTYPE exponent;
    public:
     MonomZ (eTYPE exponent_ = 1) : exponent(exponent_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
      eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const; 
  };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionDDX : public FunctionExpr<FunctionDDX<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionDDX(const Expr1& expr1_): expr1(expr1_) {}

     template <typename RTYPE, typename DTYPE>
     RTYPE 
      eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx+1,ddy,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionDDY : public FunctionExpr<FunctionDDY<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionDDY(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy+1,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionDDZ : public FunctionExpr<FunctionDDZ<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionDDZ(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz+1);
        }
   };
#ifdef SECOND_DERIVATIVES
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DX2 : public FunctionExpr<FunctionD2DX2<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DX2(const Expr1& expr1_): expr1(expr1_) {}

     template <typename RTYPE, typename DTYPE>
     RTYPE 
      eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx+2,ddy,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DXY : public FunctionExpr<FunctionD2DXY<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DXY(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx+1,ddy+1,ddz);
        }
   };
#if 0
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DYX : public FunctionExpr<FunctionD2DYX<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DYX(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx+1,ddy+1,ddz);
        }
   };
#endif
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DXZ : public FunctionExpr<FunctionD2DXZ<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DXZ(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx+1,ddy,ddz+1);
        }
   };
#if 0
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DZX : public FunctionExpr<FunctionD2DZX<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DZX(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx+1,ddy,ddz+1);
        }
   };
#endif
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DY2 : public FunctionExpr<FunctionD2DY2<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DY2(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy+2,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DYZ : public FunctionExpr<FunctionD2DYZ<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DYZ(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy+1,ddz+1);
        }
   };
#if 0
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DZY : public FunctionExpr<FunctionD2DZY<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DZY(const Expr1& expr1_): expr1(expr1_) {}
      template <typename DTYPE>
      DTYPE
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.eval(data,ddx,ddy+1,ddz+1);
        }
   };
#endif
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionD2DZ2 : public FunctionExpr<FunctionD2DZ2<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionD2DZ2(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz+2);
        }
   };
#endif
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionMinus : public FunctionExpr<FunctionMinus<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
      FunctionMinus(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return -expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, typename TYPE>
class FunctionCAdd : public FunctionExpr<FunctionCAdd<Expr1,TYPE> >
   {
     private:
      Expr1 expr1;
      TYPE constant;
     public:
      FunctionCAdd(const Expr1& expr1_, TYPE constant_)
        : expr1(expr1_) , constant(constant_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const
        {
          if ( 0 < ddx+ddy+ddz )
             {
               return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
             }
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) + Type_Manage<TYPE,RTYPE>::manage(constant);
        }
   };
//----------------------------------------------------------------------------
template <typename TYPE, class Expr1>
class FunctionCESub : public FunctionExpr<FunctionCESub<TYPE,Expr1> >
   {
     private:
      Expr1 expr1;
      TYPE constant;
     public:
      FunctionCESub(TYPE constant_, const Expr1& expr1_)
        : expr1(expr1_) , constant(constant_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const
        {
          if ( 0 < ddx+ddy+ddz )
             {
               return -expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
             }
          return Type_Manage<TYPE,RTYPE>::manage(constant) - expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, typename TYPE>
class FunctionECSub : public FunctionExpr<FunctionECSub<Expr1,TYPE> >
   {
     private:
      Expr1 expr1;
      TYPE constant;
     public:
      FunctionECSub(const Expr1& expr1_, TYPE constant_)
        : expr1(expr1_) , constant(constant_) {}
     template <typename RTYPE, typename DTYPE>
      RTYPE
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const
        {
          if ( 0 < ddx+ddy+ddz )
             {
               return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
             }
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) - Type_Manage<TYPE,RTYPE>::manage(constant);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, typename TYPE>
class FunctionCMult : public FunctionExpr<FunctionCMult<Expr1,TYPE> >
   {
     private: 
      Expr1 expr1;
      TYPE constant;
     public:
      FunctionCMult(const Expr1& expr1_, TYPE constant_)
        : expr1(expr1_) , constant(constant_) {}
     template <typename RTYPE, typename DTYPE>
      RTYPE
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) * Type_Manage<TYPE,RTYPE>::manage(constant);
        }
   };
//----------------------------------------------------------------------------
template <typename TYPE, class Expr1>
class FunctionCEDiv : public FunctionExpr<FunctionCEDiv<TYPE,Expr1> >
   {
     private:
      Expr1 expr1;
      TYPE constant;
     public:
      FunctionCEDiv( TYPE constant_, const Expr1& expr1_)
        : expr1(expr1_) , constant(constant_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return Type_Manage<TYPE,RTYPE>::manage(constant) / expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, typename TYPE>
class FunctionECDiv : public FunctionExpr<FunctionECDiv<Expr1,TYPE> >
   {
     private:
      Expr1 expr1;
      TYPE constant;
     public:
      FunctionECDiv(const Expr1& expr1_, TYPE constant_)
        : expr1(expr1_) , constant(constant_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) / Type_Manage<TYPE,RTYPE>::manage(constant);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, class Expr2>
class FunctionAdd : public FunctionExpr<FunctionAdd<Expr1,Expr2> >
   {
     private:
      Expr1 expr1;
      Expr2 expr2;
     public:
      FunctionAdd(const Expr1& expr1_, const Expr2& expr2_)
        : expr1(expr1_), expr2(expr2_) {}

     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) + expr2.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, class Expr2>
class FunctionSub : public FunctionExpr<FunctionSub<Expr1,Expr2> >
   {
     private:
      Expr1 expr1;
      Expr2 expr2;
     public:
      FunctionSub(const Expr1& expr1_, const Expr2& expr2_)
        : expr1(expr1_), expr2(expr2_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) - expr2.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, class Expr2>
class FunctionMult : public FunctionExpr<FunctionMult<Expr1,Expr2> >
   {
     private:
      Expr1 expr1;
      Expr2 expr2;
     public:
      FunctionMult(const Expr1& expr1_, const Expr2& expr2_)
        : expr1(expr1_), expr2(expr2_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
#ifdef SECOND_DERIVATEIVES
          int dx1 = ((0 < ddx)?1:0), dy1 = (((0<ddy)&&(ddx==0))?1:0), 
              dz1 = (((0<ddz)&&(ddx==0)&&(ddy==0))?1:0);
#endif
          switch( ddx+ddy+ddz )
             {
               case 0:
                 return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) * expr2.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
                 break;
               case 1:
                 return expr1.template eval<RTYPE,DTYPE>(data,0,0,0) * expr2.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) +
                        expr2.template eval<RTYPE,DTYPE>(data,0,0,0) * expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) ;
                 break;
#ifdef SECOND_DERIVATEIVES
               case 2:
                  return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) * expr2.template eval<RTYPE,DTYPE>(data,0,0,0) +
                         expr1.template eval<RTYPE,DTYPE>(data,ddx-dx1,ddy-dy1,ddz-dz1) * expr2.template eval<RTYPE,DTYPE>(data,dx1,dy1,dz1) +
                         expr1.template eval<RTYPE,DTYPE>(data,dx1,dy1,dz1) * expr2.template eval<RTYPE,DTYPE>(data,ddx-dx1,ddy-dy1,ddz-dz1) +
                         expr1.template eval<RTYPE,DTYPE>(data,0,0,0) * expr2.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) ;
#endif
                 break;
               default:
                 assert (!"This case should never happen");
                 break;
             }
          return 0.;
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, class Expr2>
class FunctionDiv : public FunctionExpr<FunctionDiv<Expr1,Expr2> >
   {
     private: 
      Expr1 expr1;
      Expr2 expr2;
     public:
      FunctionDiv(const Expr1& expr1_, const Expr2& expr2_)
        : expr1(expr1_), expr2(expr2_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
       eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          switch( ddx+ddy+ddz )
             {
               case 0:
                 return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) / expr2.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz);
                 break;
               case 1:
                 return ( expr2.template eval<RTYPE,DTYPE>(data,0,0,0) * expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) -
                          expr1.template eval<RTYPE,DTYPE>(data,0,0,0) * expr2.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz) ) /
                          std::pow(expr2.eval(data,ddx,ddy,ddz),2)  ;
                 break;
#ifdef SECOND_DERIVATEIVES
               case 2:
#endif
               default:
                 assert (!"This case should never happen");
                 break;
             }
          return 0.;
        }
   };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionExp : public FunctionExpr<FunctionExp<Expr1> >
   {
     private:
      Expr1 expr1;
     public:
     FunctionExp(const Expr1& expr1_): expr1(expr1_) {}
     template <typename RTYPE, typename DTYPE>
     RTYPE 
      eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const
        {
          if ( ddx+ddy+ddz == 0 )
             {
               return std::exp(expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz));
             }
          if ( ddx+ddy+ddz == 1 )
             {
               return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz)* std::exp(expr1.template eval<RTYPE,DTYPE>(data,0,0,0));
             }
          assert (!"This case should never happen");
          return 0.;
        }
   };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionSin : public FunctionExpr<FunctionSin<Expr1> >
   {
     private: 
      Expr1 expr1;
     public:
     FunctionSin(const Expr1& expr1_): expr1(expr1_) {}

     template <typename RTYPE, typename DTYPE>
     RTYPE 
     eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const
        {
          if ( ddx+ddy+ddz == 0 )
             {
               return std::sin(expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz));
             }
          if ( ddx+ddy+ddz == 1 )
             {
               return expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz)* std::cos(expr1.template eval<RTYPE,DTYPE>(data,0,0,0));
             }
          assert (!"This case should never happen");
          return 0.;
        }
   };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionCos : public FunctionExpr<FunctionCos<Expr1> >
   {
     private: 
      Expr1 expr1;
     public:
     FunctionCos(const Expr1& expr1_): expr1(expr1_) {}

     template <typename RTYPE, typename DTYPE>
     RTYPE 
     eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          if ( ddx+ddy+ddz == 0 )
             {
               return std::cos(expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz));
             }
          if ( ddx+ddy+ddz == 1 )
             {
               return - expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz)* std::sin(expr1.template eval<RTYPE,DTYPE>(data,0,0,0));
             }
          assert (!"This case should never happen");
          return 0.;
        }
   };
//----------------------------------------------------------------------------
template <class Expr1, typename Type>
class FunctionPow : public FunctionExpr<FunctionPow<Expr1,Type> >
   {
     private: 
      Expr1 expr1;
      Type exponent;
     public:
     FunctionPow(const Expr1& expr1_,Type expo): expr1(expr1_), exponent(expo) {}

     template <typename RTYPE, typename DTYPE>
     RTYPE 
     eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          if ( ddx+ddy+ddz == 0 )
             {
               return std::pow(expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz),(RTYPE)exponent);
             }
          if ( ddx+ddy+ddz == 1 )
             {
               return (RTYPE)(exponent)* 
                  expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz)* std::pow(expr1.template eval<RTYPE,DTYPE>(data,0,0,0),(RTYPE)(exponent-1));
             }
          assert (!"This case should never happen");
          return 0.;
        }
   };
//----------------------------------------------------------------------------
template <class Expr1>
class FunctionSqrt : public FunctionExpr<FunctionSqrt<Expr1> >
   {
     private: 
      Expr1 expr1;
     public:
     FunctionSqrt(const Expr1& expr1_): expr1(expr1_){}

     template <typename RTYPE, typename DTYPE>
     RTYPE 
     eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const 
        {
          if ( ddx+ddy+ddz == 0 )
             {
               return std::sqrt(expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz));
             }
          if ( ddx+ddy+ddz == 1 )
             {
               return 0.5* 
                  expr1.template eval<RTYPE,DTYPE>(data,ddx,ddy,ddz)/ std::sqrt(expr1.template eval<RTYPE,DTYPE>(data,0,0,0));
             }
          assert (!"This case should never happen");
          return 0.;
        }
   };
//----------------------------------------------------------------------------
template <class Expr>
struct VecFuncExpr
  {
    operator Expr()const {return *static_cast<const Expr*>(this); }
  };
//----------------------------------------------------------------------------
template <class Expr>
class FunctionExprVector1D : public VecFuncExpr<FunctionExprVector1D<Expr> >
   {
     private: 
      Expr expr;
     public:
      enum {Dimension = 1};
      FunctionExprVector1D(const Expr& expr_): expr(expr_){}

      const Expr getExpr() const
         {
           return expr;
         }
   };

//----------------------------------------------------------------------------
template <class Expr1, class Expr2>
class FunctionExprVector2D : public VecFuncExpr<FunctionExprVector2D<Expr1,Expr2> >
   {
     private: 
      Expr1 expr1;
      Expr2 expr2;
     public:
      enum {Dimension = 2};
      FunctionExprVector2D(const Expr1& expr1_, const Expr2& expr2_): expr1(expr1_), expr2(expr2_) {}

      const Expr1 getExpr1() const
         {
           return expr1;
         }
      const Expr2 getExpr2() const
         {
           return expr2;
         }
#if 0 // not needed ... 
      template <typename DTYPE>
      DTYPE eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz, SpaceDirection dir) const 
        {
          switch (dir)
             {
               case dirX:
                 return expr1.eval(data,ddx,ddy,ddz);
                 break;
               case dirY:
                 return expr2.eval(data,ddx,ddy,ddz);
                 break;
               default:
                 assert (!"This case should never happen");
                 break;
             }
          return 0.;
        }
#endif
   };

//----------------------------------------------------------------------------
template <class Expr1, class Expr2, class Expr3>
class FunctionExprVector3D : public VecFuncExpr<FunctionExprVector3D<Expr1,Expr2,Expr3> >
// : public FunctionExpr <FunctionExprVector3D<Expr1,Expr2,Expr3> >
   {
     private: 
      Expr1 expr1;
      Expr2 expr2;
      Expr3 expr3;
     public:
      enum {Dimension = 3};
      FunctionExprVector3D(const Expr1& expr1_, const Expr2& expr2_, const Expr3& expr3_)
         : expr1(expr1_), expr2(expr2_), expr3(expr3_) {}

      const Expr1 getExpr1() const
         {
           return expr1;
         }
      const Expr2 getExpr2() const
         {
           return expr2;
         }
      const Expr3 getExpr3() const
         {
           return expr3;
         }
   };
//----------------------------------------------------------------------------
template <class Expr>
inline FunctionExprVector1D<Expr> fVec(const FunctionExpr<Expr>& expr)
  {
     return FunctionExprVector1D<Expr>(expr);
  }
//----------------------------------------------------------------------------
template <class Expr1, class Expr2>
inline FunctionExprVector2D<Expr1,Expr2> fVec(const FunctionExpr<Expr1>& expr1,
                                              const FunctionExpr<Expr2>& expr2)
  {
     return FunctionExprVector2D<Expr1,Expr2>(expr1,
                                              expr2);
  }
//----------------------------------------------------------------------------
template <class Expr1, class Expr2, class Expr3>
inline FunctionExprVector3D<Expr1,Expr2,Expr3> fVec(const FunctionExpr<Expr1>& expr1,
                                                    const FunctionExpr<Expr2>& expr2,
                                                    const FunctionExpr<Expr3>& expr3)
  {
     return FunctionExprVector3D<Expr1,Expr2,Expr3>(expr1,
                                                    expr2,
                                                    expr3);
  }
//==============================================================================
