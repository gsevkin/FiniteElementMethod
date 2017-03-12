#define ET_inline inline 
//----------------------------------------------------------------------------
template <class Expression>
float
 FunctionExpr<Expression> :: operator() (float* (&d),int dx, int dy, int dz, float a) const
   {
       return static_cast<const Expression&>(*this).template eval<float,float>(d,dx,dy,dz);
   }
//----------------------------------------------------------------------------
template <class Expression>
double
 FunctionExpr<Expression> :: operator() (float* (&d),int dx, int dy, int dz, double a) const
   {
       return static_cast<const Expression&>(*this).template eval<double,float>(d,dx,dy,dz);
   }
//----------------------------------------------------------------------------
template <class Expression>
double
 FunctionExpr<Expression> :: operator() (double* (&d),int dx, int dy, int dz, double a) const
   {
       return static_cast<const Expression&>(*this).template eval<double,double>(d,dx,dy,dz);
   }
//----------------------------------------------------------------------------
template <class Expression>
std::complex<double>
 FunctionExpr<Expression> :: operator() ( float* (&d) ,int dx, int dy, int dz, std::complex<double> a) const
   {
     return static_cast<const Expression&>(*this).template eval<std::complex<double>,float>(d,dx,dy,dz);
   }
//----------------------------------------------------------------------------
template <class Expression>
std::complex<double>
 FunctionExpr<Expression> :: operator() ( double* (&d) ,int dx, int dy, int dz, std::complex<double> a) const
   {
     return static_cast<const Expression&>(*this).template eval<std::complex<double>,double>(d,dx,dy,dz);
   }
//----------------------------------------------------------------------------
template <typename eTYPE>
template <typename RTYPE, typename DTYPE>
RTYPE
 MonomX<eTYPE>::eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const
   {
     if (ddy + ddz == 0)
        {
          switch (ddx)
            {
               case 0:
                 return (RTYPE)std::pow(data[0],exponent);
                 break;
               case 1:
                 return (RTYPE)std::pow(data[0],exponent-1) * static_cast<RTYPE>(exponent);
                 break;
#ifdef SECOND_DERIVATIVES
               case 2:
                 return (RTYPE)std::pow(data[0],exponent-2) * static_cast<RTYPE>(exponent*(exponent-1));
                 break;
#endif
               default:
                 assert (!"This should never happen");
                 break;
             }
        }
     return 0;
   }
//----------------------------------------------------------------------------
template <typename eTYPE>
template <typename RTYPE, typename DTYPE>
RTYPE
 MonomY<eTYPE>::eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const
   {
     if (ddx + ddz == 0)
        {
          switch (ddy)
            {
               case 0:
                 return (RTYPE)std::pow(data[1],exponent);
                 break;
               case 1:
                 return (RTYPE)std::pow(data[1],exponent-1) * static_cast<RTYPE>(exponent);
                 break;
#ifdef SECOND_DERIVATIVES
               case 2:
                 return (RTYPE)std::pow(data[1],exponent-2) * static_cast<RTYPE>(exponent * (exponent-1));
                 break;
#endif
               default:
                 assert (!"This should never happen");
                 break;
             }
        }
     return 0;
   }
//----------------------------------------------------------------------------
template <typename eTYPE>
template <typename RTYPE, typename DTYPE>
RTYPE
 MonomZ<eTYPE>::eval(DTYPE*(&data), const int ddx, const int ddy, const int ddz) const
   {
     if (ddx + ddy == 0)
        {
          switch (ddz)
            {
               case 0:
                 return (RTYPE)std::pow(data[2],exponent);
                 break;
               case 1:
                 return (RTYPE)std::pow(data[2],exponent-1) * static_cast<RTYPE>(exponent);
                 break;
#ifdef SECOND_DERIVATIVES
               case 2:
                 return (RTYPE)std::pow(data[2],exponent-2) * static_cast<RTYPE>(exponent * (exponent-1));
                 break;
#endif
               default:
                 assert (!"This should never happen");
                 break;
             }
        }
     return 0;
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionDDX<A>
 d_dx (const FunctionExpr<A>& a)
   {
     return FunctionDDX<A>(a);
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionDDY<A>
 d_dy (const FunctionExpr<A>& a)
   {
     return FunctionDDY<A>(a);
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionDDZ<A>
 d_dz (const FunctionExpr<A>& a)
   {
     return FunctionDDZ<A>(a);
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionMinus<A>
 operator- (const FunctionExpr<A>& a)
   {
     return FunctionMinus<A>( a );
   }
//----------------------------------------------------------------------------
template <class A, class B>
ET_inline FunctionAdd<A,B>
 operator+ (const FunctionExpr<A>& a, const FunctionExpr<B>& b)
   {
     return FunctionAdd<A,B>( a, b );
   }
//----------------------------------------------------------------------------
template <class A, class B>
ET_inline FunctionSub<A,B>
 operator- (const FunctionExpr<A>& a, const FunctionExpr<B>& b)
   {
     return FunctionSub<A,B>( a, b );
   }
//----------------------------------------------------------------------------
template <class A, class B>
ET_inline FunctionMult<A,B>
 operator* (const FunctionExpr<A>& a, const FunctionExpr<B>& b)
   {
     return FunctionMult<A,B>( a, b );
   }
//----------------------------------------------------------------------------
template <class A, class B>
ET_inline FunctionDiv<A,B>
 operator/ (const FunctionExpr<A>& a, const FunctionExpr<B>& b)
   {
     return FunctionDiv<A,B>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCAdd<A,int>
 operator+ (const FunctionExpr<A>& a, int b)
   {
     return FunctionCAdd<A,int>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCAdd<A,float>
 operator+ (const FunctionExpr<A>& a, float b)
   {
     return FunctionCAdd<A,float>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCAdd<A,double>
 operator+ (const FunctionExpr<A>& a, double b)
   {
     return FunctionCAdd<A,double>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCAdd<A,std::complex<double> > 
 operator+ (const FunctionExpr<A>& a, std::complex<double> b)
   {
     return FunctionCAdd<A,std::complex<double> >( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCAdd<A,int>
 operator+ (int b, const FunctionExpr<A>& a)
   {
     return FunctionCAdd<A,int>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCAdd<A,float>
 operator+ (float b, const FunctionExpr<A>& a)
   {
     return FunctionCAdd<A,float>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCAdd<A,double>
 operator+ (double b, const FunctionExpr<A>& a)
   {
     return FunctionCAdd<A,double>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCAdd<A,std::complex<double> > 
 operator+ (std::complex<double> b, const FunctionExpr<A>& a)
   {
     return FunctionCAdd<A,std::complex<double> >( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionECSub<A,int>
 operator- (const FunctionExpr<A>& a, int b)
   {
     return FunctionECSub<A,int>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionECSub<A,float>
 operator- (const FunctionExpr<A>& a, float b)
   {
     return FunctionECSub<A,float>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionECSub<A,double>
 operator- (const FunctionExpr<A>& a, double b)
   {
     return FunctionECSub<A,double>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionECSub<A,std::complex<double> > 
 operator- (const FunctionExpr<A>& a, std::complex<double> b)
   {
     return FunctionECSub<A,std::complex<double> >( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCESub<int,A>
 operator- (int b, const FunctionExpr<A>& a)
   {
     return FunctionCESub<int,A>(b,a);
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCESub<float,A> 
 operator- (float b, const FunctionExpr<A>& a)
   {
     return FunctionCESub<float,A>(b,a);
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCESub<double,A>
 operator- (double b, const FunctionExpr<A>& a)
   {
     return FunctionCESub<double,A>(b,a);
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCESub<std::complex<double>,A> 
 operator- (std::complex<double> b, const FunctionExpr<A>& a)
   {
     return FunctionCESub<std::complex<double>,A>(b, a);
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCMult<A,int>
 operator* (const FunctionExpr<A>& a, int b)
   {
     return FunctionCMult<A,int>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCMult<A,float>
 operator* (const FunctionExpr<A>& a, float b)
   {
     return FunctionCMult<A,float>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCMult<A,double>
 operator* (const FunctionExpr<A>& a, double b)
   {
     return FunctionCMult<A,double>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCMult<A,std::complex<double> >
 operator* (const FunctionExpr<A>& a, std::complex<double> b)
   {
     return FunctionCMult<A,std::complex<double> >( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCMult<A,int>
 operator* (int b, const FunctionExpr<A>& a)
   {
     return FunctionCMult<A,int>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCMult<A,float>
 operator* (float b, const FunctionExpr<A>& a)
   {
     return FunctionCMult<A,float>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCMult<A,double>
 operator* (double b, const FunctionExpr<A>& a)
   {
     return FunctionCMult<A,double>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionCMult<A,std::complex<double> >
 operator* (std::complex<double> b, const FunctionExpr<A>& a)
   {
     return FunctionCMult<A,std::complex<double> >( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionECDiv<A,int>
 operator/ (const FunctionExpr<A>& a, int b)
   {
     return FunctionECDiv<A,int>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionECDiv<A,float>
 operator/ (const FunctionExpr<A>& a, float b)
   {
     return FunctionECDiv<A,float>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionECDiv<A,double>
 operator/ (const FunctionExpr<A>& a, double b)
   {
     return FunctionECDiv<A,double>( a, b );
   }
//----------------------------------------------------------------------------
template <class A>
ET_inline FunctionECDiv<A,std::complex<double> >
 operator/ (const FunctionExpr<A>& a, std::complex<double> b)
   {
     return FunctionECDiv<A,std::complex<double> >( a, b );
   }
//-----------------------------------------------------------------------------
template<class A>
ET_inline FunctionExp<A>
 Exp(const FunctionExpr<A>& a)
   {
     return FunctionExp<A>(a);
   }
//-----------------------------------------------------------------------------
template<class A>
ET_inline FunctionSin<A>
 Sin(const FunctionExpr<A>& a)
   {
     return FunctionSin<A>(a);
   }
//-----------------------------------------------------------------------------
template<class A>
ET_inline FunctionCos<A>
 Cos(const FunctionExpr<A>& a)
   {
     return FunctionCos<A>(a);
   }
//-----------------------------------------------------------------------------
template<class A>
ET_inline FunctionSqrt<A>
 Sqrt(const FunctionExpr<A>& a)
   {
     return FunctionSqrt<A>(a);
   }
//-----------------------------------------------------------------------------
template<class A, typename Type>
ET_inline FunctionPow<A,Type>
 Pow(const FunctionExpr<A>& a, Type exponent)
   {
     return FunctionPow<A,Type>(a,exponent);
   }
//----------------------------------------------------------------------------
template <typename eTYPE>
ET_inline MonomX<eTYPE>
 X_(eTYPE exponent)
   {
     return MonomX<eTYPE>(exponent);
   }
//----------------------------------------------------------------------------
template <typename eTYPE>
ET_inline MonomY<eTYPE>
 Y_(eTYPE exponent)
   {
     return MonomY<eTYPE>(exponent);
   }
//----------------------------------------------------------------------------
template <typename eTYPE>
ET_inline MonomZ<eTYPE>
 Z_(eTYPE exponent)
   {
     return MonomZ<eTYPE>(exponent);
   }
//----------------------------------------------------------------------------
ET_inline MonomX<int>
 X_()
   {
     return MonomX<int>(1);
   }
//----------------------------------------------------------------------------
ET_inline MonomY<int>
 Y_()
   {
     return MonomY<int>(1);
   }
//----------------------------------------------------------------------------
ET_inline MonomZ<int>
 Z_()
   {
     return MonomZ<int>(1);
   }
//==============================================================================
   
