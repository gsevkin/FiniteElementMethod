//==============================================================================
//
//  $Id: VertexComputations.C,v 1.6 2006/07/14 10:51:30 jochen Exp $
//
//==============================================================================
//----------------------------------------------------------------------------
struct ConstantForVertexExpression 
  : public PolynomialExprForVertexExpr < ConstantForVertexExpression >
    {
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints (const Element& element, const unsigned int gaussLevel) 
         { 
           return (dx+dy+dz==0)?1.:0.;
         }
      template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints (const Element& element,int dx, int dy, int dz, const unsigned int gaussLevel) 
         { 
           return (dx+dy+dz==0)?1.:0.;
         }
    };
//----------------------------------------------------------------------------
template <int p>
struct PolynomialXForVertexExpr 
  : public PolynomialExprForVertexExpr < PolynomialXForVertexExpr <p> > 
    {
      static void print (std::ostream& os)
         {
           os << "PolynomialXForVertexExpr<" << p << "> "; 
         }
      template <int dx ,int dy , int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element, const unsigned int gaussLevel )  
         {
           typename Element::bTYPE& gaussianPoint = element.getGaussianPoints()[gaussLevel][dirX];
           return ( (dy+dz==0)? 
                      Colsamm_Internal_Functions::POW<p-dx>::pow(gaussianPoint)
                              *static_cast<typename Element::bTYPE>(Colsamm_Internal_Functions::F_Fac<p,p-dx>::erg):
                      0. );
         }
      template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element,int dx ,int dy , int dz, const unsigned int gaussLevel )  
         {
           typename Element::bTYPE& gaussianPoint = element.getGaussianPoints()[gaussLevel][dirX];
           if (dy+dz==0)
              {
                switch (dx)
                  {
                    case 0:
                      return std::pow(gaussianPoint,p);
                      break;
                    case 1:
                      return std::pow(gaussianPoint,p-1) * (typename Element::TYPE)(p);
                      break;
                    case 2:
                      return std::pow(gaussianPoint,p-2)* (typename Element::TYPE)(p * (p-1));
                      break;
                    default: 
                      assert (!"This should never happen\n");
                  }
              } 
           return 0;
         }
    } ;
//----------------------------------------------------------------------------
template < int p >
struct PolynomialYForVertexExpr 
  : public PolynomialExprForVertexExpr < PolynomialYForVertexExpr <p> > 
    { 
      static void print (std::ostream& os)
         {
           os << "PolynomialYForVertexExpr<" << p << "> "; 
         }
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element, const unsigned int gaussLevel ) 
         {
           typename Element::bTYPE& gaussianPoint = element.getGaussianPoints()[gaussLevel][dirY];
           return ( 
                   ( dx + dz == 0) ? 
                     Colsamm_Internal_Functions::POW<p-dy>::pow(gaussianPoint)*
                             static_cast<typename Element::bTYPE>(Colsamm_Internal_Functions::F_Fac<p,p-dy>::erg):
                     0. 
                  ) ;
         }
      template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element,int dx ,int dy , int dz, const unsigned int gaussLevel )  
         {
           typename Element::bTYPE& gaussianPoint = element.getGaussianPoints()[gaussLevel][dirY];
           if (dx+dz==0)
              {
                switch (dy)
                  {
                    case 0:
                      return std::pow(gaussianPoint,p);
                      break;
                    case 1:
                      return std::pow(gaussianPoint,p-1) * (typename Element::TYPE)(p);
                      break;
                    case 2:
                      return std::pow(gaussianPoint,p-2)* (typename Element::TYPE)(p * (p-1));
                      break;
                    default: 
                      assert (!"This should never happen\n");
                  }
              } 
           return 0;
         }
    } ;
//----------------------------------------------------------------------------
template < int p >
struct PolynomialZForVertexExpr 
  : public PolynomialExprForVertexExpr < PolynomialZForVertexExpr<p> > 
    {
      static void print (std::ostream& os)
         {
           os << "PolynomialZForVertexExpr<" << p << "> "; 
         }
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element, const unsigned int gaussLevel ) 
         {
           typename Element::bTYPE& gaussianPoint = element.getGaussianPoints()[gaussLevel][dirZ];
           return ( 
                   ( dx + dy == 0) ? 
                     Colsamm_Internal_Functions::POW<p-dz>::pow(gaussianPoint)* 
                         static_cast<typename Element::bTYPE>(Colsamm_Internal_Functions::F_Fac<p,p-dz>::erg):
                     0. 
                   ) ;
         }
      template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element,int dx ,int dy , int dz, const unsigned int gaussLevel )  
         {
           typename Element::bTYPE& gaussianPoint = element.getGaussianPoints()[gaussLevel][dirZ];
           if (dx+dy==0)
              {
                switch (dz)
                  {
                    case 0:
                      return std::pow(gaussianPoint,p);
                      break;
                    case 1:
                      return std::pow(gaussianPoint,p-1) * (typename Element::TYPE)(p);
                      break;
                    case 2:
                      return std::pow(gaussianPoint,p-2)* (typename Element::TYPE)(p * (p-1));
                      break;
                    default: 
                      assert (!"This should never happen\n");
                  }
              } 
           return 0;
         }
    } ;
//----------------------------------------------------------------------------
template < class A , class B >
 struct PolynomialExprForVertexExprAddition 
  : public PolynomialExprForVertexExpr < PolynomialExprForVertexExprAddition<A,B> > 
    {
      static void print (std::ostream& os)
         {
           os << "PolynomialExprForVertexExprAddition<";
           A::print(os); 
           os << ",";
           B::print(os); 
           os << "> ";
         }
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element, const unsigned int gaussLevel ) 
        {
          return A::template evaluateGaussianPoints<dx,dy,dz>(element,gaussLevel) + 
                 B::template evaluateGaussianPoints<dx,dy,dz>(element,gaussLevel) ;   
        }
      template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints (const Element& element,int dx, int dy, int dz, const unsigned int gaussLevel ) 
        {
          return A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) + 
                 B::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ;   
        }
    };
//----------------------------------------------------------------------------
template < class A, class B >
inline PolynomialExprForVertexExprAddition <A,B>
operator + ( const PolynomialExprForVertexExpr <A>& a , 
             const PolynomialExprForVertexExpr <B>& b ) 
   {
     return PolynomialExprForVertexExprAddition <A,B> () ;
   } 
//----------------------------------------------------------------------------
template < class A , class B > 
struct PolynomialExprForVertexExprSubtraction 
  : public PolynomialExprForVertexExpr < PolynomialExprForVertexExprSubtraction <A,B> > 
    { 
      static void print (std::ostream& os)
         {
           os << "PolynomialExprForVertexExprSubtraction<";
           A::print(os); 
           os << ",";
           B::print(os); 
           os << "> ";
         }
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element, const unsigned int gaussLevel ) 
         {
           return A:: template evaluateGaussianPoints <dx,dy,dz>(element,gaussLevel) - 
                  B:: template evaluateGaussianPoints <dx,dy,dz>(element,gaussLevel) ;
         }

      template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element,int dx, int dy, int dz, const unsigned int gaussLevel ) 
         {
           return A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) - 
                  B::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ;
         }
    } ;
//-------------------------------------------------------------------------------------
template < class A , class B > 
inline PolynomialExprForVertexExprSubtraction <A,B>
operator - ( const PolynomialExprForVertexExpr <A>& a , 
             const PolynomialExprForVertexExpr <B>& b ) 
   {
     return PolynomialExprForVertexExprSubtraction <A,B> () ;
   }
//-------------------------------------------------------------------------------------
template < class A >
struct PolynomialExprForVertexExprUnaryMinus 
  : public PolynomialExprForVertexExpr < PolynomialExprForVertexExprUnaryMinus <A> > 
    {
      static void print (std::ostream& os)
         {
           os << "PolynomialExprForVertexExprExprUnaryMinus<";
           A::print(os); 
           os << " > ";
         }
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element, const unsigned int gaussLevel ) 
         {
           return -A:: template evaluateGaussianPoints<dx,dy,dz> (element,gaussLevel) ;
         }
      template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints ( const Element& element,int dx, int dy, int dz, const unsigned int gaussLevel ) 
         {
           return -A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ;
         }
    } ;
//-------------------------------------------------------------------------------------
template < class A >
inline PolynomialExprForVertexExprUnaryMinus <A> 
operator - ( const PolynomialExprForVertexExpr <A>& a ) 
   {
     return PolynomialExprForVertexExprUnaryMinus <A> ( ) ;
   }
//-------------------------------------------------------------------------------------
template <class A, class B, int dx, int dy, int dz, class Element>
struct Multiplication<A,B,0,dx,dy,dz,Element>
   {
      static inline typename Element::TYPE 
      derivateDependentEvaluation( const Element& element, const unsigned int gaussLevel )
         {
           return (A::template evaluateGaussianPoints<dx,dy,dz>(element,gaussLevel) * 
                   B::template evaluateGaussianPoints<dx,dy,dz>(element,gaussLevel));
         }
     template <int k> 
      static inline typename Element::TYPE 
      derivateDependentEvaluationDirection( const Element& element, const unsigned int gaussLevel )
         {
           return (A::template evaluateGaussianPoints<dx,dy,dz,k>(element,gaussLevel) * 
                   B::template evaluateGaussianPoints<dx,dy,dz>(element,gaussLevel));
         }
   };
//-------------------------------------------------------------------------------------
template <class A, class B, int dx, int dy, int dz, class Element>
struct Multiplication<A,B,1,dx,dy,dz,Element>
   {
      static inline typename Element::TYPE 
      derivateDependentEvaluationDirection( const Element& element, const unsigned int gaussLevel )
         {
           return A::template evaluateGaussianPoints<dx,dy,dz>(element,gaussLevel) * 
                  B::template evaluateGaussianPoints<0,0,0>(element,gaussLevel) +
                  A::template evaluateGaussianPoints<0,0,0>(element,gaussLevel) * 
                  B::template evaluateGaussianPoints<dx,dy,dz>(element,gaussLevel);
         }
     template <int k> 
      static inline typename Element::TYPE 
      derivateDependentEvaluationDirection( const Element& element, const unsigned int gaussLevel )
         {
           return A::template evaluateGaussianPoints<dx,dy,dz,k>(element,gaussLevel) * 
                  B::template evaluateGaussianPoints<0,0,0>(element,gaussLevel) +
                  A::template evaluateGaussianPoints<0,0,0,k>(element,gaussLevel) * 
                  B::template evaluateGaussianPoints<dx,dy,dz>(element,gaussLevel);
         }
   };
//-------------------------------------------------------------------------------------
template<class A, class B>
struct PolynomialExprForVertexExpr_Mult 
  : public PolynomialExprForVertexExpr<PolynomialExprForVertexExpr_Mult<A,B> >
    {
      static void print (std::ostream& os)
         {
           os << "PolynomialExprForVertexExprExpr_Mult<";
           A::print(os); 
           os << ",";
           B::print(os); 
           os << "> ";
         }
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints( const Element& element, const unsigned int gaussLevel )
          {
            return Multiplication<A,B,dx+dy+dz,dx,dy,dz,Element>::
                                        derivateDependentEvaluation(element,gaussLevel); 
          }
      template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints( const Element& element, int dx, int dy, int dz, const unsigned int gaussLevel )
          {
            switch (dx+dy+dz)
              {
                case 0:
                  return A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) * 
                         B::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ; 
                  break;
                case 1:
                  return A::evaluateGaussianPoints(element,0,0,0,gaussLevel) * 
                         B::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) + 
                         B::evaluateGaussianPoints(element,0,0,0,gaussLevel) * 
                         A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ; 
                  break;
                default:
                  assert (! "This should never happen\n");
                  return 0;
                  break;
              }
          }
    };
//-------------------------------------------------------------------------------------
template<class A, class B>
inline PolynomialExprForVertexExpr_Mult <A,B>
operator * ( const PolynomialExprForVertexExpr <A>& a , 
             const PolynomialExprForVertexExpr <B>& b ) 
   {
     return PolynomialExprForVertexExpr_Mult <A,B> () ;
   }
//-------------------------------------------------------------------------------------
template<class A, class D>
class PolynomialExprForVertexExpr_CMult 
  : public PolynomialExprForVertexExpr < PolynomialExprForVertexExpr_CMult <A,D> > 
    {
       static const D constantValue ;
      public :
       PolynomialExprForVertexExpr_CMult ( const D& constantValue_ ) 
          { 
            constantValue = constantValue_; 
          } 
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
        evaluateGaussianPoints( const Element& element, const unsigned int gaussLevel )
          {
            return A:: template evaluateGaussianPoints <dx,dy,dz> (element,gaussLevel) * 
                   constantValue ; 
          }
      template <class Element>
      static inline typename Element::TYPE 
        evaluateGaussianPoints( const Element& element, int dx, int dy, int dz, const unsigned int gaussLevel )
          {
            return A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) * 
                   constantValue ; 
          }
    } ;
 template<class A, class D>
 const D 
  PolynomialExprForVertexExpr_CMult <A,D > :: constantValue ;
//-------------------------------------------------------------------------------------
template < class B > 
inline PolynomialExprForVertexExpr_CMult <B,float>
 operator * (float a, const PolynomialExprForVertexExpr <B> & b ) 
    {
      return PolynomialExprForVertexExpr_CMult <B,float> (a) ;
    }
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CMult <A,float>
 operator * (const PolynomialExprForVertexExpr <A>& a , float b) 
    {
      return PolynomialExprForVertexExpr_CMult <A,float> (b);
    }
//-------------------------------------------------------------------------------------
template < class B > 
inline PolynomialExprForVertexExpr_CMult <B,double>
 operator * (double a, const PolynomialExprForVertexExpr <B> & b ) 
    {
      return PolynomialExprForVertexExpr_CMult <B,double> (a) ;
    }
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CMult <A,double>
 operator * (const PolynomialExprForVertexExpr <A>& a , double b) 
    {
      return PolynomialExprForVertexExpr_CMult <A,double> (b);
    }
//-------------------------------------------------------------------------------------
template < class B > 
inline PolynomialExprForVertexExpr_CMult < B, std::complex<double> >
 operator * (const std::complex<double>& a, const PolynomialExprForVertexExpr <B>& b ) 
    {
      return PolynomialExprForVertexExpr_CMult < B, std::complex<double> > (a);
    }
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CMult < A, std::complex<double> >
 operator * ( const PolynomialExprForVertexExpr <A>& a , const std::complex<double>& b ) 
    {
      return PolynomialExprForVertexExpr_CMult < A, std::complex<double> > (~a,b);
    }
//-------------------------------------------------------------------------------------
template<class A, class B>
struct PolynomialExprForVertexExpr_Div 
  : public PolynomialExprForVertexExpr < PolynomialExprForVertexExpr_Div <A,B> > 
    {
      static void print (std::ostream& os)
         {
           os << "PolynomialExprForVertexExprExpr_Div<";
           A::print(os); 
           os << ",";
           B::print(os); 
           os << "> ";
         }
      template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints( const Element& element, const unsigned int gaussLevel )
          {
            switch (dx + dy + dz)
               {
                 case 0:
                   return A:: template evaluateGaussianPoints <dx, dy, dz> (element,gaussLevel) / 
                          B:: template evaluateGaussianPoints <dx, dy, dz> (element,gaussLevel) ;
                   break;
                 case 1: 
                   return ( B :: template evaluateGaussianPoints <0,0,0> (element,gaussLevel) * 
                            A :: template evaluateGaussianPoints <dx,dy,dz> (element,gaussLevel) -
           	            A :: template evaluateGaussianPoints <0,0,0> (element,gaussLevel) * 
                            B :: template evaluateGaussianPoints <dx,dy,dz> (element,gaussLevel) ) / 
                          ( B :: template evaluateGaussianPoints <0,0,0> (element,gaussLevel) * 
                            B :: template evaluateGaussianPoints <0,0,0> (element,gaussLevel) );
                   break;
               }
            return 0.;
          }
      template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints( const Element& element, int dx, int dy, int dz, const unsigned int gaussLevel )
          {
            switch (dx + dy + dz)
               {
                 case 0:
                   return A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) / 
                          B::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ;
                   break;
                 case 1: 
                   return ( B ::evaluateGaussianPoints(element,0,0,0,gaussLevel) * 
                            A ::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) -
           	            A ::evaluateGaussianPoints(element,0,0,0,gaussLevel) * 
                            B ::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ) / 
                          ( B ::evaluateGaussianPoints(element,0,0,0,gaussLevel) * 
                            B ::evaluateGaussianPoints(element,0,0,0,gaussLevel) );
                   break;
               }
            return 0.;
          }
    } ;
//-------------------------------------------------------------------------------------
template<class A, class B>
inline PolynomialExprForVertexExpr_Div <A,B>
 operator / ( const PolynomialExprForVertexExpr <A>& a , 
              const PolynomialExprForVertexExpr <B>& b ) 
    {
      return PolynomialExprForVertexExpr_Div <A,B> () ;
    }
//-------------------------------------------------------------------------------------
template<class A, class D>
class PolynomialExprForVertexExpr_CDiv 
  : public PolynomialExprForVertexExpr < PolynomialExprForVertexExpr_CDiv <A,D> > 
    {
       static const D b ;
      public :
       PolynomialExprForVertexExpr_CDiv (const D& b_) 
          {  
            b = b_; 
          }
     template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints( const Element& element, const unsigned int gaussLevel )
           {
             return A:: template evaluateGaussianPoints <dx, dy, dz> (element,gaussLevel) / b ; 
           }
     template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints( const Element& element, int dx, int dy, int dz, const unsigned int gaussLevel )
           {
             return A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) / b ; 
           }
    } ;

template<class A, class D>
const D 
 PolynomialExprForVertexExpr_CDiv <A,D> :: b ;
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CDiv <A,float>
 operator / ( const PolynomialExprForVertexExpr<A> & a , float b ) 
    {
      return PolynomialExprForVertexExpr_CDiv <A,float> (b) ;
    }
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CDiv <A,double>
 operator / ( const PolynomialExprForVertexExpr<A> & a , double b ) 
    {
      return PolynomialExprForVertexExpr_CDiv <A,double> (b) ;
    }
//-------------------------------------------------------------------------------------
template < class A > 
inline PolynomialExprForVertexExpr_CDiv < A, std::complex<double> >
 operator / ( const PolynomialExprForVertexExpr <A>& a, const std::complex<double>& b ) 
    {
      return PolynomialExprForVertexExpr_CDiv < A, std::complex<double> > (b);
    }
//-----------------------------------------------------------------------------
template < class A >
struct EXP_F : public PolynomialExprForVertexExpr < EXP_F<A> >  
   {
      static void print (std::ostream& os)
         {
           os << "EXP_F<";
           A::print(os); 
           os << "> ";
         }
     template <int dx, int dy, int dz, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints( const Element& element, const unsigned int gaussLevel )
         {
           switch ( dx + dy + dz) 
              {
                case 0:
                  return std::exp( A:: template evaluateGaussianPoints <dx, dy, dz> (element,gaussLevel) ) ;
                  break;    
                case 1:
                  return exp ( A:: template evaluateGaussianPoints <0,0,0> (element,gaussLevel) ) * 
                               A:: template evaluateGaussianPoints <dx, dy, dz> (element,gaussLevel) ;
                  break;
              }
           return 0. ;
         }
     template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints( const Element& element, int dx, int dy, int dz, const unsigned int gaussLevel )
         {
           switch (dx + dy + dz) 
              {
                case 0:
                  return std::exp(A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ) ;
                  break;    
                case 1:
                  return exp (A::evaluateGaussianPoints(element,0,0,0,gaussLevel) ) * 
                              A::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel) ;
                  break;
              }
           return 0. ;
         }
   } ;
//-----------------------------------------------------------------------------
template < class A > 
inline EXP_F <A> 
 Exp ( const PolynomialExprForVertexExpr <A> & a ) 
    {
      return EXP_F <A> (~a);
    }
//-----------------------------------------------------------------------------
template < int vertexNumber >
struct Vertex_ : public VertexExpr < Vertex_ < vertexNumber > > 
   {
      static void print (std::ostream& os)
         {
           os << "Vertex_<" << vertexNumber << "> ";
         }
     template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
         {
           enum {pos = Element :: dimensionOfElement * vertexNumber + k};
           if ( dx+dy+dz == 0 ) 
              {
                return element.getCornerList()[ pos ] ; 
              }
            else 
              {
                return 0.;
              }
         }
     template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel )
         {
           const int pos = Element :: dimensionOfElement * vertexNumber + k;
           if ( dx+dy+dz == 0 ) 
              {
                return element.getCornerList()[pos] ; 
              }
           return 0.;
         }
     template <int k, class Element>
      static inline typename Element::bTYPE 
      get (const Element& element) 
         {
           enum {pos = Element :: dimensionOfElement * vertexNumber + k};
           return element.getCornerList()[ pos ] ; 
         }

   } ;
//-----------------------------------------------------------------------------
template < int vertexNumber , typename Type = double >
class Vertex_Constant : public VertexExpr < Vertex_Constant < vertexNumber, Type > > 
   {
      static Type value;
     public:
      Vertex_Constant(Type val) { value = val; }
      static void print (std::ostream& os)
         {
           os << "Vertex_Constant<" << vertexNumber << "> (" << value << ")";
         }
     template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
         {
           if ( dx+dy+dz == 0 ) 
              {
                return value;
              }
            else 
              {
                return 0.;
              }
         }
     template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel )
         {
           if ( dx+dy+dz == 0 ) 
              {
                return value;
              }
            else 
              {
                return 0.;
              }
         }
     template <int k, class Element>
      static inline typename Element::bTYPE 
      get (const Element& element) 
         {
           return value;
         }

   } ;
template < int vertexNumber , typename Type >
Type Vertex_Constant<vertexNumber,Type>::value = Type(0);
//-----------------------------------------------------------------------------
template <int number,typename Type>
inline Vertex_Constant<number,Type> 
_C(const Type t ){
  return Vertex_Constant<number,Type>(t);  
}
//-----------------------------------------------------------------------------
template <class A, class B>
struct VertexAdd : public VertexExpr < VertexAdd <A,B> > 
   {
      static void print (std::ostream& os)
         {
           os << "VertexAdd<";
           A::print(os); 
           os << ",";
           B::print(os);
           os << "> ";
         }
     VertexAdd(){}
     VertexAdd(const VertexAdd<A,B>& a){}
     template <int k, class Element>
     static inline typename Element::bTYPE 
      get (const Element& element) 
         {
           return A:: template get<k> (element) + 
                  B:: template get<k> (element) ;
         }
     template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
         {
           return A::template evaluateGaussianPoints<dx,dy,dz,k>(element,gaussLevel) + 
                  B::template evaluateGaussianPoints<dx,dy,dz,k>(element,gaussLevel);
         } 
     template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel )
         {
           return A::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel) + 
                  B::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel);
         } 
   } ;
//-----------------------------------------------------------------------------
template <class A, class B>
struct VertexSub : public VertexExpr < VertexSub <A,B> > 
   {
      static void print (std::ostream& os)
         {
           os << "VertexSub<";
           A::print(os); 
           os << ",";
           B::print(os);
           os << "> ";
         }
     template <int k, class Element>
     static inline typename Element::bTYPE 
      get (const Element& element) 
         {
           return A:: template get<k> (element) - B:: template get<k> (element);
         }
     template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
         {
           return A::template evaluateGaussianPoints<dx,dy,dz,k>(element,gaussLevel) - 
                  B::template evaluateGaussianPoints<dx,dy,dz,k>(element,gaussLevel);
         }
     template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel )
         {
           return A::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel) - 
                  B::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel);
         }
   } ;
//-----------------------------------------------------------------------------
#if 1
template <typename Type, class B>
struct ConstVertexMult : public VertexExpr < ConstVertexMult <Type,B> > 
   {
      static Type value;
      ConstVertexMult (Type val) { value=val; }
      static void print (std::ostream& os)
         {
           os << "VertexSub<";
           value;
           os << ",";
           B::print(os);
           os << "> ";
         }
     template <int k, class Element>
     static inline typename Element::bTYPE 
      get (const Element& element) 
         {
           return value * B:: template get<k> (element);
         }
     template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
         {
#if 0
           std::cout << "§§§§§§" << value << std::endl;
#endif
           return value * B::template evaluateGaussianPoints<dx,dy,dz,k>(element,gaussLevel);
         }
     template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel )
         {
           return value * B::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel);
         }
   } ;
template <typename Type, class B>
Type ConstVertexMult<Type,B>::value = Type(); 
#endif
//-----------------------------------------------------------------------------
/* This cross product (initialized by the operator *) exists only for elements
   with dimension 3. However, this is not caught yet. 
*/
template <class A, class B>
struct VertexCrossProduct : public VertexExpr < VertexCrossProduct <A,B> >  
   {
      static void print (std::ostream& os)
         {
           os << "VertexCrossProduct<";
           A::print(os); 
           os << ",";
           B::print(os);
           os << "> ";
         }
     template <int k, class Element>
     static inline typename Element::bTYPE 
      get (const Element& element) 
         {
           return A:: template get<(k+1)%Element::dimensionOfElement> (element) * 
                  B:: template get<(k+2)%Element::dimensionOfElement> (element) -
	          A:: template get<(k+2)%Element::dimensionOfElement> (element) * 
                  B:: template get<(k+1)%Element::dimensionOfElement> (element) ; 
         }
     template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
         {
            return (
             A::template evaluateGaussianPoints<dx,dy,dz,(k+1)%Element::
                      dimensionOfElement>(element,gaussLevel) * 
             B::template evaluateGaussianPoints<dx,dy,dz,(k+2)%Element::
                      dimensionOfElement>(element,gaussLevel) -
             A::template evaluateGaussianPoints<dx,dy,dz,(k+2)%Element::
                      dimensionOfElement>(element,gaussLevel) * 
             B::template evaluateGaussianPoints<dx,dy,dz,(k+1)%Element::
                      dimensionOfElement>(element,gaussLevel) 
                   );
         }
     template <class Element>
      static inline typename Element::TYPE 
       evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel )
         {
            return (
             A::evaluateGaussianPoints(element,dx,dy,dz,(k+1)%Element::dimensionOfElement,gaussLevel) * 
             B::evaluateGaussianPoints(element,dx,dy,dz,(k+2)%Element::dimensionOfElement,gaussLevel) -
             A::evaluateGaussianPoints(element,dx,dy,dz,(k+2)%Element::dimensionOfElement,gaussLevel) * 
             B::evaluateGaussianPoints(element,dx,dy,dz,(k+1)%Element::dimensionOfElement,gaussLevel) 
                   );
         }
   } ;
//-----------------------------------------------------------------------------
template <class A, class B>
struct VertexConstantProduct 
  : public VertexExpr < VertexConstantProduct <A,B> > 
    {
      static void print (std::ostream& os)
         {
           os << "VertexPolynomialProduct<";
           A::print(os); 
           os << ",";
           B::print(os);
           os << "> ";
         }
      template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
          {
            switch (dx+dy+dz)
               {
                 case 0: 
                    return (A::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel) * 
                            B::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel));
                    break;
                 case 1: 
                   return A::evaluateGaussianPoints(element,0,0,0,k,gaussLevel) * 
                          B::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel);
                    break;
               }
             return 0;
          }
      template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel )
          {
            switch (dx+dy+dz)
               {
                 case 0: 
                    return (A::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel) * 
                            B::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel));
                    break;
                 case 1: 
                   return A::evaluateGaussianPoints(element,0,0,0,k,gaussLevel) * 
                          B::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel);
                    break;
               }
             return 0;
          }
    } ;
//-----------------------------------------------------------------------------
template <class A, class B>
struct VertexPolynomialProduct 
  : public VertexExpr < VertexPolynomialProduct <A,B> > 
    {
      static void print (std::ostream& os)
         {
           os << "VertexPolynomialProduct<";
           A::print(os); 
           os << ",";
           B::print(os);
           os << "> ";
         }
      template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
          {
            return Multiplication<A,B,dx+dy+dz,dx,dy,dz,Element>::
                           template derivateDependentEvaluationDirection<k>(element,gaussLevel); 
          }
      template <class Element>
      static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel )
          {
            switch (dx+dy+dz)
               {
                 case 0: 
                    return (A::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel) * 
                            B::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel));
                    break;
                 case 1: 
                   return A::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel) * 
                          B::evaluateGaussianPoints(element,0,0,0,gaussLevel) +
                          A::evaluateGaussianPoints(element,0,0,0,k,gaussLevel) * 
                          B::evaluateGaussianPoints(element,dx,dy,dz,gaussLevel);
                    break;
               }
             return 0;
          }
    } ;
//-----------------------------------------------------------------------------
template  < class A >
struct VertexMin : public VertexExpr < VertexMin < A > > 
   {
      static void print (std::ostream& os)
         {
           os << "VertexMin<";
           A::print(os); 
           os << "> ";
         }
    template <int k, class Element>
     static inline typename Element::bTYPE
      get (const Element& element) 
         {
           return - A :: template get<k>(element) ; 
         }
     template <int dx, int dy, int dz, int k, class Element>
      static inline typename Element::TYPE 
     evaluateGaussianPoints(const Element& element,const unsigned int gaussLevel )
         {
           return -A::template evaluateGaussianPoints<dx,dy,dz,k>(element,gaussLevel);
         }
     template <class Element>
      static inline typename Element::TYPE 
     evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k, const unsigned int gaussLevel)
         {
           return -A::evaluateGaussianPoints(element,dx,dy,dz,k,gaussLevel);
         }
   } ;
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexPolynomialProduct <A,B>
 operator * (const VertexExpr <A>& a , const PolynomialExprForVertexExpr <B>& b) 
    {
      return VertexPolynomialProduct <A,B> (); 
    }
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexPolynomialProduct <A,B>
 operator * ( const PolynomialExprForVertexExpr <B>& b , const VertexExpr <A>& a) 
    {
      return VertexPolynomialProduct <A,B> () ; 
    }
//-----------------------------------------------------------------------------
template <int index, typename Type, class B>
inline VertexConstantProduct <Vertex_Constant<index,Type>,B>
 operator * ( const Vertex_Constant<index,Type>& b , const VertexExpr <B>& a) 
    {
      return VertexConstantProduct <Vertex_Constant<index,Type>,B>() ; 
    }
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexAdd<A,B>
 operator + (const VertexExpr <A>& a , const VertexExpr <B>& b) 
    {
      return VertexAdd <A,B> () ; 
    }
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexSub <A,B>
 operator - (const VertexExpr <A>& a, const VertexExpr <B>& b) 
    {
      return VertexSub <A,B>() ;
    }
//-----------------------------------------------------------------------------
template < class A > 
inline VertexMin <A> 
 operator - (const VertexExpr <A>& a) 
    {
      return VertexMin <A> () ; 
    }
//-----------------------------------------------------------------------------
template <class Element, class A, class B>
inline typename Element::bTYPE
 vec (const Element& element, const VertexExpr <A>& a, const VertexExpr <B>& b) 
    {
      typename Element::TYPE result = 
         A:: template get<0> (element) * B:: template get<0> (element) +
         A:: template get<1> (element) * B:: template get<1> (element) ;
      if ((DIMENSION)(Element::dimensionOfElement) == D3) 
         {
           result += A:: template get<2> (element) * B:: template get <2> (element);
         }
      return result;
    }
//-----------------------------------------------------------------------------
template <class Element, class A>
inline typename Element::bTYPE
 length ( const Element& element, const VertexExpr <A>& a ) 
    {
      if ( (DIMENSION)(Element::dimensionOfElement) == D2 ) 
         {
           return A:: template get<0> (element) * A:: template get<0> (element) +
 	          A:: template get<1> (element) * A:: template get<1> (element) ;
         }
      return A:: template get<0> (element) * A:: template get<0> (element) +
             A:: template get<1> (element) * A:: template get<1> (element) +
             A:: template get<2> (element) * A:: template get<2> (element) ;
    }
//-----------------------------------------------------------------------------
template <class Element, class A, class B>
inline typename Element::bTYPE
 det (const Element& element, const VertexExpr <A>& a, const VertexExpr <B>& b ) 
    {
      if ( (DIMENSION)(Element::dimensionOfElement) == D2 ) 
         {
           return A:: template get<0> (element) * B:: template get<1> (element) +
 	          A:: template get<1> (element) * B:: template get<0> (element) ;
         }
      return 0.;
    }
//-----------------------------------------------------------------------------
template <class A, class B>
inline VertexCrossProduct <A,B>
 operator * (const VertexExpr <A>& a , const VertexExpr <B>& b ) 
   {
     return VertexCrossProduct <A,B> () ;
   }
//-----------------------------------------------------------------------------
template < class A >
struct D_DX_Vertex: public VertexExpr < D_DX_Vertex<A> > 
   {
     template <int dx, int dy, int dz, int k, class Element>
     static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element, const unsigned int gaussLevel )
         {
           return A:: template evaluateGaussianPoints <dx+1, dy, dz, k> (element,gaussLevel);
         }
     template <class Element>
     static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k,  const unsigned int gaussLevel )
         {
           return A::evaluateGaussianPoints(element,dx+1, dy, dz, k,gaussLevel);
         }
   };
//-----------------------------------------------------------------------------
template < class A >
struct D_DY_Vertex : public VertexExpr < D_DY_Vertex<A> > 
   {
     template <int dx, int dy, int dz, int k, class Element>
     static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element, const unsigned int gaussLevel )
         {
           return A:: template evaluateGaussianPoints <dx, dy+1, dz, k> (element,gaussLevel);
         }
     template <class Element>
     static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element,int dx, int dy, int dz, int k,  const unsigned int gaussLevel )
         {
           return A::evaluateGaussianPoints(element,dx, dy+1, dz, k, gaussLevel);
         }
   } ;
//-----------------------------------------------------------------------------
template < class A >
struct D_DZ_Vertex : public VertexExpr < D_DZ_Vertex<A> > 
   {
     template <int dx, int dy, int dz, int k, class Element>
     static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element, const unsigned int gaussLevel )
         {
           return A::template evaluateGaussianPoints <dx, dy, dz+1, k> (element,gaussLevel);
         }
     template <class Element>
     static inline typename Element::TYPE 
      evaluateGaussianPoints(const Element& element, int dx, int dy, int dz, int k, const unsigned int gaussLevel )
         {
           return A::evaluateGaussianPoints(element,dx, dy, dz+1, k, gaussLevel);
         }
   };
//-----------------------------------------------------------------------------
template < class A > 
inline D_DX_Vertex<A> 
 d_dx ( const VertexExpr<A>& a ) 
    {
      return D_DX_Vertex<A>() ;
    }
//-----------------------------------------------------------------------------
template < class A > 
inline D_DY_Vertex<A> 
 d_dy ( const VertexExpr<A>& a ) 
    {
      return D_DY_Vertex<A>() ;
    }
//-----------------------------------------------------------------------------
template < class A > 
inline D_DZ_Vertex<A> 
 d_dz ( const VertexExpr<A>& a ) 
    {
      return D_DZ_Vertex<A>();
    }
//-----------------------------------------------------------------------------
template <int p1, int p2, int p3, int p4, int p5>
template <class Element>
inline typename Element::bTYPE
PyramidVolume<p1,p2,p3,p4,p5>::compute(const Element& element) 
   {
#if 0
     std::cout << vec (element, Vertex_<Element::numberOfCorners>()-Vertex_<p1>(), 
                        (Vertex_<p2>()-Vertex_<p4>())*(Vertex_<p3>() - Vertex_<p1>()) ) << " ";
     std::cout << vec (element, Vertex_<p5>()-Vertex_<p1>(), 
                        (Vertex_<p2>()-Vertex_<p4>())*(Vertex_<p3>() - Vertex_<p1>()) ) << "      " ;
     std::cout << fabs(
           vec (element, Vertex_<Element::numberOfCorners>()-Vertex_<p1>(), 
                        (Vertex_<p2>()-Vertex_<p4>())*(Vertex_<p3>() - Vertex_<p1>()) ) /
           vec (element, Vertex_<p5>()-Vertex_<p1>(), 
                         (Vertex_<p2>()-Vertex_<p4>()) * (Vertex_<p3>()-Vertex_<p1>())) ) << std::endl;
     return vec (element, Vertex_<Element::numberOfCorners>()-Vertex_<p1>(), 
                        (Vertex_<p2>()-Vertex_<p4>())*(Vertex_<p3>() - Vertex_<p1>()) ) /
           vec (element, Vertex_<p5>()-Vertex_<p1>(), 
                         (Vertex_<p2>()-Vertex_<p4>()) * (Vertex_<p3>()-Vertex_<p1>())) ;
#else
     return fabs(
           vec (element, Vertex_<Element::numberOfCorners>()-Vertex_<p1>(), 
                        (Vertex_<p2>()-Vertex_<p4>())*(Vertex_<p3>() - Vertex_<p1>()) ) /
           vec (element, Vertex_<p5>()-Vertex_<p1>(), 
                         (Vertex_<p2>()-Vertex_<p4>()) * (Vertex_<p3>()-Vertex_<p1>())) );
#endif
   }
//-----------------------------------------------------------------------------
template <int p1, int p2, int p3, int p4>
template <class Element>
inline typename Element::bTYPE
TetrahedronVolume<p1,p2,p3,p4>::compute(const Element& element)   
         {
           return fabs (
                  vec(element, Vertex_<p2>()-Vertex_<p1>(), 
                               (Vertex_<p3>()-Vertex_<p1>()) * (Vertex_<Element::numberOfCorners>()-Vertex_<p1>()) ) /
	          vec(element, Vertex_<p2>()-Vertex_<p1>(), 
                               (Vertex_<p3>()-Vertex_<p1>()) * (Vertex_<p4>()-Vertex_<p1>())) );
         } 
//-----------------------------------------------------------------------------
template <int p1, int p2, int p3>
template <class Element>
inline typename Element::bTYPE
TriangleVolume<p1,p2,p3>::compute(const Element& element)
         {
           if ((DIMENSION)(Element::dimensionOfElement) == D2)
              {
                 return fabs (
                         det(element,Vertex_<p2>() - Vertex_<p1>() , Vertex_<Element::numberOfCorners>() - Vertex_<p1>() ) /
                         det(element,Vertex_<p2>() - Vertex_<p1>() , Vertex_<p3>() - Vertex_<p1>() ) );
              }
#if 0
           return fabs(length(element,vec(element,Vertex_<p2>() - Vertex_<p1>() , 
                                                    Vertex_<Element::numberOfCorners>() - Vertex_<p1>() ) ) / 
                        length(element,vec(element,Vertex_<p2>() - Vertex_<p1>() , 
                                                               Vertex_<p3>() - Vertex_<p1>() ) ) );
#endif
           return fabs(length(element,(Vertex_<p2>() - Vertex_<p1>())* Vertex_<Element::numberOfCorners>() - Vertex_<p1>() ) / 
                        length(element,(Vertex_<p2>() - Vertex_<p1>() )*(Vertex_<p3>() - Vertex_<p1>() ) ) );
        }
//==============================================================================
