//==============================================================================
//
//  $Id: Loops.C,v 1.12 2006/07/19 11:30:48 jochen Exp $
//
//==============================================================================
// Some arbitrary enums, in order to improve the readability of the 
// implementation!
//-----------------------------------------------------------------------------
template <typename Type, class Integrand>
struct ResultingType<Type***,Integrand,0>
   {
      typedef Type DimType;
   };
//-----------------------------------------------------------------------------
template <typename Type, class Integrand>
struct ResultingType<Type***,Integrand,1>
   {
      typedef Type* DimType;
   };
//-----------------------------------------------------------------------------
template <typename Type, class Integrand>
struct ResultingType<Type***,Integrand,2>
   {
      typedef Type** DimType;
   };
//-----------------------------------------------------------------------------
template <typename Type, class Integrand>
struct ResultingType<Type***,Integrand,3>
   {
      typedef Type*** DimType;
   };
//-----------------------------------------------------------------------------
template <typename Type, class Integrand>
struct ResultingType<std::vector<std::vector<std::vector<Type> > >,Integrand ,0>
   {
      typedef Type DimType;
   };
//-----------------------------------------------------------------------------
template <typename Type, class Integrand>
struct ResultingType<std::vector<std::vector<std::vector<Type> > >,Integrand ,1>
   {
      typedef typename std::vector<Type> DimType;
   };
//-----------------------------------------------------------------------------
template <typename Type, class Integrand>
struct ResultingType<std::vector<std::vector<std::vector<Type> > >,Integrand ,2>
   {
      typedef typename std::vector<std::vector<Type> > DimType;
   };
//-----------------------------------------------------------------------------
template <typename Type, class Integrand>
struct ResultingType<std::vector<std::vector<std::vector<Type> > >,Integrand ,3>
   {
      typedef typename std::vector<std::vector<std::vector<Type> > > DimType;
   };
//-----------------------------------------------------------------------------
StencilManagement::StencilManagement()
   {
     recompute = true;
   }
//-----------------------------------------------------------------------------
template <class PointArrayTYPE>
inline void
StencilManagement:: printStencil(const PointArrayTYPE& array) const
    {
       std::cout << "***** Stencil **** " << std::endl;
       std::cout << "\033[33m" ;
       std::cout.width(12);
       std::cout << array << std::endl;
       std::cout << "\033[m" ;
    }
// ---------------------------------------------------------------------------
template <class PointArrayTYPE>
inline void 
StencilManagement:: printStencil(const PointArrayTYPE& array, const int size) const
    {
       std::cout << "***** Stencil **** " << std::endl;
       std::cout << "\033[33m" ;
       for (int iter = 0; iter < size ; ++iter)
          {
            std::cout.width(12);
            std::cout << array[iter] << " ";
          }
       std::cout << std::endl;
       std::cout << "\033[m" ;
    }
// ---------------------------------------------------------------------------
template <class PointArrayTYPE>
inline void 
StencilManagement:: printStencil(const PointArrayTYPE& array, const int size1, const int size2) const
    {
       std::cout << "***** Stencil **** " << std::endl;
       std::cout << "\033[33m" ;
       for (int iter1 = 0; iter1 < size1 ; ++iter1)
          {
            for (int iter2 = 0; iter2 < size2 ; ++iter2)
               {
                 std::cout.width(12);
                 std::cout << array[iter1][iter2] << " ";
               }
            std::cout << std::endl;
          }
       std::cout << std::endl;
       std::cout << "\033[m" ;
    }
// ---------------------------------------------------------------------------
template <class PointArrayTYPE>
inline void
StencilManagement:: printStencil(const PointArrayTYPE& array, const int size1, const int size2, const int size3) const
    {
       std::cout << "***** Stencil **** " << std::endl;
       std::cout << "\033[33m" ;
       for (int iter1 = 0; iter1 < size1 ; ++iter1)
          {
            for (int iter2 = 0; iter2 < size2 ; ++iter2)
               {
                  for (int iter3 = 0; iter3 < size3 ; ++iter3)
                     {
                       std::cout.width(12);
                       std::cout << array[iter1][iter2][iter3]<< " ";
                     }
                  std::cout << std::endl;
               }
            std::cout << std::endl;
          }
        std::cout << std::endl;
       std::cout << "\033[m" ;
    }
//-----------------------------------------------------------------------------
template <typename TYPE, int basisFn1, int basisFn2>
Stencil_Init<TYPE, STLVector, basisFn1, basisFn2> ::Stencil_Init()
   {
     stencil = std::vector < std::vector < std::vector < TYPE > > >  
         ( basisFn2+1 , std::vector < std::vector < TYPE > > 
         ( basisFn1 , std::vector < TYPE > (basisFn1, 0. ) ) ) ;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int basisFn1, int basisFn2>
Stencil_Init<TYPE, STLVector, basisFn1, basisFn2> ::~Stencil_Init()
   {
     for (int i=0; i < basisFn2+1; ++i)
        {
          for (int j=0; j < basisFn1; ++j)
             {
               stencil[i][j].clear();
             }
          stencil[i].clear();
        }
     stencil.clear();

   }
//-----------------------------------------------------------------------------
template <typename TYPE, int basisFn1, int basisFn2>
typename Stencil_Init<TYPE, STLVector, basisFn1, basisFn2>::StencilType&
Stencil_Init<TYPE, STLVector, basisFn1, basisFn2> ::getElementStencil()
   {
     return stencil;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int basisFn1, int basisFn2>
Stencil_Init<TYPE, Array, basisFn1, basisFn2> ::Stencil_Init()
   {
     stencil = new TYPE** [basisFn2+1];
     for (int i=0; i < basisFn2+1; ++i)
        {
          stencil[i] = new TYPE*[basisFn1];
          for (int j=0; j < basisFn1; ++j)
             {
               stencil[i][j] = new TYPE [basisFn1];
#ifdef INITIALIZE_ARRAYS
               for (int k=0; k < basisFn1; ++k)
                  {
                    stencil[i][j][k] = 0.;
                  }
#endif
             }
         }
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int basisFn1, int basisFn2>
Stencil_Init<TYPE, Array, basisFn1, basisFn2> ::~Stencil_Init()
   {
     for (int i=0; i < basisFn2+1; ++i)
        {
          for (int j=0; j < basisFn1; ++j)
             {
               delete [] stencil[i][j];
             }
          delete [] stencil[i];
        }
     delete [] stencil;
   }
//-----------------------------------------------------------------------------
template <typename TYPE, int basisFn1, int basisFn2>
typename Stencil_Init<TYPE, Array, basisFn1, basisFn2>::StencilType&
Stencil_Init<TYPE, Array, basisFn1, basisFn2> ::getElementStencil()
   {
     return (*this);
   }
//-----------------------------------------------------------------------------
template <class Element, class Integrand>
typename Element::TYPE 
Compile_Time_Loop <0>::gaussianQuadrature ( Element& element , const BasisFunctionExpr<Integrand>& expr) 
   {
     int iterator[1] = {0};
     typename Element::TYPE value = 0.;
     for (iterator[0]=0; iterator[0]< Element::numberGaussianPoints; ++iterator[0]) 
        {
          value += expr.eval(element,iterator) * element.getTransformationFactorTimesWeight(iterator[0]);
        }
     return value;
   }
//-----------------------------------------------------------------------------
template <class Element, class Integrand>
void
Compile_Time_Loop <0>::gaussianQuadrature ( Element& element , typename Element::TYPE & stencil, const BasisFunctionExpr<Integrand>& expr) 
   {
     int iterator[1] = {0};
     stencil = 0.;
     for (iterator[0]=0; iterator[0]< Element::numberGaussianPoints; ++iterator[0]) 
        {
          stencil += expr.eval(element,iterator) * element.getTransformationFactorTimesWeight(iterator[0]);
        }
   }
//-----------------------------------------------------------------------------
template <class Element,class Integrand>
typename Element::StencilType1Dim& 
Compile_Time_Loop <1>::gaussianQuadrature (Element& element , const BasisFunctionExpr<Integrand>& expr) 
   {
    const int loopSize = element.getNumberOfBasisFunctions((0<Integrand::hasP)?BasisSet_2:BasisSet_1);
    int iterator[2] = {0,0};
    typename Element::StencilType1Dim& data = element.getElementStencil()[0][0];
    for (iterator[0]=0; iterator[0]<loopSize; ++iterator[0]) 
        {
          typename Element::TYPE value = 0.;
          for (iterator[1]=0; iterator[1]<Element::numberGaussianPoints; ++iterator[1]) 
             {
               value += expr.eval(element,iterator) * element.getTransformationFactorTimesWeight(iterator[1]);
             }
          data[iterator[0]] = value;
        }
     return data ;
   }
//-----------------------------------------------------------------------------
template <class Element, class Stencil, class Integrand>
void
Compile_Time_Loop <1>::gaussianQuadrature (Element& element , Stencil& stencil, const BasisFunctionExpr<Integrand>& expr) 
   {
    const int loopSize = element.getNumberOfBasisFunctions((0<Integrand::hasP)?BasisSet_2:BasisSet_1);
    int iterator[2] = {0,0};
     for (iterator[0]=0; iterator[0]<loopSize; ++iterator[0]) 
        {
          typename Element::TYPE value = 0.;
          for (iterator[1]=0; iterator[1]<Element::numberGaussianPoints; ++iterator[1]) 
             {
               value += expr.eval(element,iterator) * element.getTransformationFactorTimesWeight(iterator[1]);
             }
          stencil[iterator[0]] = value;
        }
   }
//-----------------------------------------------------------------------------
template <class Element,class Integrand>
typename Element::StencilType2Dim& 
Compile_Time_Loop <2>::gaussianQuadrature (Element& element , const BasisFunctionExpr<Integrand>& expr) 
   { 
     const int loopSize1 = element.getNumberOfBasisFunctions(BasisSet_1);
     const int loopSize2 = element.getNumberOfBasisFunctions((0<Integrand::hasP)?BasisSet_2:BasisSet_1);
#if 0
     std::cout << "loopSize1 " << loopSize1 << "!      loopSize2 " << loopSize2 << std::endl;
#endif
     int iterator[3] = {0,0,0};
     typename Element::StencilType2Dim& data = element.getElementStencil()[0];
     for (iterator[0]=0; iterator[0]<loopSize1; ++iterator[0]) 
        {
          for (iterator[1]=0; iterator[1]<loopSize2; ++iterator[1]) 
             {
               typename Element::TYPE value = 0.;
               for (iterator[2]=0; iterator[2]<Element::numberGaussianPoints; ++iterator[2]) 
                  {
#if 0
                    std::cout << "iterator ["<< iterator[0] << "," <<  iterator[1] << "," <<  iterator[2] << "]" << std::endl;
#endif
                    value += expr.eval(element,iterator) * element.getTransformationFactorTimesWeight(iterator[2]);
                  }
               data[iterator[1]][iterator[0]] = value;
             }
        }
     return data;
   }
//-----------------------------------------------------------------------------
template <class Element, class Stencil,class Integrand>
void
Compile_Time_Loop <2>::gaussianQuadrature (const Element& element , Stencil& stencil, const BasisFunctionExpr<Integrand>& expr) 
   { 
     const int loopSize1 = element.getNumberOfBasisFunctions(BasisSet_1);
     const int loopSize2 = element.getNumberOfBasisFunctions((0<Integrand::hasP)?BasisSet_2:BasisSet_1);
#if 0
     std::cout << "loopSize1 " << loopSize1 << "!      loopSize2 " << loopSize2 << std::endl;
#endif
     int iterator[3] = {0,0,0};
     for (iterator[0]=0; iterator[0]<loopSize1; ++iterator[0]) 
        {
          for (iterator[1]=0; iterator[1]<loopSize2; ++iterator[1]) 
             {
               typename Element::TYPE value = 0.;
               for (iterator[2]=0; iterator[2]<Element::numberGaussianPoints; ++iterator[2]) 
                  {
#if 0
                    std::cout << "iterator ["<< iterator[0] << "," <<  iterator[1] << "," <<  iterator[2] << "]" << std::endl;
#endif
                    value += expr.eval(element,iterator) * element.getTransformationFactorTimesWeight(iterator[2]);
                  }
               stencil[iterator[1]][iterator[0]] = value;
             }
        }
   }
//-----------------------------------------------------------------------------
template <class Element,class Integrand>
typename Element::StencilType3Dim& 
Compile_Time_Loop <3>::gaussianQuadrature ( Element& element , const BasisFunctionExpr<Integrand>& expr) 
   {
    const int loopSize1 = element.getNumberOfBasisFunctions(BasisSet_1);
    const int loopSize2 = element.getNumberOfBasisFunctions(BasisSet_2);
    int iterator[4] = {0,0,0,0};
#if 0
    std:: cout << "loopSize1 " << loopSize1 << std::endl;
    std:: cout << "loopSize2 " << loopSize2 << std::endl;
#endif
    typename Element::StencilType3Dim& data = element.getElementStencil();
     for (iterator[0]=0; iterator[0]<loopSize1; ++iterator[0])
        {
          for (iterator[1]=0; iterator[1]<loopSize1; ++iterator[1]) 
             {
               for (iterator[2]=0; iterator[2]<loopSize2; ++iterator[2]) 
                  {
                    typename Element::TYPE value = 0.;
                    for (iterator[3]=0; iterator[3]<Element::numberGaussianPoints; ++iterator[3]) 
                       {
#if 0
                         std::cout << "It(" << iterator[0] <<","<< iterator[1] <<","<< iterator[2] <<","<< iterator[3] <<")" << std::endl;
#endif
                         value += expr.eval(element,iterator) * element.getTransformationFactorTimesWeight(iterator[3]);
                       }
                    data[iterator[2]][iterator[1]][iterator[0]] = value;
                  }
             }
        }
     return data;
   }
//----------------------------------------------------------------------------
inline float 
Colsamm_Internal_Functions::ABS<float> :: abs_(const float& c) 
   { 
     return fabs(c);
   } 
//----------------------------------------------------------------------------
inline double 
Colsamm_Internal_Functions::ABS <double> :: abs_(const double& c) 
   { 
     return fabs(c);
   } 
//----------------------------------------------------------------------------
inline double 
Colsamm_Internal_Functions::ABS <std::complex<double> > :: abs_ (const std::complex<double>&  c ) 
   {
     return abs(c);
   }
//----------------------------------------------------------------------------
namespace Colsamm_Internal_Functions
   {
//----------------------------------------------------------------------------
     template <typename TYPE> 
     struct DECIDE <std::complex<double>,TYPE >
        {
          typedef std::complex<double> I; 
        };
//----------------------------------------------------------------------------
     template <typename TYPE> 
     struct DECIDE <TYPE,std::complex<double> >
        {
          typedef std::complex<double> I; 
        };
//----------------------------------------------------------------------------
     template <> 
     struct DECIDE<double,double>
        {
          typedef double I; 
        };
//----------------------------------------------------------------------------
     template <> 
     struct DECIDE <double,float>
        {
          typedef double I; 
        };
//----------------------------------------------------------------------------
     template <> 
     struct DECIDE<float,double>
        {
          typedef double I; 
        };
//----------------------------------------------------------------------------
     template <> 
     struct DECIDE<float,float>
        {
          typedef float I; 
        };
//----------------------------------------------------------------------------
     template <int a, int b> 
     struct F_Fac 
        {
          enum { erg = (F_Fac <a,b+1>::erg*(b+1) ) }; 
        };
//----------------------------------------------------------------------------
     template <int a> 
     struct F_Fac<a,a>
        {  
          enum{erg=1};
      };
//----------------------------------------------------------------------------
     template <int a > 
     struct F_Fac<a,-1>
        { 
          enum{erg=0}; 
        };
    }
//----------------------------------------------------------------------------
template <unsigned int p>
template <class T>
inline T 
Colsamm_Internal_Functions::POW<p> :: pow(const T & d) 
   {
     return POW < p-1 > :: pow(d) * d;
   } 
//----------------------------------------------------------------------------
template <>
template <class T>
inline T 
Colsamm_Internal_Functions::POW<1> :: pow(const T & d) 
   {
     return d;
   } 
//----------------------------------------------------------------------------
template <>
template <class T>
inline T 
Colsamm_Internal_Functions::POW<0> :: pow(const T & d) 
   {
     return 1.;
   } 
//==============================================================================

