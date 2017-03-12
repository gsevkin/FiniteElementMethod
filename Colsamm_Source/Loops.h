//==============================================================================
//
//  $Id: Loops.h,v 1.34 2006/07/19 11:30:48 jochen Exp $
//
//==============================================================================
// Some arbitrary enums, in order to improve the readability of the 
// implementation!
#define _MAX_(a,b) ((static_cast<int>(a)>static_cast<int>(b)))?(static_cast<int>(a)):(static_cast<int>(b))
//----------------------------------------------------------------------------
enum IntegrationAccuracy 
   {
     Gauss1 = 1 , 
     Gauss2 = 2 , 
     Gauss3 = 3
   } ;
//----------------------------------------------------------------------------
enum ElementType 
   {
     cuboid , 
     tetrahedron , 
     pyramid , 
     prism , 
     quadrangle , 
     triangle , 
     interval
   } ;
//----------------------------------------------------------------------------
enum what 
   {
     fn = 0 ,
     dx = 1 ,
     dy = 2 ,
     dz = 3 ,
     dxx = 4 ,
     dxy = 5 ,
     dxz = 6 ,
     dyy = 7 ,
     dyz = 8 ,
     dzz = 9
   };
//----------------------------------------------------------------------------
enum SpaceDirection
   {
     dirX = 0 ,
     dirY = 1 ,
     dirZ = 2 ,
     weight = 3
   };
//----------------------------------------------------------------------------
enum 
   {
     _functionValues = 1,
     _functionAndFirstDerivativeValues = 4,
     _functionAndFirstAndSecondDerivativeValues = 10
   };
//----------------------------------------------------------------------------
enum ExtensionValues
   {
     determinant        = 0 ,
     normalLength       = 1 ,
     trafoFactor        = 2 ,
     _extensionSize = _functionAndFirstDerivativeValues
   };
//----------------------------------------------------------------------------
enum stencilType
   {
     STLVector = 0 ,
     Array = 1 
   };
//----------------------------------------------------------------------------
enum domainSecification
   {
     interior = 1 , 
     boundary = 2 ,
     boundaryPlain = 3
   } ;
//----------------------------------------------------------------------------
enum DIMENSION
   { 
     D1 = 1 , 
     D2 = 2 , 
     D3 = 3 
   } ;
//----------------------------------------------------------------------------
enum BasisFunctionSet
   {
     BasisSet_1 = 0,
     BasisSet_2 = 1, 
     BasisSet_3 = 2
   };
//----------------------------------------------------------------------------
enum IntegrandType
   {
     plain = 0,
     containV = 1, 
     containW = 2, 
     containVandW = 3, 
     containP = 4, 
     containVandP = 5, 
     containWandP = 6,
     containVandWandP = 7
   };
//-----------------------------------------------------------------------------
enum BasisFunctionType
   {
      oneDimensional = 0,
      vectorBasisFunctions = 1,
      derivativeBasisFunctions = 2
   };
//-----------------------------------------------------------------------------
enum TypeDOF{
   functionDOF = 0,
   first_derivativeDOF = 1
};
//-----------------------------------------------------------------------------
template <class Expr>
struct BasisFunctionExpr;
//-----------------------------------------------------------------------------
template <class Expr, int l>
struct BasisFunctionVectorExpr;
//-----------------------------------------------------------------------------
template <typename StencilType, class Integrand, int dimension>
struct ResultingType{};
//-----------------------------------------------------------------------------
template <typename StencilType, class Integrand>
struct ResultType
   { 
      enum {dimension = (Integrand::hasV>0?1:0)+(Integrand::hasW>0?1:0)+(Integrand::hasP>0?1:0)};
      typedef typename ResultingType<StencilType,Integrand,dimension>::DimType DimType;
   //   typedef typename ResultingType<StencilType,Integrand,dimension+1>::DimType Vec_DimType;
   };
//-----------------------------------------------------------------------------
class StencilManagement
   {
     protected: 
      bool recompute ;
     public:
      inline StencilManagement(); 
      template <class PointArrayTYPE>
      inline void printStencil(const PointArrayTYPE& array) const;
      template <class PointArrayTYPE>
      inline void printStencil(const PointArrayTYPE& array, const int size) const;
      template <class PointArrayTYPE>
      inline void printStencil(const PointArrayTYPE& array, const int size1, const int size2) const;
      template <class PointArrayTYPE>
      inline void printStencil(const PointArrayTYPE& array, const int size1, const int size2, const int size3) const;
   };
//-----------------------------------------------------------------------------
template <typename Type, stencilType stencType, int basisFn1, int basisFn2>
class Stencil_Init;
//-----------------------------------------------------------------------------
template <typename Type, int basisFn1, int basisFn2>
class Stencil_Init<Type, STLVector,basisFn1,basisFn2> : public StencilManagement
   {
     public:
       typedef typename std::vector<std::vector<std::vector<Type> > > StencilType3Dim; 
       typedef typename std::vector<std::vector<Type> > StencilType2Dim; 
       typedef typename std::vector<Type> StencilType1Dim; 
       typedef Type StencilType0Dim; 
       typedef typename std::vector<std::vector<std::vector<Type> > > StencilType; 
     protected: 
       StencilType stencil;
     public:
       inline Stencil_Init();
       inline ~Stencil_Init();
       inline StencilType& getElementStencil();
   };
//-----------------------------------------------------------------------------
template <typename Type, int basisFn1, int basisFn2>
class Stencil_Init<Type, Array,basisFn1,basisFn2> : public StencilManagement
   {
     public:
       typedef Type*** StencilType3Dim; 
       typedef Type** StencilType2Dim; 
       typedef Type* StencilType1Dim; 
       typedef Type StencilType0Dim; 
       typedef Type*** StencilType; 
     protected: 
       StencilType stencil;
     public:
       inline Stencil_Init();
       inline ~Stencil_Init();
       inline StencilType& getElementStencil();
   };
//-----------------------------------------------------------------------------
template <int LoopSize>
struct Compile_Time_Loop; 
//-----------------------------------------------------------------------------
template <>
struct Compile_Time_Loop <0> 
   {
     template <class Element, class Integrand>
     static inline typename Element::TYPE
      gaussianQuadrature ( Element& element , const BasisFunctionExpr<Integrand>& expr); 
     template <class Element, class Integrand>
     static inline void
      gaussianQuadrature ( Element& element , typename Element::TYPE & stencil , const BasisFunctionExpr<Integrand>& expr); 
   };
//-----------------------------------------------------------------------------
template <>
struct Compile_Time_Loop <1> 
   {
     template <class Element, class Integrand>
     static inline typename Element::StencilType1Dim&
      gaussianQuadrature ( Element& element , const BasisFunctionExpr<Integrand>& expr); 
     template <class Element, class Stencil, class Integrand>
     static inline void
      gaussianQuadrature ( Element& element , Stencil& stencil, const BasisFunctionExpr<Integrand>& expr); 
   };
//-----------------------------------------------------------------------------
template <>
struct Compile_Time_Loop <2> 
   {
     template <class Element, class Integrand>
     static inline typename Element::StencilType2Dim&
      gaussianQuadrature ( Element& element , const BasisFunctionExpr<Integrand>& expr); 
     template <class Element, class Stencil, class Integrand>
     static inline void
      gaussianQuadrature (const Element& element , Stencil& stencil, const BasisFunctionExpr<Integrand>& expr); 
   };
//-----------------------------------------------------------------------------
template <>
struct Compile_Time_Loop <3> 
   {
     template <class Element, class Integrand>
     static inline typename Element::StencilType3Dim&
      gaussianQuadrature ( Element& element , const BasisFunctionExpr<Integrand>& expr); 
    };
namespace Colsamm_Internal_Functions
  {
//----------------------------------------------------------------------------
   template <typename TYPE>
   struct ABS;
//----------------------------------------------------------------------------
   template <>
   struct ABS<float>
      {
        inline static float abs_(const float& c); 
      };
//----------------------------------------------------------------------------
   template <>
   struct ABS<double>
      {
        inline static double abs_(const double& c); 
       };
//----------------------------------------------------------------------------
   template <>
   struct ABS<std::complex<double> >
      {
        inline static double abs_(const std::complex<double>& c); 
      };
//----------------------------------------------------------------------------
   template <class A, class B> 
   struct DECIDE;
//----------------------------------------------------------------------------
   template <int a, int b> 
   struct F_Fac; 
//----------------------------------------------------------------------------
   template <unsigned int p>
   struct POW
      {
        template <class T>
        inline static T pow(const T &d); 
      };
   }
//==============================================================================

