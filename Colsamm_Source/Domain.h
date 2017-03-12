//==============================================================================
//
//  $Id: Domain.h,v 1.58 2006/07/19 09:05:18 jochen Exp $
//
//==============================================================================

#define TRANSFORMATION( B ) typedef __typeof__(B)  TRAFO
#define Define_Transformation( B, A ) namespace A##FOMEL{ __typeof__(B) formel(B); } typedef __typeof__(B) A
#define Define_Element_Transformation( B ) typedef __typeof__(B) 
//-----------------------------------------------------------------------------
#define TEMPLATES_OF_CLASS_DOMAIN_DECL \
          template < unsigned int cornersOfElement, \
                     DIMENSION init_dimension,  \
                     class TRAFO, \
                     domainSecification interiorIboundary, \
                     unsigned int numberOfAnsatzFunctions, \
                     unsigned int numberOfSecondAnsatzFunctions, \
                     ElementType elType, \
                     IntegrationAccuracy intAcc, \
                     typename ETYPE = double, \
                     stencilType STENCIL = STLVector,\
                     typename eTYPE = double, \
                     BasisFunctionType BFType = oneDimensional >
//-----------------------------------------------------------------------------
#define TEMPLATES_OF_CLASS_DOMAIN \
          template < unsigned int cornersOfElement, \
                     DIMENSION init_dimension,  \
                     class TRAFO, \
                     domainSecification interiorIboundary, \
                     unsigned int numberOfAnsatzFunctions, \
                     unsigned int numberOfSecondAnsatzFunctions, \
                     ElementType elType, \
                     IntegrationAccuracy intAcc, \
                     typename ETYPE, \
                     stencilType STENCIL,\
                     typename eTYPE, \
                     BasisFunctionType BFType >
//-----------------------------------------------------------------------------
#define DOMAIN_WITH_TEMPLATES \
          _Domain_ < cornersOfElement,\
                     init_dimension,\
                     TRAFO,\
                     interiorIboundary,\
                     numberOfAnsatzFunctions, \
                     numberOfSecondAnsatzFunctions,\
                     elType,\
                     intAcc,\
                     ETYPE,\
                     STENCIL,\
                     eTYPE,\
                     BFType>
//-----------------------------------------------------------------------------
#define STENCIL_TYPE \
    typename ResultType< \
       typename Stencil_Init<ETYPE,STENCIL,numberOfAnsatzFunctions,numberOfSecondAnsatzFunctions>::StencilType,\
        Integrand>::DimType
//-----------------------------------------------------------------------------
#define VEC_STENCIL_TYPE \
    typename ResultType< \
       typename Stencil_Init<ETYPE,STENCIL,numberOfAnsatzFunctions,numberOfSecondAnsatzFunctions>::StencilType,\
        Integrand>::Vec_DimType
//=============================================================================
// Introducing the additional template type, in order to steer the type of the
// resulting stencil, i.e. nested stl vectors or an array
//=============================================================================
/**
 The class Domain serves as center within the Colsamm coding! However, only 
 inherited classes of Domain are used and get allocated. 
*/
TEMPLATES_OF_CLASS_DOMAIN_DECL
class _Domain_
  : public Gaussian_Points<eTYPE,intAcc,elType>,
    public Basis_Functions<ETYPE,init_dimension,Gaussian_Points<eTYPE,intAcc,elType>::numberGaussianPoints+1,
                          BFType,numberOfAnsatzFunctions,numberOfSecondAnsatzFunctions>,
    public Mappings<ETYPE,init_dimension,interiorIboundary,TRAFO,Gaussian_Points<eTYPE,intAcc,elType>::numberGaussianPoints>,
    public Corner_Classes<eTYPE,init_dimension,cornersOfElement,elType>,
    public Stencil_Init<ETYPE,STENCIL,numberOfAnsatzFunctions,numberOfSecondAnsatzFunctions>
   {
  // -----------------------------------------------------------------------
     public:
      typedef ETYPE TYPE;
      typedef eTYPE bTYPE;
      typedef Basis_Functions<ETYPE,init_dimension,Gaussian_Points<eTYPE,intAcc,elType>::numberGaussianPoints+1,
                          BFType,numberOfAnsatzFunctions,numberOfSecondAnsatzFunctions>  BasisFunctions;
      using Corner_Classes<eTYPE,init_dimension,cornersOfElement,elType>::numberOfCorners;
  // -----------------------------------------------------------------------
    /**
      Method for evaluating the Dirac functions at particular points  
     */
      template < class Integrand >
      inline const STENCIL_TYPE
       value ( const BasisFunctionExpr<Integrand> & a_ );
    /**
      Method for evaluating the Dirac functions at particular points  
     */
      template < class Integrand , int length>
      inline const VEC_STENCIL_TYPE
       value ( const BasisFunctionVectorExpr<Integrand,length> & a_ );

    /**
      Method for integrating an integrand of a weak formulation,
      returning the local stencil
     */
      template < class Integrand >
      inline STENCIL_TYPE
       integrate ( const BasisFunctionExpr < Integrand > & a_ );
    /**
      Method for integrating an integrand of a weak formulation,
      computing the local stencil in to the first component
     */
      template < typename PointArray, class Integrand >
      inline void
       integrate ( PointArray& array, const BasisFunctionExpr < Integrand > & a_ );

    /**
     Method for setting the basis functions into set 1
     */
      template < class FunctionExpression>
      inline void
       Set ( const FunctionExpression & a_ );
    /** 
     Method for setting the basis functions into any other set 
     by specifying the set Number
     */
      template < class FunctionExpression >
      inline void
       Set (BasisFunctionSet setNumber, const FunctionExpression & a_ );

    /**
     Method for setting the basis functions into set 1
     */
      template < class FunctionExpression>
      inline void
       Set ( const FunctionExpression & a_, const unsigned int cnr_1, const unsigned int cnr_2);


    /**
     Method for setting a special point for evaluation
     */
      template < class PointArrayTYPE >
      void
       Point_ ( const PointArrayTYPE & point );

    /**
     Method for setting a special point for evaluation
     */
      void
       Point_ ();

    /**
     Method for setting a special point for evaluation including
     the test whether the point lies inside or outside the actual
     element.
     */ 
      template < class PointArrayTYPE >
      bool
       Point_Test ( const PointArrayTYPE & point );

    /**
     Method for passing the vertices of the actual element
     */
      template < class PointArrayTYPE >
      inline DOMAIN_WITH_TEMPLATES&
       operator ( ) ( const PointArrayTYPE & cc , const int numberPoints = cornersOfElement);
   };
 //-------------------------------------------------------------------------------
  template<int numberCorners, DIMENSION dim, class TRAFO, int numberOfFunctions, ElementType elType,
           IntegrationAccuracy intAcc, typename ETYPE = double, stencilType STENCIL = STLVector, typename eTYPE = double>
  class Simple_Element : 
    public _Domain_ <numberCorners, dim, TRAFO, interior, numberOfFunctions,0,elType,intAcc,ETYPE,STENCIL,eTYPE,oneDimensional>
   {};
 //-------------------------------------------------------------------------------
  template<int numberCorners, DIMENSION dim, class TRAFO, int numberOfFunctions, ElementType elType,
           IntegrationAccuracy intAcc, typename ETYPE = double, stencilType STENCIL = STLVector, typename eTYPE = double>
  class Simple_Boundary_Element : 
    public _Domain_ <numberCorners, dim, TRAFO, boundary, numberOfFunctions,0,elType,intAcc,ETYPE,STENCIL,eTYPE,oneDimensional>
   {};
 //-------------------------------------------------------------------------------
  template<int numberCorners, DIMENSION dim, class TRAFO, int numberOfFunctions1, int numberOfFunctions2, ElementType elType,
           IntegrationAccuracy intAcc, domainSecification interiorIboundary = interior, 
           typename ETYPE = double,stencilType STENCIL = STLVector, typename eTYPE = double>
  class Mixed_Element : 
    public _Domain_<numberCorners,dim,TRAFO,interiorIboundary,numberOfFunctions1,numberOfFunctions2,elType,intAcc,ETYPE,STENCIL,eTYPE,oneDimensional>
   {};
 //-------------------------------------------------------------------------------
  template<int numberCorners, DIMENSION dim, class TRAFO, int numberOfFunctions, ElementType elType,
           IntegrationAccuracy intAcc, domainSecification interiorIboundary = interior, 
           typename ETYPE = double,stencilType STENCIL = STLVector, typename eTYPE = double>
  class Vectorial_Element : 
    public _Domain_ <numberCorners,dim,TRAFO,interiorIboundary,numberOfFunctions,0,elType,intAcc,ETYPE,STENCIL,eTYPE,vectorBasisFunctions>
   {};

//==============================================================================

