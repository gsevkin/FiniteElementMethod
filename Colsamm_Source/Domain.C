// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template < class Integrand >
inline const STENCIL_TYPE
DOMAIN_WITH_TEMPLATES :: value ( const BasisFunctionExpr < Integrand > & expr ) 
   {
     int iterator[3] = {0,0,this->numberGaussianPoints};
     int &lauf = iterator[0];
     this->precomputeDerivatives(*this,BasisSet_1);
     for ( lauf = 0 ; lauf < this->numberSet[0] ; ++lauf ) 
        {
          iterator[1] = lauf ;
          this->stencil[0] [0] [lauf] = expr.eval(*this,iterator);
        }
     return this->stencil[0] [0] ;
   }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template < class Integrand , int length>
inline const VEC_STENCIL_TYPE
DOMAIN_WITH_TEMPLATES :: value ( const BasisFunctionVectorExpr < Integrand, length > & expr ) 
   {
     int iterator[3] = {0,0,this->numberGaussianPoints};
     int &lauf = iterator[0];
     this->point(*this,this->numberGaussianPoints);
     this->precomputeDerivativesAdditional(*this,BasisSet_1);
     const Integrand& integrand(expr);
     switch (Integrand::length)
        {
           case 2:
             for ( lauf = 0 ; lauf < this->numberSet[0] ; ++lauf ) 
                {
                  iterator[1] = lauf ;
                  this->stencil[0] [0] [lauf] = integrand.getHeadOfVector().eval(*this,iterator);
                  this->stencil[0] [1] [lauf] = integrand.evalLastComponent(*this,iterator);
                }
              break;
           case 3:
           default: 
             assert (!"case not implemented .... \n" );
        }
     return this->stencil[0] ;
   }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template < class Integrand >
inline STENCIL_TYPE
DOMAIN_WITH_TEMPLATES :: integrate ( const BasisFunctionExpr < Integrand > & expr ) 
    {
      if (this->recompute == true)
         {
           this->precomputeDerivatives(*this,BasisSet_1);
           if (0 < numberOfSecondAnsatzFunctions)
              {
                this->precomputeDerivatives(*this,BasisSet_2);
              }
           this->recompute = false;
         }
      enum {LoopSize=Integrand::hasV+Integrand::hasW+Integrand::hasP};
      return Compile_Time_Loop <LoopSize>::gaussianQuadrature(*this,expr);
    }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template <typename PointArray, class Integrand>
inline void
DOMAIN_WITH_TEMPLATES :: integrate (PointArray& array, const BasisFunctionExpr < Integrand > & expr)
    {
      if (this->recompute == true)
         {
           this->precomputeDerivatives(*this,BasisSet_1);
#if 1
           if (0 < numberOfSecondAnsatzFunctions)
              {
                this->precomputeDerivatives(*this,BasisSet_2);
              }
#endif
           this->recompute = false;
         }
      enum {LoopSize=Integrand::hasV+Integrand::hasW+Integrand::hasP};
      Compile_Time_Loop<LoopSize>::gaussianQuadrature(*this,array,expr);
    }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template <class FunctionExpression>
inline void
DOMAIN_WITH_TEMPLATES :: Set (const FunctionExpression & expr)
    {
      BasisFunctions::Set_BFN(BasisSet_1,expr, this->gaussianPoints);
    }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template <class FunctionExpression>
inline void
 DOMAIN_WITH_TEMPLATES:: Set (BasisFunctionSet setNumber, const FunctionExpression& expr)
    {
      BasisFunctions::Set_BFN (setNumber,expr,this->gaussianPoints);
    } 
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template <class FunctionExpression>
inline void
DOMAIN_WITH_TEMPLATES :: Set (const FunctionExpression & expr, const unsigned int cnr_1, const unsigned int cnr_2 )
    {
      BasisFunctions::Set_BFN(BasisSet_1,expr, cnr_1, cnr_2, this->gaussianPoints);
    }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
inline void
DOMAIN_WITH_TEMPLATES:: Point_ ()
    {
     
      for (int pointIter = 0; pointIter<init_dimension; ++ pointIter)
         {
           this->coordinatesOfVertices[cornersOfElement*init_dimension + pointIter] = this->stress[pointIter];
         }
      this->gaussianPoints[this->numberGaussianPoints] = this->stress;
      this->eval (this->gaussianPoints);
      this->point(*this, this->numberGaussianPoints);
      this->precomputeDerivativesAdditional(*this,BasisSet_1);
      if (1 < this->numberOfBasisFunctionSets_)
         {
           this->precomputeDerivativesAdditional(*this,BasisSet_2);
         }
    }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template <class PointArrayTYPE>
inline void
DOMAIN_WITH_TEMPLATES:: Point_ (const PointArrayTYPE& pointArray)
    {
      for (int pointIter = 0; pointIter<init_dimension; ++ pointIter)
         {
           this->coordinatesOfVertices[cornersOfElement*init_dimension + pointIter] = pointArray[pointIter];
         }
      this->compute_Barycentric_Coordinates(*this);
      this->gaussianPoints[this->numberGaussianPoints] = this->barys;
      this->eval (this->gaussianPoints);
      this->point(*this, this->numberGaussianPoints);
      this->precomputeDerivativesAdditional(*this,BasisSet_1);
      if (1 < this->numberOfBasisFunctionSets_)
         {
           this->precomputeDerivativesAdditional(*this,BasisSet_2);
         }
    }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template <class PointArrayTYPE>
inline bool
 DOMAIN_WITH_TEMPLATES:: Point_Test(const PointArrayTYPE& pointArray)
    {
      eTYPE test_bary = 0.;
      bool is_in_element = true;
      Point_ ( pointArray );
      for (unsigned int cornerIter = 0; cornerIter<cornersOfElement; ++cornerIter)
         {
           test_bary += this->barys[cornerIter+4];
           std::cout << "  " << this->barys[cornerIter+4] ;
         }
       std::cout << std::endl;
      if (1.e-20 <= test_bary - 1. )
         {
           is_in_element = false;
         }
#if 1
      std::cout << "point lies ";
      if ( is_in_element == false )
         {
           std::cout << "NOT ";
         }
      std::cout << "in the element!!" << std::endl;
#endif     
      return is_in_element;
    }
// ---------------------------------------------------------------------------
TEMPLATES_OF_CLASS_DOMAIN
template <class PointArrayTYPE>
inline DOMAIN_WITH_TEMPLATES& 
 DOMAIN_WITH_TEMPLATES::operator()(const PointArrayTYPE & arrayOfPointCoordinates, const int numberPoints)
    {
      const int pointCopySize = init_dimension * numberPoints;
   // Copying the point   
      for (int pointIter = 0 ; pointIter < pointCopySize ; ++pointIter ) 
         {
           this->coordinatesOfVertices[pointIter] = arrayOfPointCoordinates[pointIter];
         }
   // Applying the transformation
      for (int gaussIter = 0; gaussIter<this->numberGaussianPoints; ++gaussIter)
         {
            this->point(*this,gaussIter);
         }
      for (int gaussIter = 0; gaussIter<this->numberGaussianPoints; ++gaussIter)
         {
           this->extensionValues[gaussIter][trafoFactor] *= this->gaussianPoints[gaussIter][weight];
         }
      this->recompute = true;
      return (*this);
    }
//==============================================================================

