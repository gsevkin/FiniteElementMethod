//==============================================================================
//
//  $Id: Mapping.h,v 1.45 2006/07/25 08:11:31 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension_, int _numberGaussianPoints>
class MappingManagement
   {
      public:
        enum{dimensionOfElement = dimension_};
      protected:
        enum{dimMap = dimension_ +1 , 
             dimVal = dimensionOfElement*dimensionOfElement, 
             dimExt = _extensionSize};
        TYPE ***mappingValues,
             *jacobianMatrixEntries,
             **extensionValues;
      public:
        MappingManagement();
        ~MappingManagement();
        inline int dimension() const {return (int)dimensionOfElement;}
        inline TYPE** getMappingValues(int gaussLevel) const ;
        inline TYPE* getValuesOfJacobianMatrix(int gaussLevel) const ;
        inline const TYPE& getDeterminant(int gaussLevel) const ;
        inline const TYPE& getNormalLength(int gaussLevel) const ;
        inline const TYPE& getTransformationFactorTimesWeight(int gaussLevel) const ;
   };
//-----------------------------------------------------------------------------
template <typename _TYPE, DIMENSION d, domainSecification int_bound, class TRAFO, int _numberGaussianPoints>
class Mappings ;
//-----------------------------------------------------------------------------
template <typename TYPE,class TRAFO, int _numberGaussianPoints>
class Mappings <TYPE,D3,interior,TRAFO,_numberGaussianPoints> : public MappingManagement<TYPE,D3,_numberGaussianPoints>
    {
      private: 
       const TRAFO trafo;
      public:
       Mappings() : trafo(TRAFO()) {} 
       template <class Element>
        inline void point ( const Element& element, const unsigned int gaussLevel );
   };
//-----------------------------------------------------------------------------
template <typename TYPE,class TRAFO, int _numberGaussianPoints>
class Mappings<TYPE,D2,interior,TRAFO,_numberGaussianPoints>  : public MappingManagement<TYPE,D2,_numberGaussianPoints> 
     {
      private:
       const TRAFO trafo;
      public:
       Mappings() : trafo(TRAFO()) {} 
       template <class Element>
        inline void point(const Element& element, const unsigned int gaussLevel );
   };
//-----------------------------------------------------------------------------
template <typename TYPE,class TRAFO, int _numberGaussianPoints>
class Mappings<TYPE,D1,interior,TRAFO,_numberGaussianPoints> : public MappingManagement<TYPE,D1,_numberGaussianPoints>
     {
      private:
       const TRAFO trafo;
      public:
       Mappings() : trafo(TRAFO()) {} 
       template <class Element>
        inline void point(const Element& element, const unsigned int gaussLevel );
     };
// -----------------------------------------------------------------------------
template <typename TYPE,class TRAFO, int _numberGaussianPoints>
class Mappings<TYPE, D2,boundary,TRAFO,_numberGaussianPoints> : public MappingManagement<TYPE,D2,_numberGaussianPoints> 
     {
      private: 
       const TRAFO trafo;
      public:
       Mappings() : trafo(TRAFO()) {} 
       template <class Element>
        inline void point ( const Element& element, const unsigned int gaussLevel );
     };
// -----------------------------------------------------------------------------
template <typename TYPE,class TRAFO, int _numberGaussianPoints>
class Mappings<TYPE,D3,boundary,TRAFO,_numberGaussianPoints> : public MappingManagement<TYPE,D3,_numberGaussianPoints> 
     {
      private: 
       const TRAFO trafo;
      public:
       Mappings() : trafo(TRAFO()) {} 
       template <class Element>
        inline void point ( const Element& element, const unsigned int gaussLevel );
      };
// -----------------------------------------------------------------------------
template <typename TYPE,class TRAFO, int _numberGaussianPoints>
class Mappings<TYPE,D3,boundaryPlain,TRAFO,_numberGaussianPoints> : public MappingManagement<TYPE,D3,_numberGaussianPoints> 
    {
      private: 
       const TRAFO trafo;
      public:
       Mappings() : trafo(TRAFO()) {} 
       template <class Element>
        inline void point(const Element& element, const unsigned int gaussLevel);
    };
// -----------------------------------------------------------------------------
template <typename TYPE,class TRAFO, int _numberGaussianPoints>
class Mappings<TYPE,D2,boundaryPlain,TRAFO,_numberGaussianPoints> : public MappingManagement<TYPE,D2,_numberGaussianPoints> 
    {
      private: 
       const TRAFO trafo;
      public:
       Mappings() : trafo(TRAFO()) {} 
       template <class Element>
        inline void point(const Element& element, const unsigned int gaussLevel);
    };
//==============================================================================
