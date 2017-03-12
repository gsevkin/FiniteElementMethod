//==============================================================================
//
//  $Id: Gauss_Points.h,v 1.24 2006/07/14 10:51:30 jochen Exp $
//
//==============================================================================
template <int p1, int p2, int p3, int p4, int p5>
struct PyramidVolume
   {
     template <class Element>
     static inline typename Element::bTYPE
      compute(const Element& element);
   };
//-----------------------------------------------------------------------------
template <int p1, int p2, int p3, int p4>
struct TetrahedronVolume
   {
     template <class Element>
     static inline typename Element::bTYPE
      compute(const Element& element);
   };
//-----------------------------------------------------------------------------
template <int p1, int p2, int p3>
struct TriangleVolume
   {
     template <class Element>
     static inline typename Element::bTYPE
      compute(const Element& element);
   };
//-----------------------------------------------------------------------------
template <typename TYPE, IntegrationAccuracy t_, ElementType et>
class Gaussian_Points;
//-----------------------------------------------------------------------------
template <typename TYPE, int numberOfGaussianPoints>
struct Gauss_Init
   {
     protected:
       TYPE **gaussianPoints,
             *stress;
     public: 
     enum {numberGaussianPoints = numberOfGaussianPoints} ;
     Gauss_Init();
     ~Gauss_Init();
     inline TYPE** getGaussianPoints() const;
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss3, cuboid>
  : public Gauss_Init <TYPE,27>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss2, cuboid>
  : public Gauss_Init <TYPE,8>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss1, cuboid>
  : public Gauss_Init <TYPE,1>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss2, pyramid>
  : public Gauss_Init <TYPE,8>
   {
     Gaussian_Points();
   };
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss1, pyramid>
  : public Gauss_Init <TYPE,2>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss3, prism>
  : public Gauss_Init <TYPE,21>
   {
     Gaussian_Points();
   };
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss2, prism>
  : public Gauss_Init <TYPE,6>
   {
     Gaussian_Points();
   };
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss1, prism >
  : public Gauss_Init <TYPE,1>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss2, tetrahedron>
  : public Gauss_Init <TYPE,4>
   {
     Gaussian_Points();
   };
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss1, tetrahedron>
  : public Gauss_Init <TYPE,1>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss3, quadrangle>
  : public Gauss_Init <TYPE,9>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss2, quadrangle>
  : public Gauss_Init <TYPE,4>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss1, quadrangle>
  : public Gauss_Init <TYPE,1>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss3, triangle>
  : public Gauss_Init <TYPE,7>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss2, triangle>
  : public Gauss_Init <TYPE,3>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss1, triangle>
  : public Gauss_Init <TYPE,1>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss3, interval>
  : public Gauss_Init <TYPE,3>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss2 , interval>
  : public Gauss_Init <TYPE,2>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE>
struct Gaussian_Points<TYPE, Gauss1, interval>
  : public Gauss_Init <TYPE,1>
   {
     Gaussian_Points();
   } ;
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
class CornerManagement 
   {
     protected:
       TYPE *barys,
            *coordinatesOfVertices;
     public:
       CornerManagement();
       ~CornerManagement(); 
       int getNumberOfCorners() const {return numberCorners;}
       inline TYPE* getCornerList() const;
       inline TYPE* getArrayOfBarycentricCoordinates() const;
       inline TYPE edge_length (const unsigned int cnr_1, const unsigned int cnr_2) const;
       inline TYPE edge_(SpaceDirection direction) const;
   };
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners, ElementType elType>
struct Corner_Classes;
//-----------------------------------------------------------------------------
template <typename TYPE, int numberCorners>
struct Corner_Classes <TYPE,D3,numberCorners,cuboid> : CornerManagement<TYPE,D3,numberCorners>
   {
     enum {numberOfCorners = numberCorners};
     template <class Element>
     inline void
     compute_Barycentric_Coordinates (const Element& element);
   };
//-----------------------------------------------------------------------------
template <typename TYPE, int numberCorners>
struct Corner_Classes <TYPE,D3,numberCorners,pyramid> : CornerManagement<TYPE,D3,numberCorners>
   {
     enum {numberOfCorners = numberCorners};
     template <class Element>
     inline void
     compute_Barycentric_Coordinates (const Element& element);
   };
//-----------------------------------------------------------------------------
template <typename TYPE, int numberCorners>
struct Corner_Classes <TYPE,D3,numberCorners,prism> : CornerManagement<TYPE,D3,numberCorners>
   {
     enum {numberOfCorners = numberCorners};
     template <class Element>
     inline void
     compute_Barycentric_Coordinates (const Element& element);
   };
//-----------------------------------------------------------------------------
template <typename TYPE, int numberCorners>
struct Corner_Classes <TYPE,D3,numberCorners,tetrahedron> : CornerManagement<TYPE,D3,numberCorners>
   {
     enum {numberOfCorners = numberCorners};
     template <class Element>
     inline void
     compute_Barycentric_Coordinates (const Element& element);
   };
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
struct Corner_Classes <TYPE,dimension,numberCorners,quadrangle> : CornerManagement<TYPE,dimension,numberCorners>
   {
     enum {numberOfCorners = numberCorners};
     template <class Element>
     inline void
     compute_Barycentric_Coordinates (const Element& element);
   };
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
struct Corner_Classes <TYPE,dimension,numberCorners,triangle> : CornerManagement<TYPE,dimension,numberCorners>
   {
     enum {numberOfCorners = numberCorners};
     template <class Element>
     inline void
     compute_Barycentric_Coordinates (const Element& element);
   };
//-----------------------------------------------------------------------------
template <typename TYPE, DIMENSION dimension, int numberCorners>
struct Corner_Classes <TYPE,dimension,numberCorners,interval> : CornerManagement<TYPE,dimension,numberCorners>
   {
     enum {numberOfCorners = numberCorners};
     template <class Element>
     inline void
     compute_Barycentric_Coordinates (const Element& element);
   };
//-----------------------------------------------------------------------------
//==============================================================================
