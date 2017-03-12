//==============================================================================
//
//  $Id: Elements.h,v 1.62 2006/07/25 13:52:27 jochen Exp $
//
//==============================================================================
//-----------------------------------------------------------------------------
namespace ELEMENTS
   {
    //-----------------------------------------------------------------------------
    // Defining the mappings ...
    //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (   P_0() +
              ( P_1() - P_0() ) * _U() +
              ( P_3() - P_0() + (P_0() - P_1() + P_2() - P_3() )*_U() ) *_V() +
              ( ( P_4() - P_0()) + (P_0() - P_1() - P_4() + P_5() )*_U() +
                ( P_0()- P_3() - P_4() + P_7() +
                  ( P_1()-P_0()-P_2()+P_3()+P_4()-P_5()+P_6()-P_7() )*
                                                _U() )*_V() )*_W()
            )
              Hexahedron_Transformation;
    //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (   P_0() +
              ( P_1() - P_0() ) * _U() +
              ( P_2() - P_0() + (P_0() - P_1() + P_3() - P_2() )*_U() ) *_V() +
              ( ( P_4() - P_0()) + (P_0() - P_1() - P_4() + P_5() )*_U() +
                ( P_0()- P_2() - P_4() + P_6() +
                  ( P_1()-P_0()-P_3()+P_2()+P_4()-P_5()+P_7()-P_6() )*
                                                _U() )*_V() )*_W()
            )
              Hexahedron_Transformation_2;
   //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (
              P_0() +
              ( P_1() - P_0()) * _U() +
              ( P_3() - P_0()) * _V() +
              ( P_4() - P_0()) * _W()
            )
              Cuboid_Transformation;
   //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (
              P_0() +
              ( P_1() - P_0()) * _U() +
              ( P_2() - P_0()) * _V() +
              ( P_3() - P_0()) * _W()
            )
              Cuboid_Transformation_2;
   //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (
              P_0() +
              ( P_1() - P_0()) * _U() +
              ( P_2() - P_0()) * _V() +
              ( P_3() - P_0()) * _W()
            )
              Tetrahedron_Transformation;
   //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (
              P_0() +
              ( P_1() - P_0() ) * _U() +
              ( P_3() - P_0() +
               ( P_0() - P_1() + P_2() - P_3() ) * _U()) * _V() +
               ( P_4() - P_0() ) * _W()
            )
              Pyramid_Transformation;
   //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (
              P_0() +
              ( P_1() - P_0() ) * _U() +
              ( P_2() - P_0() ) * _V() +
              ( P_3() - P_0() + (P_4() + P_0() - P_1() - P_3() )* _U() +
                                (P_5() + P_0() - P_2() - P_3() )* _V() ) * _W()
            )
              Prism_Transformation;
   //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (
                P_0() +
              ( P_1() - P_0() ) * _U() +
              ( P_2() - P_0() ) * _V()
            )
              Triangle_Transformation;
     //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (
              P_0() +
              ( P_1() - P_0() ) * _U() +
              ( P_3() - P_0() + ( P_0() - P_1() + P_2() - P_3() )* _U() ) * _V()
            )
              Quadrangle_Transformation;
   //-----------------------------------------------------------------------------
     Define_Element_Transformation
            (
                P_0() +
              ( P_1() - P_0() ) * _U()
            )
              Interval_Transformation;
//--------------------------------------------------------------------------------
/*--------------------------------------------------------------------------------
   _Domain_ < \
   	number of vertices,
  	dimension D1 / D2 / D3 ,
   	mapping from reference to a general element,
  	interior / boundary,
  	number of basis functions of set 1,
  	number of basis functions of set 2,
        geometry of reference element,
   	Gaussian quadrature formula, Gauss1 / Gauss2 / Gauss3
   	double / complex<double>,
        STLVector / Array
           >
*/
//--------------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Hexahedron_
       : public _Domain_ <8,D3,Hexahedron_Transformation,interior,8,0,cuboid,accuracy,dTYPE,sTYPE>
        {
          _Hexahedron_()
             {
               this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) * ( 1. - Z_(1) ) ) ;
               this->Set (        X_(1)   * ( 1. - Y_(1) ) * ( 1. - Z_(1) ) ) ;
               this->Set (        X_(1)   *        Y_(1)   * ( 1. - Z_(1) ) ) ;
               this->Set ( ( 1. - X_(1) ) *        Y_(1)   * ( 1. - Z_(1) ) ) ;
               this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) *        Z_(1)   ) ;
               this->Set (        X_(1)   * ( 1. - Y_(1) ) *        Z_(1)   ) ;
               this->Set (        X_(1)   *        Y_(1)   *        Z_(1)   ) ;
               this->Set ( ( 1. - X_(1) ) *        Y_(1)   *        Z_(1)   ) ;
             }
        } ;
   // Simple type for the standard element
     typedef _Hexahedron_<> Hexahedron;
//--------------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Hexahedron_2
       : public _Domain_ <8,D3,Hexahedron_Transformation_2,interior,8,0,cuboid,accuracy,dTYPE,sTYPE>
        {
          _Hexahedron_2()
             {
               this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) * ( 1. - Z_(1) ) ) ;
               this->Set (        X_(1)   * ( 1. - Y_(1) ) * ( 1. - Z_(1) ) ) ;
               this->Set ( ( 1. - X_(1) ) *        Y_(1)   * ( 1. - Z_(1) ) ) ;
               this->Set (        X_(1)   *        Y_(1)   * ( 1. - Z_(1) ) ) ;
               this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) *        Z_(1)   ) ;
               this->Set (        X_(1)   * ( 1. - Y_(1) ) *        Z_(1)   ) ;
               this->Set ( ( 1. - X_(1) ) *        Y_(1)   *        Z_(1)   ) ;
               this->Set (        X_(1)   *        Y_(1)   *        Z_(1)   ) ;
             }
        } ;
   // Simple type for the standard element
     typedef _Hexahedron_2<> Hexahedron_2;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Hexahedron_Boundary_
       : public _Domain_ <8,D3,Hexahedron_Transformation,boundary,8,0,quadrangle,accuracy,dTYPE,sTYPE>
        {
          _Hexahedron_Boundary_()
             {
#if 0
               this->Set ( ( 1 - X_() ) * ( 1 - Y_() ) * ( 1 - Z_() ) ) ;
               this->Set (       X_()   * ( 1 - Y_() ) * ( 1 - Z_() ) ) ;
               this->Set (       X_()   *       Y_()   * ( 1 - Z_() ) ) ;
               this->Set ( ( 1 - X_() ) *       Y_()   * ( 1 - Z_() ) ) ;
#else
               this->Set ( ( 1 - X_() ) * ( 1 - Y_() ) * ( 1 - Z_() ) ) ;
               this->Set (       X_()   * ( 1 - Y_() ) * ( 1 - Z_() ) ) ;
               this->Set (       X_()   *       Y_()   * ( 1 - Z_() ) ) ;
               this->Set ( ( 1 - X_() ) *       Y_()   * ( 1 - Z_() ) ) ;
               this->Set ( ( 1 - X_() ) * ( 1 - Y_() ) *       Z_()   ) ;
               this->Set (       X_()   * ( 1 - Y_() ) *       Z_()   ) ;
               this->Set (       X_()   *       Y_()   *       Z_()   ) ;
               this->Set ( ( 1 - X_() ) *       Y_()   *       Z_()   ) ;
#endif
             }
        } ;
   // Simple type for the standard element
     typedef _Hexahedron_Boundary_<> Hexahedron_Boundary;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Cuboid_
       : public _Domain_<8,D3,Cuboid_Transformation,interior,8,0,cuboid,accuracy,dTYPE,sTYPE>
        {
          _Cuboid_()
             {
               this->Set ( ( 1 - X_() ) * ( 1 - Y_() ) * ( 1 - Z_() ) ) ;
               this->Set (       X_()   * ( 1 - Y_() ) * ( 1 - Z_() ) ) ;
               this->Set (       X_()   *       Y_()   * ( 1 - Z_() ) ) ;
               this->Set ( ( 1 - X_() ) *       Y_()   * ( 1 - Z_() ) ) ;
               this->Set ( ( 1 - X_() ) * ( 1 - Y_() ) *       Z_()   ) ;
               this->Set (       X_()   * ( 1 - Y_() ) *       Z_()   ) ;
               this->Set (       X_()   *       Y_()   *       Z_()   ) ;
               this->Set ( ( 1 - X_() ) *       Y_()   *       Z_()   ) ;
             }
        } ;
   // Simple type for the standard element
     typedef _Cuboid_<> Cuboid;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Cuboid_Quadratic_
       : public _Domain_<8,D3,Cuboid_Transformation,interior,27,0,cuboid,accuracy,dTYPE,sTYPE>
        {
          _Cuboid_Quadratic_()
             {
    		this->Set((1-X_())*(1-2*X_()) * (1-Y_())*(1-2*Y_()) * (1-Z_())*(1-2*Z_()));
    		this->Set( 4*(1-X_()) * X_()  * (1-Y_())*(1-2*Y_()) * (1-Z_())*(1-2*Z_()));
    		this->Set(   X_() *(2*X_()-1) * (1-Y_())*(1-2*Y_()) * (1-Z_())*(1-2*Z_()));

    		this->Set((1-X_())*(1-2*X_()) * 4*(1-Y_()) * Y_()   * (1-Z_())*(1-2*Z_()));
    		this->Set( 4*(1-X_()) * X_()  * 4*(1-Y_()) * Y_()   * (1-Z_())*(1-2*Z_()));
    		this->Set(   X_() *(2*X_()-1) * 4*(1-Y_()) * Y_()   * (1-Z_())*(1-2*Z_()));

    		this->Set((1-X_())*(1-2*X_()) *    Y_() *(2*Y_()-1) * (1-Z_())*(1-2*Z_()));
    		this->Set( 4*(1-X_()) * X_()  *    Y_() *(2*Y_()-1) * (1-Z_())*(1-2*Z_()));
    		this->Set(   X_() *(2*X_()-1) *    Y_() *(2*Y_()-1) * (1-Z_())*(1-2*Z_()));

    		this->Set((1-X_())*(1-2*X_()) * (1-Y_())*(1-2*Y_()) * 4.*(1-Z_()) * Z_() );
    		this->Set( 4*(1-X_()) * X_()  * (1-Y_())*(1-2*Y_()) * 4.*(1-Z_()) * Z_() );
    		this->Set(   X_() *(2*X_()-1) * (1-Y_())*(1-2*Y_()) * 4.*(1-Z_()) * Z_() );

    		this->Set((1-X_())*(1-2*X_()) * 4*(1-Y_()) * Y_()   * 4.*(1-Z_()) * Z_() );
    		this->Set( 4*(1-X_()) * X_()  * 4*(1-Y_()) * Y_()   * 4.*(1-Z_()) * Z_() );
    		this->Set(   X_() *(2*X_()-1) * 4*(1-Y_()) * Y_()   * 4.*(1-Z_()) * Z_() );

    		this->Set((1-X_())*(1-2*X_()) *    Y_() *(2*Y_()-1) * 4.*(1-Z_()) * Z_() );
    		this->Set( 4*(1-X_()) * X_()  *    Y_() *(2*Y_()-1) * 4.*(1-Z_()) * Z_() );
    		this->Set(   X_() *(2*X_()-1) *    Y_() *(2*Y_()-1) * 4.*(1-Z_()) * Z_() );

    		this->Set((1-X_())*(1-2*X_()) * (1-Y_())*(1-2*Y_()) *    Z_() *(2*Z_()-1));
    		this->Set( 4*(1-X_()) * X_()  * (1-Y_())*(1-2*Y_()) *    Z_() *(2*Z_()-1));
    		this->Set(   X_() *(2*X_()-1) * (1-Y_())*(1-2*Y_()) *    Z_() *(2*Z_()-1));

    		this->Set((1-X_())*(1-2*X_()) * 4*(1-Y_()) * Y_()   *    Z_() *(2*Z_()-1));
    		this->Set( 4*(1-X_()) * X_()  * 4*(1-Y_()) * Y_()   *    Z_() *(2*Z_()-1));
    		this->Set(   X_() *(2*X_()-1) * 4*(1-Y_()) * Y_()   *    Z_() *(2*Z_()-1));

    		this->Set((1-X_())*(1-2*X_()) *    Y_() *(2*Y_()-1) *    Z_() *(2*Z_()-1));
    		this->Set( 4*(1-X_()) * X_()  *    Y_() *(2*Y_()-1) *    Z_() *(2*Z_()-1));
    		this->Set(   X_() *(2*X_()-1) *    Y_() *(2*Y_()-1) *    Z_() *(2*Z_()-1));

             }
        } ;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct Vector_Cuboid_
       : public Vectorial_Element<8,D3,Cuboid_Transformation,12,cuboid,accuracy,interior,dTYPE,sTYPE>
        {
          Vector_Cuboid_()
             {
               this->Set ( fVec(  (1 - Y_()) * (1 - Z_()), Null_(), Null_()), 0,1 ) ;
               this->Set ( fVec(  Y_() * (1 - Z_())      , Null_(), Null_()), 3,2 ) ;
               this->Set ( fVec(  (1 - Y_()) * Z_()      , Null_(), Null_()), 4,5 ) ;
               this->Set ( fVec(   Y_() * Z_()           , Null_(), Null_()), 7,6 ) ;

               this->Set ( fVec(  Null_(), (1 - X_()) * (1 - Z_()), Null_()), 0,3 ) ;
               this->Set ( fVec(  Null_(),  X_() * (1 - Z_())     , Null_()), 1,2 ) ;
               this->Set ( fVec(  Null_(),  (1 - X_()) * Z_()     , Null_()), 4,7 ) ;
               this->Set ( fVec(  Null_(),   X_() * Z_()          , Null_()), 5,6 ) ;

               this->Set ( fVec(  Null_(), Null_(), (1 - X_()) * (1 - Y_())), 0,4 ) ;
               this->Set ( fVec(  Null_(), Null_(),  X_() * (1 - Y_()))     , 1,5 ) ;
               this->Set ( fVec(  Null_(), Null_(), (1 - X_()) * Y_())      , 3,7 ) ;
               this->Set ( fVec(  Null_(), Null_(),  X_() * Y_())           , 2,6 ) ;
             }
        } ;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector >
     struct _Cuboid_Boundary_
       : public _Domain_ <8,D3,Cuboid_Transformation,boundary,4,0,quadrangle,accuracy,double,STLVector>
        {
          _Cuboid_Boundary_()
             {
               this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) ) ;
               this->Set (        X_(1)   * ( 1. - Y_(1) ) ) ;
               this->Set (        X_(1)   *        Y_(1)   ) ;
               this->Set ( ( 1. - X_(1) ) *        Y_(1)   ) ;
             }
        } ;
     //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Cuboid_Mixed_
       : public _Domain_<8,D3,Cuboid_Transformation,interior,2,9,cuboid,accuracy,dTYPE,sTYPE>
         {
           _Cuboid_Mixed_()
              {
                this->Set ( BasisSet_1,X_(1)) ;
                this->Set ( BasisSet_1,Y_(1)) ;
#if 0
                this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) * ( 1. - Z_(1) ) ) ;
                this->Set (        X_(1)   * ( 1. - Y_(1) ) * ( 1. - Z_(1) ) ) ;
                this->Set (        X_(1)   *        Y_(1)   * ( 1. - Z_(1) ) ) ;
                this->Set ( ( 1. - X_(1) ) *        Y_(1)   * ( 1. - Z_(1) ) ) ;
                this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) *        Z_(1)   ) ;
                this->Set (        X_(1)   * ( 1. - Y_(1) ) *        Z_(1)   ) ;
                this->Set (        X_(1)   *        Y_(1)   *        Z_(1)   ) ;
                this->Set ( ( 1. - X_(1) ) *        Y_(1)   *        Z_(1)   ) ;
#endif     
                this->Set ( BasisSet_2,( 1. - X_(1) ) * ( 1. - Y_(1) ) * ( 1. - Z_(1) ) + 1. ) ;
                this->Set ( BasisSet_2,       X_(1)   * ( 1. - Y_(1) ) * ( 1. - Z_(1) ) + 1. ) ;
                this->Set ( BasisSet_2,       X_(1)   *        Y_(1)   * ( 1. - Z_(1) ) + 1. ) ;
                this->Set ( BasisSet_2,( 1. - X_(1) ) *        Y_(1)   * ( 1. - Z_(1) ) + 1. ) ;
                this->Set ( BasisSet_2,( 1. - X_(1) ) * ( 1. - Y_(1) ) *        Z_(1)   + 1. ) ;
                this->Set ( BasisSet_2,       X_(1)   * ( 1. - Y_(1) ) *        Z_(1)   + 1. ) ;
                this->Set ( BasisSet_2,1. + Y_(1) + Z_(1) + 1. ) ;
                this->Set ( BasisSet_2,1. + Y_(1) + Z_(1) + 1. ) ;
                this->Set ( BasisSet_2,1. + Y_(1) + Z_(1) + 1. ) ;
              }
        };
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Tetrahedron_
       : public Simple_Element<4,D3,Tetrahedron_Transformation,4,tetrahedron,accuracy,dTYPE,sTYPE>
        {
          _Tetrahedron_()
            {
              this->Set ( 1. - X_(1) - Y_(1) - Z_(1) ) ;
              this->Set ( X_(1)  ) ;
              this->Set ( Y_(1)  ) ;
              this->Set ( Z_(1)  ) ;
            }
        } ;
   // Simple type for the standard element
      typedef _Tetrahedron_ <> Tetrahedron;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct Vector_Tetrahedron_
       : public Vectorial_Element<4,D3,Tetrahedron_Transformation,6,tetrahedron,accuracy,interior,dTYPE,sTYPE>
        {
          Vector_Tetrahedron_()
             {
               this->Set ( fVec( ( 1 - Y_()  - Z_() ) , X_()  , X_() ), 0,1);
               this->Set ( fVec(  Y_()  , 1 - X_() - Z_(), Y_() ), 0, 2); 
               this->Set ( fVec(  Z_(), Z_() , 1 - X_() - Y_() ), 0, 3);
               this->Set ( fVec( Y_(),  - X_()  , Null_() ), 1, 2);
               this->Set ( fVec( Z_(), Null_(), - X_() ), 1, 3); 
               this->Set ( fVec( Null_(), Z_(), - Y_() ), 2, 3);
             }
        } ;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Tetrahedron_Boundary_
       : public _Domain_ <4,D3,Tetrahedron_Transformation,boundary,3,0,triangle,accuracy,dTYPE,sTYPE>
        {
          _Tetrahedron_Boundary_()
            {
              this->Set ( 1. - X_(1) - Y_(1) ) ;
              this->Set ( X_(1) ) ;
              this->Set ( Y_(1) ) ;
            }
        } ;
   // Simple type for the standard element
      typedef _Tetrahedron_Boundary_ <> Tetrahedron_Boundary;
    //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Pyramid_
       : public _Domain_<5,D3,Pyramid_Transformation,interior,5,0,pyramid,accuracy,dTYPE,sTYPE>
        {
          _Pyramid_()
             {
               this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) * ( 1. -  Z_(1) ) ) ;
               this->Set (        X_(1)   * ( 1. - Y_(1) ) * ( 1. -  Z_(1) ) ) ;
               this->Set (        X_(1)   *        Y_(1)   * ( 1. -  Z_(1) ) ) ;
               this->Set ( ( 1. - X_(1) ) *        Y_(1)   * ( 1. -  Z_(1) ) ) ;
               this->Set (                                       Z_(1)   ) ;
             }
        } ;
   // Simple type for the standard element
      typedef _Pyramid_ <> Pyramid;
    //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Prism_
       : public _Domain_<6,D3,Prism_Transformation,interior,6,0,prism,accuracy,dTYPE,sTYPE>
        {
          _Prism_()
             {
               this->Set ( ( 1. - X_(1)  - Y_(1) ) * ( 1. -  Z_(1) ) ) ;
               this->Set (        X_(1)   *  ( 1. -  Z_(1) ) ) ;
               this->Set (        Y_(1)   *  ( 1. -  Z_(1) ) ) ;
               this->Set ( ( 1. - X_(1)  - Y_(1) ) *  Z_(1)  ) ;
               this->Set (        X_(1)   *  Z_(1) ) ;
               this->Set (        Y_(1)   *  Z_(1) ) ;
             }
        } ;
   // Simple type for the standard element
      typedef _Prism_ <> Prism;
     //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Triangle_
       : public Simple_Element<3,D2,Triangle_Transformation,3,triangle,accuracy,dTYPE,sTYPE>
        {
          _Triangle_()
             {
               this->Set ( 1. - X_(1) - Y_(1) ) ;
               this->Set ( X_(1) ) ;
               this->Set ( Y_(1) ) ;
             }
        } ;
   // Simple type for the standard element
      typedef _Triangle_<> Triangle;
     //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Triangle_Quadratic_
       : public _Domain_<3,D2,Triangle_Transformation,interior,6,0,triangle,accuracy,dTYPE,sTYPE>
        {
          _Triangle_Quadratic_()
             {
               this->Set ( 2 * (1 - X_(1) - Y_(1) )  * ( 0.5 - X_(1) - Y_(1) ) ) ;
               this->Set ( 2 * ( X_(1)  * (  X_(1) - 0.5 ) ) ) ;
               this->Set ( 2 * ( Y_(1)  * (  Y_(1) - 0.5 ) ) ) ;
               this->Set ( 4 * (1 - X_(1) - Y_(1) ) * X_(1) ) ;
               this->Set ( 4 * ( X_(1) * Y_(1) ) ) ;
               this->Set ( 4 * (1 - X_(1) - Y_(1) ) * Y_(1) ) ;
             }
        } ;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct Vector_Triangle_
       : public Vectorial_Element<3,D2,Triangle_Transformation,3,triangle,accuracy,interior,dTYPE,sTYPE>
        {
          Vector_Triangle_()
             {
               this->Set ( fVec( 1 - Y_() ,  X_()     ), 0, 1);
               this->Set ( fVec(   - Y_() ,  X_()     ), 1, 2);
               this->Set ( fVec(   - Y_() ,  X_() - 1 ), 2, 0);
             }
        } ;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Quadrangle_
       : public _Domain_<4,D2,Quadrangle_Transformation,interior,4,0,quadrangle,accuracy,dTYPE,sTYPE>
        {
          _Quadrangle_()
             {
               this->Set ( ( 1 - X_(1) ) * ( 1 - Y_(1) ) ) ;
               this->Set ( X_(1) * ( 1 - Y_(1) ) ) ;
               this->Set ( X_(1) * Y_(1) ) ;
               this->Set ( ( 1 - X_(1) ) * Y_(1) ) ;
             }
        } ;
   // Simple type for the standard element
      typedef _Quadrangle_ <> Quadrangle;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct Vector_Quadrangle_
       : public Vectorial_Element<4,D2,Quadrangle_Transformation,4,quadrangle,accuracy,interior,dTYPE,sTYPE>
        {
          Vector_Quadrangle_()
             {
               this->Set ( fVec( ( 1 - Y_() ), Null_() ), 0, 1);
               this->Set ( fVec(       Y_()  , Null_() ), 3, 2); 
               this->Set ( fVec(    Null_()  , 1 -X_() ), 0, 3);
               this->Set ( fVec(    Null_()  ,    X_() ), 1, 2);
             }
        } ;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Triangle3D_
       : public _Domain_<3,D3,Triangle_Transformation,boundary,3,0,triangle,accuracy,dTYPE,sTYPE>
        {
          _Triangle3D_ ()
             {
               this->Set ( 1. - X_(1) - Y_(1) ) ;
               this->Set ( X_(1)  ) ;
               this->Set ( Y_(1)  ) ;
             }
        } ;
     //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Quadrangle3D_
       : public _Domain_<4,D3,Quadrangle_Transformation,boundary,4,0,quadrangle,accuracy,dTYPE,sTYPE>
        {
          _Quadrangle3D_()
             {
               this->Set ( ( 1. - X_(1) ) * ( 1. - Y_(1) ) ) ;
               this->Set ( X_(1)* ( 1. - Y_(1) ) ) ;
               this->Set ( X_(1)* Y_(1) ) ;
               this->Set ( ( 1. - X_(1) ) * Y_(1)) ;
             }
        } ;
   //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct Vector_Edge_2D_
       : public _Domain_<2,D2,Interval_Transformation,boundary,1,0,interval,accuracy,dTYPE,sTYPE,double,vectorBasisFunctions>
        {
          Vector_Edge_2D_()
             {
                 this->Set ( fVec(c_(1)), 0,1 );
             }
        } ;
     //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Quadrangle_Boundary_
       : public _Domain_<4,D2,Quadrangle_Transformation,boundary,4,0,interval,accuracy,dTYPE,sTYPE>
        {
          _Quadrangle_Boundary_()
             {
               this->Set ( ( 1 - X_() ) * ( 1 - Y_() ) ) ;
               this->Set (       X_()   * ( 1 - Y_() ) ) ;
               this->Set (       X_()   *       Y_()   ) ;
               this->Set ( ( 1 - X_() ) *       Y_()   ) ;
             }
        } ;
     //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Interval_
       : public _Domain_<2,D1,Interval_Transformation,interior,2,0,interval,accuracy,dTYPE,sTYPE>
        {
          _Interval_()
             {
               this->Set ( 1. - X_() ) ;
               this->Set ( X_() ) ;
             }
        } ;
      typedef _Interval_ <> Interval;
     //-----------------------------------------------------------------------------
     template < IntegrationAccuracy accuracy = Gauss2, typename dTYPE = double, stencilType sTYPE = STLVector>
     struct _Edge2D_
       : public _Domain_<2,D2,Triangle_Transformation,boundary,2,0,interval,accuracy,dTYPE,sTYPE>
        {
          _Edge2D_()
             {
               this->Set ( 1. - X_(1) ) ;
               this->Set ( X_(1) ) ;
             }
        } ;
     //-----------------------------------------------------------------------------
   }
//==============================================================================
