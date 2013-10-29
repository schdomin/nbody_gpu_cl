#ifndef CVECTOR_CUH_
#define CVECTOR_CUH_

#include <fstream> //ds << operator
#include <math.h>  //ds sqrt



namespace NBody
{

//ds basic structure
class CVector
{

//ds not accessible from outside
private:

    //ds three elements
    double m_dElement0;
    double m_dElement1;
    double m_dElement2;

//ds ctor/dtor
public:

    //ds custom constructor
    CVector( const double& p_dElement0, const double& p_dElement1, const double& p_dElement2 );

    //ds default constructor - we need to have this since operator overloading for pointers is not allowed..
    CVector( );

    //ds default destructor
    ~CVector( );


//ds operators
public:

    //ds r/w indexing - ( ) is used instead of [ ] to mark the difference that this is not readonly - CAREFUL this allows manipulation of the data
    double& operator( )( const unsigned int& p_uIndex );

    //ds readonly indexing
    double operator[ ]( const unsigned int& p_uIndex ) const;

    //ds setting
    void operator=( const CVector& p_cRightHandSide );

    //ds adding and assign
    CVector& operator+=( const CVector& p_cRightHandSide );

    //ds dividing and assign
    CVector& operator/=( const double& p_dRightHandSide );

    //ds multiplication and assign
    CVector& operator*=( const double& p_dRightHandSide );

    //ds dividing
    CVector& operator/( const double& p_dRightHandSide );

    //ds simple addition
    friend const CVector operator+( const CVector& p_cLeftHandSide, const CVector& p_cRightHandSide );

    //ds simple subtraction (needed for distance calculations)
    friend const CVector operator-( const CVector& p_cLeftHandSide, const CVector& p_cRightHandSide );

    //ds simple multiplication
    friend const CVector operator*( const double& p_dLeftHandSide, const CVector& p_cRightHandSide );

    //ds printing
    friend std::ostream& operator<<( std::ostream& p_cLeftHandSide, const CVector& p_cRightHandSide );

//ds static functions
public:

    static double absoluteValue( const CVector& p_cVector );
    static double absoluteValue( const double p_vecVector[3] );
    static const CVector crossProduct( const CVector& p_cVector1, const CVector& p_cVector2 );
    static const CVector crossProduct( const double p_vecVector1[3], const double p_vecVector2[3] );

};

} //ds namespace NBody



#endif //ds CVECTOR_CUH_
