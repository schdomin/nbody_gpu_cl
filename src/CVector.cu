#include "CVector.cuh"



namespace NBody
{

//ds ctor/dtor
//ds default constructor
CVector::CVector( const double& p_dElement0, const double& p_dElement1, const double& p_dElement2 ): m_dElement0( p_dElement0 ),
                                                                                                     m_dElement1( p_dElement1 ),
                                                                                                     m_dElement2( p_dElement2 )
{
    //ds nothing to do
}

//ds default constructor - we need to have this since operator overloading for pointers is not allowed..
CVector::CVector( ): m_dElement0( 0.0 ),
                     m_dElement1( 0.0 ),
                     m_dElement2( 0.0 )
{
    //ds nothing to do
}

//ds default destructor
CVector::~CVector( )
{
    //ds nothing to do
}

//ds operators
//ds r/w indexing - ( ) is used instead of [ ] to mark the difference that this is not a "real" array - CAREFUL this allows manipulation of the data
double& CVector::operator( )( const unsigned int& p_uIndex )
{
    //ds map the index operator to the element
    if     ( 0 == p_uIndex ){ return m_dElement0; }
    else if( 1 == p_uIndex ){ return m_dElement1; }
    else if( 2 == p_uIndex ){ return m_dElement2; }

    //ds if an index greater 2 is required throw an exception
    else
    {
        //ds TODO implement exceptions
        throw std::exception( );
    }
}

//ds readonly indexing - ( ) is used instead of [ ] to mark the difference that this is not a "real" array
double CVector::operator[ ]( const unsigned int& p_uIndex ) const
{
    //ds map the index operator to the element
    if     ( 0 == p_uIndex ){ return m_dElement0; }
    else if( 1 == p_uIndex ){ return m_dElement1; }
    else if( 2 == p_uIndex ){ return m_dElement2; }

    //ds if an index greater 2 is required throw an exception
    else
    {
        //ds TODO implement exceptions
        throw std::exception( );
    }
}

//ds setting
void CVector::operator=( const CVector& p_cRightHandSide )
{
    //ds get all the elements
    m_dElement0 = p_cRightHandSide.m_dElement0;
    m_dElement1 = p_cRightHandSide.m_dElement1;
    m_dElement2 = p_cRightHandSide.m_dElement2;
}

//ds adding and assign
CVector& CVector::operator+=( const CVector& p_cRightHandSide )
{
    //ds add all the elements
    m_dElement0 += p_cRightHandSide.m_dElement0;
    m_dElement1 += p_cRightHandSide.m_dElement1;
    m_dElement2 += p_cRightHandSide.m_dElement2;

    return *this;
}

//ds dividing and assign
CVector& CVector::operator/=( const double& p_dRightHandSide )
{
    //ds add all the elements
    m_dElement0 /= p_dRightHandSide;
    m_dElement1 /= p_dRightHandSide;
    m_dElement2 /= p_dRightHandSide;

    return *this;
}

//ds multiplication and assign
CVector& CVector::operator*=( const double& p_dRightHandSide )
{
    //ds add all the elements
    m_dElement0 *= p_dRightHandSide;
    m_dElement1 *= p_dRightHandSide;
    m_dElement2 *= p_dRightHandSide;

    return *this;
}

//ds dividing
CVector& CVector::operator/( const double& p_dRightHandSide )
{
    //ds add all the elements
    m_dElement0 /= p_dRightHandSide;
    m_dElement1 /= p_dRightHandSide;
    m_dElement2 /= p_dRightHandSide;

    return *this;
}

//ds simple addition
const CVector operator+( const CVector& p_cLeftHandSide, const CVector& p_cRightHandSide )
{
    return CVector( p_cLeftHandSide.m_dElement0 + p_cRightHandSide.m_dElement0,
                    p_cLeftHandSide.m_dElement1 + p_cRightHandSide.m_dElement1,
                    p_cLeftHandSide.m_dElement2 + p_cRightHandSide.m_dElement2 );
}

//ds simple subtraction (needed for distance calculations)
const CVector operator-( const CVector& p_cLeftHandSide, const CVector& p_cRightHandSide )
{
    return CVector( p_cLeftHandSide.m_dElement0 - p_cRightHandSide.m_dElement0,
                    p_cLeftHandSide.m_dElement1 - p_cRightHandSide.m_dElement1,
                    p_cLeftHandSide.m_dElement2 - p_cRightHandSide.m_dElement2 );
}

//ds simple multiplication
const CVector operator*( const double& p_dLeftHandSide, const CVector& p_cRightHandSide )
{
    //ds add all the elements
    return CVector( p_dLeftHandSide*p_cRightHandSide.m_dElement0,
                    p_dLeftHandSide*p_cRightHandSide.m_dElement1,
                    p_dLeftHandSide*p_cRightHandSide.m_dElement2 );
}

//ds printing
std::ostream& operator<<( std::ostream& p_cLeftHandSide, const CVector& p_cRightHandSide )
{
    //ds build the string
    p_cLeftHandSide << "0: " << p_cRightHandSide.m_dElement0 << "\n"
                    << "1: " << p_cRightHandSide.m_dElement1 << "\n"
                    << "2: " << p_cRightHandSide.m_dElement2;

    return p_cLeftHandSide;
}

//ds static functions
double CVector::absoluteValue( const CVector& p_cVector )
{
    return sqrt( pow( p_cVector.m_dElement0, 2 ) + pow( p_cVector.m_dElement1, 2 ) + pow( p_cVector.m_dElement2, 2 ) );
}

double CVector::absoluteValue( const double p_vecVector[3] )
{
    return sqrt( pow( p_vecVector[0], 2 ) + pow( p_vecVector[1], 2 ) + pow( p_vecVector[2], 2 ) );
}

const CVector CVector::crossProduct( const CVector& p_cVector1, const CVector& p_cVector2 )
{
    return CVector( p_cVector1.m_dElement1*p_cVector2.m_dElement2-p_cVector1.m_dElement2*p_cVector2.m_dElement1,
                    p_cVector1.m_dElement2*p_cVector2.m_dElement0-p_cVector1.m_dElement0*p_cVector2.m_dElement2,
                    p_cVector1.m_dElement0*p_cVector2.m_dElement1-p_cVector1.m_dElement1*p_cVector2.m_dElement0 );
}

const CVector CVector::crossProduct( const double p_vecVector1[3], const double p_vecVector2[3] )
{
    return CVector( p_vecVector1[1]*p_vecVector2[2]-p_vecVector1[2]*p_vecVector2[1],
                    p_vecVector1[2]*p_vecVector2[0]-p_vecVector1[0]*p_vecVector2[2],
                    p_vecVector1[0]*p_vecVector2[1]-p_vecVector1[1]*p_vecVector2[0] );
}

} //ds namespace NBody
