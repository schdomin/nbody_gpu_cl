#ifndef CCUBICDOMAIN_CUH_
#define CCUBICDOMAIN_CUH_

#include "CVector.cuh"            //ds basic structure for 3d information
#include <stdlib.h>               //ds labs, drand48
#include "thrust/sequence.h"      //ds thrust
#include "thrust/sort.h"          //ds thrust
#include "thrust/copy.h"          //ds thrust
#include "thrust/device_ptr.h"    //ds thrust
#include "thrust/device_vector.h" //ds thrust



namespace NBody
{

//ds inherent structure for particle information
struct CParticle
{
    //ds fix the particle index at construction
    CParticle( const unsigned int& p_uIndexParticle ): m_uIndexParticle( p_uIndexParticle ){ };

    //ds constant properties
    double m_dMass;

    //ds position and velocity and acceleration information - thrust::device doesn't allow nested structs
    double m_cPosition[3];
    double m_cVelocity[3];
    double m_cAcceleration[3];
    double m_cNewAcceleration[3];

    //ds indexi
    unsigned int m_uIndexParticle; //<- should be const, but because we're using vectors it has to remain changeable (copy ctor of vector, etc)
    unsigned int m_uIndexCell;
};

//ds domain to solve n-body problems
class CCubicDomain
{

//ds ctor/dtor
public:

    //ds default constructor requires environmental parameters: N number of bodies, dT time step, T number of time steps
    CCubicDomain( const std::pair< double, double >& p_pairBoundaries,
                  const unsigned int& p_uNumberOfParticles,
                  const double& p_dMinimumDistance,
                  const unsigned int& p_uNumberOfCells,
                  const unsigned int& p_uMaximumCellIndex );

    //ds default destructor
    ~CCubicDomain( );

//ds attributes
private:

    //ds internal particle storage (reference to host vector in main)
    std::vector< NBody::CParticle > m_vecParticles;

    //ds support structure (does not grow dynamically, is of size N) and returns the start and end index of the current cell
    unsigned int *m_arrCellIndexRange;

    //ds domain properties
    const std::pair< double, double > m_pairBoundaries;
    const double m_dDomainWidth;
    const unsigned int m_uNumberOfParticles;

    //ds cell list properties
    const unsigned int m_uNumberOfCells;
    const unsigned int m_uMaximumCellIndex;
    const unsigned int m_uMaximumNeighborCellIndexRange;

    //ds streams for offline data - needed for the movie and graphs
    std::string m_strParticleInformation;
    std::string m_strIntegralsInformation;

//ds accessors
public:

    thrust::host_vector< NBody::CParticle > createParticlesUniformFromNormalDistribution( const double& p_dTargetKineticEnergy, const double& p_dParticleMass = 1.0 );
    void saveParticlesToStream( );
    void saveIntegralsToStream( const double& p_dTotalEnergy );
    void writeParticlesToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps );
    void writeIntegralsToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps, const double& p_dTimeStepSize );

    //ds getter/setter to communicate with CUDA kernels
    unsigned int* getCellIndexRange( );
    unsigned int getMaximumNeighborCellIndexRange( ) const;

//ds helpers
private:

    CVector _getCenterOfMass( ) const;
    CVector _getAngularMomentum( ) const;
    CVector _getLinearMomentum( ) const;
    double _getLennardJonesPotential( const CParticle& p_CParticle1,  const CParticle& p_CParticle2, const double& p_dMinimumDistance, const double& p_dPotentialDepth ) const;
    CVector _getLennardJonesForce( const CParticle& p_CParticle1,  const CParticle& p_CParticle2, const double& p_dMinimumDistance, const double& p_dPotentialDepth ) const;
    double _getUniformlyDistributedNumber( ) const;
    double _getNormallyDistributedNumber( ) const;

public:

    //ds cell lists
    void updateCellList( thrust::host_vector< NBody::CParticle >& p_vecParticles );
    void updateCellIndexRange( );
    unsigned int getCellIndex( const double& p_dX, const double& p_dY, const double& p_dZ ) const;

//ds internal statics
private:

    //ds needed for sorting
    static bool _isCurrentCellIndexSmaller( const CParticle& p_pcParticle1, const CParticle& p_pcParticle2 );

};

} //ds namespace NBody



#endif //ds CCUBICDOMAIN_CUH_
