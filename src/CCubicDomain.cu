#include "CCubicDomain.cuh"
#include <algorithm>      //ds std::sort



namespace NBody
{

//ds default constructor requires environmental parameters: N number of bodies, dT time step, T number of time steps inclusive calculations (ugly but required for the const attributes)
CCubicDomain::CCubicDomain( const std::pair< double, double >& p_pairBoundaries,
                            const unsigned int& p_uNumberOfParticles,
                            const double& p_dMinimumDistance,
                            const unsigned int& p_uNumberOfCells,
                            const unsigned int& p_uMaximumCellIndex ): m_pairBoundaries( p_pairBoundaries ),
                                                                       m_dDomainWidth( fabs( m_pairBoundaries.first ) + fabs( m_pairBoundaries.second ) ),
                                                                       m_uNumberOfParticles( p_uNumberOfParticles ),
                                                                       m_uNumberOfCells( p_uNumberOfCells ),
                                                                       m_uMaximumCellIndex( p_uMaximumCellIndex ),
                                                                       m_uMaximumNeighborCellIndexRange( getCellIndex( m_pairBoundaries.second, m_pairBoundaries.second, m_pairBoundaries.second )
                                                                                                       - getCellIndex( ( m_pairBoundaries.second - 2.5*p_dMinimumDistance/3 ),
                                                                                                                       ( m_pairBoundaries.second - 2.5*p_dMinimumDistance/3 ),
                                                                                                                       ( m_pairBoundaries.second - 2.5*p_dMinimumDistance/3 ) ) ),
                                                                       m_strParticleInformation( "" ),
                                                                       m_strIntegralsInformation( "" )
{
    //ds initialize particles vector
    m_vecParticles.reserve( p_uNumberOfParticles );

    //ds initialize support structure - pairs are initialized with 0,0 by default
    m_arrCellIndexRange = new unsigned int[2*m_uMaximumCellIndex];

    //ds initialize
    for( unsigned int u = 0; u < m_uMaximumCellIndex; ++u )
    {
        m_arrCellIndexRange[2*u+0] = 0;
        m_arrCellIndexRange[2*u+1] = 0;
    }
}

//ds default destructor
CCubicDomain::~CCubicDomain( )
{
    //ds clear internal structure
    m_vecParticles.clear( );

    //ds deallocated memory
    delete[] m_arrCellIndexRange;
}

//ds accessors
thrust::host_vector< NBody::CParticle > CCubicDomain::createParticlesUniformFromNormalDistribution( const double& p_dTargetKineticEnergy, const double& p_dParticleMass )
{
    //ds host vector to return to main scope
    thrust::host_vector< NBody::CParticle > vecParticles;

    //ds kinetic energy to derive from initial situation
    double dKineticEnergy( 0.0 );

    //ds set particle information for each
    for( unsigned int u = 0; u < m_uNumberOfParticles; ++u )
    {
        //ds get a particle instance
        CParticle cParticle( u );

        //ds set the particle mass (same for all particles in this case)
        cParticle.m_dMass = p_dParticleMass;

        //ds set the position: uniformly distributed between boundaries in this case
        cParticle.m_cPosition[0] = _getUniformlyDistributedNumber( );
        cParticle.m_cPosition[1] = _getUniformlyDistributedNumber( );
        cParticle.m_cPosition[2] = _getUniformlyDistributedNumber( );

        //ds set velocity values: from normal distribution
        cParticle.m_cVelocity[0] = _getNormallyDistributedNumber( );
        cParticle.m_cVelocity[1] = _getNormallyDistributedNumber( );
        cParticle.m_cVelocity[2] = _getNormallyDistributedNumber( );

        //ds set accelerations
        cParticle.m_cAcceleration[0]    = 0.0;
        cParticle.m_cAcceleration[1]    = 0.0;
        cParticle.m_cAcceleration[2]    = 0.0;
        cParticle.m_cNewAcceleration[0] = 0.0;
        cParticle.m_cNewAcceleration[1] = 0.0;
        cParticle.m_cNewAcceleration[2] = 0.0;

        //ds add the resulting kinetic component (needed below)
        dKineticEnergy += cParticle.m_dMass/2*pow( NBody::CVector::absoluteValue( cParticle.m_cVelocity ), 2 );

        //ds add it to the vectors
        vecParticles.push_back( cParticle );
        m_vecParticles.push_back( cParticle );
    }

    //ds calculate the rescaling factor
    const double dRescalingFactor( sqrt( p_dTargetKineticEnergy/dKineticEnergy ) );

    //ds for each particle
    for( unsigned int u = 0; u < m_uNumberOfParticles; ++u )
    {
        //ds rescale the velocity component
        vecParticles[u].m_cVelocity[0] *= dRescalingFactor;
        vecParticles[u].m_cVelocity[1] *= dRescalingFactor;
        vecParticles[u].m_cVelocity[2] *= dRescalingFactor;
    }

    return vecParticles;
}

/*void CCubicDomain::updateParticlesVelocityVerlet( const double& p_dTimeStep, const double& p_dMinimumDistance, const double& p_dPotentialDepth )
{
    //ds allocate a temporary array to hold the accelerations
    CVector *arrNewAccelerations = new CVector[m_uNumberOfParticles];

    //ds for each cell id (dangerous games with possible unsigned integer overflows here to avoid tons of if cases and worse readability)
    for( unsigned int u = 0; u < m_uMaximumCellIndex; ++u )
    {
        //ds handle the current cell particles - this only works because the data is ordered also in cell order (indexing is shifted by 1, [0,0] means there is no particle in the cell)
        for( unsigned int v = m_arrCellIndexRange[u].first; v < m_arrCellIndexRange[u].second; ++v )
        {
            //ds for each particle in this cell
            for( int i = static_cast< int >( u-m_uMaximumNeighborCellIndexRange ); i < static_cast< int >( u+m_uMaximumNeighborCellIndexRange+1 ); ++i )
            {
                //ds force to sum up
                CVector cTotalForce;

                //ds internal "real" index
                unsigned int uNeighborCellIndex( 0 );

                //ds check for periodic boundary conditions
                if( i < 0 )
                {
                    //ds shift up by maximum cell index
                    uNeighborCellIndex += m_uMaximumCellIndex;
                }
                else if( i > static_cast< int >( m_uMaximumCellIndex ) )
                {
                    //ds shift down by maximum cell index
                    uNeighborCellIndex -= m_uMaximumCellIndex;
                }
                else
                {
                    //ds regular case (i must be positive)
                    uNeighborCellIndex = i;
                }

                //ds interact with all neighbor particles
                for( unsigned int w = m_arrCellIndexRange[uNeighborCellIndex].first; w < m_arrCellIndexRange[uNeighborCellIndex].second; ++w )
                {
                    cTotalForce += _getLennardJonesForce( m_vecParticles[v-1], m_vecParticles[w-1], p_dMinimumDistance, p_dPotentialDepth );
                }

                //ds if we got the total force calculate the resulting acceleration and save it to our array
                arrNewAccelerations[v] = cTotalForce/m_vecParticles[v].m_dMass;
            }
        }
    }

    //ds for each particle we have to calculate the effect of the acceleration now
    for( unsigned int u = 0; u < m_uNumberOfParticles; ++u )
    {
        //ds velocity-verlet for position
        m_vecParticles[u].m_cPosition = m_vecParticles[u].m_cPosition + p_dTimeStep*m_vecParticles[u].m_cVelocity + 1/2*pow( p_dTimeStep, 2 )*m_vecParticles[u].m_cAcceleration;

        //ds produce periodic boundary shifting - check each element: x,y,z
        for( unsigned int v = 0; v < 3; ++v )
        {
            //ds check if we are below the boundary
            while( m_pairBoundaries.first > m_vecParticles[u].m_cPosition( v ) )
            {
                //ds map the particle to the other boundary by shifting it up to the boundary
                m_vecParticles[u].m_cPosition( v ) += m_dDomainWidth;
            }

            //ds check if we are above the boundary
            while( m_pairBoundaries.second < m_vecParticles[u].m_cPosition ( v ) )
            {
                //ds map the particle to the other boundary by shifting it back to the boundary
                m_vecParticles[u].m_cPosition( v ) -= m_dDomainWidth;
            }
        }

        //ds velocity-verlet for velocity
        m_vecParticles[u].m_cVelocity = m_vecParticles[u].m_cVelocity + p_dTimeStep/2*( arrNewAccelerations[u] + m_vecParticles[u].m_cAcceleration );

        //ds update the acceleration
        m_vecParticles[u].m_cAcceleration = arrNewAccelerations[u];
    }

    //ds deallocate the temporary accelerations array
    delete[] arrNewAccelerations;

    //ds update our cell list structures
    _updateCellList( );
}*/

void CCubicDomain::saveParticlesToStream( )
{
    //ds format: X Y Z U V W

    //ds for each particle
    for( unsigned int u = 0; u < m_uNumberOfParticles; ++u )
    {
        //ds get a buffer for snprintf (max possible buffer size checked)
        char chBuffer[256];

        //ds get the particle stream
        std::snprintf( chBuffer, 100, "%f %f %f %f %f %f", m_vecParticles[u].m_cPosition[0], m_vecParticles[u].m_cPosition[1], m_vecParticles[u].m_cPosition[2],
                                                           m_vecParticles[u].m_cVelocity[0], m_vecParticles[u].m_cVelocity[1], m_vecParticles[u].m_cVelocity[2] );

        //ds append the buffer to our string
        m_strParticleInformation += chBuffer;
        m_strParticleInformation += "\n";
    }
}

void CCubicDomain::saveIntegralsToStream( const double& p_dTotalEnergy )
{
    //ds format: E X Y Z X Y Z X Y Z

    //ds get information - caution, memory gets allocated
    const double dTotalEnergy               ( p_dTotalEnergy );
    const NBody::CVector vecCenterOfMass    ( _getCenterOfMass( ) );
    const NBody::CVector vecAngularMomentum ( _getAngularMomentum( ) );
    const NBody::CVector vecLinearMomentum  ( _getLinearMomentum( ) );

    //ds buffer for snprintf
    char chBuffer[256];

    //ds get the integrals stream
    std::snprintf( chBuffer, 256, "%f %f %f %f %f %f %f %f %f %f", dTotalEnergy,
                                                                   vecCenterOfMass[0]   , vecCenterOfMass[1]   , vecCenterOfMass[2],
                                                                   vecAngularMomentum[0], vecAngularMomentum[1], vecAngularMomentum[2],
                                                                   vecLinearMomentum[0] , vecLinearMomentum[1] , vecLinearMomentum[2] );

    //ds append the buffer to our string
    m_strIntegralsInformation += chBuffer;
    m_strIntegralsInformation += "\n";
}

void CCubicDomain::writeParticlesToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps )
{
    //ds ofstream object
    std::ofstream ofsFile;

    //ds open the file for writing
    ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

    //ds if it worked
    if( ofsFile.is_open( ) )
    {
        //ds first dump setup information number of particles and timesteps
        ofsFile << m_uNumberOfParticles << " " << p_uNumberOfTimeSteps << "\n" << m_strParticleInformation;
    }

    //ds close the file
    ofsFile.close( );
}

void CCubicDomain::writeIntegralsToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps, const double& p_dTimeStepSize )
{
    //ds ofstream object
    std::ofstream ofsFile;

    //ds open the file for writing
    ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

    //ds if it worked
    if( ofsFile.is_open( ) )
    {
        //ds dump first integrals information
        ofsFile << p_uNumberOfTimeSteps << " " << p_dTimeStepSize << "\n" << m_strIntegralsInformation;
    }

    //ds close the file
    ofsFile.close( );
}

unsigned int* CCubicDomain::getCellIndexRange( )
{
    return m_arrCellIndexRange;
}

unsigned int CCubicDomain::getMaximumNeighborCellIndexRange( ) const
{
    return m_uMaximumNeighborCellIndexRange;
}

CVector CCubicDomain::_getCenterOfMass( ) const
{
    //ds center to find
    CVector cCenter;

    //ds total mass
    double dMassTotal( 0.0 );

    //ds for each particle
    for( unsigned int u = 0; u < m_uNumberOfParticles; ++u )
    {
        //ds add the current relative mass
        cCenter( 0 ) += m_vecParticles[u].m_dMass*m_vecParticles[u].m_cPosition[0];
        cCenter( 1 ) += m_vecParticles[u].m_dMass*m_vecParticles[u].m_cPosition[1];
        cCenter( 2 ) += m_vecParticles[u].m_dMass*m_vecParticles[u].m_cPosition[2];

        //ds add the current mass
        dMassTotal += m_vecParticles[u].m_dMass;
    }

    //ds divide by total mass
    cCenter /= dMassTotal;

    return cCenter;
}

CVector CCubicDomain::_getAngularMomentum( ) const
{
    //ds momentum
    CVector cMomentum;

    //ds for each particle
    for( unsigned int u = 0; u < m_uNumberOfParticles; ++u )
    {
        //ds temporary array for conversion
        double vecMassTimesVelocity[3];

        //ds set the values
        vecMassTimesVelocity[0] = m_vecParticles[u].m_dMass*m_vecParticles[u].m_cVelocity[0];
        vecMassTimesVelocity[1] = m_vecParticles[u].m_dMass*m_vecParticles[u].m_cVelocity[1];
        vecMassTimesVelocity[2] = m_vecParticles[u].m_dMass*m_vecParticles[u].m_cVelocity[2];

        //ds add the current momentum
        cMomentum += NBody::CVector::crossProduct( m_vecParticles[u].m_cPosition, vecMassTimesVelocity );
    }

    return cMomentum;
}

CVector CCubicDomain::_getLinearMomentum( ) const
{
    //ds momentum
    CVector cMomentum;

    //ds for each particle
    for( unsigned int u = 0; u < m_uNumberOfParticles; ++u )
    {
        //ds add the current momentum
        cMomentum( 0 ) += m_vecParticles[u].m_dMass*m_vecParticles[u].m_cVelocity[0];
        cMomentum( 1 ) += m_vecParticles[u].m_dMass*m_vecParticles[u].m_cVelocity[1];
        cMomentum( 2 ) += m_vecParticles[u].m_dMass*m_vecParticles[u].m_cVelocity[2];
    }

    return cMomentum;
}

//ds helpers
double CCubicDomain::_getLennardJonesPotential( const CParticle& p_CParticle1,  const CParticle& p_CParticle2, const double& p_dMinimumDistance, const double& p_dPotentialDepth ) const
{
    //ds cutoff distance
    const double dDistanceCutoff( 2.5*p_dMinimumDistance );

    //ds potential to calculate - default 0
    double dPotential( 0.0 );

    //ds get the distance between the particles
    const CVector cRadius( p_CParticle2.m_cPosition[0] - p_CParticle1.m_cPosition[0],
                           p_CParticle2.m_cPosition[1] - p_CParticle1.m_cPosition[1],
                           p_CParticle2.m_cPosition[2] - p_CParticle1.m_cPosition[2] );

    //ds get the current distance between 2 and 1
    const double dDistanceAbsolute( NBody::CVector::absoluteValue( cRadius ) );

    //ds if we are between the minimum distance and the cutoff range
    if( p_dMinimumDistance < dDistanceAbsolute && dDistanceCutoff > dDistanceAbsolute )
    {
        //ds add the potential
        dPotential = 4*p_dPotentialDepth*( pow( p_dMinimumDistance/dDistanceAbsolute, 12 )
                                         - pow( p_dMinimumDistance/dDistanceAbsolute, 6 ) );
    }

    return dPotential;
}

CVector CCubicDomain::_getLennardJonesForce( const CParticle& p_CParticle1,  const CParticle& p_CParticle2, const double& p_dMinimumDistance, const double& p_dPotentialDepth ) const
{
    //ds get the distance between the particles
    const CVector cRadius( p_CParticle2.m_cPosition[0] - p_CParticle1.m_cPosition[0],
                           p_CParticle2.m_cPosition[1] - p_CParticle1.m_cPosition[1],
                           p_CParticle2.m_cPosition[2] - p_CParticle1.m_cPosition[2] );

    //ds get the current distance between 2 and 1
    const double dDistanceAbsolute( NBody::CVector::absoluteValue( cRadius ) );

    //ds return the resulting force - no cutoff check here and no periodic boundaries calculation, our cell list takes care of that
    return -24*p_dPotentialDepth*( 2*pow( p_dMinimumDistance/dDistanceAbsolute, 12 ) - pow( p_dMinimumDistance/dDistanceAbsolute, 6  ) )
                                *1/pow( dDistanceAbsolute, 2 )*cRadius;
}

double CCubicDomain::_getUniformlyDistributedNumber( ) const
{
    //ds drand48 returns [0,1], we need [-1,1] -> therefore 2x[0,1] -> [0,2] -> -1 ->[0-1,2-1] = [-1,1] (for domain size 2)
    return ( m_dDomainWidth*drand48( ) - m_dDomainWidth/2 );
}

double CCubicDomain::_getNormallyDistributedNumber( ) const
{
    //ds calculate the uniform number first [0,1]
    const double dUniformNumber( drand48( ) );

    //ds return the normal one (formula)
    return sqrt( -2*log( dUniformNumber ) )*cos( 2*M_PI*dUniformNumber );
}

void CCubicDomain::updateCellList( thrust::host_vector< NBody::CParticle>& p_vecParticles )
{
    //ds construct new cell data
    for( unsigned int u = 0; u < m_uNumberOfParticles; ++u )
    {
        //ds set the cell index of the current particle
        p_vecParticles[u].m_uIndexCell = getCellIndex( p_vecParticles[u].m_cPosition[0], p_vecParticles[u].m_cPosition[1], p_vecParticles[u].m_cPosition[2] );
    }

    //ds sort the particles according to cell index
    thrust::sort( p_vecParticles.begin( ), p_vecParticles.end( ), _isCurrentCellIndexSmaller );

    //ds copy vectors to our particles handle
    thrust::copy( p_vecParticles.begin( ), p_vecParticles.end( ), m_vecParticles.begin( ) );

    //ds update the index range
    updateCellIndexRange( );
}

void CCubicDomain::updateCellIndexRange( )
{
    //ds check if index range is allocated
    if( 0 != m_arrCellIndexRange )
    {
        //ds loop parameters: current cell index
        unsigned int uCurrentCellIndex( m_vecParticles[0].m_uIndexCell );

        //ds already save the first element (has to be the first since vecParticles is sorted by cell indexi)
        m_arrCellIndexRange[2*uCurrentCellIndex+0] = 1;

        //ds loop over all particles/cells - shifted, index 0 means there is no particle
        for( unsigned int u = 1; u < m_uNumberOfParticles+1; ++u )
        {
            //ds if we get the first mismatch or we're at the end
            if( uCurrentCellIndex < m_vecParticles[u-1].m_uIndexCell || m_uNumberOfParticles == u )
            {
                //ds add the index before if its not the last
                if( m_uNumberOfParticles-1 > u )
                {
                    //ds take the one before
                    m_arrCellIndexRange[2*uCurrentCellIndex+1] = u-1;
                }
                else
                {
                    //ds add the last
                    m_arrCellIndexRange[2*uCurrentCellIndex+1] = u;
                }

                //ds update cell indexi
                uCurrentCellIndex  = m_vecParticles[u-1].m_uIndexCell;

                //ds check limits
                if( m_uMaximumCellIndex < uCurrentCellIndex )
                {
                    //ds TODO implement exceptions
                    throw std::exception( );
                }

                //ds save the current index and go on
                m_arrCellIndexRange[2*uCurrentCellIndex+0] = u;
            }
        }
    }
}

unsigned int CCubicDomain::getCellIndex( const double& p_dX, const double& p_dY, const double& p_dZ ) const
{
    //ds calculate the cell index and return it
    return ( floor( ( p_dX + m_dDomainWidth/2 )/2*m_uNumberOfCells )
           + floor( ( p_dY + m_dDomainWidth/2 )/2*m_uNumberOfCells )*m_uNumberOfCells
           + floor( ( p_dZ + m_dDomainWidth/2 )/2*m_uNumberOfCells )*m_uNumberOfCells*m_uNumberOfCells );
}

bool CCubicDomain::_isCurrentCellIndexSmaller( const CParticle& p_pcParticle1, const CParticle& p_pcParticle2 )
{
    return ( p_pcParticle1.m_uIndexCell < p_pcParticle2.m_uIndexCell );
}

} //ds namespace NBody
