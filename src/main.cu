#include "CCubicDomain.cuh"           //ds domain structure
#include "Timer.cuh"                  //ds time measurement
#include <iostream>                   //ds cout
#include <cuda.h>                     //ds needed for eclipse indexer only (not for compilation)
#include <cuda_runtime.h>             //ds needed for eclipse indexer only (not for compilation)
#include <device_launch_parameters.h> //ds needed for eclipse indexer only (not for compilation)



//ds cuda kernels - optimized for one block
__global__ void computeAccelerationsLennardJones( NBody::CParticle* p_vecParticles,
                                                  const unsigned int* p_arrCellIndexRange,
                                                  const double p_dMinimumDistance,
                                                  const double p_dPotentialDepth,
                                                  const unsigned int p_uMaximumCellIndex,
                                                  const unsigned int p_uMaximumNeighborCellIndexRange )
{
    //ds get the current cell index
    const unsigned int uCurrentCellIndex1D( threadIdx.x );

    //ds handle the current cell particles - this only works because the data is ordered also in cell order (indexing is shifted by 1, [0,0] means there is no particle in the cell)
    for( unsigned int v = p_arrCellIndexRange[2*uCurrentCellIndex1D+0]; v < p_arrCellIndexRange[2*uCurrentCellIndex1D+1]; ++v )
    {
        //ds get the particle index
        const unsigned int uCurrentParticleIndex1D( ( v-1 ) );

        //ds for each particle in this cell
        for( int i = static_cast< int >( uCurrentCellIndex1D-p_uMaximumNeighborCellIndexRange ); i < static_cast< int >( uCurrentCellIndex1D+p_uMaximumNeighborCellIndexRange+1 ); ++i )
        {
            //ds force to sum up
            double vecTotalForce[3];

            //ds initialize the vector
            vecTotalForce[0] = 0.0;
            vecTotalForce[1] = 0.0;
            vecTotalForce[2] = 0.0;

            //ds internal "real" index
            unsigned int uNeighborCellIndex( 0 );

            //ds check for periodic boundary conditions
            if( i < 0 )
            {
                //ds shift up by maximum cell index
                uNeighborCellIndex += p_uMaximumCellIndex;
            }
            else if( i > static_cast< int >( p_uMaximumCellIndex ) )
            {
                //ds shift down by maximum cell index
                uNeighborCellIndex -= p_uMaximumCellIndex;
            }
            else
            {
                //ds regular case (i must be positive)
                uNeighborCellIndex = i;
            }

            //ds interact with all neighbor particles
            for( unsigned int w = p_arrCellIndexRange[uNeighborCellIndex+0]; w < p_arrCellIndexRange[uNeighborCellIndex+1]; ++w )
            {
                //ds get the radial vector between the particles
                double vecRadius[3];

                //ds calculate the distance: domain + particle2 - particle1
                vecRadius[0] = p_vecParticles[w-1].m_cPosition[0] - p_vecParticles[uCurrentParticleIndex1D].m_cPosition[0];
                vecRadius[1] = p_vecParticles[w-1].m_cPosition[1] - p_vecParticles[uCurrentParticleIndex1D].m_cPosition[1];
                vecRadius[2] = p_vecParticles[w-1].m_cPosition[2] - p_vecParticles[uCurrentParticleIndex1D].m_cPosition[2];

                //ds get the absolute distance
                const double dDistanceAbsolute( sqrt( pow( vecRadius[0], 2 ) + pow( vecRadius[1], 2 ) + pow( vecRadius[2], 2 ) ) );

                //ds check if above minimum distance only (using cell lists)
                if( p_dMinimumDistance < dDistanceAbsolute )
                {
                    //ds calculate the lennard jones force prefix
                    const double dLJFPrefix( -24*p_dPotentialDepth*( 2*pow( p_dMinimumDistance/dDistanceAbsolute, 12 ) - pow( p_dMinimumDistance/dDistanceAbsolute, 6  ) )
                                                                  *1/pow( dDistanceAbsolute, 2 ) );

                    //ds add the information to the force including the radial component
                    vecTotalForce[0] += dLJFPrefix*vecRadius[0];
                    vecTotalForce[1] += dLJFPrefix*vecRadius[1];
                    vecTotalForce[2] += dLJFPrefix*vecRadius[2];
                }
            }

            //ds get particle mass
            const double dParticleMass( p_vecParticles[uCurrentParticleIndex1D].m_dMass );

            //ds if we got the total force calculate the resulting acceleration and save it to our array
            p_vecParticles[uCurrentParticleIndex1D].m_cNewAcceleration[0] = vecTotalForce[0]/dParticleMass;
            p_vecParticles[uCurrentParticleIndex1D].m_cNewAcceleration[1] = vecTotalForce[1]/dParticleMass;
            p_vecParticles[uCurrentParticleIndex1D].m_cNewAcceleration[2] = vecTotalForce[2]/dParticleMass;
        }
    }
}

__global__ void updateParticlesVelocityVerlet( NBody::CParticle* p_vecParticles,
                                               const double p_dLowerBoundary,
                                               const double p_dUpperBoundary,
                                               const double p_dTimeStepSize )
{
    //ds particle index
    const unsigned int uIndex1D( threadIdx.x );

    //ds calculate domain size
    const double dDomainSize( fabs( p_dLowerBoundary ) + fabs( p_dUpperBoundary ) );

    //ds velocity-verlet for position
    p_vecParticles[uIndex1D].m_cPosition[0] = p_vecParticles[uIndex1D].m_cPosition[0] + p_dTimeStepSize*p_vecParticles[uIndex1D].m_cVelocity[0] + 1.0/2*pow( p_dTimeStepSize, 2 )*p_vecParticles[uIndex1D].m_cAcceleration[0];
    p_vecParticles[uIndex1D].m_cPosition[1] = p_vecParticles[uIndex1D].m_cPosition[1] + p_dTimeStepSize*p_vecParticles[uIndex1D].m_cVelocity[1] + 1.0/2*pow( p_dTimeStepSize, 2 )*p_vecParticles[uIndex1D].m_cAcceleration[1];
    p_vecParticles[uIndex1D].m_cPosition[2] = p_vecParticles[uIndex1D].m_cPosition[2] + p_dTimeStepSize*p_vecParticles[uIndex1D].m_cVelocity[2] + 1.0/2*pow( p_dTimeStepSize, 2 )*p_vecParticles[uIndex1D].m_cAcceleration[2];

    //ds produce periodic boundary shifting - check each element: x,y,z
    for( unsigned int v = 0; v < 3; ++v )
    {
        //ds check if we are below the boundary
        while( p_dLowerBoundary > p_vecParticles[uIndex1D].m_cPosition[v] )
        {
            //ds map the particle to the other boundary by shifting it up to the boundary
            p_vecParticles[uIndex1D].m_cPosition[v] += dDomainSize;
        }

        //ds check if we are above the boundary
        while( p_dUpperBoundary < p_vecParticles[uIndex1D].m_cPosition[v] )
        {
            //ds map the particle to the other boundary by shifting it back to the boundary
            p_vecParticles[uIndex1D].m_cPosition[v] -= dDomainSize;
        }
    }

    //ds velocity-verlet for velocity
    p_vecParticles[uIndex1D].m_cVelocity[0] = p_vecParticles[uIndex1D].m_cVelocity[0] + ( p_dTimeStepSize/2 )*( p_vecParticles[uIndex1D].m_cNewAcceleration[0] + p_vecParticles[uIndex1D].m_cAcceleration[0] );
    p_vecParticles[uIndex1D].m_cVelocity[1] = p_vecParticles[uIndex1D].m_cVelocity[1] + ( p_dTimeStepSize/2 )*( p_vecParticles[uIndex1D].m_cNewAcceleration[1] + p_vecParticles[uIndex1D].m_cAcceleration[1] );
    p_vecParticles[uIndex1D].m_cVelocity[2] = p_vecParticles[uIndex1D].m_cVelocity[2] + ( p_dTimeStepSize/2 )*( p_vecParticles[uIndex1D].m_cNewAcceleration[2] + p_vecParticles[uIndex1D].m_cAcceleration[2] );

    //ds update the old accelerations
    p_vecParticles[uIndex1D].m_cAcceleration[0] = p_vecParticles[uIndex1D].m_cNewAcceleration[0];
    p_vecParticles[uIndex1D].m_cAcceleration[1] = p_vecParticles[uIndex1D].m_cNewAcceleration[1];
    p_vecParticles[uIndex1D].m_cAcceleration[2] = p_vecParticles[uIndex1D].m_cNewAcceleration[2];
}

__global__ void getTotalEnergy( NBody::CParticle* p_vecParticles,
                                const double p_fMinimumDistance,
                                const double p_fPotentialDepth,
                                const unsigned int p_uNumberOfParticles,
                                double* p_dTotalEnergy );

int main( int argc, char** argv )
{
    //ds check simple input arguments - CAUTION: the implementation expects real numbers, the simulation will be corrupted if invalid values are entered
    if( 4 != argc )
    {
        //ds inform
        std::cout << "usage: nbody_cpu_cl [Number of particles] [Number of time steps] [Target energy]" << std::endl;
        return 0;
    }

    //ds start timing
    Timer tmTimer; tmTimer.start( );

    //ds domain configuration
    const std::pair< double, double > pairBoundaries( -1.0, 1.0 );
    const double dDomainWidth( fabs( pairBoundaries.first ) + fabs( pairBoundaries.second ) );
    const unsigned int uNumberOfParticles( atoi( argv[1] ) );

    //ds current simulation configuration
    const double dTimeStepSize( 0.0001 );
    const unsigned int uNumberOfTimeSteps( atoi( argv[2] ) );
    const double dMinimumDistance( pow( 1.0/uNumberOfParticles, 1.0/3 ) );
    const double dPotentialDepth( 1.0 );

    //ds target kinetic energy
    const double dTargetKineticEnergy( atol( argv[3] ) );

    //ds cell list information
    const unsigned int uNumberOfCells1D( floor( dDomainWidth/( 2.5*dMinimumDistance ) ) );
    const unsigned int uMaximumCellIndex( uNumberOfCells1D + pow( uNumberOfCells1D, 2 ) + pow( uNumberOfCells1D, 3 ) + 1 );

    std::cout << "------- CPU SETUP -----------------------------------------------------------" << std::endl;
    std::cout << "  Number of particles: " << uNumberOfParticles << std::endl;
    std::cout << "        Boundary (3D): [" << pairBoundaries.first << ", " << pairBoundaries.second << "]" << std::endl;
    std::cout << "         Domain Width: " << dDomainWidth << std::endl;
    std::cout << "     Minimum distance: " << dMinimumDistance << std::endl;
    std::cout << "      Cutoff distance: " << 2.5*dMinimumDistance << std::endl;
    std::cout << "      Potential depth: " << dPotentialDepth << std::endl;
    std::cout << "Target kinetic energy: " << dTargetKineticEnergy << std::endl;
    std::cout << " Number of time steps: " << uNumberOfTimeSteps << std::endl;
    std::cout << "       Time step size: " << dTimeStepSize << std::endl;
    std::cout << "------- CELL LISTS ----------------------------------------------------------" << std::endl;
    std::cout << " Number of cells 1D M: " << uNumberOfCells1D << std::endl;
    std::cout << "   Maximum cell index: " << uMaximumCellIndex << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;

    //ds allocate a domain to work with specifying number of particles and timing
    NBody::CCubicDomain cDomain( pairBoundaries, uNumberOfParticles, dMinimumDistance, uNumberOfCells1D, uMaximumCellIndex );

    //ds get particles for host and device
    thrust::host_vector< NBody::CParticle >   h_vecParticles( cDomain.createParticlesUniformFromNormalDistribution( dTargetKineticEnergy ) );
    thrust::device_vector< NBody::CParticle > d_vecParticles( h_vecParticles );

    //ds get a raw pointers for kernel usage
    NBody::CParticle *d_vecParticlesRaw( thrust::raw_pointer_cast( &d_vecParticles[0] ) );

    //ds get cell list specific parameters
    const unsigned int uMaximumNeighborCellIndexRange( cDomain.getMaximumNeighborCellIndexRange( ) );

    //ds support structure for host and device
    unsigned int *h_arrCellIndexRange( cDomain.getCellIndexRange( ) );
    unsigned int *d_arrCellIndexRange( 0 );

    //ds total energy to be calculated with CUDA
    double h_dTotalEnergy ( 0.0 );
    double* d_dTotalEnergy( 0 );

    //ds allocate memory on device
    cudaMalloc( (void **)&d_arrCellIndexRange  , uMaximumCellIndex*2*sizeof( unsigned int ) );
    cudaMalloc( (void **)&d_dTotalEnergy       , sizeof( double ) ) ;

    //ds information
    std::cout << "               Status:  0% done - current step: 0";

    //ds start simulation
    for( unsigned int uCurrentTimeStep = 1; uCurrentTimeStep < uNumberOfTimeSteps+1; ++uCurrentTimeStep )
    {
        //ds calculate percentage done
        const double dPercentageDone( 100.0*uCurrentTimeStep/uNumberOfTimeSteps );

        //ds get a formatted string -> 100% -> 3 digits
        char chBuffer[4];

        //ds fill the buffer
        std::snprintf( chBuffer, 4, "%3.0f", dPercentageDone );

        //ds print info
        std::cout << '\xd';
        std::cout << "               Status: " << chBuffer << "% done - current step: " << uCurrentTimeStep;

        //ds copy support structure memory to gpu
        cudaMemcpy( d_arrCellIndexRange, h_arrCellIndexRange, uMaximumCellIndex*2*sizeof( unsigned int ), cudaMemcpyHostToDevice );

        //ds compute accelerations for all cells - launch as many threads as we have cells
        computeAccelerationsLennardJones<<< 1, uMaximumCellIndex >>>( d_vecParticlesRaw,
                                                                      d_arrCellIndexRange,
                                                                      dMinimumDistance,
                                                                      dPotentialDepth,
                                                                      uMaximumCellIndex,
                                                                      uMaximumNeighborCellIndexRange );
        //ds update particle properties
        updateParticlesVelocityVerlet<<< 1, uNumberOfParticles >>>( d_vecParticlesRaw,
                                                                    pairBoundaries.first,
                                                                    pairBoundaries.second,
                                                                    dTimeStepSize );

        //ds compute total energy
        getTotalEnergy<<< 1, uNumberOfParticles, uNumberOfParticles*sizeof( double )  >>>( d_vecParticlesRaw,
                                                                                           dMinimumDistance,
                                                                                           dPotentialDepth,
                                                                                           uNumberOfParticles,
                                                                                           d_dTotalEnergy );

        //ds copy device vector to host
        thrust::copy( d_vecParticles.begin( ), d_vecParticles.end( ), h_vecParticles.begin( ) );

        //ds update cell lists
        cDomain.updateCellList( h_vecParticles );

        //ds get the integrals information from gpu to cpu
        cudaMemcpy( &h_dTotalEnergy , d_dTotalEnergy, sizeof( double ), cudaMemcpyDeviceToHost );

        //ds record situation (we will write the stream to the file in one operation afterwards )
        cDomain.saveParticlesToStream( );
        cDomain.saveIntegralsToStream( h_dTotalEnergy );
    }

    //ds deallocate memory
    cudaFree( d_arrCellIndexRange );
    cudaFree( d_dTotalEnergy );

    //ds save the streams to a file
    cDomain.writeParticlesToFile( "bin/simulation.txt", uNumberOfTimeSteps );
    cDomain.writeIntegralsToFile( "bin/integrals.txt", uNumberOfTimeSteps, dTimeStepSize );

    //ds stop timing
    const double dDurationSeconds( tmTimer.stop( ) );

    //ds cause an output ostream
    std::cout << std::endl;
    std::cout << "     Computation time: " << dDurationSeconds << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;

    return 0;
}

__global__ void getTotalEnergy( NBody::CParticle* p_vecParticles,
                                const double p_fMinimumDistance,
                                const double p_fPotentialDepth,
                                const unsigned int p_uNumberOfParticles,
                                double* p_dTotalEnergy )
{
    //ds dynamic shared total energy to sum up by first thread
    extern __shared__ double s_arrTotalEnergy[];

    //ds regular index
    const unsigned int uIndex1D( threadIdx.x );

    //ds make sure the shared memory is empty (each thread does this)
    s_arrTotalEnergy[uIndex1D] = 0.0;

    //ds wait until all threads are done
    __syncthreads( );

    //ds add the kinetic component of the current particle
    s_arrTotalEnergy[uIndex1D] += p_vecParticles[uIndex1D].m_dMass/2*pow( sqrt( pow( p_vecParticles[uIndex1D].m_cVelocity[0], 2 )
                                                                              + pow( p_vecParticles[uIndex1D].m_cVelocity[1], 2 )
                                                                              + pow( p_vecParticles[uIndex1D].m_cVelocity[2], 2 ) ), 2 );

    //ds cutoff
    const float fDistanceCutoff( 2.5*p_fMinimumDistance );

    //ds calculate the total energy of the new configuration - loop over all other particles (dont do the same particles twice)
    for( unsigned int u = uIndex1D+1; u < p_uNumberOfParticles; ++u )
    {
        //ds get the absolute distance
        const float fDistanceAbsolute( sqrt( pow( p_vecParticles[u].m_cPosition[0] - p_vecParticles[uIndex1D].m_cPosition[0], 2 )
                                           + pow( p_vecParticles[u].m_cPosition[1] - p_vecParticles[uIndex1D].m_cPosition[1], 2 )
                                           + pow( p_vecParticles[u].m_cPosition[2] - p_vecParticles[uIndex1D].m_cPosition[2], 2 ) ) );

        //ds if we are between the minimum distance and the cutoff range
        if( p_fMinimumDistance < fDistanceAbsolute && fDistanceCutoff > fDistanceAbsolute )
        {
            //ds add the potential component
            s_arrTotalEnergy[uIndex1D] += 4*p_fPotentialDepth*( pow( p_fMinimumDistance/fDistanceAbsolute, 12 ) - pow( p_fMinimumDistance/fDistanceAbsolute, 6 ) );
        }
    }

    //ds wait until all threads are done
    __syncthreads( );

    //ds thread 0 calculates the total energy
    if( 0 == uIndex1D )
    {
        //ds total energy to sum up
        float dTotalEnergy( 0.0 );

        for( unsigned int u = 0; u < p_uNumberOfParticles; ++u )
        {
            dTotalEnergy += s_arrTotalEnergy[u];
        }

        //ds set the return value
        *p_dTotalEnergy = dTotalEnergy;
    }
}
