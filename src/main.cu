#include "CCubicDomain.cuh"           //ds domain structure
#include "Timer.cuh"                  //ds time measurement
#include <iostream>                   //ds cout
#include <cuda.h>                     //ds needed for eclipse indexer only (not for compilation)
#include <cuda_runtime.h>             //ds needed for eclipse indexer only (not for compilation)
#include <device_launch_parameters.h> //ds needed for eclipse indexer only (not for compilation)



//ds cuda kernels - optimized for one block
__global__ void computeAccelerationsLennardJones( const unsigned int p_uNumberOfParticles,
                                                  const NBody::CParticle* p_vecParticles,
                                                  const std::pair< unsigned int, unsigned int >* p_arrCellIndexRange,
                                                  const double p_dMinimumDistance,
                                                  const double p_dPotentialDepth,
                                                  const unsigned int p_uMaximumCellIndex,
                                                  const unsigned int p_uMaximumNeighborCellIndexRange,
                                                  double* p_arrNewAccelerations )
{
    //ds get the current cell index
    const unsigned int uCurrentCellIndex( threadIdx.x );

    //ds handle the current cell particles - this only works because the data is ordered also in cell order (indexing is shifted by 1, [0,0] means there is no particle in the cell)
    for( unsigned int v = p_arrCellIndexRange[uCurrentCellIndex].first; v < p_arrCellIndexRange[uCurrentCellIndex].second; ++v )
    {
        //ds get the particle index
        const unsigned int uCurrentParticleIndex1D( ( v-1 ) );

        //ds for each particle in this cell
        for( int i = static_cast< int >( uCurrentCellIndex-p_uMaximumNeighborCellIndexRange ); i < static_cast< int >( uCurrentCellIndex+p_uMaximumNeighborCellIndexRange+1 ); ++i )
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
            for( unsigned int w = p_arrCellIndexRange[uNeighborCellIndex].first; w < p_arrCellIndexRange[uNeighborCellIndex].second; ++w )
            {
                //ds get the radial vector between the particles
                double vecRadius[3];

                //ds calculate the distance: domain + particle2 - particle1
                vecRadius[0] = p_vecParticles[w-1].m_cPosition[0] - p_vecParticles[uCurrentParticleIndex1D].m_cPosition[0];
                vecRadius[1] = p_vecParticles[w-1].m_cPosition[1] - p_vecParticles[uCurrentParticleIndex1D].m_cPosition[1];
                vecRadius[2] = p_vecParticles[w-1].m_cPosition[2] - p_vecParticles[uCurrentParticleIndex1D].m_cPosition[2];

                //ds get the absolute distance
                const double dDistanceAbsolute( sqrt( pow( vecRadius[0], 2 ) + pow( vecRadius[1], 2 ) + pow( vecRadius[2], 2 ) ) );

                //ds calculate the lennard jones force prefix
                const double dLJFPrefix( -24*p_dPotentialDepth*( 2*pow( p_dMinimumDistance/dDistanceAbsolute, 12 ) - pow( p_dMinimumDistance/dDistanceAbsolute, 6  ) )
                                                              *1/pow( dDistanceAbsolute, 2 ) );

                //ds add the information to the force including the radial component
                vecTotalForce[0] += dLJFPrefix*vecRadius[0];
                vecTotalForce[1] += dLJFPrefix*vecRadius[1];
                vecTotalForce[2] += dLJFPrefix*vecRadius[2];
            }

            //ds get particle mass
            const double dParticleMass( p_vecParticles[uCurrentParticleIndex1D].m_dMass );

            //ds if we got the total force calculate the resulting acceleration and save it to our array
            p_arrNewAccelerations[3*uCurrentParticleIndex1D+0] = vecTotalForce[0]/dParticleMass;
            p_arrNewAccelerations[3*uCurrentParticleIndex1D+1] = vecTotalForce[1]/dParticleMass;
            p_arrNewAccelerations[3*uCurrentParticleIndex1D+2] = vecTotalForce[2]/dParticleMass;
        }
    }
}

__global__ void updateParticlesVelocityVerlet( const unsigned int p_uNumberOfParticles,
                                               NBody::CParticle* p_vecParticles,
                                               const float p_fTimeStepSize )
{

}

int main( int argc, char** argv )
{
    //ds start timing
    Timer tmTimer; tmTimer.start( );

    //ds domain configuration
    const std::pair< double, double > pairBoundaries( -1.0, 1.0 );
    const double dDomainWidth( labs( pairBoundaries.first ) + labs( pairBoundaries.second ) );
    const unsigned int uNumberOfParticles( 100 );

    //ds current simulation configuration
    const double dTimeStepSize( 0.0001 );
    const unsigned int uNumberOfTimeSteps( 5000 );
    const double dMinimumDistance( pow( 1.0/uNumberOfParticles, 1.0/3 ) );
    const double dPotentialDepth( 1.0 );

    //ds target kinetic energy
    const double dTargetKineticEnergy( 1000.0 );

    //ds cell list information
    const unsigned int uNumberOfCells1D( floor( dDomainWidth/( 2.5*dMinimumDistance ) ) );
    const unsigned int uMaximumCellIndex( uNumberOfCells1D + pow( uNumberOfCells1D, 2 ) + pow( uNumberOfCells1D, 3 ) + 1 );

    std::cout << "--------CPU SETUP------------------------------------------------------------" << std::endl;
    std::cout << "  Number of particles: " << uNumberOfParticles << std::endl;
    std::cout << "        Boundary (3D): [" << pairBoundaries.first << ", " << pairBoundaries.second << "]" << std::endl;
    std::cout << "         Domain Width: " << dDomainWidth << std::endl;
    std::cout << "     Minimum distance: " << dMinimumDistance << std::endl;
    std::cout << "      Cutoff distance: " << 2.5*dMinimumDistance << std::endl;
    std::cout << "      Potential depth: " << dPotentialDepth << std::endl;
    std::cout << "Target kinetic energy: " << dTargetKineticEnergy << std::endl;
    std::cout << " Number of time steps: " << uNumberOfTimeSteps << std::endl;
    std::cout << "       Time step size: " << dTimeStepSize << std::endl;
    std::cout << "--------CELL LISTS-----------------------------------------------------------" << std::endl;
    std::cout << " Number of cells 1D M: " << uNumberOfCells1D << std::endl;
    std::cout << "   Maximum cell index: " << uMaximumCellIndex << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;

    //ds allocate a domain to work with specifying number of particles and timing
    NBody::CCubicDomain cDomain( pairBoundaries, uNumberOfParticles, dMinimumDistance, uNumberOfCells1D, uMaximumCellIndex );

    //ds create particles uniformly from a normal distribution
    cDomain.createParticlesUniformFromNormalDistribution( dTargetKineticEnergy );

    //ds support structure for host and device
    std::pair< unsigned int, unsigned int > *h_arrCellIndexRange = cDomain.getCellIndexRange( );
    std::pair< unsigned int, unsigned int > *d_arrCellIndexRange = 0;

    //ds accelerations buffer on the GPU
    double* d_arrNewAccelerations( 0 ); //Nx3

    //ds allocate memory on device
    cudaMalloc( (void **)&d_arrCellIndexRange, uMaximumCellIndex*sizeof( std::pair< unsigned int, unsigned int > ) );
    cudaMalloc( (void **)&d_arrNewAccelerations, uNumberOfParticles*3*sizeof( double ) ) ;

    //ds get particles for the device
    thrust::device_vector< NBody::CParticle > d_vecParticles( 100 ); //cDomain.getParticles( ) );

    //ds get a raw pointer for kernel usage
    NBody::CParticle *vecParticles( thrust::raw_pointer_cast( &d_vecParticles[0] ) );

    //ds get cell list specific parameters
    const unsigned int uMaximumNeighborCellIndexRange( cDomain.getMaximumNeighborCellIndexRange( ) );

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
        cudaMemcpy( d_arrCellIndexRange, h_arrCellIndexRange, uMaximumCellIndex*sizeof( std::pair< unsigned int, unsigned int > ), cudaMemcpyHostToDevice );

        //ds compute accelerations for all cells - launch as many threads as we have cells
        computeAccelerationsLennardJones<<< 1, uMaximumCellIndex >>>( uNumberOfParticles,
                                                                      vecParticles,
                                                                      d_arrCellIndexRange,
                                                                      dMinimumDistance,
                                                                      dPotentialDepth,
                                                                      uMaximumCellIndex,
                                                                      uMaximumNeighborCellIndexRange,
                                                                      d_arrNewAccelerations );

        //ds copy particles back to the domain - this call also updates the cell lists and changes the support structure
        //cDomain.setParticles( d_vecParticles );

        //ds record situation (we will write the stream to the file in one operation afterwards )
        cDomain.saveParticlesToStream( );
        cDomain.saveIntegralsToStream( dMinimumDistance, dPotentialDepth );
    }

    //ds deallocate memory
    cudaFree( d_arrCellIndexRange );
    cudaFree( d_arrNewAccelerations );

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
