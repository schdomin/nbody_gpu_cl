function createVideoFromSimulationFile( )

%ds filepath
strFilepath = '../bin/simulation.txt';

%ds open the file
fileID = fopen( strFilepath );

%ds get the first line
cCell = textscan( fileID, '%u %u', 1 );

%ds get number of particles N and timesteps T
uNumberOfParticles = cCell{1};
uNumberOfTimesteps = cCell{2};

%ds informative
disp( ['[N] Number of particles: ', num2str( uNumberOfParticles ) ] );
disp( ['[T] Number of timesteps: ', num2str( uNumberOfTimesteps ) ] );

%ds data structures NxT (each timestep is a column)
matX = zeros( uNumberOfParticles, uNumberOfTimesteps );
matY = zeros( uNumberOfParticles, uNumberOfTimesteps );
matZ = zeros( uNumberOfParticles, uNumberOfTimesteps );
matU = zeros( uNumberOfParticles, uNumberOfTimesteps );
matV = zeros( uNumberOfParticles, uNumberOfTimesteps );
matW = zeros( uNumberOfParticles, uNumberOfTimesteps );

disp( [ 'starting data import from: ', strFilepath ] ); 
tic;

%ds for each timestep
for uCurrentTimestep = 1:1:uNumberOfTimesteps
   
	%ds get the current line
    cCell = textscan( fileID, '%f %f %f %f %f %f', uNumberOfParticles );

    %ds save the information in the respective vectors
    matX( :, uCurrentTimestep ) = cCell{1};
    matY( :, uCurrentTimestep ) = cCell{2};
    matZ( :, uCurrentTimestep ) = cCell{3};
    matU( :, uCurrentTimestep ) = cCell{4};
    matV( :, uCurrentTimestep ) = cCell{5};
    matW( :, uCurrentTimestep ) = cCell{6};   
    
end;

disp( [ 'finished data import - time: ', num2str( toc ) ] );

%ds calculate the absolute velocity of each particle for each timestep
matAbsoluteVelocity = sqrt( matU.^2 + matV.^2 + matW.^2 );

%ds get the maximum of them - required to colorize particles
dMaximumAbsoluteVelocity = max( max( matAbsoluteVelocity ) );
disp( [ 'maximum absolute velocity = ', num2str( dMaximumAbsoluteVelocity ) ] );

%ds scale the velocity matrix to [0,1] for colorization
matAbsoluteVelocity = matAbsoluteVelocity./dMaximumAbsoluteVelocity;

%ds prepare video writer
%ds information and code snippets from:
%   http://www.mathworks.ch/ch/help/matlab/ref/videowriterclass.html
writerObj = VideoWriter('simulation.avi');
writerObj.FrameRate = 25;
open( writerObj );
disp( [ 'timestep: ', num2str( 1 ) ] );

%ds get the color vector
vecColor = [ matAbsoluteVelocity( :, 1 ), zeros( uNumberOfParticles, 1 ), 1-matAbsoluteVelocity( :, 1 ) ];

%ds create initial data with the first timestep
scatter3( matX( :, 1 ), matY( :, 1 ), matZ( :, 1 ), 25, vecColor, 'fill', 'MarkerEdgeColor', 'black' );
axis( [ -1, 1, -1, 1, -1, 1 ] );
set( gca, 'nextplot' ,'replacechildren' );
set( gcf, 'Renderer' ,'zbuffer' );

% %ds for each remaining timestep
for uCurrentTimestep = 2:1:uNumberOfTimesteps
    
    %ds only save every 20th frame instead 500 -> 25 frames per sec
    if mod( uCurrentTimestep, 20 ) == 0
        
        %ds color vector
        vecColor = [ matAbsoluteVelocity( :, uCurrentTimestep ), zeros( uNumberOfParticles, 1 ), 1-matAbsoluteVelocity( :, uCurrentTimestep ) ];
    
        %ds create a figure
        scatter3( matX( :, uCurrentTimestep ), matY( :, uCurrentTimestep ), matZ( :, uCurrentTimestep ), 25, vecColor, 'fill', 'MarkerEdgeColor', 'black' );
        frame = getframe( gcf );
        writeVideo( writerObj, frame );
    
    end

    disp( [ 'timestep: ', num2str( uCurrentTimestep ) ] );
    
end

%ds write video file
close( writerObj );
