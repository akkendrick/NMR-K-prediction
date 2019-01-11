function plotKwithDepth(K,z,Kest)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    figure
    box on
    grid on
    
    hold on
    plot(K,z,'+','MarkerSize',3)
    plot(Kest,z,'*','MarkerSize',3)

    xlabel('Hydraulic Conductivity (m/s)')
    ylabel('Depth (m)')
    set(gca, 'YDir','reverse')



end

