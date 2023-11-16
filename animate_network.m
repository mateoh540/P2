%% animate_network.m
%{
The purpose of this script is to allow the user to visualize how the
network evolves over the course of the simulation.  

                            *** NOTE *** 
Either run_lloyds.m or run_consensus.m must have been successfully 
executed BEFORE animate_network.m is run.  This script uses the t, P and G 
arrays to produce the visualization and so the data must already exist in
the Workspace. 
%}

%% PARAMETERS - can change w/o compromising script execution
MAX_EDGE_WEIGHT = 5;
UPDATE_LAG = 0.1;       % in seconds

%% SIMULATION
n=length(t)-1;
for i = 1:n 
    Gt = graph(G(:,:,i),'omitselfloops');
    wts = MAX_EDGE_WEIGHT * Gt.Edges.Weight/max(Gt.Edges.Weight);
    p = plot(Gt,'LineWidth',wts);
    p.XData = P(i,:,1);
    p.YData = P(i,:,2);
    title(strcat('t = ',num2str(t(i),2)));
    xlabel('x');
    ylabel('y');
    drawnow;
    pause(UPDATE_LAG);
end