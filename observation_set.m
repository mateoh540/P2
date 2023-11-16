function obset_data = observation_set(p, x, y, d, ro)
%{ 
Let "n" be the number of agents in the observation set.
Let "m' be the length of the arena
Let "w" be the number of data points in the observation set

--- Inputs ---
p: n by 2 array holding position data for agents in obs set
   e.g., p(i,2) is the y coordinate for agent-i

x, y, d: m by m meshgrids holding arena data 

ro: radius of observation (scalar) used by all agents

--- Outputs ---
obset_data: w by 2 array holding all weighted data in observation set.  
   e.g., z(:,2) are all of the y data in the obs set weighted by density

--- Hints ---
1:  Make sure that you are not accidentally including duplicate xy
points in your output.  Depending on how your design your code, this is 
easy to miss when the agents' area of observation overlap.  Consider 
using the unique() function to resolve. 

2: Running x.*d will weight the x meshgrid by the density meshgrid 
%}



end