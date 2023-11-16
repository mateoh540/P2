function [P1, E1] = move_agents(P0, Pdot, E0, dt)
%{  
Let "n" be the number of agents.

--- Inputs ---
P0: 1 by n by 2 array holding positional data for all agents at timestep i
    e.g, P0(1,j,2) is the y position data for agent-j at timestep i 

Pdot:  1 by n by 2 array holding positional data for all agents at 
          timestep i+1 BEFORE it has been constrained. Same structure as
          P0.

E0:  1 by n array holding energy data for all agents at timestep i
    e.g., E0(j) is the energy of agent-j at timestep i 

dt: timestep (scalar)

--- Outputs ---
P1:  1 by n by 2 array holding positional data for all agents at 
         timestep i+1 AFTER it has been constrained. Same structure as
          P0.

E1: 1 by n array holding energy data for all agents at timestep i+1.
    Same structure as E0.
%}



end