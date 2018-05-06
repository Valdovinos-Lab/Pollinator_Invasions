%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code created for simulating invasion of pollinator species into
% plant-pollinator networks, using Valdovinos et al. 2013, 2016 model
% Authors: Fernanda S. Valdovinos & Pablo Moisset de Espanes
% Last Modification: May 5, 2018
% Published: Nature Communications 2018 (see paper for more in detail
% explanation of parameters, variables, and simulation design)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs (change manually below to run in your PC, use
% 'Run_Invasions_cluster' to vary them automatically in a cluster
% frG_native: 1 (all natives are adaptive foragers)
% frG_Inv: 0 & 1 (adaptive foraging of invaders: 0 without, 1 with)
% k_level: 0 & 1 (generalism level of invaders: 0 specialist, 1 generalist)
% opW: 0, 1 & 2 (invader connect to: 0 random, 1 least, 2 most connected plant
% ftauI: 1 & 2 (invader foraging efficiency: 1 low, 2 high)
% muAP: 1, 3 (animal mortality scenario: 1 high, 3 low)
% indxP: Was left there for technical reasons, just leave it empty: []
% indxA: Indicate the pollinator species that is introduced, which in this
% code is always the firs pollinator in the matrix, so leave it as 1.
% vectG: defines which pollinator species exhibit adaptive foraging.
% In: Incidence matrix, i.e. who pollinates whom.
% Defines the mortality scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VP, VA, A]=Run_Invasions_PC(r_i)

global J_pattern

frG_native=1;
frG_Inv=1;
k_level=1;
opW=0;
ftauI=2;
muAP=1;
sem=0;
indxA=1;
indxP=[];

rand('seed',sem+r_i);

M=load(sprintf('1200_networks/m%04d.txt',r_i));% change here the network to be studied
[sA, iA]=sort(sum(M));% sorting pollinators by their ascending degree
[sP, iP]=sort(sum(M,2));% sorting plants by their ascending degree
M=M(iP,iA);% Network is sorted by ascending degree of plant and animal species
           % That is: Most speacialized plants in the first rows (top of the matrix)
           % and most generalized plants in the last rows (bottom of the matrix)

[In, indx]=invasionA(M,k_level,opW);
    
[rows, cols]= size(In);
J_pattern = J_zero_pattern(In) ;

vectG_native=frG_native*ones(1,cols-1);
vectG_Inv=frG_Inv;
vectG=[vectG_Inv vectG_native];
    
[VP, VA, A]=Val_InvA(indxP,indxA,vectG,In,muAP,ftauI);

out_dir='results/'; 
save(sprintf('%sVP%d_vGinv%d_K%d_T%d_muAP%d_m%04d.mat',out_dir,frG_native,frG_Inv,k_level,ftauI,muAP,r_i),'VP')
save(sprintf('%sVA%d_vGinv%d_K%d_T%d_muAP%d_m%04d.mat',out_dir,frG_native,frG_Inv,k_level,ftauI,muAP,r_i),'VA')
save(sprintf('%sAlpha%d_vGinv%d_K%d_muAP%d_T%d_m%04d.mat',out_dir,frG_native,frG_Inv,k_level,ftauI,muAP,r_i),'A')
    
end
