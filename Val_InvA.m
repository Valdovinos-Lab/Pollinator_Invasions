%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code created for simulating invasion of pollinator species into
% plant-pollinator networks, using Valdovinos et al. 2013, 2016 model
% Authors: Fernanda S. Valdovinos & Pablo Moisset de Espanes
% Last Modification: May 5, 2018
% Published: Nature Communications 2018 (see paper for more in detail
% explanation of parameters, variables, and simulation design)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function (Val_InvA.m) generates parameters of the dynamic model and 
% safe them in the global structure 'network_metadata' (using
% 'create_metadata.m') so they can be called repeatedly, and integrates the
% ordinary differential equations of the model specified in
% 'Valdovinos2013_rhs.m'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% indxP: Was left there for technical reasons, just leave it empty: []
% indxA: Indicate the pollinator species that is introduced, which in this
% code is always the firs pollinator in the matrix, so leave it as 1.
% vectG: defines which pollinator species exhibit adaptive foraging.
% In: Incidence matrix, i.e. who pollinates whom.
% Defines the mortality scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VP, VA, A]=Val_InvA(indxP,indxA,vectG,In,muAP,ftauI)

global network_metadata indRemP indRemA J_pattern

tmax=10000;
EUp=2e-2;
EUa=1e-3;

[m, n]=size(In);
B=sparse(In);
VP=cell(1,2);
VA=cell(1,2);
A=cell(1,2);

% Parameters of the uniform distribution from where the parameters of the
% dynamic model are drawn:
varp=1e-1;% variance of plant parameters
vara=0;% variance of animal parameters
mC=0.2; vC=vara;
mE=0.8; vE=varp;
mb=0.4; vb=vara;
mU=0.06; vU=varp;
mw=1.2; vw=varp;
mB=0.2; vB=varp;
mG=2; vG=vara;
mg=0.4; vg=varp;
mphi=0.04; vphi=varp;
mtau=1; vtau=vara;
mepsilon=1; vepsilon=varp;
vmA=vara; vmP=varp;

if muAP==1
    mmA=0.05; mmP=0.001; % high pollinator mortality 
elseif muAP==2
    mmA=0.001; mmP=0.02; % high plant mortality
elseif muAP==3
    mmA=0.001; mmP=0.001; % low plant and animal mortality
elseif muAP==4
    mmA=0.03; mmP=0.005; % high plant and animal mortality
end

% Parameters are drawn from uniform distribution (see Valdovinos et al.
% 2018, Nature Communications for complete description of the model and
% parameter definition)

% (10%meanP)-meanP+(10%meanP); (0.01%meanA)-meanA+(0.01%meanA)
c=uniform_rand(mC,vC,m,n).*B;
e=uniform_rand(mE,vE,m,n).*B;
b=uniform_rand(mb,vb,m,n).*B;

u=uniform_rand(mU,vU,m,1);
Beta=uniform_rand(mB,vB,m,1);
G=uniform_rand(mG,vG,n,1).*vectG';
g=uniform_rand(mg,vg,m,1);
mu_a=uniform_rand(mmA,vmA,n,1);
mu_p=uniform_rand(mmP,vmP,m,1);
w=uniform_rand(mw,vw,m,1);
phi=uniform_rand(mphi,vphi,m,1);
epsilon=uniform_rand(mepsilon,vepsilon,m,1);

tau=uniform_rand(mtau,vtau,n,1);
tau(indxA)=ftauI*tau(indxA);

%Create structure 
network_metadata = create_metadata(B, e, mu_p, mu_a, c, b, u, w, Beta, G, g, phi, tau, epsilon) ;

%Give initial state
mz=0.5; vz=1e-1;

yzero=uniform_rand(mz,vz,2*m+n,1);

initial_plants=yzero(1:m);
initial_nectar=yzero(m+1:2*m);

initial_animals=yzero(2*m+1:2*m+n);
initial_animals(indxA)=0;
initial_alphas=B;

%Normalization and packing.
initial_alphas=initial_alphas*diag(sum(initial_alphas).^(-1));
saving_initial_alpha_InvA=initial_alphas(:,indxA);% to use it when invasion starts
initial_alphas(:,indxA)=0;

initial_alphas=initial_alphas(network_metadata.nz_pos) ;

% Combining all initial variables
initial_state=full([initial_plants;initial_nectar;initial_animals;initial_alphas]);
tspan = [0 tmax];

indRemP=indxP; % I need this global variable for the sistem of equations (polinizacion rhs) to keep the species in 0
indRemA=indxA;

%% Integrating the dynamic model for the first tmax=10,000 time steps
options = odeset('JPattern', J_pattern,'NonNegative',1:2*m+n) ;
[t, y]=ode15s(@Valdovinos2013_rhs,tspan,initial_state, options) ;

yf = y(end,:)';

% Retriving final densities and foraging efforts for t=10,000
[plantsf, nectarf, animalsf, alphasf] = unpack(yf, network_metadata);

tmp = alphasf * diag(sparse(animalsf.*tau)) ;
sigma = diag(sparse(plantsf.*epsilon)) * alphasf ; %sigma nxm sparse
sigma = sigma * diag(sparse(1./(sum(sigma)+realmin))) ;
pol_event= sum(sigma .* tmp, 2);% per-capita summed pollination events

tmp = alphasf * diag(sparse(tau)) ;
tmp = (diag(sparse(nectarf)) * tmp) .* b;
N_extract = sum(tmp)'; % sum of the resources that each individual extracts

% Plant and animal species that went extinct
vzp=plantsf<EUp;
vza=animalsf<EUa;

% Saving results for the first tmax=10,000 time steps
A{1}=alphasf;
VP{1}=[vzp plantsf nectarf pol_event];
VA{1}=[vza animalsf N_extract];

% Checking everything is correct
assert(all(abs(plantsf(indxP))<1e-4));
assert(all(abs(nectarf(indxP))<1e-4));
assert(all(abs(animalsf(indxA))<1e-4));
assert(all(all(alphasf>-1e-4 & alphasf<1.0001)));
assert( all( abs(sum(alphasf(:,indxA+1:end)) - 1.0)<1e-3 ) ) ;

% Keeping extinct species in zero for the next dynamic run
indxA2 = find(vza);
indxA2 = indxA2(2:end);

indxP2=indxP;
indRemP=indxP2; % I need this global variable for the sistem of equations (Valdovinos2013_rhs) to keep the species in 0
indRemA=indxA2;

initial_plants=plantsf;
initial_plants(indxP2)=0;% keeps extinct plants in zero

initial_nectar=nectarf;
initial_nectar(indxP2)=0;% keeps the floral rewards of extinct plants in zero

initial_animals=animalsf;

%% Introducing the alien pollinator
initial_animals(indxA)=EUa+0.5e-3;% the initial density of the introduced animal is slightly above the extinction treshold for animals
initial_animals(indxA2)=0;% keeps extinct pollinators in zero

alphasf(:,indxA)=full(saving_initial_alpha_InvA);
initial_alphas=alphasf(network_metadata.nz_pos) ;

%Combine initial variables to run the second 10,000 time steps of the model
%after introducing the alien pollinator
initial_state=full([initial_plants;initial_nectar;initial_animals;initial_alphas]);

[t2, y2]=ode15s(@Valdovinos2013_rhs,tspan,initial_state,options) ;

yf2 = y2(end,:)';

% Retriving final densities and foraging efforts for t=20,000
[p, N, a, Alpha] = unpack(yf2, network_metadata ) ;

assert(all(all(Alpha>-1e-4 & Alpha<1.0001)));

%% Mechanistic variables at t=20,000

% Plants
tmp = Alpha * diag(sparse(a.*tau)) ;
sigma = diag(sparse(p.*epsilon)) * Alpha ; %sigma nxm sparse
sigma = sigma * diag(sparse(1./(sum(sigma)+realmin))) ;
pol_event= sum(sigma .* tmp, 2);

% Animals
tmp = Alpha * diag(sparse(tau)) ;
tmp = (diag(sparse(N)) * tmp) .* b;
N_extract = sum(tmp)'; % sum of the resources that each individual extracts


%% Number of extinct plant and animal species (below EU)
vzp=p<EUp;
vza=a<EUa;

%% Saving results

A{2}=Alpha;
VP{2}=[vzp p N pol_event];
VA{2}=[vza a N_extract];

end
