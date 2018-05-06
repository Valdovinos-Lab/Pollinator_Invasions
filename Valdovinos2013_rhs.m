%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code created for simulating invasion of pollinator species into
% plant-pollinator networks, using Valdovinos et al. 2013, 2016 model
% Authors: Fernanda S. Valdovinos & Pablo Moisset de Espanes
% Last Modification: May 5, 2018
% Published: Nature Communications 2018 (see paper for more in detail
% explanation of parameters, variables, and simulation design)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Valdovinos et al. 2013 ODE's RHS. 
% New version. Avoids 0/0 cases to simplify simulations where species 
% are removed or added in the middle of the run

function dx = Valdovinos2013_rhs(t,x)

global network_metadata indRemP indRemA

plant_qty  = network_metadata.plant_qty ;
animal_qty = network_metadata.animal_qty ;
nz_pos  = network_metadata.nz_pos ;
e       = network_metadata.e ;
mu_p    = network_metadata.mu_p ;
mu_a    = network_metadata.mu_a ;
c       = network_metadata.c ;
b       = network_metadata.b ;
u       = network_metadata.u ;
w       = network_metadata.w ;
Beta    = network_metadata.Beta ;
G       = network_metadata.G ;
g       = network_metadata.g ;
phi     = network_metadata.phi ;
tau     = network_metadata.tau ;
epsilon = network_metadata.epsilon ;
In      = network_metadata.In ;

[p, N, a, Alpha] = unpack(x, network_metadata ) ;
p(indRemP)=0;% We force extinct species to zero to avoid resucitation due to stiff integrator
N(indRemP)=0;
a(indRemA)=0;

% Model's specific computation begins here

sigma = diag(sparse(p.*epsilon)) * Alpha ; %sigma nxm sparse
sigma = sigma * diag(sparse(1./(sum(sigma)+realmin))) ;

Gamma = g .* (1 - u'*p - w.*p + u.*p) ;% non sparse

tmp = Alpha * diag(sparse(a.*tau)) ; %tmp nxm sparse

dp = ( (Gamma .* sum(e .* sigma .* tmp, 2)) - mu_p) .* p;
tmp = (diag(sparse(N)) * tmp) .* b ;
da = sum(c .* tmp, 1)' - mu_a .* a ;
dN = Beta .* p - phi.*N - sum(tmp, 2) ;

% Adaptive dynamics starts here

% Fitness function
DH = diag(sparse(N)) * sparse(c.*b) ; %nxm sparse
DH(Alpha<0)=-DH(Alpha<0) ;

wavg = sum(Alpha.*DH) ; %Weights for average. nxm sparse

%This is the replicator equation
dAlpha = Alpha.*DH - Alpha*diag(sparse(wavg)) ;
dAlpha = dAlpha*diag(sparse(G)) ;


% Now pack the answer
dx = full([dp; dN; da; dAlpha(nz_pos)]) ;
