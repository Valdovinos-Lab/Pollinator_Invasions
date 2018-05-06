%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code created for simulating invasion of pollinator species into
% plant-pollinator networks, using Valdovinos et al. 2013, 2016 model
% Authors: Fernanda S. Valdovinos & Pablo Moisset de Espanes
% Last Modification: May 5, 2018
% Published: Nature Communications 2018 (see paper for more in detail
% explanation of parameters, variables, and simulation design)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function introduces alien pollinators into the native network
% Introduced alien can be generalist (k_level=1) or specialists (k_level=0),
% and they can be linked more likely to generalist natives (opW=1),
% to specialists (opW=2) or randomly to any native (opW=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [In, indx]=invasionA(M,k_level,opW)

[m, n]=size(M);
    
% Defining the degree of native species
knative=sum(M);% same guild (animals)
knative2=sum(M,2)';% the other guild (plants)
InvM=zeros(m,1); % preparing the sub-matrix for the alien
fk=round(length(knative)*0.3);
    
% Defining the degree of the introduced species    

    ascend_knative=sort(knative);
    
    if k_level==0
        k=round(mean(ascend_knative(1:fk)));% degree of alien is equal to
                                            % the mean of 30% most
                                            % specialist natives
    elseif k_level==1
        k=round(mean(ascend_knative(end-fk:end)));% degree of alien is equal
                                            % to the mean of 30% most
                                            % generalist natives
    end
    
% Defining the species to which the introduced pollinator is linked to    
    native_label2=1:length(knative2);
    knative2_label=[native_label2; knative2];    

  if opW==0
     R=randperm(length(knative2));
     Links_Inv(1,:)=R(1:k);
    
  elseif opW==1
         
     Links_Inv=zeros(nInv,k);
     knative2_labelB=knative2_label;
         
     for l=1:k
         Pm=knative2_labelB(2,:)./sum(knative2_labelB(2,:));% p porportional to the degree, normalized
         Pm2=suma_one(Pm);

         R = mnrnd(1,Pm2');   % rand numbers from multinomial distribution
         fW=find(R);          % find the chosen index ofthe linked native 
            
         Links_Inv(i_nInv,l)=knative2_labelB(1,fW);
         knative2_labelB(:,fW)=[];
                        
     end

    elseif opW==2
        [degree, sp_label]=sort(knative2,'descend');
        kmin=min(degree(1:k));
        generalists=sp_label(degree>=kmin);
        perm_generalists=generalists(randperm(length(generalists)));
        Links_Inv=perm_generalists(1:k);
    
   end
            
% Introducing the alien pollinator
  InvM(Links_Inv(1,:),1)=1;
  In=[InvM M];
  indx=1;
  
end

    
    

