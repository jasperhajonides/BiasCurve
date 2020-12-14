function [biases angmus precs]=  Calculate_BiasCurve(response,target,other,itrl,config);
    
%This function calculates the a bias curve where the relative to another feature

% .  -Responses : the responses a participant makes on a 2*pi scale
% .  -Target    : the target feature, yet again on a 2*pi scale
% .  -Other     : the other feature you would expect a bias towards or away
%                     from. These should be on a 2*pi scale and should be
%                     independent from the target orientation.
%(c) Nick Myers 2018



%% Calculate BIAS
%get a few useful variables 

angdif = circ_dist(other,target)';%(tm{isub}(:,12))/180*pi*2;
rdif = circ_dist(response,target);

nbin = config.nbin; %
pbin  = config.pbin; %1/4;
i    = circ_bini(angdif(itrl),nbin,pbin);
for ibin = 1:nbin
    biases(ibin) = circ_mean(rdif(itrl(i(:,ibin))));
    angmus(ibin) = circ_mean(angdif(itrl(i(:,ibin)))');
    [~,sd] = circ_std(rdif(itrl(i(:,ibin))));
    precs(ibin) = 1./sd;
end
fprintf(' - done!\n');


