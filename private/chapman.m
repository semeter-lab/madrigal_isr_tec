%*****************************************************
% Return a Chapman-layer ionosphere given parameters z0,Nmax,H
%*****************************************************
function N=chapman(z,z0,Nmax,H)
arguments
    z (1,:) {mustBeNumeric, mustBePositive}
    z0 (1,1) {mustBeNumeric, mustBePositive}
    Nmax (1,1) {mustBeNumeric, mustBePositive}
    H (1,1) {mustBeNumeric, mustBePositive}
end

% INPUTS:
% independent variable (altitudes to be evaluated)
% Altitude of ionospheric peak (km)
% Peak density (1/m^3)
% characteristic width of the layer (km)
%
% OUTPUTS:
% N(z)
%
% Nominal values for F-layer, z0=300, Nmax=1e12, H=65
%

N=Nmax*exp(.5*(1-(z-z0)/H - exp((z0-z)/H)));

return
