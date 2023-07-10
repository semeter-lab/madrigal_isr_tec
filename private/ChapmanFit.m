function [estimated_guess,N] = ChapmanFit(v,x,z,z0,Nmax,H)
arguments
    v (1,:) {mustBePositive}
    x (1,:) {mustBePositive}
    z (1,:) {mustBePositive}
    z0 (1,1) {mustBePositive}
    Nmax (1,1) {mustBePositive}
    H (1,1) {mustBePositive}
end

% plot(v,x,'ko');
% hold on;

start = [z0 Nmax H];



%ylim([0 6])
% set(gca, 'XMinorTick','on', 'XMinorGrid','on')
% set(gca, 'YMinorTick','on', 'YMinorGrid','on')



%start = [-2; 1; 4; 1; 2];
%options = optimset('TolX',0.1);
% opt = optimset(MaxFunEvals=1e6, MaxIter=1e3);
opt = optimset();
estimated_guess = fminsearch(@(guess) ChapmanErr(guess,x,v),start, opt);


N = chapman(z,estimated_guess(1),estimated_guess(2),estimated_guess(3));


% plot(N,z,'b-');
% hold off;
%
% xlabel('Electron Density')
% ylabel('Altitude (km)')
% legend('Data','Chapman',...
%     'Location','NorthEast','box','off')
% title('Chapman Least Squares Fitting')

end
