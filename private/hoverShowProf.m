function hoverShowProf (~, ~, sideProfile, sideProfPts, xData, date, T)

C = get (gca, 'CurrentPoint');

ind = find(abs(floor(C(1,1)) - xData) == min(abs(floor(C(1,1)) - xData)));


% disp("index: ")
% disp(ind)

NE = sortrows(T(T.time == profTimes(ind),:), "gdalt");

% TEC = zeros(length(numTimes), 1);
% chapParams = zeros(length(numTimes),3);
z = 100:10:500;
H = 50;

% average values at very close/same altitudes
% then interpolate
perAlt = 1;
jBeg = 1;       %need these for indexing for Mvec
jEnd = 1;
x = [];     %vector of non-repeating altitudes
v = [];     %vector of averaged Ne's
for j = 2:length(NE)
    if NE(j,1) - NE(j-1,1) <= 1 %group similar altitudes
        perAlt = perAlt + 1;
        jEnd = jEnd + 1;
    else
        Mvec = NE(jBeg:jEnd, 2);  %average Ne of similar altitudes
        NeAvg = mean(Mvec);
        perAlt = 1;
        jEnd = jEnd + 1;
        jBeg = jEnd;
        x = cat(1,x,NE(j-1,1));
        v = cat(1,v,NeAvg);
    end
end

Mvec = NE(jBeg:jEnd, 2);  %does same as above else for last segment
NeAvg = mean(Mvec);
x = cat(1,x,NE(length(NE),1));
v = cat(1,v,NeAvg);

[Nmax,I] = max(v);
z0 = x(I);

[estimated_guess,N] = ChapmanFit(v,x,z,z0,Nmax,H);

% chapParams(i,:) = estimated_guess;


% TEC(i) = trapz(z,N);

% date(i) = datestr(times(i)/86400 + datenum(1970,1,1), 'dd-mmm-yyyy HH:MM:SS');
% figure
% plot(N,z)
sideProfile.XData = N;
sideProfile.YData = z;
sideProfPts.XData = NE(:,2);
sideProfPts.YData = NE(:,1);
xlabel('Electron density (Ne in m-3)'), ylabel('Altitude (km)')
timeStr = date(ind);%datestr(pt(1), 'dd-mmm-yyyy HH:MM:SS');
timeStrtxt = string(timeStr);
title('Profile w/ Chapman Fit at ' + timeStrtxt)


end % <--- optional if this is embedded into a function
