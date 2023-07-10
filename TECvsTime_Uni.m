% Universial TECvsTime function

%Takes parameter dataset, char vector of name of text file to analyze.
%Dataset must have the following columns, matching the column number:
%Col 1: NEL, Col 3: GDALT, Col 4: UT1_UNIX
%Make sure dataset is in directory


function TECvsTime_Uni(dataset)

unsortedT = readtable(dataset);
T = sortrows(unsortedT, 4); %sort based on time
clear unsortedT
Trows = height(T);

% count number of valid rows of NEL and create array of
% each time value
NELrows = 0;
profTimes = [];
timesCnt = 0;

for i = 1:Trows
    if ~isequaln(T.Var1(i), NaN) && ~isequaln(T.Var1(i), 0) && ~isequaln(T.Var4(i), 0) && ~isequaln(T.Var4(i), NaN)  %if valid NEL and time value
        NELrows = NELrows + 1;

        if isempty(profTimes) || (~isequaln(T.Var4(i), profTimes(timesCnt)))
            profTimes = cat(1,profTimes,T.Var4(i));
            timesCnt = timesCnt + 1;

        end
    end
end

% create array of number of NEL values at each
% specific time, numTimes
numTimes = [];
numNELatTime = 0;
thisTime = 1;
for i = 1:length(T.Var4)
    if T.Var4(i) == profTimes(thisTime)  %%equal to current time
       if ~isequaln(T.Var1(i), NaN) && ~isequaln(T.Var1(i), 0)
           numNELatTime = numNELatTime + 1;
       end
    elseif ~isequaln(T.Var4(i), 0) && ~isequaln(T.Var4(i), NaN)  %%equal to next time
        numTimes = cat(1,numTimes,numNELatTime);
        thisTime = thisTime + 1;
        if ~isequaln(T.Var1(i), NaN) && ~isequaln(T.Var1(i), 0)
            numNELatTime = 1;
        else
            numNELatTime = 0;
        end
    end
end
numTimes = cat(1,numTimes,numNELatTime);
numTimes = [numTimes, profTimes];


NEL = zeros(NELrows, 1);
GDALT = zeros(NELrows, 1);
% fill arrays of zeros with valid NEL values
% and corresponding GDALT values
NELcnt = 1;

for i = 1:Trows
    if ~isequaln(T.Var1(i), NaN) && ~isequaln(T.Var1(i), 0) && ~isequaln(T.Var3(i), NaN) && ~isequaln(T.Var3(i), 0)
        NEL(NELcnt) = T.Var1(i);
        GDALT(NELcnt) = T.Var3(i);
        NELcnt = NELcnt + 1;
    end
end

% Reverse the log_10 for NE
NE = 10 .^ NEL;

TEC = zeros(length(numTimes), 1);
chapParams = zeros(length(numTimes),3);
z = 100:10:500;
H = 50;
segBegin = 1;
profSegs = zeros(length(numTimes),2);
profSegs(1,1) = segBegin;
% For EACH PROFILE, integrate GDALT vs. NE for TEC values
for i = 1:length(numTimes)
    if numTimes(i) ~= 1
        segEnd = segBegin + numTimes(i) - 1;
        profSegs(i,2) = segEnd;

%         if i==22101 || i==12547 || i==439% || i==2148
        unsortedNE = [GDALT(segBegin:segEnd) NE(segBegin:segEnd)];
        sortedNE = sortrows(unsortedNE, 1);         %sort based on GDALT in order

        % average values at very close/same altitudes
        % then interpolate
        perAlt = 1;
        jBeg = 1;       %need these for indexing for Mvec
        jEnd = 1;
        x = [];     %vector of non-repeating altitudes
        v = [];     %vector of averaged Ne's
        for j = 2:length(sortedNE)
            if sortedNE(j,1) - sortedNE(j-1,1) <= 1 %group similar altitudes
                perAlt = perAlt + 1;
                jEnd = jEnd + 1;
            else
                Mvec = sortedNE(jBeg:jEnd, 2);  %average Ne of similar altitudes
                NeAvg = mean(Mvec);
                perAlt = 1;
                jEnd = jEnd + 1;
                jBeg = jEnd;
                x = cat(1,x,sortedNE(j-1,1));
                v = cat(1,v,NeAvg);
            end
        end
        Mvec = sortedNE(jBeg:jEnd, 2);  %does same as above else for last segment
        NeAvg = mean(Mvec);
        x = cat(1,x,sortedNE(length(sortedNE),1));
        v = cat(1,v,NeAvg);

        [Nmax,I] = max(v);
        z0 = x(I);

        [estimated_guess,N] = ChapmanFit(v,x,z,z0,Nmax,H);

        chapParams(i,:) = estimated_guess;


        TEC(i) = trapz(z,N);

%         GMTdate = datestr(times(i)/86400 + datenum(1970,1,1), 'dd-mmm-yyyy HH:MM:SS');
%
%         figure
%         plot(N,z,':.',v,x,'o');
%         xlabel('Electron Density (Ne in m-3)'), ylabel('Altitude (km)');
%         title(strcat('Ne vs. Alt w/ Chapman Interpolation (', GMTdate, ')'));

        % TEC(i) = trapz(xq,vq);

%         if (TEC(i) > 6e15)
%             disp("outlier at ");
%             disp(times(i));
%         end
    end
    segBegin = segBegin + numTimes(i);
    profSegs(i+1,1) = segBegin;
end

%find Chapman fit outliers in the profiles
[outliers,whereOutliers] = chapOutliers(chapParams);
sond_mmTEC = movmean(TEC, 250);

% convert unix time to universal time
date = strings(length(profTimes), 1);
for i = 1:length(date)
     date(i) = datestr(profTimes(i)/86400 + datenum(1970,1,1), 'dd-mmm-yyyy HH:MM:SS');
end
% convert date to numerical value for matlab plotting
sond_plotDate = datenum(date);

%Side-by-side display of profile
figure
sideProfile = plot(0,0);
hold on
sideProfPts = scatter(0,0);

%showProfileFcn generates profile of point clicked on TEC graph
figure
sondFig = gcf;
sondrestromPlot = plot(sond_plotDate, sond_mmTEC);

set(sondrestromPlot,'HitTest','off')
set(gca,'ButtonDownFcn',{@showProfileFcn, profTimes, profSegs, GDALT,...
    NE})
datetick('x', 'yyyy-mm-dd');
xlabel('Time (UTC)'), ylabel('TEC (m^-^2)');
title('Total Electron Count vs Time');

%Shows x-coordinate of point being hovered over
t = annotation('textbox');
t.Position = [0.72 0.78 0.15 0.09];
hold on
pnt = scatter(sond_plotDate(1),sond_mmTEC(1),'HitTest','off');
vert = xline(sond_plotDate(1),'LineStyle',"--",'HitTest','off');
% set(sondFig, 'WindowButtonMotionFcn', {@hoverShowCoord, t, pnt, vert,...
%     date, sond_plotDate, sond_mmTEC})

%set(sondFig, 'WindowButtonMotionFcn', {@hoverShowProf, sideProfile, sond_plotDate, date, profSegs,...
%    GDALT, NE})
set(sondFig, 'WindowButtonMotionFcn', {@hover, sideProfile, sideProfPts, profSegs, GDALT, NE,...
    t, pnt, vert, date, sond_plotDate, sond_mmTEC})




% hold on

% %profile outlier arrays, for plotting
% outlierTEC = [];
% outlierTime = [];
% for i = 1:length(whereOutliers)
%       outlierTEC = cat(1,outlierTEC,TEC(whereOutliers(i)));
%       outlierTime = cat(1,outlierTime,plotDate(whereOutliers(i)));
%
% end
%
% scatter(outlierTime,outlierTEC,8,'filled')
