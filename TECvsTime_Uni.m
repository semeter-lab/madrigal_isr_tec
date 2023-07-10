% Universial TECvsTime function

function TECvsTime_Uni(data_file)
arguments
  data_file (1,1) string {mustBeFile}
end

%% load entire data file, can take large amount of RAM
% a priori knownledge of how Madrigal lays out HDF5 groups / datasets
compound_ds = "/Data/Table Layout";

% Matlab h5read only reads entire compound dataset.
% Python h5py trivially reads invididual fields of compound dataset.
tic
raw = h5read(data_file, compound_ds);

T.time = datetime(1970, 1, 1, 0, 0, 0) + seconds(raw.ut1_unix);
T.gdalt = raw.gdalt;
T.nel = raw.nel;

clear('raw')
T = struct2table(T);
%% sorting
T = sortrows(T, "time");
disp("Loaded and sorted data in " + num2str(toc, "%.1f") + " seconds.")

%% count number of valid rows of NEL and create array of each time value
goodRows = ~isnan(T.gdalt) & ~isnan(T.nel) & T.nel > 0 & ~isnat(T.time);
T = T(goodRows,:);

profTimes = unique(T.time);

% Reverse the log_10 for NE
T.nel = 10 .^ T.nel;

Nt = length(profTimes);

TEC = zeros(Nt, 1);
chapParams = zeros(Nt, 3);
z = 100:10:500;
H = 50;

%% For EACH PROFILE, integrate GDALT vs. NE for TEC values
tic
for i = 1:Nt
    NE = sortrows(T(T.time == profTimes(i),:), "gdalt");

    [x, ia] = unique(NE.gdalt);
    NE = NE(ia,:);
    % NOTE: only takes first element of unique altitude data.
    % previous code averaged same altitude data.

    [Nmax, I] = max(NE.nel);
    z0 = x(I);

    [chapParams(i,:), N] = ChapmanFit(NE.nel, x, z, z0, Nmax, H);

    TEC(i) = trapz(z, N);

    if ~mod(i,100)
        disp("Integrating " + num2str(i/Nt*100, "%.1f") + " %")
    end
end

disp("Integration took " + num2str(toc, "%.1f") + " seconds.")

%find Chapman fit outliers in the profiles
% [outliers,whereOutliers] = chapOutliers(chapParams);
sond_mmTEC = movmean(TEC, 250);

%Side-by-side display of profile
figure
sideProfile = plot(0,0);
hold on
sideProfPts = scatter(0,0);

%showProfileFcn generates profile of point clicked on TEC graph
figure
sondFig = gcf;
sondrestromPlot = plot(profTimes, sond_mmTEC);

set(sondrestromPlot,'HitTest','off')
set(gca,'ButtonDownFcn',{@showProfileFcn, profTimes, T})
datetick('x', 'yyyy-mm-dd');
xlabel('Time (UTC)'), ylabel('TEC (m^-^2)');
title('Total Electron Count vs Time');

%Shows x-coordinate of point being hovered over
t = annotation('textbox');
t.Position = [0.72 0.78 0.15 0.09];
hold on
pnt = scatter(T.time(1), sond_mmTEC(1), HitTest='off');
vert = xline(T.time(1), LineStyle="--", HitTest='off');
% set(sondFig, 'WindowButtonMotionFcn', {@hoverShowCoord, t, pnt, vert,...
%     date, sond_plotDate, sond_mmTEC})

%set(sondFig, 'WindowButtonMotionFcn', {@hoverShowProf, sideProfile, sond_plotDate, date, profSegs,...
%    GDALT, NE})
set(sondFig, WindowButtonMotionFcn={@hover, sideProfile, sideProfPts, T,...
    t, pnt, vert, profTimes, T.time(1), sond_mmTEC})


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

end
