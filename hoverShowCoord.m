function hoverShowCoord (~, ~, t, pnt, vert, date, xData, yData)

C = get (gca, 'CurrentPoint');

% title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);

ind = find(abs(floor(C(1,1)) - xData) == min(abs(floor(C(1,1)) - xData)));

pnt.XData = xData(ind);
pnt.YData = yData(ind);
vert.Value = xData(ind);


t.String = ["Time: " + num2str(date(ind)), "TEC: " + sprintf('%.4e',yData(ind))];

% annotation('textbox', [0.75, 0.1, 0.1, 0.1], 'String', "pi value is " + pi)

end
