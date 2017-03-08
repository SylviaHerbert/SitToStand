function plotInitialStates(g,data,schemeData)
max_v = (pi/8);       % allowing for some sway
standing_min = [-pi/15, -max_v, -pi/15, -max_v];
standing_max = [pi/15, max_v, 0.15, max_v];
data0 = shapeRectangleByCorners(g, standing_min, standing_max);


[gPos,dataPos]=proj(g,data,[0 1 0 1],'min');
[~,dataPos0]=proj(g,data0,[0 1 0 1],'min');

figure(5)
clf
[~,hL] = contour(gPos.xs{1},gPos.xs{2},dataPos0,[0 0],'Color','r');
hold on
[~,hV] = contour(gPos.xs{1},gPos.xs{2},dataPos,[0 0]),'Color','b';

figure(6)
clf
subplot(1,2,1)
hCost = surfc(gPos.xs{1},gPos.xs{2},dataPos0);
subplot(1,2,2)
hValue = surfc(gPos.xs{1},gPos.xs{2},dataPos);
 

figure(7)
clf
BRT = dataPos<=0;
BRT0 = dataPos0 <=0;

ind = find(BRT == 1);
for i = 1:length(ind)
  [indK, indH]=ind2sub([21,21],ind(i));
  z(1) = gPos.vs{1}(indK);
  z(2) = gPos.vs{2}(indH);
  [x,y]=transferCoordinates(z,schemeData);
  h{i}=plot(x,y,'Linewidth',3);
  hold on
  axis([-1 1 0 1.5])
end

end
function [x,y]=transferCoordinates(z,schemeData)
x=zeros(1,4);
y=zeros(1,4);
L = [0 schemeData.dynSys.L0 schemeData.dynSys.L1 schemeData.dynSys.R2]; %length between each point
y(2)=L(2); %draw a straight line up from angle to knee
x(3) = x(2)-L(3)*sin(z(1));
y(3) = y(2) + L(3)*cos(z(1));
x(4) = x(3)-L(4)*sin(z(1)+z(2));
y(4) = y(3)+L(4)*cos(z(1)+z(2));
end