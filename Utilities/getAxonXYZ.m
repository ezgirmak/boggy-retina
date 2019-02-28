function [X,Y,Z,AISInd,PL] = getAxonXYZ(tree,d_z,mode)

% Read in reconstruction data (swapping Y and Z)
X = tree.X;
Z = tree.Y;
Y = tree.Z;
R = tree.R;
type = typeN_tree(tree);

% Consider only axon points
X = X(R == 2);
Y = Y(R == 2);
Z = Z(R == 2);
R = R(R == 2);
type = type(R == 2);

X = [0; X];
Y = [0; Y];
Z = [0; Z];

% Find distance between adjacent points and terminate axon if > 100 um
adjDist = sqrt((X(2:end)-X(1:end-1)).^2 + ...
    (Y(2:end)-Y(1:end-1)).^2 + ...
    (Z(2:end)-Z(1:end-1)).^2);
adjDist(end+1) = 151;

axonTerm = find(adjDist > 150);
X = X(1:axonTerm(1));
Y = Y(1:axonTerm(1));
Z = Z(1:axonTerm(1));

% Get array of cumulative distance and find max distance
cumDist = zeros(1,length(X));

for j = 2:length(X)
    cumDist(j) = cumDist(j-1) + sqrt((X(j)-X(j-1))^2 + (Y(j)-Y(j-1))^2 + (Z(j)-Z(j-1))^2);
end

maxDist = max(cumDist);

% Find vector of last 50 um of axon
[~, ind_last] = min(abs(maxDist-cumDist-50));
last_proj = [X(end)-X(ind_last); 0; Z(end)-Z(ind_last)];

% Align the last axon segment with the z-axis (parallel fibres)
angle_proj = -atan2(last_proj(1),-last_proj(3));

rot_mat = [cos(angle_proj) 0 -sin(angle_proj);...
    0 1 0;...
    sin(angle_proj) 0 cos(angle_proj)];

XYZ_rot = zeros(length(X),3);

for i = 1:length(X)
    XYZ_rot(i,:) = rot_mat*[X(i); Y(i); Z(i)];
end

X = XYZ_rot(:,1);
Y = XYZ_rot(:,2);
Z = XYZ_rot(:,3);

% Flip the y-axis if necessary
if Y(1) > Y(end)
    Y = -Y;
end
if Z(1) > Z(end)
    Z = -Z;
end

% Ensure path length is a multiple of 10
frac10 = roundn(maxDist,1)/maxDist;
X = X*frac10;
Y = Y*frac10;
Z = Z*frac10;
maxDist = maxDist*frac10;

% Find AIS coordinates
[~, ind_AIS] = min(abs(cumDist-50));
AIS = [X(ind_AIS); Y(ind_AIS); Z(ind_AIS)];

% Sample axon at fixed arc lengths
XYZ = interparc(round(maxDist/d_z)+1,X,Y,Z,'linear');
X = XYZ(:,1)*1e-6;
Y = XYZ(:,2)*1e-6;
Z = XYZ(:,3)*1e-6;

if strcmpi(mode,'extend')
    
    % Append mirror image of axon to make it symmetrical
    X = [flipud(-X(2:end)); X];
    Y = [flipud(-Y(2:end)); Y];
    Z = [flipud(-Z(2:end)); Z];
    
    % Get path length
    PL = 2*maxDist;
    
    AISInd = knnsearch([X Y Z], AIS'*1e-6);
    
elseif strcmpi(mode,'centre')
    
    % Translate axon so the midpoint is at the origin
    X = X-X(round(length(X)/2));
    Y = Y-Y(round(length(X)/2));
    Z = Z-Z(round(length(X)/2));
    
    % Get path length
    PL = maxDist;
    
    AISInd = knnsearch([X-X(1) Y-Y(1) Z-Z(1)], AIS'*1e-6);
    
end