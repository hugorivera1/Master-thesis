%% Import data
clear all

%Setting crystal symmetry (cs) and specimen symmetry (ssO)
cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ssO = specimenSymmetry('orthorhombic');

% Input file
datadir = 'D:\Master_thesis\EBSD\Hugo\J4';
file = fullfile(datadir, 'pattern.osc');

% Loading orientation data
ebsd1 = EBSD.load(file, cs, 'convertEuler2SpatialReferenceFrame','setting 2');

% Set reference frame
setMTEXpref('xAxisDirection', 'south');
setMTEXpref('zAxisDirection', 'outOfThePlane');
setMTEXpref('figSize','large')
% X || north        || TD
% Y || west         || RD
% Z || outOfPlane   || ND
rot = rotation.byAxisAngle(xvector,0*degree);
ebsd = rotate(ebsd1,rot,'keepEuler');

%% Filter away low IQ-measurements (Second phase particle data)

iqfilterVal = 950; %Insert appropriate filterval
ebsd(ebsd.imagequality < iqfilterVal).phase = -1;

%% Calculate ODF

%Setting crystal symmetry (cs) and specimen symmetry (ssO)
cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ssO = specimenSymmetry('orthorhombic');

% Calculate ODF
h = [Miller(1, 1, 1, cs{2}), Miller(1, 0, 0, cs{2}), Miller(1, 1, 0, cs{2})];
odf = calcDensity(ebsd('indexed').orientations, 'resolution', 7.5*degree,...
    'halfwidth', 5*degree);

% Rotate ODF so that these axes are parallel
% X || north        || RD: Old Y to new X: [0 1 0]
% Y || west         || TD: Old Z to new Y: [0 0 1]
% Z || OutOfPlane   || ND: Old X to new Z: [1 0 0]
rotODF = rotation.byMatrix([1 0 0; 0 0 1; 0 1 0]);
odf2 = rotate(odf, rotODF);

% Final rotation correcting for a slight misalignment of the sample in the
% microscope
[odf3, rot_inv3] = centerSpecimen(odf2, yvector);

% Setting correct specimensymmetry for ODF
odf3.SS = specimenSymmetry('orthorhombic');
%% Plotting phi2 sections

% Defining texture components
br = orientation.byEuler(35*degree, 45*degree, 90*degree, cs, ssO);
cu = orientation.byEuler(90*degree, 35*degree, 45*degree, cs, ssO);
cube = orientation.byEuler(0.001*degree, 0.001*degree, 0.001*degree, cs, ssO);
cubeND = orientation.byEuler(22*degree,0.001*degree,0.001*degree, cs, ssO);
cubeND45 = orientation.byEuler(45*degree, 0.001*degree, 0.001*degree, cs, ssO);
goss = orientation.byEuler(0, 45*degree, 0, cs, ssO);
p = orientation.byMiller([0 1 1], [1 2 2], cs, ssO);
q = orientation.byMiller([0 1 3], [2 3 1], cs, ssO);
s = orientation.byEuler(59*degree, 37*degree, 63*degree, cs, ssO);

% Setting contour levels for the phi2 sections
levelsODF = [0, 2, 4, 6, 8, 10, 15, 20, 25, 30];
phi2_colorlim = [0 30];

%Plotting phi2 sections and annotate texture components
figure
plot(odf3, 'phi2', [0 45 65]*degree, 'contourf', levelsODF, 'minmax')
CLim(gcm,phi2_colorlim);

hold on
annotate(br.symmetrise, 'marker', 'd', 'markerfacecolor', 'g','DisplayName','Br')
annotate(cu.symmetrise, 'marker', '^', 'markerfacecolor', 'b','DisplayName','Cu')
annotate(cube.symmetrise, 'marker', 's', 'markerfacecolor', 'r','DisplayName','Cube')
annotate(cubeND.symmetrise, 'marker', 's', 'markerfacecolor', 'c','DisplayName','CubeND')
annotate(goss.symmetrise, 'marker', 'o', 'markerfacecolor', 'k','DisplayName','Goss')
annotate(p.symmetrise, 'marker', '>', 'markerfacecolor',[1 0.6 0] ,'DisplayName','P')
annotate(s.symmetrise, 'marker', 'd', 'markerfacecolor', 'm','DisplayName','S')
hold on
[h,icons] = legend('Location','southoutside','Orientation','Horizontal');
n = ceil(numel(icons)/2);
newicons= icons(n+1:end);
for k=1:(length(newicons))
    newicons(k).Children.MarkerSize = 12;  
end
hold off
