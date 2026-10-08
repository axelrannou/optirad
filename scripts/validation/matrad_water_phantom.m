function matrad_water_phantom(outFile, matRadRoot)
% Water-phantom reference: same geometry as optirad_phantom (WaterPhantom.cpp), computed with matRad.
% Writes the dose as an ODOSE1 file (see src/validation/DoseRawIO.hpp).
%   matlab -batch "matrad_water_phantom('/tmp/matrad_water.odose', '/path/to/matRad')"

if nargin < 2
    matRadRoot = fullfile(getenv('HOME'), 'Documents', 'Workspace', 'matrad_env', 'matRad');
end
run(fullfile(matRadRoot, 'matRad_rc.m'));

% Keep in sync with WaterPhantomConfig. The local matRad copy has debug code indexing dose voxel
% (148, 79, 112), so the cube must be at least that large.
res = 2.5; halfLat = 200; depth = 400; isoDepth = 100; field = 100; bixel = 5;
nx = round(2 * halfLat / res); ny = round(depth / res); nz = nx;

ct = struct();
ct.cubeDim = [ny nx nz];                   % matRad cube order: [y, x, z]
ct.resolution = struct('x', res, 'y', res, 'z', res);
ct.x = -halfLat + (0:nx-1) * res;
ct.y = -isoDepth + (0:ny-1) * res;
ct.z = -halfLat + (0:nz-1) * res;
ct.numOfCtScen = 1;
ct.cubeHU = {zeros(ct.cubeDim)};

[X, Y, Z] = meshgrid(ct.x, ct.y, ct.z);    % meshgrid output is [y, x, z]
target = abs(X) <= field/2 & abs(Z) <= field/2 & abs(Y) <= 2*res;

cst = cell(2, 6);
cst(1,:) = {1, 'PTV',  'TARGET', {find(target(:))}, struct('Priority', 1, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', true), []};
cst(2,:) = {2, 'BODY', 'OAR',    {(1:numel(target))'}, struct('Priority', 2, 'alphaX', 0.1, 'betaX', 0.05, 'Visible', true), []};

pln.radiationMode = 'photons';
pln.machine = 'Generic';
pln.bioModel = 'none';
pln.multScen = 'nomScen';
pln.numOfFractions = 1;
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = bixel;
pln.propStf.numOfBeams = 1;
pln.propStf.isoCenter = [0 0 0];
pln.propDoseCalc.doseGrid.resolution = struct('x', res, 'y', res, 'z', res);

stf = matRad_generateStf(ct, cst, pln);
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
result = matRad_calcCubes(ones(dij.totalNumOfBixels, 1), dij);
dose = double(result.physicalDose);

dg = dij.doseGrid;
f = fopen(outFile, 'w');
assert(f > 0, 'cannot open %s', outFile);
fwrite(f, [uint8('ODOSE1') 0 0], 'uint8');
fwrite(f, int32(size(dose)), 'int32');
fwrite(f, [dg.resolution.y dg.resolution.x dg.resolution.z], 'double');
fwrite(f, [dg.x(1) dg.y(1) dg.z(1)], 'double');
fwrite(f, dose(:), 'double');
fclose(f);
fprintf('Written %s (max dose %.4f, bixels %d)\n', outFile, max(dose(:)), dij.totalNumOfBixels);
end
