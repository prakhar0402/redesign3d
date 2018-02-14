clear all
close all
clc

addpath(genpath('../../../alecjacobson-gptoolbox-00c124c'));

% inputFile = '../data/cactus/cactus';
% inputFile = '../data/part2_thins_uni/part2_thins_uni';
% inputFile = '../data/part2_thins_uni/part2_thins_uni_f001';
% inputFile = '../data/part2_thins_uni/part2_thins_uni_f001_f001';
% inputFile = '../data/puzzlepart/puzzlepart_uni';
% inputFile = '../data/fidgetflyer/fidgetflyer';
% inputFile = '../data/horse/horse2';
inputFile = '../data/rabbit/rabbit';
% inputFile = '../data/chess/chess1';
% inputFile = '../data/table/table';
% inputFile = '../data/chair/chair4';
% inputFile = '../data/batman/batman_5';
% inputFile = '../data/horse/horse2';

[V, F] = readOFF([inputFile, '.off']);

outputIdentifier = 'f002';
filename = [inputFile, '_', outputIdentifier];

[V, F] = readOFF([filename, '.off']);
fid = fopen([filename, '.out']);

line = fgetl(fid);
fileloc = fgetl(fid);
line = fgetl(fid);
identifier = fgetl(fid);
line = fgetl(fid);
T1 = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
K_SDF = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
K_S = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
K_D = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
VERTEX_MASS = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
TIME_STEP = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
TOTAL_TIME = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
STEPS = fscanf(fid, '%d\n', 1);
line = fgetl(fid);
min_sdf = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
max_sdf = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
nMove = fscanf(fid, '%d\n', 1);
line = fgetl(fid);
nV = fscanf(fid, '%d\n', 1);
line = fgetl(fid);
norm_sdf = fscanf(fid, '%f\n', nV);
line = fgetl(fid);
movable_index = fscanf(fid, '%f\n', nV);
line = fgetl(fid);
normals = fscanf(fid, '%f %f %f\n', [3 nV])';
line = fgetl(fid);
max_change = fscanf(fid, '%f\n', STEPS);

fclose(fid);

sdf = norm_sdf * (max_sdf - min_sdf) + min_sdf;

figure
patch('Faces',F,'Vertices',V,'FaceVertexCData',sdf,'FaceColor','interp');
colorbar
axis equal
xlabel({inputFile, outputIdentifier}, 'fontweight', 'bold')
set(gca, 'XTick', '', 'YTick', '')
box on

saveas(gcf, [filename, '.jpg'])

% T1 = 0.08;

remain_frac = sum(sdf < T1)/nV;
fprintf('Fraction of remaining vertices below threshold diameter = %f\n', remain_frac);

figure
patch('Faces',F,'Vertices',V,'FaceVertexCData',1*(sdf < T1) + 0*(movable_index > 0),'FaceColor','interp');
colorbar
axis equal
xlabel({inputFile, outputIdentifier}, 'fontweight', 'bold')
% set(gca, 'XTick', '', 'YTick', '')
box on
hold on
% quiver3(V(:, 1), V(:, 2), V(:, 3), normals(:, 1), normals(:, 2), normals(:, 3));

figure
plot(TIME_STEP:TIME_STEP:TIME_STEP*STEPS, max_change, 'LineWidth', 2)
set(gca, 'FontSize', 16);
xlabel('Time Steps', 'FontWeight', 'bold', 'FontSize', 20);
ylabel({'Maxima of','Vertex Displacements'}, 'FontWeight', 'bold', 'FontSize', 20)
axis square