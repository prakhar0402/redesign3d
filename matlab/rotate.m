clear all
close all
clc

addpath(genpath('../../../alecjacobson-gptoolbox-00c124c'));

tx = 0;
ty = 0;
tz = 0;

%% part2_thin
% inputFile = '../data/part2_thins_uni/part2_thins_uni';
% 
% [V, F] = readOFF([inputFile, '.off']);
% 
% tx = -20;
% ty = 30;
% tz = 30;

% V = V*rotx(tx)*roty(ty)*rotz(tz);

%% cactus
inputFile = '../data/cactus/cactus_f001';

[V, F] = readOFF([inputFile, '.off']);

tx = 30;
ty = 30;
tz = 0;

V = V*roty(ty)*rotz(tz)*rotx(tx);

%% rabbit
inputFile = '../data/rabbit/rabbit_f002';

[V, F] = readOFF([inputFile, '.off']);

tx = 0;
ty = -30;
tz = 0;

V = V*rotx(tx)*roty(ty)*rotz(tz);

%% chair
inputFile = '../data/chair/chair4_f001';

[V, F] = readOFF([inputFile, '.off']);

tx = 60;
ty = 60;
tz = 0;

V = V*roty(ty)*rotz(tz)*rotx(tx);

%% horse
inputFile = '../data/horse/horse2_f001';

[V, F] = readOFF([inputFile, '.off']);

tx = 0;
ty = -30;
tz = 0;

V = V*roty(ty)*rotz(tz)*rotx(tx);

%%


writeOFF([inputFile, '_rot.off'], V, F);

