
% NOTE: Before running this script, first provide the path to the compiled Java
% code for the model, as in the example below:
% MODEL_PATH='C:\Users\bortk\Documents\NetBeansProjects\Sim\dist\Sim.jar';

MODEL_PATH='C:\Users\bortk\Documents\NetBeansProjects\Sim\dist\Sim.jar';
javaaddpath(MODEL_PATH);

% NOTE: Also provide paths to the mat files "LGNToCortex.mat" and
% "fbInputFileUnif.mat" which came with the model code, as in the example
% below:
% LGNInputFile = 'C:\Users\bortk\Downloads\JNSModel\LGNToCortex.mat';
% fbInputFile = 'C:\Users\bortk\Downloads\JNSModel\fbInputFileUnif.mat';

LGNInputFile = 'C:\Users\bortk\Downloads\JNSModel\LGNToCortex.mat';
fbInputFile = 'C:\Users\bortk\Downloads\JNSModel\fbInputFileUnif.mat';

numHCRows=3;
numExcitatoryPerHC=3000;
numInhibitoryPerHC=1000;
numLGNBlockRows=6;
numLGNPerRowPerBlock=2;
numLGNPerColPerBlock=2.5;
connectivityEE=.15;
connectivityIE=.6;
connectivityEI=.6;
connectivityII=.6;
egalEHelper=(.05 : (-.05 - .05) / 10 : -.05)';
egalIHelper=0.0;
%LGNInputFile = 'LGNToCortex16Even3INoise.mat';


disp('setting main network parameters...')
a = sim.Network(numHCRows,numExcitatoryPerHC,numInhibitoryPerHC,numLGNBlockRows,numLGNPerRowPerBlock,numLGNPerColPerBlock);

highegal = @(n) (-0.14/4)*(n-4) + 0.15;
lowegal  = @(n) (-0.08/8)*(n-4) + 0.15;

%egal = [.2 .2 .2  .19 .18 .14 .11 .1 .1];
%egal = [.22 .22 .22 .20 .18 .14 .11 .1 .1];
%egal = [.24 .24 .24 .21 .18 .15 .12 .1 .1];
%egal = [.2 .2 .2 .18 .16 .14 .12 .1 .1];
%egal = [1 1 1 1 1 1 1 1 1] * 0.15;

%egal = [.28 .26 .24 .21 .18 .15 .12 .1 .1];

%egal = linspace(.3, 0, 9);
%egal  = [.24 .23 .22 .20 .18 .15 .12 .1 .1];
%egal  = [.22 .20 .18  .16 .16 .13 .12 .1 .1];
%18.8000   17.6250   16.4500   15.2750   14.5700   13.5830   12.6900   11.7030   10.7160
%egal  = linspace(18.5,11.5,9)*.01;
egal  = linspace(14+5,14-5,9)*.01;
%addlow
M = [   3/2   3/2   1   1/2 0  ;...
        2     2     4/3 2/3 0  ;...
        3     3     2   1   1/4;...
        4     4     8/3 4/3 1/2;...
        5     5     7/2 2   3/4;...
        6     6     4   2   1   ] / 100;
egal(1:5)=egal(1:5) + M(1,:);

%350 36
%320 38
%320 38
%300 42
%270 44
%240 46
%

%%addhigh
%egal(1:2)=egal(1:2)+.015;egal(3)=egal(3)+.01;egal(4)=egal(4)+.005;
%egal  = linspace(19,11,9) * .01;egal(1:4) = egal(1:4) - 0.005;
%egal  = [0.1980    0.1863    0.1835    0.1707    0.1507    0.1358    0.1269    0.1170    0.1072];
%egal  = [19.8000   18.9250   18.5500   17.0750   15.0700   13.5830  12.6900   11.7030   10.7160] * .01;%up3
%egal  = [19.0000   18.0000   16.9500   16.2750   15.5700   13.5830   12.6900   11.7030   10.7160] * .01;
%linspace(.20, .10, 9);
%egal([5,6,7,8,9]) = egal([5,6,7,8,9])+ 1*[.0025, .0035, .0050, .0060, .0070];% egal = egal * .94;
%egal = [.19 .18 .17  .16 .15 .14 .13 .12 .12];

%egal = [.2 .2 .2  .18 .16 .14 .12 .1 .1];
%egal = [.24 .24 .24 .21 .18 .16 .14 .12 .12];
%egal =  [.30 .26 .23  .20 .17 .13 .11 .09 .09];

a.fbEgalMod = linspace(.225,.075,9) / .15;

a.plainsEE=egal; %[.11 .2];%[.01 .29]; %
a.plainsEI=.6*ones(1,9);
%a.plainsIE=.6*ones(1,9);
a.plainsIE= [.6   .6   .6   .6   .6   .53  .46  .40  .40];
a.bndryIE = [1.67 1.67 1.67 1.67 1.67 1.48 1.28 1.11 1.11];
a.plainsII=.6*ones(1,9);
a.bndryEE=2*egal; %[.22 .4];% [.02 .56]; %
a.bndryEI=1.67*ones(1,9);
a.bndryIE=1.67*ones(1,9);
a.bndryII=1.67*ones(1,9);

%turn off connectivity
%a.plainsEE = [0 0 0 0 0 0 0 0 0];
%a.plainsIE = [0 0 0 0 0 0 0 0 0];
%a.plainsEI = [0 0 0 0 0 0 0 0 0];
%a.plainsII = [0 0 0 0 0 0 0 0 0];
%a.bndryEE  = [0 0 0 0 0 0 0 0 0];
%a.bndryIE  = [0 0 0 0 0 0 0 0 0];
%a.bndryEI  = [0 0 0 0 0 0 0 0 0];
%a.bndryII  = [0 0 0 0 0 0 0 0 0];

disp('allocating connectivity matrix...')
a.allocateConnectivity(connectivityEE,connectivityIE,connectivityEI,connectivityII,egalEHelper,egalIHelper,LGNInputFile,fbInputFile);

disp('building network...')
a.buildNetwork();