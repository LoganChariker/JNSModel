%%

%====================%
%temporary code below%
%====================%


duration = 30000;
    

fbEgalMod = { [(egal / .15 - 1)*1.5+1] ,...
              [.3 .25 .2 .16 .12 .08 .04 .04 .04] / .15 ,...
              [linspace(.20,.10,9) / .15] ,...
              [linspace(.21,.09,9) / .15] ,...
              [linspace(.24,.06,9) / .15] ,...
              [linspace(.25,.07,9) / .15] ,...
              [linspace(.225,.075,9) / .15],...
              [linspace(.2025,.0975,9) / .15],...
              [linspace(.15+.045*1.5,.15-.045*1.5,9) / .15],...
              };
          
fbISD = { [80] };

fbIPeak = { [.60], [.55] };
          
fbExtraProb = { [0] };

pLow = { 0 };
pHigh = { 2 };

gHighACh = [ .350, 1.05, 2.5, 1.0 ];
gLowACh  = [ .300, 1.00, 1.8, 0.6 ];
g = gHighACh;
gEEFac = 1;%(1+(g(2)-1)*1/2);

tauEFac = { [3.5]/3 };
tauIFac = { [5.625] /5 };

ambientStE = { .200 };
ambientStI = { .420 };
ambientWkE = { .225 };
ambientWkI = { g(1) };
sEE        = { .028 };
justSEE    = { .028 * gEEFac };
sEI        = { 2.0 };
sIE        = { .36 };
eps        = { 0 };
threshI    = { [.9 1.1] };
threshE    = { [.95 1.05] };
sIIMod     = { [.65 .85] };
LGNFileInd = { 1 };
conn       = { [.15 .6 .6 .6] };
egalInd    = { 1 };
bgRate     = { .00 };
wavelengths= { 800*2.5./[1.25;  2;  2.5;  3;  4] };%arrayfun(@(p) 500 * 1.10^p,0 : 9)';
angles     = { [1;2;3;4;5;6;7;8] };
sLGN       = { 2.1 / gEEFac * g(2) };
sLGNInh    = { 3.0 };
shotNoise  = { 0.075 };
vthresh    = { 1.2 };
lgndelayst = { 20 };
egal      = egal; %linspace(15+4.5,15-4.5,9)*.01;% egal([5 6 7 8 9]) = egal([5 6 7 8 9]) - [1/4 1/2 3/4 1 1]*.01; %egal([1 2 3 4]) = egal([1 2 3 4]) + [.001 .001 .0025 .0025]; %
%egal1      = [0.1980    0.1863    0.1835    0.1707    0.1507    0.1358    0.1269    0.1170    0.1072];
%egal1      = [18.8000   17.6250   16.9500   16.2750   15.5700   13.5830   12.6900   11.7030   10.7160] * .01;
%linspace(.20, .10, 9); egal1([5,6,7,8,9]) = egal1([5,6,7,8,9])+ 1*[.0025, .0035, .0050, .0060, .0070];% egal1 = egal1 * .94;
egbd1      = 2*egal;
fb         = { g(3:4) };
plainsEE  = { egal };
bndryEE   = { egbd1 };

%==matrix of parameters we will choose from==%
%==parameters' columns: drive, see, sei/see, sie/see, eps, lowerIthresh higherIthresh ==%
parameters = [...    
    rowCombineMany({ambientStE{1} ambientWkE{1} sEE{1} sEI{1} sIE{1}  eps{1} threshI{1} sIIMod{1} LGNFileInd{1} conn{1} egalInd{1} bgRate{1}  angles{1} wavelengths{1} sLGN{1} sLGNInh{1} shotNoise{1} vthresh{1} lgndelayst{1} ambientWkI{1} plainsEE{1} bndryEE{1} justSEE{1} threshE{1} ambientStI{1} fb{1}(1,:) fbEgalMod{1} tauEFac{1} tauIFac{1} fbISD{1} fbExtraProb{1} fbIPeak{1} pLow{1} pHigh{1} });...
    ];
    %amst amwk see     sei sie eps  threshI siimod  LGNfile     connectivity
    %Egal RB

paramSelector = 17;

%==choose parameters==%
    SF = 1/(parameters(paramSelector, 19) / 800) * 2.5;
    fac = exp(-0.5 * ((SF-2.5)/1.5).^2);
    fac2 = fac;
    if fac < .26
        fac2 = .3;
    end
    if fac < .01
        fac2 = .1;
    end
    a.pValues = [1.3, 1.0, 1.0, 0.7, 0.7, 0.4, 0.4, 0.7, 0.7, 1.0, 1.0, 1.3, 1.3] * fac2;
    a.pLow = parameters(paramSelector, 64);
    a.pHigh = parameters(paramSelector, 65);
    a.fbEgalMod = parameters(paramSelector, 50:58);
    a.fbStart = 8000;
    a.fbEnd = 8800;
    a.fbClock = 8000;
    a.currentFBSpike = 0;
    a.fbIPeak = parameters(paramSelector, 63);
    a.fbISD = parameters(paramSelector,61);
    %a.fbExtraProb = parameters(paramSelector,62);
    a.fbPref = parameters(paramSelector, 48);
    a.fbOrth = parameters(paramSelector, 49);
    a.setAmbient(parameters(paramSelector, 1), parameters(paramSelector, 2), parameters(paramSelector,47), parameters(paramSelector, 25));
    a.setSEE(parameters(paramSelector, 3));
    a.setJustSEE(parameters(paramSelector, 44));
    a.setSEIOverSEE(parameters(paramSelector, 4));
    a.setSIEOverSEE(parameters(paramSelector, 5));
    a.setContrast(parameters(paramSelector, 6));
    a.setThreshIRange(parameters(paramSelector, 7), parameters(paramSelector, 8));
    a.setThreshERange(parameters(paramSelector, 45), parameters(paramSelector, 46));
    a.setModIIRange(parameters(paramSelector, 9), parameters(paramSelector, 10));
    %a.plainsEE = parameters(paramSelector, 26:34); %next few lines are hacks
    %a.bndryEE = parameters(paramSelector, 35:43);   
    %a.LGNInputFile = '';
    %a.setLGNInputsIfNew(LGNInputFiles(parameters(paramSelector, 11)));
    %a.loadFB();
    
    a.tauAMPA = 1.0 * parameters(paramSelector, 59);
    a.tauGABA = 1.67  * parameters(paramSelector, 60);
    
    %disp('setting egalitarian...');a.setEgalIfNew(parameters(paramSelector, 12),...
    %    parameters(paramSelector, 13),...
    %    parameters(paramSelector, 14),...
    %    parameters(paramSelector, 15),...
    %    egalType{parameters(paramSelector,16)}, [0.0]);
    
    %disp('setting background rate...'); a.RB = parameters(paramSelector, 17);
    %disp('setting up feedback...');
    %a.expEPreExc(1) = mean(a.numEPreExc);
    %a.expEPreInh(1) = mean(a.numEPreInh);
    %a.spontaneousRate = parameters(paramSelector, 18);
    %a.feedbackFractionExc = parameters(paramSelector, 19);
    %a.feedbackFractionInh = parameters(paramSelector, 20);
    a.setKOrient( sim.MyPoint(cos((parameters(paramSelector, 18)-1)*pi/8), sin((parameters(paramSelector, 18)-1)*pi/8)) );
    a.loadFB();
    a.setKFreq( 2*pi / parameters(paramSelector, 19) );    
    
        %set fast I
        tau = 1.67;
        dt = 1.5/parameters(paramSelector,60);
        accum = (1/(1-0.0335616))*exp(-dt/tau) / (6.0 * tau);
        g3 = accum;
        accum = accum * dt/tau;
        g2 = accum;
        accum = accum * dt/tau;
        g1 = accum;
        accum = accum * dt/tau;
        g0 = accum;
        a.fastSpikeInitial = [g3 g2 g1 g0];
  
    
    %===============%
    %========THE FOLLOWING IS A TEMPORARY TEST====%
    a.sLGNOverSEEFactor = parameters(paramSelector,20);
    a.sLGNInh = parameters(paramSelector,21) * a.sEE;
    a.shotNoise= parameters(paramSelector,22); 
    a.threshLGN= parameters(paramSelector,23); 
    a.lgnSpikeDelaySteps = parameters(paramSelector,24);
    a.w = 4.0 / 1000.0;
    %=============================================%
    a.resetNetState();
    disp('net state and parameters reset')
    disp('parameters:')
    disp(parameters(paramSelector,:));
    a.setPhases();   
    %a.saveConnectivity(  '/home/charikar/mfiles/Sim/connections.dat');
    %a.saveFBConnectivity('/home/charikar/mfiles/Sim/fbconnections.dat');
    a.setFBPValues();
    %==chose parameters==%

%% graph past-50ms-fr for 2 sec
durationInMs = 2000;
skipInMs = 50;

skipInSteps = skipInMs / a.dt;
a.statSkipInSteps = skipInSteps;
durationInSteps = durationInMs / a.dt;

%collect data from neurons within vertical-preferring patch (receiving strong drive) and
%horizontal-preferring patch (receiving weak drive)
drivePatchCenter=[500,750];
noDrivePatchCenter=[1000,750];
patchRadius=100;
locExc=zeros(a.numExcitatory,2);
locInh=zeros(a.numExcitatory,2);
drivePatchExc=false(a.numExcitatory,1);
noDrivePatchExc=false(a.numExcitatory,1);
drivePatchInh=false(a.numInhibitory,1);
noDrivePatchInh=false(a.numInhibitory,1);
for i=1:a.numExcitatory
    loc=[a.locExc(i).x, a.locExc(i).y];
    if norm(loc-drivePatchCenter) < patchRadius
        drivePatchExc(i)=true;
    end
    if norm(loc-noDrivePatchCenter) < patchRadius
        noDrivePatchExc(i)=true;
    end
end
for i=1:a.numInhibitory
    loc=[a.locInh(i).x, a.locInh(i).y];
    if norm(loc-drivePatchCenter) < patchRadius
        drivePatchInh(i)=true;
    end
    if norm(loc-noDrivePatchCenter) < patchRadius
        noDrivePatchInh(i)=true;
    end
end

numInDrivePatchExc=sum(drivePatchExc);
numInDrivePatchInh=sum(drivePatchInh);
numInNoDrivePatchExc=sum(noDrivePatchExc);
numInNoDrivePatchInh=sum(noDrivePatchInh);

a.membershipExc(drivePatchExc)=1;
a.membershipExc(noDrivePatchExc)=2;
a.membershipInh(drivePatchInh)=1;
a.membershipInh(noDrivePatchInh)=2;


maxylim=20;

dataCollectTimes = 0 : skipInMs : durationInMs;
%frdata = zeros(size(dataCollectTimes));
fr1data= zeros(size(dataCollectTimes));
fr2data= zeros(size(dataCollectTimes));
fr1Idata= zeros(size(dataCollectTimes));
fr2Idata= zeros(size(dataCollectTimes));


f = figure('position', [400, 100, 1600, 900]);
ax = axes('position',[.4,.1,.5,.8]);
hold on; grid on;
%pl   = plot(ax,dataCollectTimes,frdata,'r-');
plp1 = plot(ax,dataCollectTimes,fr1data,'r--');
plp2 = plot(ax,dataCollectTimes,fr2data,'r:');
plp1I = plot(ax,dataCollectTimes,fr1Idata,'b--');
plp2I = plot(ax,dataCollectTimes,fr2Idata,'b:');
lin = line([0 0],[0 maxylim],'color','g');
xlim([0 durationInMs]);
ylim([0 maxylim]);
xlabel('time(ms)')
ylabel('past-50ms-firing-rate')

textAxes = axes('position',[.05, .1, .15, .9],...
                'visible','off');    
            
uiPanel = uipanel('position',get(textAxes,'position') + [.1 0 0 0]);
        
sietextbase = 'SIE/SEE = ';
drivetextbase = 'Contrast = ';
seitextbase = 'SEI/SEE = ';
ambetextbase = 'Amb->E = ';
ambitextbase = 'Amb->I = ';
amballtextbase = 'WkAmb->E =';
bgratetextbase = 'bgrate = ';
fbetextbase = 'fb->E = ';
fbitextbase = 'fb->I = ';
inctextbase = 'increment = ';
angletextbase = 'angle = ';
freqtextbase = 'wavelength = ';
phasetextbase = 'sp phase (rad/2\pi)= ';
ambwkItextbase= 'WkAmb->I ';
LGNEtextbase = 'LGN->E ';
LGNItextbase = 'LGN->I ';

changeamt = .01;
angle = 0;
wavelen = wavelengths{1};
phase = 0;

textPosition = .9:-.05 : 0.05;
sietext =    text(.05,textPosition(1),[sietextbase num2str(a.sIE/a.sEE)]);
drivetext =  text(.05,textPosition(2),[drivetextbase num2str(a.eps)]);
seitext =    text(.05,textPosition(3),[seitextbase num2str(a.sEI/a.sEE)]);
ambetext =   text(.05,textPosition(4),[ambetextbase num2str(a.rateStrongAmbE)]);
amballtext = text(.05,textPosition(5),[amballtextbase num2str(a.rateWeakAmbE)]);
%bgratetext = text(.05,textPosition(6),[bgratetextbase num2str(a.spontaneousRate)]);
%fbetext =    text(.05,textPosition(7),[fbetextbase num2str(a.feedbackFractionExc)]);
%fbitext =    text(.05,textPosition(8),[fbitextbase num2str(a.feedbackFractionInh)]);
inctext =    text(.05,textPosition(6),[inctextbase num2str(double(changeamt))]);
angletext =  text(.05,textPosition(7),[angletextbase num2str(angle)]);
freqtext =   text(.05,textPosition(8),[freqtextbase num2str(2*pi/a.kFreq)]);
phasetext =  text(.05,textPosition(9),[phasetextbase num2str(phase)]);
ambwkItext =  text(.05,textPosition(10),[ambwkItextbase num2str(a.rateWeakAmbI)]);
LGNEtext =   text(.05,textPosition(11),[LGNEtextbase num2str(a.sLGNOverSEEFactor)]);
LGNItext =   text(.05,textPosition(12),[LGNItextbase num2str(a.sLGNInh/a.sEE)]);
ambitext =   text(.05,textPosition(13),[ambitextbase num2str(a.rateStrongAmbI)]);
ns = 0;

leftCallback = {...
    'a.setSIEOverSEE(a.sIE/a.sEE - changeamt);settext(sietext,sietextbase,a.sIE/a.sEE);';...
    'a.eps = a.eps - .1;a.setContrast(a.eps);settext(drivetext,drivetextbase,a.eps);';...
    'a.setSEIOverSEE(a.sEI/a.sEE - changeamt);settext(seitext,seitextbase,a.sEI/a.sEE);';...
    'a.setAmbient(a.rateStrongAmbE - changeamt,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI);settext(ambetext,ambetextbase,a.rateStrongAmbE);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE - changeamt,a.rateStrongAmbI,a.rateWeakAmbI);settext(amballtext,amballtextbase,a.rateWeakAmbE);';...
%    'a.spontaneousRate = a.spontaneousRate - changeamt;settext(bgratetext,bgratetextbase,a.spontaneousRate);';...
%    'a.feedbackFractionExc = a.feedbackFractionExc - changeamt;settext(fbetext,fbetextbase,a.feedbackFractionExc);';...
%    'a.feedbackFractionInh = a.feedbackFractionInh - changeamt;settext(fbitext,fbitextbase,a.feedbackFractionInh);';...
    'changeamt = changeamt / 10;settext(inctext,inctextbase,changeamt);';...
    'angle = angle - changeamt; a.kOrient.x = cos(angle); a.kOrient.y = sin(angle);settext(angletext,angletextbase,angle);';...
    'wavelen = wavelen - changeamt; a.setKFreq(2*pi/wavelen); settext(freqtext,freqtextbase,2*pi/a.kFreq);';...
    'phase = phase - changeamt; a.kPhase = 2*pi*phase; settext(phasetext,phasetextbase,a.kPhase / (2*pi));';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI - changeamt);settext(ambwkItext,ambwkItextbase,a.rateWeakAmbI);';...
    'a.sLGNOverSEEFactor = a.sLGNOverSEEFactor - .25; settext(LGNEtext,LGNEtextbase,a.sLGNOverSEEFactor);';...
    'a.sLGNInh = (a.sLGNInh/a.sEE - .25)*a.sEE; settext(LGNItext,LGNItextbase,a.sLGNInh/a.sEE);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE,a.rateStrongAmbI - changeamt,a.rateWeakAmbI);settext(ambitext,ambitextbase,a.rateStrongAmbI);';...
};
rightCallback = {...
    'a.setSIEOverSEE(a.sIE/a.sEE + changeamt);settext(sietext,sietextbase,a.sIE/a.sEE);';...
    'a.eps = a.eps + .1;a.setContrast(a.eps);settext(drivetext,drivetextbase,a.eps);';...
    'a.setSEIOverSEE(a.sEI/a.sEE + changeamt);settext(seitext,seitextbase,a.sEI/a.sEE);';...
    'a.setAmbient(a.rateStrongAmbE + changeamt,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI);settext(ambetext,ambetextbase,a.rateStrongAmbE);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE + changeamt,a.rateStrongAmbI,a.rateWeakAmbI);settext(amballtext,amballtextbase,a.rateWeakAmbE);';...
%    'a.spontaneousRate = a.spontaneousRate + changeamt;settext(bgratetext,bgratetextbase,a.spontaneousRate);';...
%    'a.feedbackFractionExc = a.feedbackFractionExc + changeamt;settext(fbetext,fbetextbase,a.feedbackFractionExc);';...
%    'a.feedbackFractionInh = a.feedbackFractionInh + changeamt;settext(fbitext,fbitextbase,a.feedbackFractionInh);';...
    'changeamt = changeamt * 10;settext(inctext,inctextbase,changeamt);';...
    'angle = angle + changeamt; a.kOrient.x = cos(angle); a.kOrient.y = sin(angle);settext(angletext,angletextbase,angle);';...
    'wavelen = wavelen + changeamt; a.setKFreq(2*pi/wavelen); settext(freqtext,freqtextbase,2*pi/a.kFreq);';...
    'phase = phase + changeamt; a.kPhase = phase; settext(phasetext,phasetextbase,a.kPhase);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI + changeamt);settext(ambwkItext,ambwkItextbase,a.rateWeakAmbI);';...
    'a.sLGNOverSEEFactor = a.sLGNOverSEEFactor + .25; settext(LGNEtext,LGNEtextbase,a.sLGNOverSEEFactor);';...
    'a.sLGNInh = (a.sLGNInh/a.sEE + .25)*a.sEE; settext(LGNItext,LGNItextbase,a.sLGNInh/a.sEE);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE,a.rateStrongAmbI + changeamt,a.rateWeakAmbI);settext(ambitext,ambitextbase,a.rateStrongAmbI);';...
};

for i = 1 : length(leftCallback)
    uicontrol(uiPanel,'Style', 'pushbutton', 'String', '<',...
    'units','normalized',...
    'position', [.1 textPosition(i)-.02 .05 .05],...
    'callback', leftCallback{i});
uicontrol(uiPanel,'Style', 'pushbutton', 'String', '>',...
    'units','normalized',...
    'position', [.15 textPosition(i)-.02 .05 .05],...
    'callback', rightCallback{i});
end

prefvis = 1;
orthvis = 1;
fullvis = 1;
visstring = 'visible';
onstring = 'on';
offstring = 'off';
uicontrol(f,'Style', 'pushbutton', 'String', '',...
    'position', [50 50 20 20],...
    'callback', 'if fullvis==1; set(pl,visstring,offstring);fullvis=0;else;set(pl,visstring,onstring);fullvis=1;end;')

keepgoing = 1;
while keepgoing
    
%frdata(1) = a.timeLocalFR(a.currentCounter+1)/50/27;
fr1data(1) = a.patchOneFRExc/skipInMs/numInDrivePatchExc * 1000;
fr2data(1) = a.patchTwoFRExc/skipInMs/numInNoDrivePatchExc * 1000;
fr1Idata(1) = a.patchOneFRInh/skipInMs/numInDrivePatchInh * 1000;
fr2Idata(1) = a.patchTwoFRInh/skipInMs/numInNoDrivePatchInh * 1000;
set(lin,'xdata',dataCollectTimes(1)*[0 0])
for i = 2 : length(dataCollectTimes)
    disp(['i=' num2str(i)])
    for j = 1 : skipInSteps
        a.step();
    end
    title(ax,i)
    
    %frdata(i) = a.timeLocalFR(a.currentCounter+1)/50/27;
    fr1data(i) = a.patchOneFRExc/skipInMs/numInDrivePatchExc * 1000;
    fr2data(i) = a.patchTwoFRExc/skipInMs/numInNoDrivePatchExc * 1000;
    fr1Idata(i) = a.patchOneFRInh/skipInMs/numInDrivePatchInh * 1000;
    fr2Idata(i) = a.patchTwoFRInh/skipInMs/numInNoDrivePatchInh * 1000;
    
    %set(pl,'YData',frdata);
    set(plp1, 'YData',fr1data);
    set(plp2, 'YData',fr2data);
    set(plp1I, 'YData',fr1Idata);
    set(plp2I, 'YData',fr2Idata);
    set(lin, 'XData', dataCollectTimes(i)*[1 1])
    
    title({[num2str(a.t) ' ms']...
        ['S^{IE}/S^{EE}=' num2str(a.sIE/a.sEE)]...
        ['contrast=' num2str(a.eps)]});% num2str(a.currentCounter) num2str(a.numStepsMod10) num2str(a.eSpikes / 12 / a.t,2) num2str(ns / 12 / a.t,2)});
   
    drawnow();
end
end

%% view individual neurons in real time
durationInMs = 1000;
skipInMs = 2;

skipInSteps = skipInMs / a.dt;
a.statSkipInSteps = skipInSteps;
durationInSteps = durationInMs / a.dt;

maxylim=1;
minylim=-.1;

dataCollectTimes = 0 : skipInMs : durationInMs;
%frdata = zeros(size(dataCollectTimes));
vdata= zeros(size(dataCollectTimes));

f = figure('position', [400, 100, 1600, 900]);
ax = axes('position',[.4,.1,.5,.8]);
hold on; grid on;
%pl   = plot(ax,dataCollectTimes,frdata,'r-');
plp1 = plot(ax,dataCollectTimes,vdata,'r');
lin = line([0 0],[0 maxylim],'color','g');
xlim([0 durationInMs]);
ylim([minylim maxylim]);
xlabel('time(ms)')
ylabel('v')

textAxes = axes('position',[.05, .1, .15, .9],...
                'visible','off');    
            
uiPanel = uipanel('position',get(textAxes,'position') + [.1 0 0 0]);
        
sietextbase = 'SIE/SEE = ';
drivetextbase = 'Contrast = ';
seitextbase = 'SEI/SEE = ';
ambetextbase = 'Amb->E = ';
amballtextbase = 'Amb->All =';
bgratetextbase = 'bgrate = ';
fbetextbase = 'fb->E = ';
fbitextbase = 'fb->I = ';
inctextbase = 'increment = ';
angletextbase = 'angle = ';
freqtextbase = 'wavelength = ';
phasetextbase = 'sp phase (rad/2\pi)= ';
indextextbase = 'neuron index = ';

changeamt = .01;
angle = 0;
wavelen = wavelengths{1};
phase = 0;
whichNeuron = 100; title(ax,a.numInputsExc(inDrivePatchExc(whichNeuron)));

textPosition = .9:-.05 : 0.05;
sietext =    text(.05,textPosition(1),[sietextbase num2str(a.sIE/a.sEE)]);
drivetext =  text(.05,textPosition(2),[drivetextbase num2str(a.eps)]);
seitext =    text(.05,textPosition(3),[seitextbase num2str(a.sEI/a.sEE)]);
ambetext =   text(.05,textPosition(4),[ambetextbase num2str(a.rateStrongAmbE)]);
amballtext = text(.05,textPosition(5),[amballtextbase num2str(a.rateWeakAmbE)]);
%bgratetext = text(.05,textPosition(6),[bgratetextbase num2str(a.spontaneousRate)]);
%fbetext =    text(.05,textPosition(7),[fbetextbase num2str(a.feedbackFractionExc)]);
%fbitext =    text(.05,textPosition(8),[fbitextbase num2str(a.feedbackFractionInh)]);
inctext =    text(.05,textPosition(6),[inctextbase num2str(double(changeamt))]);
angletext =  text(.05,textPosition(7),[angletextbase num2str(angle)]);
freqtext =   text(.05,textPosition(8),[freqtextbase num2str(wavelen)]);
phasetext =  text(.05,textPosition(9),[phasetextbase num2str(phase)]);
indextext =  text(.05,textPosition(10),[indextextbase num2str(whichNeuron)]);

ns = 0;

leftCallback = {...
    'a.setSIEOverSEE(a.sIE/a.sEE - changeamt);settext(sietext,sietextbase,a.sIE/a.sEE);';...
    'a.eps = a.eps - 1;a.setContrast(a.eps);settext(drivetext,drivetextbase,a.eps);';...
    'a.setSEIOverSEE(a.sEI/a.sEE - changeamt);settext(seitext,seitextbase,a.sEI/a.sEE);';...
    'a.setAmbient(a.rateStrongAmbE - changeamt,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI);settext(ambetext,ambetextbase,a.rateStrongAmbE);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE - changeamt,a.rateStrongAmbI,a.rateWeakAmbI);settext(amballtext,amballtextbase,a.rateWeakAmbE);';...
%    'a.spontaneousRate = a.spontaneousRate - changeamt;settext(bgratetext,bgratetextbase,a.spontaneousRate);';...
%    'a.feedbackFractionExc = a.feedbackFractionExc - changeamt;settext(fbetext,fbetextbase,a.feedbackFractionExc);';...
%    'a.feedbackFractionInh = a.feedbackFractionInh - changeamt;settext(fbitext,fbitextbase,a.feedbackFractionInh);';...
    'changeamt = changeamt / 10;settext(inctext,inctextbase,changeamt);';...
    'angle = angle - changeamt; a.kOrient.x = cos(angle); a.kOrient.y = sin(angle);settext(angletext,angletextbase,angle);';...
    'wavelen = wavelen - changeamt; a.setKFreq(2*pi/wavelen); settext(freqtext,freqtextbase,2*pi/a.kFreq);';...
    'phase = phase - changeamt; a.kPhase = 2*pi*phase; settext(phasetext,phasetextbase,a.kPhase / (2*pi));';...
    'whichNeuron = whichNeuron - 1;title(ax,a.numInputsExc(inDrivePatchExc(whichNeuron)));settext(indextext,indextextbase,whichNeuron);';...
};
rightCallback = {...
    'a.setSIEOverSEE(a.sIE/a.sEE + changeamt);settext(sietext,sietextbase,a.sIE/a.sEE);';...
    'a.eps = a.eps + 1;a.setContrast(a.eps);settext(drivetext,drivetextbase,a.eps);';...
    'a.setSEIOverSEE(a.sEI/a.sEE + changeamt);settext(seitext,seitextbase,a.sEI/a.sEE);';...
    'a.setAmbient(a.rateStrongAmbE + changeamt,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI);settext(ambetext,ambetextbase,a.rateStrongAmbE);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE + changeamt,a.rateStrongAmbI,a.rateWeakAmbI);settext(amballtext,amballtextbase,a.rateWeakAmbE);';...
%    'a.spontaneousRate = a.spontaneousRate + changeamt;settext(bgratetext,bgratetextbase,a.spontaneousRate);';...
%    'a.feedbackFractionExc = a.feedbackFractionExc + changeamt;settext(fbetext,fbetextbase,a.feedbackFractionExc);';...
%    'a.feedbackFractionInh = a.feedbackFractionInh + changeamt;settext(fbitext,fbitextbase,a.feedbackFractionInh);';...
    'changeamt = changeamt * 10;settext(inctext,inctextbase,changeamt);';...
    'angle = angle + changeamt; a.kOrient.x = cos(angle); a.kOrient.y = sin(angle);settext(angletext,angletextbase,angle);';...
    'wavelen = wavelen + changeamt; a.setKFreq(2*pi/wavelen); settext(freqtext,freqtextbase,2*pi/a.kFreq);';...
    'phase = phase + changeamt; a.kPhase = phase; settext(phasetext,phasetextbase,a.kPhase);';...
    'whichNeuron = whichNeuron + 1;title(ax,a.numInputsExc(inDrivePatchExc(whichNeuron)));settext(indextext,indextextbase,whichNeuron);';...
};

for i = 1 : length(leftCallback)
    uicontrol(uiPanel,'Style', 'pushbutton', 'String', '<',...
    'units','normalized',...
    'position', [.1 textPosition(i)-.02 .05 .05],...
    'callback', leftCallback{i});
uicontrol(uiPanel,'Style', 'pushbutton', 'String', '>',...
    'units','normalized',...
    'position', [.15 textPosition(i)-.02 .05 .05],...
    'callback', rightCallback{i});
end

prefvis = 1;
orthvis = 1;
fullvis = 1;
visstring = 'visible';
onstring = 'on';
offstring = 'off';
uicontrol(f,'Style', 'pushbutton', 'String', '',...
    'position', [50 50 20 20],...
    'callback', 'if fullvis==1; set(pl,visstring,offstring);fullvis=0;else;set(pl,visstring,onstring);fullvis=1;end;')

keepgoing = 1;
while keepgoing
    
%frdata(1) = a.timeLocalFR(a.currentCounter+1)/50/27;
vdata(1) = a.vExc(inNoDrivePatchExc(whichNeuron));
set(lin,'xdata',dataCollectTimes(1)*[0 0])
for i = 2 : length(dataCollectTimes)
    disp(['i=' num2str(i)])
    for j = 1 : skipInSteps
        a.step();
    end
    title(ax,a.numInputsExc(inNoDrivePatchExc(whichNeuron)));
    
    %frdata(i) = a.timeLocalFR(a.currentCounter+1)/50/27;
    vdata(i) = a.vExc(inNoDrivePatchExc(whichNeuron));
    
    %set(pl,'YData',frdata);
    set(plp1, 'YData',vdata);
    set(lin, 'XData', dataCollectTimes(i)*[1 1])
            
    drawnow();
end
end

%% view individual neuron's conductances in real time
durationInMs = 1000;
skipInMs = 2;

skipInSteps = skipInMs / a.dt;
a.statSkipInSteps = skipInSteps;
durationInSteps = durationInMs / a.dt;

maxylim=1;
minylim=-.1;

dataCollectTimes = 0 : skipInMs : durationInMs;
%frdata = zeros(size(dataCollectTimes));
gLGNdata= zeros(size(dataCollectTimes));
gAMPAdata=zeros(size(dataCollectTimes));
gGABAdata=zeros(size(dataCollectTimes));
%vdata = zeros(size(dataCollectTimes));

f = figure('position', [400, 100, 1600, 900]);
ax = axes('position',[.4,.1,.5,.8]);
hold on; grid on;
%pl   = plot(ax,dataCollectTimes,frdata,'r-');
%plp4 = plot(ax,dataCollectTimes,vdata,'k');
plp1 = plot(ax,dataCollectTimes,gLGNdata,'g');
plp2 = plot(ax,dataCollectTimes,gAMPAdata,'r');
plp3 = plot(ax,dataCollectTimes,gGABAdata,'b');

lin = line([0 0],[0 maxylim],'color','g');
xlim([0 durationInMs]);
ylim([minylim maxylim]);
xlabel('time(ms)')
ylabel('v')

textAxes = axes('position',[.05, .1, .15, .9],...
                'visible','off');    
            
uiPanel = uipanel('position',get(textAxes,'position') + [.1 0 0 0]);
        
sietextbase = 'SIE/SEE = ';
drivetextbase = 'Contrast = ';
seitextbase = 'SEI/SEE = ';
ambetextbase = 'Amb->E = ';
amballtextbase = 'Amb->All =';
bgratetextbase = 'bgrate = ';
fbetextbase = 'fb->E = ';
fbitextbase = 'fb->I = ';
inctextbase = 'increment = ';
angletextbase = 'angle = ';
freqtextbase = 'wavelength = ';
phasetextbase = 'sp phase (rad/2\pi)= ';
indextextbase = 'neuron index = ';

changeamt = .01;
angle = 0;
wavelen = wavelengths{1}(1);
phase = 0;
whichNeuron = 100; title(ax,a.numInputsExc(inDrivePatchExc(whichNeuron)));

textPosition = .9:-.05 : 0.05;
sietext =    text(.05,textPosition(1),[sietextbase num2str(a.sIE/a.sEE)]);
drivetext =  text(.05,textPosition(2),[drivetextbase num2str(a.eps)]);
seitext =    text(.05,textPosition(3),[seitextbase num2str(a.sEI/a.sEE)]);
ambetext =   text(.05,textPosition(4),[ambetextbase num2str(a.rateStrongAmbE)]);
amballtext = text(.05,textPosition(5),[amballtextbase num2str(a.rateWeakAmbE)]);
%bgratetext = text(.05,textPosition(6),[bgratetextbase num2str(a.spontaneousRate)]);
%fbetext =    text(.05,textPosition(7),[fbetextbase num2str(a.feedbackFractionExc)]);
%fbitext =    text(.05,textPosition(8),[fbitextbase num2str(a.feedbackFractionInh)]);
inctext =    text(.05,textPosition(6),[inctextbase num2str(double(changeamt))]);
angletext =  text(.05,textPosition(7),[angletextbase num2str(angle)]);
freqtext =   text(.05,textPosition(8),[freqtextbase num2str(wavelen)]);
phasetext =  text(.05,textPosition(9),[phasetextbase num2str(phase)]);
indextext =  text(.05,textPosition(10),[indextextbase num2str(whichNeuron)]);

ns = 0;

leftCallback = {...
    'a.setSIEOverSEE(a.sIE/a.sEE - changeamt);settext(sietext,sietextbase,a.sIE/a.sEE);';...
    'a.eps = a.eps - .1;a.setContrast(a.eps);settext(drivetext,drivetextbase,a.eps);';...
    'a.setSEIOverSEE(a.sEI/a.sEE - changeamt);settext(seitext,seitextbase,a.sEI/a.sEE);';...
    'a.setAmbient(a.rateStrongAmbE - changeamt,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI);settext(ambetext,ambetextbase,a.rateStrongAmbE);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE - changeamt,a.rateStrongAmbI,a.rateWeakAmbI);settext(amballtext,amballtextbase,a.rateWeakAmbE);';...
%    'a.spontaneousRate = a.spontaneousRate - changeamt;settext(bgratetext,bgratetextbase,a.spontaneousRate);';...
%    'a.feedbackFractionExc = a.feedbackFractionExc - changeamt;settext(fbetext,fbetextbase,a.feedbackFractionExc);';...
%    'a.feedbackFractionInh = a.feedbackFractionInh - changeamt;settext(fbitext,fbitextbase,a.feedbackFractionInh);';...
    'changeamt = changeamt / 10;settext(inctext,inctextbase,changeamt);';...
    'angle = angle - changeamt; a.kOrient.x = cos(angle); a.kOrient.y = sin(angle);settext(angletext,angletextbase,angle);';...
    'wavelen = wavelen - changeamt; a.setKFreq(2*pi/wavelen); settext(freqtext,freqtextbase,2*pi/a.kFreq);';...
    'phase = phase - changeamt; a.kPhase = 2*pi*phase; settext(phasetext,phasetextbase,a.kPhase / (2*pi));';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI - changeamt);settext(ambwkItext,ambwkItextbase,a.rateWeakAmbI);';...
};
rightCallback = {...
    'a.setSIEOverSEE(a.sIE/a.sEE + changeamt);settext(sietext,sietextbase,a.sIE/a.sEE);';...
    'a.eps = a.eps + .1;a.setContrast(a.eps);settext(drivetext,drivetextbase,a.eps);';...
    'a.setSEIOverSEE(a.sEI/a.sEE + changeamt);settext(seitext,seitextbase,a.sEI/a.sEE);';...
    'a.setAmbient(a.rateStrongAmbE + changeamt,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI);settext(ambetext,ambetextbase,a.rateStrongAmbE);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE + changeamt,a.rateStrongAmbI,a.rateWeakAmbI);settext(amballtext,amballtextbase,a.rateWeakAmbE);';...
%    'a.spontaneousRate = a.spontaneousRate + changeamt;settext(bgratetext,bgratetextbase,a.spontaneousRate);';...
%    'a.feedbackFractionExc = a.feedbackFractionExc + changeamt;settext(fbetext,fbetextbase,a.feedbackFractionExc);';...
%    'a.feedbackFractionInh = a.feedbackFractionInh + changeamt;settext(fbitext,fbitextbase,a.feedbackFractionInh);';...
    'changeamt = changeamt * 10;settext(inctext,inctextbase,changeamt);';...
    'angle = angle + changeamt; a.kOrient.x = cos(angle); a.kOrient.y = sin(angle);settext(angletext,angletextbase,angle);';...
    'wavelen = wavelen + changeamt; a.setKFreq(2*pi/wavelen); settext(freqtext,freqtextbase,2*pi/a.kFreq);';...
    'phase = phase + changeamt; a.kPhase = phase; settext(phasetext,phasetextbase,a.kPhase);';...
    'a.setAmbient(a.rateStrongAmbE,a.rateWeakAmbE,a.rateStrongAmbI,a.rateWeakAmbI + changeamt);settext(ambwkItext,ambwkItextbase,a.rateWeakAmbI);';...
};

for i = 1 : length(leftCallback)
    uicontrol(uiPanel,'Style', 'pushbutton', 'String', '<',...
    'units','normalized',...
    'position', [.1 textPosition(i)-.02 .05 .05],...
    'callback', leftCallback{i});
uicontrol(uiPanel,'Style', 'pushbutton', 'String', '>',...
    'units','normalized',...
    'position', [.15 textPosition(i)-.02 .05 .05],...
    'callback', rightCallback{i});
end

prefvis = 1;
orthvis = 1;
fullvis = 1;
visstring = 'visible';
onstring = 'on';
offstring = 'off';
uicontrol(f,'Style', 'pushbutton', 'String', '',...
    'position', [50 50 20 20],...
    'callback', 'if fullvis==1; set(pl,visstring,offstring);fullvis=0;else;set(pl,visstring,onstring);fullvis=1;end;')

keepgoing = 1;
while keepgoing
    
%frdata(1) = a.timeLocalFR(a.currentCounter+1)/50/27;
gLGNdata(1) = a.gLGNExc(inDrivePatchExc(whichNeuron),1);
gAMPAdata(1)= a.gAMPAExc(inDrivePatchExc(whichNeuron),1);
gGABAdata(1)= a.gGABAExc(inDrivePatchExc(whichNeuron),1);
%vdata(1) =    a.vExc(inDrivePatchExc(whichNeuron),1);
set(lin,'xdata',dataCollectTimes(1)*[0 0])
for i = 2 : length(dataCollectTimes)
    disp(['i=' num2str(i)])
    for j = 1 : skipInSteps
        a.step();
    end
    title(ax,a.numInputsExc(inDrivePatchExc(whichNeuron)));
    
    %frdata(i) = a.timeLocalFR(a.currentCounter+1)/50/27;
    gLGNdata(i) = a.gLGNExc(inDrivePatchExc(whichNeuron),1);
    gAMPAdata(i)= a.gAMPAExc(inDrivePatchExc(whichNeuron),1);
    gGABAdata(i)= a.gGABAExc(inDrivePatchExc(whichNeuron),1);
    %vdata(i) =    a.vExc(inDrivePatchExc(whichNeuron),1);
    
    %set(pl,'YData',frdata);
    set(plp1, 'YData',gLGNdata);
    set(plp2, 'YData',gAMPAdata);
    set(plp3, 'YData',gGABAdata);
    %set(plp4, 'YData',vdata);
    set(lin, 'XData', dataCollectTimes(i)*[1 1])
            
    drawnow();
end
end

%% animate E and I spikes
f = figure('position', [400, 400, 850, 850]);
duration = 300000;

%record patch frs
optiPatchNS = 0;
orthPatchNS = 0;
tStart = a.t;

locsX = zeros(1,a.numExcitatory);
locsY = zeros(1,a.numExcitatory);
for i=1:a.numExcitatory;locsX(i)=a.locExc(i).x;locsY(i)=a.locExc(i).y;end;

%I locations
iLocsX = zeros(1,a.numExcitatory);
iLocsY = zeros(1,a.numExcitatory);
for i=1:a.numInhibitory;iLocsX(i)=a.locInh(i).x;iLocsY(i)=a.locInh(i).y;end;

sleeping = 1 : a.numExcitatory;
threshold = ones(1,a.numExcitatory);
inRefrac = 1;%sleeping(a.sleepExc < 0);
hold on;
p2 = plot(locsX(inRefrac),locsY(inRefrac),'.','color','red','markersize',20);
hold off;

iInRefrac = 1;%sleeping(a.sleepInh < 0);
hold on;
p3 = plot(iLocsX(iInRefrac),iLocsY(iInRefrac),'.','color','blue','markersize',20);
hold off;

xlim([0 a.patchWidth]);ylim([0 a.patchWidth]);title([num2str(i*a.dt) ' ms']);%view([-90, 0])
%figure;
ns = 0;
changeamt = .01;
uicontrol(f,'Style', 'pushbutton', 'String', 'dec sie',...
    'position', [20 20 100 20],...
    'callback', 'a.setSIEOverSEE(a.sIE/a.sEE - changeamt);');
uicontrol(f,'Style', 'pushbutton', 'String', 'inc sie',...
    'position', [140 20 100 20],...
    'callback', 'a.setSIEOverSEE(a.sIE/a.sEE + changeamt);');
uicontrol(f,'Style', 'pushbutton', 'String', 'dec changeamt',...
    'position', [260 20 100 20],...
    'callback', 'changeamt = changeamt / 10;');
uicontrol(f,'Style', 'pushbutton', 'String', 'inc changeamt',...
    'position', [380 20 100 20],...
    'callback', 'changeamt = changeamt * 10;');
uicontrol(f,'Style', 'pushbutton', 'String', 'inc drive',...
    'position', [500 20 100 20],...
    'callback', 'a.eps = a.eps + 1;');
uicontrol(f,'Style', 'pushbutton', 'String', 'dec drive',...
    'position', [620 20 100 20],...
    'callback', 'a.eps = a.eps - 1;');
uicontrol(f,'Style', 'pushbutton', 'String', 'reset fr',...
    'position', [740 20 100 20],...
    'callback', 'tStart=a.t;optiPatchNS=0;orthPatchNS=0;');

%[~, I] = sort(a.numInputsExc);
%[~, IInv] = sort(I);
%[~, IInh] = sort(a.numInputsInh);
%[~, IInvInh] = sort(IInh);

for i = 1 : duration
    %plot3(1:a.numExcitatory, a.vExc,(a.gLGNExc(:,1) + a.gAMPAExc(:,1) + a.gNMDAExc(:,1))*a.eRev + a.gGABAExc(:,1) *a.iRev,'.','color','r');
    
    %find neuron in refractory
    inRefrac = sleeping(a.sleepExc < 0);
    set(p2,'xdata',locsX(inRefrac),...
        'ydata',locsY(inRefrac)); 
    
    iInRefrac = sleeping(a.sleepInh < 0);
    set(p3,'xdata',iLocsX(iInRefrac),...
        'ydata',iLocsY(iInRefrac));
    ns = ns + a.numSpikesNowExc;
    title({[num2str(i*a.dt) ' ms']...
        ['fr=' num2str(a.timeLocalFR(a.currentCounter+1)/50/12,2)]...
        ['opt patch fr=' num2str(1000*optiPatchNS/(numInDrivePatchExc*(a.t-tStart)))]...
        ['orth. patch fr=' num2str(1000*orthPatchNS/(numInDrivePatchExc*(a.t-tStart)))]...
        ['S^{IE}/S^{EE}=' num2str(a.sIE/a.sEE)]...
        ['contrast=' num2str(a.eps)]});% num2str(a.currentCounter) num2str(a.numStepsMod10) num2str(a.eSpikes / 12 / a.t,2) num2str(ns / 12 / a.t,2)});
    
    %set the fr counters    
    optiPatchNS = optiPatchNS + sum(ismember(a.spikesNowExc(1:a.numSpikesNowExc)+1, inDrivePatchExc(1:numInDrivePatchExc)));
    orthPatchNS = orthPatchNS + sum(ismember(a.spikesNowExc(1:a.numSpikesNowExc)+1, inNoDrivePatchExc(1:numInNoDrivePatchExc)));
    
    %pause(.1)
    drawnow();%
    a.step();
end

%% animate LGN conductance to E
f = figure('position', [400, 400, 750, 750]);
duration = 300000;

disp('Turning off cortical synapses, do not forget to reset!')
a.sEE = 0;
a.sIE = 0;
a.sEI = 0;
a.sII = 0;

%convert java locations to MATLAB
locsX = zeros(1,a.numExcitatory);
locsY = zeros(1,a.numExcitatory);
for i=1:a.numExcitatory;locsX(i)=a.locExc(i).x;locsY(i)=a.locExc(i).y;end;

im = vec2mat(a.gLGNExc(:,1),110) * 10;

%initially plot locations of E cells
%facecolors = [ones(a.numExcitatory,1),zeros(a.numExcitatory,2)];
p = image(im); colormap(gray);% truesize([500,500],);
%truesize;
%p2 = plot3(locsX,locsY,a.gLGNExc(:,1),'.','color','red');
xlim([0 a.patchWidth]);ylim([0 a.patchWidth]);title([num2str(i*a.dt) ' ms']);%view([-90, 0])

%pushbuttons
uicontrol(f,'Style', 'pushbutton', 'String', 'inc drive',...
    'position', [500 20 100 20],...
    'callback', 'a.eps = a.eps + 1;');
uicontrol(f,'Style', 'pushbutton', 'String', 'dec drive',...
    'position', [620 20 100 20],...
    'callback', 'a.eps = a.eps - 1;');

for i = 1 : duration
    %set the color to be the LGN conductance
    %facecolors(:,1) = a.gLGNExc(:,1);
    %set(p,'cdata',vec2mat(a.gLGNExc(:,1),110) * 10);
    %truesize([1500,1500]);
    p = image(vec2mat(a.gLGNExc(:,1),110) * 1000);% truesize([500,500],);
    %truesize;
    title({[num2str(i*a.dt) ' ms'] num2str(a.timeLocalFR(a.currentCounter+1)/50/12,2) ['S^{IE}/S^{EE}=' num2str(a.sIE/a.sEE)] ['contrast=' num2str(a.eps)] ['sample LGN cell fr=' num2str(a.RON(a.locLGNON(1),a.t)*1000)] ['sample LGN cell fr=' num2str(a.RON(a.locLGNON(2),a.t)*1000)]});% num2str(a.currentCounter) num2str(a.numStepsMod10) num2str(a.eSpikes / 12 / a.t,2) num2str(ns / 12 / a.t,2)});
    
    %pause(.1)
    drawnow();%
    a.step();
end

%% animate LGN conductance to I
f = figure('position', [400, 400, 750, 750]);
duration = 300000;

disp('Turning off cortical synapses, do not forget to reset!')
a.sEE = 0;
a.sIE = 0;
a.sEI = 0;
a.sII = 0;

%convert java locations to MATLAB
locsX = zeros(1,a.numInhibitory);
locsY = zeros(1,a.numInhibitory);
for i=1:a.numInhibitory;locsX(i)=a.locInh(i).x;locsY(i)=a.locInh(i).y;end;

%im = vec2mat(a.gLGNExc(:,1),110) * 10;

%initially plot locations of E cells
%facecolors = [ones(a.numExcitatory,1),zeros(a.numExcitatory,2)];
%p = image(im); colormap(gray);% truesize([500,500],);
%truesize;
p2 = plot3(locsX,locsY,a.gLGNInh(:,1),'.','color','blue');
xlim([0 a.patchWidth]);ylim([0 a.patchWidth]);title([num2str(i*a.dt) ' ms']);%view([-90, 0])

%pushbuttons
uicontrol(f,'Style', 'pushbutton', 'String', 'inc drive',...
    'position', [500 20 100 20],...
    'callback', 'a.eps = a.eps + 1;');
uicontrol(f,'Style', 'pushbutton', 'String', 'dec drive',...
    'position', [620 20 100 20],...
    'callback', 'a.eps = a.eps - 1;');

for i = 1 : duration
    %set the color to be the LGN conductance
    %facecolors(:,1) = a.gLGNExc(:,1);
    set(p2,'zdata',a.gLGNInh(:,1));
    %set(p,'cdata',vec2mat(a.gLGNExc(:,1),110) * 10);
    %truesize([1500,1500]);
    %p = image(vec2mat(a.gLGNExc(:,1),110) * 1000);% truesize([500,500],);
    %truesize;
    title({[num2str(i*a.dt) ' ms'] num2str(a.timeLocalFR(a.currentCounter+1)/50/12,2) ['S^{IE}/S^{EE}=' num2str(a.sIE/a.sEE)] ['contrast=' num2str(a.eps)] ['sample LGN cell fr=' num2str(a.RON(a.locLGNON(1),a.t)*1000)] ['sample LGN cell fr=' num2str(a.RON(a.locLGNON(2),a.t)*1000)]});% num2str(a.currentCounter) num2str(a.numStepsMod10) num2str(a.eSpikes / 12 / a.t,2) num2str(ns / 12 / a.t,2)});
    
    zlim([0 .1]);
    %pause(.1)
    drawnow();%
    a.step();
end

%%
figure;
%a.resetNetState();
numNeurons = a.numExcitatory + a.numInhibitory;

locsX = zeros(1,a.numExcitatory);
locsY = zeros(1,a.numExcitatory);
for i=1:a.numExcitatory;locsX(i)=a.locExc(i).x;locsY(i)=a.locExc(i).y;end;

%I locations
iLocsX = zeros(1,a.numExcitatory);
iLocsY = zeros(1,a.numExcitatory);
for i=1:a.numInhibitory;iLocsX(i)=a.locInh(i).x;iLocsY(i)=a.locInh(i).y;end;


p1 = plot3(locsX,locsY,a.vExc,'.','color','r');

sleeping = 1 : a.numExcitatory;
threshold = ones(1,a.numExcitatory);
inRefrac = 1;%sleeping(a.sleepExc < 0);
hold on;
p2 = plot3(locsX(inRefrac),locsY(inRefrac), threshold(inRefrac),'.','color',[0 .5 0],'markersize',20);
hold off;

iInRefrac = 1;%sleeping(a.sleepInh < 0);
hold on;
p3 = plot3(iLocsX(iInRefrac),iLocsY(iInRefrac), threshold(iInRefrac),'.','color','blue','markersize',20);
hold off;

xlim([0 a.patchWidth]);ylim([0 a.patchWidth]);zlim([0 1]);title([num2str(i*a.dt) ' ms']);%view([-90, 0])
%figure;
for i = 1 : duration
    %plot3(1:a.numExcitatory, a.vExc,(a.gLGNExc(:,1) + a.gAMPAExc(:,1) + a.gNMDAExc(:,1))*a.eRev + a.gGABAExc(:,1) *a.iRev,'.','color','r');
    set(p1,'zdata',a.vExc);
    %find neuron in refractory
    inRefrac = sleeping(a.sleepExc < 0);
    set(p2,'xdata',locsX(inRefrac),...
        'ydata',locsY(inRefrac),...
        'zdata',threshold(inRefrac));
    
    iInRefrac = sleeping(a.sleepInh < 0);
    set(p3,'xdata',iLocsX(iInRefrac),...
        'ydata',iLocsY(iInRefrac),...
        'zdata',threshold(iInRefrac));
    
    title([num2str(i*a.dt) ' ms']); box on;
    
    drawnow();%
    a.stepRecordingSpikes();
    
end


%% animate time-local firing rates for LGN ON cells
figure;
numNeurons = a.numExcitatory + a.numInhibitory;

locsX = zeros(1,a.numLGN);
locsY = zeros(1,a.numLGN);
for i=1:a.numLGN;locsX(i)=a.locLGNON(i).x;locsY(i)=a.locLGNON(i).y;end;

LGNRate = zeros(1,a.numLGN);
for i=1:a.numLGN;LGNRate(i)=a.RONDrift(i,a.t);end;
p1 = plot3(locsX,locsY,LGNRate,'.','color','r');

xlim([0 a.patchWidth]);ylim([0 a.patchWidth]);zlim([0 1]);title([num2str(i*a.dt) ' ms']);%view([-90, 0])
%figure;
for i = 1 : duration
    %plot3(1:a.numExcitatory, a.vExc,(a.gLGNExc(:,1) + a.gAMPAExc(:,1) + a.gNMDAExc(:,1))*a.eRev + a.gGABAExc(:,1) *a.iRev,'.','color','r');
    for j=1:a.numLGN;LGNRate(j)=a.RONDrift(j,a.t);end;
    set(p1,'zdata',LGNRate);
    %find neuron in refractory
        
    title([num2str(i*a.dt) ' ms']); box on;
    
    drawnow();%
    a.stepRecordingSpikes();
    
end

%% animate type 2
duration = 3000;

a.changeSLGN(.04);
a.changeSIE(0.003);

numNeurons = a.numExcitatory + a.numInhibitory;

p1 = plot3(1:a.numExcitatory, a.vExc,a.gNMDAExc(:,1),'.','color','r');
hold on;
p2 = plot3(a.numExcitatory+1:numNeurons, a.vInh,a.gNMDAInh(:,1),'.','color','b');
hold off;
xlim([0 numNeurons]);ylim([-.66 1]);zlim([-.01 .05]);title([num2str(i*a.dt) ' ms']);view([-90, 0])
%figure;
for i = 1 : duration
    %plot3(1:a.numExcitatory, a.vExc,(a.gLGNExc(:,1) + a.gAMPAExc(:,1) + a.gNMDAExc(:,1))*a.eRev + a.gGABAExc(:,1) *a.iRev,'.','color','r');
    set(p1,'ydata', a.vExc);
    set(p1,'zdata',a.gNMDAExc(:,1));
    %plot3(a.numExcitatory+1:numNeurons, a.vInh,(a.gLGNInh(:,1) + a.gAMPAInh(:,1) + a.gNMDAInh(:,1))*a.eRev + a.gGABAInh(:,1) *a.iRev,'.','color','b');
    set(p2,'ydata', a.vInh);
    set(p2,'zdata', a.gNMDAInh(:,1));
    title([num2str(i*a.dt) ' ms']);
    
    %F(i) = getframe(gcf);
    drawnow();%
    a.step();
    
end

%% animate type 3
duration = 3000;

a.changeSLGN(.04);
a.changeSIE(0.003);

numNeurons = a.numExcitatory + a.numInhibitory;

p1 = plot3(a.gGABAExc(:,1),a.vExc,a.gAMPAExc(:,1),'.','color','r');
hold on;
p2 = plot3(a.gGABAInh(:,1), a.vInh,a.gAMPAInh(:,1),'.','color','b');
hold off;
xlim([0 .07]);ylim([-.66 1]);zlim([0 .01]);title([num2str(i*a.dt) ' ms']);%view([-90, 0])
%figure;
for i = 1 : duration
    %plot3(1:a.numExcitatory, a.vExc,(a.gLGNExc(:,1) + a.gAMPAExc(:,1) + a.gNMDAExc(:,1))*a.eRev + a.gGABAExc(:,1) *a.iRev,'.','color','r');
    set(p1,'xdata', a.gGABAExc(:,1));
    set(p1,'ydata', a.vExc);
    set(p1,'zdata',a.gAMPAExc(:,1));
    %plot3(a.numExcitatory+1:numNeurons, a.vInh,(a.gLGNInh(:,1) + a.gAMPAInh(:,1) + a.gNMDAInh(:,1))*a.eRev + a.gGABAInh(:,1) *a.iRev,'.','color','b');
    set(p2,'xdata', a.gGABAInh(:,1));
    set(p2,'ydata', a.vInh);
    set(p2,'zdata', a.gAMPAInh(:,1));
    title([num2str(i*a.dt) ' ms']);
    
    %F(i) = getframe(gcf);
    drawnow();%
    a.step();
    
end

%% animate time-local LGN input rates cortical cells
figure;
numNeurons = a.numExcitatory + a.numInhibitory;
a.eps=10;
locsX = zeros(1,a.numExcitatory);
locsY = zeros(1,a.numExcitatory);
for i=1:a.numExcitatory;locsX(i)=a.locExc(i).x;locsY(i)=a.locExc(i).y;end;

cELGNON = a.cELGNON + 1;
nELGNON = a.nELGNON;

cELGNOFF = a.cELGNOFF + 1;
nELGNOFF = a.nELGNOFF;

t=a.t;
dt=1;
inputRate = zeros(1,a.numExcitatory);

for i=1:a.numLGN
    %go through LGN ON cells and distribute firing rates to the cortical cells
    outRateON = a.RON(a.locLGNON(i),t);
    outRateOFF = a.ROFF(a.locLGNOFF(i),t);
    inputRate(cELGNON(i,1:nELGNON(i))) = inputRate(cELGNON(i,1:nELGNON(i))) + outRateON;
    inputRate(cELGNOFF(i,1:nELGNOFF(i))) = inputRate(cELGNOFF(i,1:nELGNOFF(i))) + outRateOFF;
end
convMatrix = [ones(1,10) zeros(1,100) ones(1,10) zeros(1,100) ones(1,10)]/30;
avgRate = conv(inputRate, convMatrix,'same');

p1 = plot3(locsX,locsY,avgRate,'.','color','r');

xlim([0 a.patchWidth]);ylim([0 a.patchWidth]);zlim([0 1]);title([num2str(i*a.dt) ' ms']);%view([-90, 0])
%figure;
for t = 0:dt:1000
    %plot3(1:a.numExcitatory, a.vExc,(a.gLGNExc(:,1) + a.gAMPAExc(:,1) + a.gNMDAExc(:,1))*a.eRev + a.gGABAExc(:,1) *a.iRev,'.','color','r');
    inputRate = zeros(1,a.numExcitatory);
    
    for i=1:a.numLGN
        %go through LGN ON cells and distribute firing rates to the cortical cells
        outRateON = a.RON(a.locLGNON(i),t);
        outRateOFF = a.ROFF(a.locLGNOFF(i),t);
        inputRate(cELGNON(i,1:nELGNON(i))) = inputRate(cELGNON(i,1:nELGNON(i))) + outRateON;
        inputRate(cELGNOFF(i,1:nELGNOFF(i))) = inputRate(cELGNOFF(i,1:nELGNOFF(i))) + outRateOFF;
    end
    avgRate = conv(inputRate, convMatrix,'same');
    set(p1,'zdata',avgRate);
    %find neuron in refractory
        
    title([num2str(t) ' ms']); box on;
    
    drawnow();%
    %a.stepRecordingSpikes();
    
end

%% wait
duration = 1000;
a.changeSLGN(.00);
a.changeSIE(0.0);
for i = 1 : duration; a.step();end

%% look at any type of presynaptic connections
preType = 'E'
postType= 'I'
switch([preType '->' postType])
    case 'E->E'
        numPre = a.numExcitatory;
        numPost= a.numExcitatory;
        getLocPostX = @(i) a.locExc(i).x;
        getLocPostY = @(i) a.locExc(i).y;
        getAllPre = @(i) a.cPreEExc(i,1:a.numEPreExc(i)) + 1;
        getNumPre = @(i) a.numEPreExc(i);
        getLocPreX = @(j) a.locExc(j).x;
        getLocPreY = @(j) a.locExc(j).y;
        titletext = ' E cells (blue) presynaptic to sample E cell (red)';
    case 'I->E'
        numPre = a.numInhibitory;
        numPost= a.numExcitatory;
        getLocPostX = @(i) a.locExc(i).x;
        getLocPostY = @(i) a.locExc(i).y;
        getAllPre = @(i) a.cPreIExc(i,1:a.numIPreExc(i)) + 1;
        getNumPre = @(i) a.numIPreExc(i);
        getLocPreX = @(j) a.locInh(j).x;
        getLocPreY = @(j) a.locInh(j).y;
        titletext = ' I cells (blue) presynaptic to sample E cell (red)';
    case 'E->I'
        numPre = a.numExcitatory;
        numPost= a.numInhibitory;
        getLocPostX = @(i) a.locInh(i).x;
        getLocPostY = @(i) a.locInh(i).y;
        getAllPre = @(i) a.cPreEInh(i,1:a.numEPreInh(i)) + 1;
        getNumPre = @(i) a.numEPreInh(i);
        getLocPreX = @(j) a.locExc(j).x;
        getLocPreY = @(j) a.locExc(j).y;
        titletext = ' E cells (blue) presynaptic to sample I cell (red)';
    case 'I->I'
        numPre = a.numInhibitory;
        numPost= a.numInhibitory;
        getLocPostX = @(i) a.locInh(i).x;
        getLocPostY = @(i) a.locInh(i).y;
        getAllPre = @(i) a.cPreIInh(i,1:a.numIPreInh(i)) + 1;
        getNumPre = @(i) a.numIPreInh(i);
        getLocPreX = @(j) a.locInh(j).x;
        getLocPreY = @(j) a.locInh(j).y;
        titletext = ' I cells (blue) presynaptic to sample I cell (red)';
end

figure; 
preX = zeros(1,numPre);
preY = zeros(1,numPre);
for k = 1 : 1
    i = randi(numPost,1);
    postX = getLocPostX(i);
    postY = getLocPostY(i);
    allPreSyn = getAllPre(i);
    
    for j = 1 : getNumPre(i)
        preX(j) = getLocPreX(allPreSyn(j));
        preY(j) = getLocPreY(allPreSyn(j));
    end
    
    scatter(preX(1:getNumPre(i)),preY(1:getNumPre(i)),'b'); 
      
    hold on;
    scatter(postX,postY,'r','facecolor','r');
    hold off;
    
    xlim([0 a.patchWidth]);ylim([0 a.patchWidth])
    xlabel('micrometers');
    ylabel('micrometers');
    title([ num2str(getNumPre(i)) titletext])
    drawnow();
    
    pause(3)
end

%% make 1 sec runs
runDurationInMs = 1000;
runDurationInSteps = runDurationInMs / a.dt;

numExcitatory = a.numExcitatory;
numInhibitory = a.numInhibitory;

%==stats collecting vectors==%
MAXSPIKESExc = .04 * runDurationInMs * numExcitatory;
MAXSPIKESInh = .08 * runDurationInMs * numInhibitory;
spikesExc = zeros(1, MAXSPIKESExc + numExcitatory);
spikeTimesExc = zeros(1,MAXSPIKESExc + numExcitatory);
spikesInh = zeros(1, MAXSPIKESInh + numInhibitory);
spikeTimesInh = zeros(1,MAXSPIKESInh + numInhibitory);

disp('running');
t1 = clock;
numSpikesExc = 0; numSpikesInh = 0; recordLengthSteps = 0;
for step = 1 : runDurationInSteps
    if mod(step, 1000) == 1
        disp(step);
    end
    a.step();

    %record this step
    numSpikesNowExc = a.numSpikesNowExc;
    spikesExc(numSpikesExc+1 : numSpikesExc+numSpikesNowExc) = a.spikesNowExc(1 : numSpikesNowExc);
    spikeTimesExc(numSpikesExc+1 : numSpikesExc+numSpikesNowExc) = a.spikeTimesNowExc(1 : numSpikesNowExc);

    numSpikesExc = numSpikesExc + numSpikesNowExc;
    if numSpikesExc > MAXSPIKESExc; break; end;

    numSpikesNowInh = a.numSpikesNowInh;
    spikesInh(numSpikesInh+1 : numSpikesInh+numSpikesNowInh) = a.spikesNowInh(1 : numSpikesNowInh);
    spikeTimesInh(numSpikesInh+1 : numSpikesInh+numSpikesNowInh) = a.spikeTimesNowInh(1 : numSpikesNowInh);

    numSpikesInh = numSpikesInh + numSpikesNowInh;
    if numSpikesInh > MAXSPIKESInh; break; end;

    %end recording this step
end


figure('position',[500,100,1200,700]);
        %==============show raster================%
        startTime = a.t-runDurationInMs;
        timeWindow = [startTime startTime+1000];
        
        excitatoryDrivePatchEnumerator = zeros(1,numExcitatory);
        inhibitoryDrivePatchEnumerator = zeros(1,numInhibitory);
        excitatoryDrivePatchEnumerator(inDrivePatchExc(1:numInDrivePatchExc)) = 1:numInDrivePatchExc;
        inhibitoryDrivePatchEnumerator(inDrivePatchInh(1:numInDrivePatchInh)) = numInDrivePatchExc + 1:numInDrivePatchInh + numInDrivePatchExc;
        
        whichInDrivePatchExc = ismember(spikesExc(1:numSpikesExc)+1, inDrivePatchExc(1:numInDrivePatchExc));
        spikesInDrivePatchExc = spikesExc(whichInDrivePatchExc)+1;
        tempx = spikeTimesExc(whichInDrivePatchExc);
        %[~,sortedIndexes] = sort(numInputsExc(inDrivePatchExc(1:numInDrivePatchExc))); %sort by number of inputs        
        sortBy = double(a.numEPreExc(inDrivePatchExc(1:numInDrivePatchExc))) ...
              ./ double(a.numIPreExc(inDrivePatchExc(1:numInDrivePatchExc))) ...
              +  double(1e5 * a.numInputsExc(inDrivePatchExc(1:numInDrivePatchExc)));
        [sorted,sortedIndexes] = sort(sortBy); %sort by number of inputs        
        newIndexes = zeros(numInDrivePatchExc,1);
        newIndexes(sortedIndexes) = 1 : numInDrivePatchExc;
        tempy = newIndexes(excitatoryDrivePatchEnumerator(spikesInDrivePatchExc));
        %cutoff = find(tempx > 20000);
        scrollsubplot(3,1,1); scatter(tempx,tempy,'r','.');        
        hold on;
        for i = 2 : numInDrivePatchExc
            if sorted(i) - sorted(i-1) > 1e4
                line([0,10000],(i-.5)*[1 1]);
            end
        end

        whichInDrivePatchInh = ismember(spikesInh(1:numSpikesInh)+1, inDrivePatchInh(1:numInDrivePatchInh));
        spikesInDrivePatchInh = spikesInh(whichInDrivePatchInh)+1;
        tempx = spikeTimesInh(whichInDrivePatchInh);
        tempy = inhibitoryDrivePatchEnumerator(spikesInDrivePatchInh);
        %cutoff = find(tempx > 20000);
        scatter(tempx,tempy,'b','.');
        hold off;

        title(['raster (optimal patch); S^{EI}/S^{EE} = ' num2str(parameters(paramSelector,4)) '; drive=+' num2str(parameters(paramSelector,6)) ' contrast; S^{IE}/S^{EE}=' num2str(parameters(paramSelector,5)) '; rates (E, I)=(' num2str(length(spikesInDrivePatchExc)/(runDurationInMs/1000)/numInDrivePatchExc) ', ' num2str(length(spikesInDrivePatchInh)/(runDurationInMs/1000)/numInDrivePatchInh) ')']);
        ylim([0 numInDrivePatchExc+numInDrivePatchInh])
        xlim(timeWindow)

        %=====Undriven patch=====%
        excitatoryNoDrivePatchEnumerator = zeros(1,numExcitatory);
        inhibitoryNoDrivePatchEnumerator = zeros(1,numInhibitory);
        excitatoryNoDrivePatchEnumerator(inNoDrivePatchExc(1:numInNoDrivePatchExc)) = 1:numInNoDrivePatchExc;
        inhibitoryNoDrivePatchEnumerator(inNoDrivePatchInh(1:numInNoDrivePatchInh)) = numInNoDrivePatchExc + 1:numInNoDrivePatchInh + numInNoDrivePatchExc;


        whichInNoDrivePatchExc = ismember(spikesExc(1:numSpikesExc)+1, inNoDrivePatchExc(1:numInNoDrivePatchExc));
        spikesInNoDrivePatchExc = spikesExc(whichInNoDrivePatchExc)+1;
        [~,sortedIndexes] = sort(numInputsExc(inNoDrivePatchExc(1:numInNoDrivePatchExc))); %sort by number of inputs        
        newIndexes = zeros(numInNoDrivePatchExc,1);
        newIndexes(sortedIndexes) = 1 : numInNoDrivePatchExc;
        scrollsubplot(3,1,2); scatter(spikeTimesExc(whichInNoDrivePatchExc),newIndexes(excitatoryNoDrivePatchEnumerator(spikesInNoDrivePatchExc)),'r','.');
        hold on;

        whichInNoDrivePatchInh = ismember(spikesInh(1:numSpikesInh)+1, inNoDrivePatchInh(1:numInNoDrivePatchInh));
        spikesInNoDrivePatchInh = spikesInh(whichInNoDrivePatchInh)+1;
        scatter(spikeTimesInh(whichInNoDrivePatchInh),inhibitoryNoDrivePatchEnumerator(spikesInNoDrivePatchInh),'b','.');
        hold off;

        title(['raster (perpendicular patch); rates (E, I)=(' num2str(length(spikesInNoDrivePatchExc)/(runDurationInMs/1000)/numInNoDrivePatchExc) ', ' num2str(length(spikesInNoDrivePatchInh)/(runDurationInMs/1000)/numInNoDrivePatchInh) ')']);
        ylim([0 numInNoDrivePatchExc+numInNoDrivePatchInh])
        xlim(timeWindow)