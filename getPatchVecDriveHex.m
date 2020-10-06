% get a vector which picks out only the neurons in a small patch of cortex
% warning: do not use in an i-loop
orientation = 0*pi/8;
orientationVec = [cos(orientation) sin(orientation)]; %perp to stripes
HCCenter = 1.5*a.HCWidth * [1, 1];
m = max(abs(cos(pi-2*orientation)), abs(sin(pi-2*orientation)));
prefPatchCenter = HCCenter + a.HCWidth/2*[cos(pi-2*orientation) sin(pi-2*orientation)]/m;
orthPatchCenter = HCCenter - a.HCWidth/2*[cos(pi-2*orientation) sin(pi-2*orientation)]/m;
radiusFraction = 1/16;
recordingPatchRadius = a.patchWidth * radiusFraction;
recordingPatchCenterDrive = sim.MyPoint(prefPatchCenter(1),prefPatchCenter(2));
%a.patchWidth*(1/3.0+1/12),a.patchWidth*(1/2.0+1/9)

numInDrivePatchExc = 0; inDrivePatchExc = zeros(1,ceil(a.numExcitatory * radiusFraction^2) + 500);
numInDrivePatchInh = 0; inDrivePatchInh = zeros(1,ceil(a.numInhibitory * radiusFraction^2) + 500);
for i = 1 : a.numExcitatory
    if a.locExc(i).distance(recordingPatchCenterDrive) < recordingPatchRadius
        %then include neuron i
        numInDrivePatchExc = numInDrivePatchExc + 1;
        inDrivePatchExc(numInDrivePatchExc) = i;
    end
end
for i = 1 : a.numInhibitory
    if a.locInh(i).distance(recordingPatchCenterDrive) < recordingPatchRadius
        %then include neuron i
        numInDrivePatchInh = numInDrivePatchInh + 1;
        inDrivePatchInh(numInDrivePatchInh) = i;
    end
end
%==driven part done, do undriven part==%
recordingPatchCenterNoDrive = sim.MyPoint(orthPatchCenter(1),orthPatchCenter(2));
%a.patchWidth*(2/3.0-1/12),a.patchWidth*(1/2.0+1/9)
numInNoDrivePatchExc = 0; inNoDrivePatchExc = zeros(1,ceil(a.numExcitatory * radiusFraction^2) + 500);
numInNoDrivePatchInh = 0; inNoDrivePatchInh = zeros(1,ceil(a.numInhibitory * radiusFraction^2) + 500);
for i = 1 : a.numExcitatory
    if a.locExc(i).distance(recordingPatchCenterNoDrive) < recordingPatchRadius
        %then include neuron i
        numInNoDrivePatchExc = numInNoDrivePatchExc + 1;
        inNoDrivePatchExc(numInNoDrivePatchExc) = i;
    end
end
for i = 1 : a.numInhibitory
    if a.locInh(i).distance(recordingPatchCenterNoDrive) < recordingPatchRadius
        %then include neuron i
        numInNoDrivePatchInh = numInNoDrivePatchInh + 1;
        inNoDrivePatchInh(numInNoDrivePatchInh) = i;
    end
end
%==combine them==%
numInPatchExc = numInDrivePatchExc + numInNoDrivePatchExc; inPatchExc = zeros(1,numInPatchExc);
inPatchExc = [inDrivePatchExc(1:numInDrivePatchExc), inNoDrivePatchExc(1:numInNoDrivePatchExc)];
numInPatchInh = numInDrivePatchInh + numInNoDrivePatchInh; inPatchInh = zeros(1,numInPatchInh);
inPatchInh = [inDrivePatchInh(1:numInDrivePatchInh), inNoDrivePatchInh(1:numInNoDrivePatchInh)];

%==setup the patch stats in java==%
for i=1:numInDrivePatchExc; a.membershipExc(inDrivePatchExc(i)) = 1; end;
for i=1:numInNoDrivePatchExc; a.membershipExc(inNoDrivePatchExc(i)) = 2; end;
for i=1:numInDrivePatchInh; a.membershipInh(inDrivePatchInh(i)) = 1; end;
for i=1:numInNoDrivePatchInh; a.membershipInh(inNoDrivePatchInh(i)) = 2; end;

%% view them
for i=1:numInPatchExc; locsX(i) = a.locExc(inPatchExc(i)).x;locsY(i) = a.locExc(inPatchExc(i)).y;end;
figure;scatter(locsX,locsY);xlim([0 1500]);ylim([0 1500])
