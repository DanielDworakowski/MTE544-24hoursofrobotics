% Set up environment
rng('default')
posMinBound = [0 0];
posMaxBound = [50 40];
numObsts = 70;
endPos = [1 1];
startPos = [48.5 38.5];

minLen.a = 1;
maxLen.a = 3;
minLen.b = 2;
maxLen.b = 6;

obstBuffer = 0.5;
maxCount = 10000;

[aObsts,bObsts,obsPtsStore] = polygonal_world(posMinBound, posMaxBound, minLen, maxLen, numObsts, startPos, endPos, obstBuffer, maxCount);

plotEnvironment(obsPtsStore,posMinBound, posMaxBound, startPos, endPos, 1);
