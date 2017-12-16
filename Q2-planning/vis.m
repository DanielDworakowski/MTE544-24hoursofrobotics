clc; clear; close all;
run('planningEnv.m');


% Create visibility graph
tic;
allGraphPts = [];
obsEdges = [];
for i=1:numObsts
    obsPts{i} = [obsPtsStore(:,2*(i-1)+1:2*i)];
    allGraphPts = [allGraphPts;obsPts{i}];
    obsEdges = [obsEdges; [obsPts{i}(1:end,:) obsPts{i}([2:end 1],:)]];
end
allGraphPts = [allGraphPts; startPos; endPos];
n = length(allGraphPts(:,1));
% Set initial link possibilities
A = ones(n,n) - eye(n);
% Remove links in obstacles, but not on edges
ptCount = 0;
for i=1:numObsts
    numOPts = length(obsPts{i});
    A(ptCount+1:ptCount+numOPts,ptCount+1:ptCount+numOPts) = zeros(numOPts,numOPts);
    for j = 1:numOPts-1
        A(ptCount+j,ptCount+j+1) = 1;
        A(ptCount+j+1,ptCount+j) = 1;
    end
    A(ptCount+1,ptCount+numOPts) = 1;
    A(ptCount+numOPts,ptCount+1) = 1;
    ptCount=ptCount+numOPts;
end


% Check for collisions among remaining links
for i=1:n
    for j=i:n
        inColl = CheckCollision(allGraphPts(i,:),allGraphPts(j,:),obsEdges);
        if (inColl)
            A(i,j) = 0;
            A(j,i) = 0;
        end
        if (A(i,j)~=0)
            D(i,j) = norm(allGraphPts(i,:)-allGraphPts(j,:));
            D(j,i) = D(i,j);
        end
    end
end

for i=1:n
    for j=i:n
        if (A(i,j))
            hold on;
            plot([allGraphPts(i,1) allGraphPts(j,1)],[allGraphPts(i,2) allGraphPts(j,2)],'k');

        end
    end
end
toc;