

% load data
data = readmatrix('WarehouseExample.txt');

dims = size(data);
m = dims(1);
n = dims(2);
depot = [0,50];

xMat = data(:,1:2:end);
yMat = data(:,2:2:end);

%VARIABLES
mu = 100; % Number of parents per generation
lambda = 300; % Number of children per generation
comProb = 1; % Prob of combination
mutProb = 1; % Prob of mutation
noNewCount = 20; % Number of generations without children until end of program



%%% INIT:
totalDists = zeros(1,mu); % List of distances
items_order = zeros(mu,m); % Order of items selected
piles = zeros(mu,m); % Piles the items are selected from

for k = 1:mu
    items_order(k,:) = randsample(1:m,m); % randomly generate initial order of items to collect
    piles(k,:) = randi([1,5],1,m); % randomly generate the initial olocations to visit for each item
    index = sub2ind(size(xMat),items_order(k,:),piles(k,:));
    x_pos = [depot(1),xMat(index),depot(1)];
    y_pos = [depot(2),yMat(index),depot(2)];
    
    for i = 2:length(x_pos)
        dist = sqrt((x_pos(i-1)- x_pos(i))^2 + (y_pos(i-1) - y_pos(i))^2); 
        totalDists(k) = totalDists(k) + dist; 
    end
end

parent_io_mat = items_order;
parent_p_mat = piles;


%%% MAIN LOOP:

while noNewCount > 0
    
    % Create distribution for selecting the parents
    [~, ~, ranking] = unique(totalDists);
    weights = ranking/sum(ranking);
    
    % child item order and pile matrices
    child_io_mat = repelem(0,lambda,m);
    child_p_mat = repelem(0,lambda,m);
    totalDists_child = repelem(0,1,lambda);
    
    % Select 2*lambda parents randomly and create child lambda times
    parents = repelem(0,lambda,2);
    
    for i = 1:lambda
        % Select 2 parents randomly (weighted dist)
        parentsRank = randsample(ranking,2,true,weights); 
        parents(i,:) = parentsRank; % Crate list of parents used for each child
        parents_io = parent_io_mat(parentsRank,:); % Get item order of parents
        parents_p = parent_p_mat(parentsRank,:); % Get selected piles for parents
        
        % COMBINATION
        
        % select 2 points to form a subset to cut from parent 1
        cut = datasample(1:m,2,"Replace",false);
        minCut = min(cut);
        maxCut = max(cut);
        % part of child from parent 1
        child_io_p1 = parents_io(1,minCut:maxCut);
        % part of child from parent 2 (all items not taken by parent 1)
        child_io_p2 = parents_io(2,boolean(abs(ismember(parents_io(2,:),child_io_p1)-1)));
        
        % combine the child parts into the new child
        if minCut == 1 && maxCut == m
            child_io = child_io_p1;
            child_p = parents_p(1,:);
        elseif minCut == 1
            child_io = [child_io_p1,child_io_p2];
            child_p = [parents_p(1,minCut:maxCut),parents_p(2,maxCut+1:end)];
        elseif maxCut == m
            child_io = [child_io_p2,child_io_p1];
            child_p = [parents_p(2,1:minCut-1),parents_p(1,minCut:maxCut)];
        else
            child_io = [child_io_p2(1:minCut),child_io_p1,child_io_p2(minCut+1:end)];
            child_p = [parents_p(2,1:minCut-1),parents_p(1,minCut:maxCut),parents_p(2,maxCut+1:end)];
        end
        
            
        % MUTATION 1
        
        mut1 = true;
        denom = 1;
        % perform a random number of random mutations
        while mut1 == true
            prob = rand(1);
            
            if prob < mutProb/denom % prob of new mutation decreases for every mutation
                % switch two random elements around
                mutation = randsample(1:m,2); % choose 2 random points to switch
                mut_io_1 = child_io(mutation(1)); mut_p_1 = child_p(mutation(1));
                mut_io_2 = child_io(mutation(2)); mut_p_2 = child_p(mutation(2));
                child_io(mutation(1)) = mut_io_2; child_p(mutation(1)) = mut_p_2;
                child_io(mutation(2)) = mut_io_1; child_p(mutation(2)) = mut_p_1;
                denom = denom + 0.5;
            else
                mut1 = false;
            end
        end
        
        
        % MUTATION 2
        
        mut2 = true;
        denom = 1;
        % perform a random number of random mutations
        while mut2 == true
            prob = rand(1);
            
            if prob < mutProb/denom % prob of new mutation decreases for every mutation
                % switch two random elements around
                mutPos = randi([1,m],1,1);
                mutValue = randi([1,n/2],1,1); % choose 2 random points to switch
                child_p(mutPos) = mutValue;
                denom = denom + 1;
            else
                mut2 = false;
            end
        end
        
        % Find x and y value of child at each point and add depot
        index = sub2ind(size(xMat),child_io,child_p);
        % Add depot coords to start and end of child
        x_pos = [depot(1),xMat(index),depot(1)];
        y_pos = [depot(2),yMat(index),depot(2)];
        
        % EVALUATE NEW CHILD
        for j = 2:(m+2)
            dist = sqrt((x_pos(j-1)- x_pos(j))^2 + (y_pos(j-1) - y_pos(j))^2); 
            totalDists_child(i) = totalDists_child(i) + dist; 
        end
        
        % Add new child to the matrix
        child_io_mat(i,:) = child_io;
        child_p_mat(i,:) = child_p;
    end
    
    % Combine best parents and best children to form new generation mu+lambda
    totalDists = [totalDists,totalDists_child]; %#ok<AGROW>
    [totalDists,minIndex] = mink(totalDists,mu);
    parent_minIndex = minIndex(minIndex <= mu);
    child_minIndex = minIndex(minIndex > mu) - mu;
    parent_io_mat = [parent_io_mat(parent_minIndex,:);child_io_mat(child_minIndex,:)];
    parent_p_mat = [parent_p_mat(parent_minIndex,:);child_p_mat(child_minIndex,:)];
    
    % End program when no new children paths are surviving
    new = mu - length(parent_minIndex);
    if new == 0
        noNewCount = noNewCount - 1;
        disp(noNewCount)
    end
end


% find the optimal solution
[minDist,minInd] = min(totalDists);
items_order = parent_io_mat(minInd,:);
piles = parent_p_mat(minInd,:);
index = sub2ind(size(xMat),items_order,piles);
x_pos = [depot(1),xMat(index),depot(1)];
y_pos = [depot(2),yMat(index),depot(2)];

%%% PLOT SOLUTION:
title('Evolutionnary Path')
hold on

for i = 1:m
    scatter(xMat(i,:),yMat(i,:),"bla",".") % plot each row of the data as a scatter plot
    b = num2str(i); c = cellstr(b); % change number i to string
    dx = 0.5; dy = 0.5; % displacement so the text does not overlap with the data points
    text(xMat(i,:)+dx, yMat(i,:)+dy, c); % Assign a different number for each item
end

plot(x_pos,y_pos,"r") % plot the path
plot(depot(1),depot(2),"*","color","blue") % put a point for the depot
hold off

% Display the optimal distance
disp("The length of the best found path is:")
disp(minDist)
