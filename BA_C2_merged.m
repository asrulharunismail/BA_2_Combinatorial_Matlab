clc;
clear;
close all;
tic
%% Problem Definition
[typeOfFunction] = 'Eil51'; 
Instance=Tsplib(typeOfFunction);
Dims=Instance.dim;
ObjFunction=@(x) Instance.evaluation( x );             % Objective Function
VarSize=[1 Dims];                                     

%% Bees Algorithm Parameters
MaxEval = 1000000;
n=30;                               % Number of Scout Bees
nep=30;                             
MaxIt=2000;
%% Initialization
Empty_Bees.Position=[];
Empty_Bees.Cost=[];
Empty_Bees.counter=[];
Bees=repmat(Empty_Bees,n,1);
counter=0;
% Generate Initial Solutions
for i=1:n
    Bees(i).Position=randperm(Dims);
    Bees(i).Cost=ObjFunction(Bees(i).Position);
    counter=counter+1;
    Bees(i).counter= counter;
end
recruitment = round(linspace(nep,1,n));
size = linspace(0,1,n);
ColonySize=sum(recruitment);        % total number of foragers

%% Sites Selection 
[~, RankOrder]=sort([Bees.Cost]);
Bees=Bees(RankOrder);

%% Bees Algorithm Local and Global Search
for it=1:MaxIt
    
    if counter >= MaxEval
        break;
    end

    % All Sites (Exploitation and Exploration)
    for i=1:n

        bestnewbee.Cost=inf;

        assigntment= linspace(0,1,recruitment(i));
        

        for j=1:recruitment(i)
            newbee.Position= Foraging_Combi(Bees(i).Position,assigntment(j)* Dims);
            newbee.Cost=ObjFunction(newbee.Position);
            counter=counter+1;
            newbee.counter= counter;
            if newbee.Cost<bestnewbee.Cost
                bestnewbee=newbee;
            end
        end

        if bestnewbee.Cost<Bees(i).Cost
            Bees(i)=bestnewbee;
        end

    end

    % SORTING
    [~, RankOrder]=sort([Bees.Cost]);
    Bees=Bees(RankOrder);

    % Update Best Solution Ever Found
    OptSol=Bees(1);

    % taking of result
    OptCost(it)=OptSol.Cost;
    Counter(it)=counter;
    Time(it)=toc;
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(OptCost(it)) ' --> Time = ' num2str(Time(it)) ' seconds' '; Fittness Evaluations = ' num2str(Counter(it))]);
    
    %figure(1);
    %PlotSolution(OptSol.Position,Instance);
end


%% Results
%figure;
%semilogy(OptCost,'LineWidth',2);
%xlabel('Iteration');
%ylabel('Best Cost');


function qnew=Foraging_Combi(q,sz)
    nVar=numel(q);
    ngh= ceil(triangular(0,sz,1)*nVar);
    m=randi([1 4]);
    
    switch m
        case 1
            % Do Swap
            qnew=Swap(q, ngh);
            
        case 2
            % Do Reversion
            qnew=Reversion(q, ngh);

        case 3
            % Do Insertion
            qnew=Insertion(q, ngh);

        case 4
            % Do Reversion
            qnew=Reversion(q, ngh);
    end

end

function qnew=Swap(q, ngh)

    n=numel(q);
    i1=randi([1 n]);
    i2=i1 + randi([1 ngh]);
    i2(i2>n)=n;
    qnew=q;
    qnew([i1 i2])=q([i2 i1]);
    
end

function qnew=Reversion(q, ngh)

    n=numel(q);
    
    i1=randi([1 n]);
    i2=i1 + randi([1 ngh]);
    i2(i2>n)=n;
    
    qnew=q;
    qnew(i1:i2)=q(i2:-1:i1);

end

function qnew=Insertion(q, ngh)

    n=numel(q);
    a=randi(2);
    switch a
        case 1
            i1=randi([1 n]);
        case 2
            i1=randi([1 n-1]);
            i1=[i1 i1+1];
    end
    i2=i1(end) + randi([-ngh ngh]);
    i2(i2>n)=n;
    i2(i2<1)=1;
    
    
    if i1<i2
        qnew=[q(1:i1-1) q(i1+1:i2) q(i1) q(i2+1:end)];
    else
        qnew=[q(1:i2) q(i1) q(i2+1:i1-1) q(i1+1:end)];
    end

end


function random_number = triangular(min, mode, max)
    r = rand();
    prop_mode = (mode - min) / (max - min);

    if r <= prop_mode
        random_number = min + sqrt(r * (max - min) * (mode - min));
    else
        random_number = max - sqrt((1 - r) * (max - min) * (max - mode));
    end
end

function PlotSolution(tour,model)
    tour=[tour tour(1)];
    plot(model.x(tour),model.y(tour),'k-o',...
    'MarkerSize',2.5,...
    'MarkerFaceColor','y',...
    'LineWidth',.5);
    xlabel('x');
    ylabel('y');
    axis equal;
    grid off;
    alpha = 0.1;
    xmin = min(model.x);
    xmax = max(model.x);
    dx = xmax - xmin;
    xmin = floor((xmin - alpha*dx)/10)*10;
    xmax = ceil((xmax + alpha*dx)/10)*10;
    xlim([xmin xmax]);
    ymin = min(model.y);
    ymax = max(model.y);
    dy = ymax - ymin;
    ymin = floor((ymin - alpha*dy)/10)*10;
    ymax = ceil((ymax + alpha*dy)/10)*10;
    ylim([ymin ymax]);
end


function L=TourLength(tour,model)

    n=numel(tour);

    tour=[tour tour(1)];
    
    L=0;
    for i=1:n
        L=L+model.D(tour(i),tour(i+1));
    end

end

