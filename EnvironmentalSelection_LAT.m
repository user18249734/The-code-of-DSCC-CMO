function [return_pop,return_Fitness] = EnvironmentalSelection_LAT(Population,N,VAR)
% Multi-objective-based CHT is used to sort the Population of local
% auxiliary task

input_cons = Population.cons;
input_cons(input_cons<0) = 0;
input_cons = sum(input_cons,2);

findex = find(input_cons<=VAR);
fPopulation = Population(findex);

if isempty(fPopulation)
    fPopulation = [];
    fFitness = [];
elseif length(fPopulation) <= N
    cons = fPopulation.cons;
    cons(cons<0)=0;
    cons = sum(cons,2);
    fFitness = CalFitness([fPopulation.objs,cons]);

    % Sort the population
    [fFitness,rank] = sort(fFitness);
    fPopulation = fPopulation(rank);
    fFitness = fFitness(rank);
elseif length(fPopulation) > N
    cons = fPopulation.cons;
    cons(cons<0)=0;
    cons = sum(cons,2);
    fFitness = CalFitness([fPopulation.objs,cons]);
    Next = fFitness < 1;
    if sum(Next) <= N
        [~,Rank] = sort(fFitness);
        Next(Rank(1:N )) = true;
    elseif sum(Next) > N
        Del  = Truncation(fPopulation(Next).objs, sum(Next)-N );
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end

    fPopulation = fPopulation(Next);
    fFitness   = fFitness(Next);
    % Sort the population
    [fFitness,rank] = sort(fFitness);
    fPopulation = fPopulation(rank);

end

return_pop = [fPopulation];
return_Fitness = [fFitness];
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

%% Truncation
Distance = pdist2(PopObj,PopObj);
Distance(logical(eye(length(Distance)))) = inf;
Del = false(1,size(PopObj,1));
while sum(Del) < K
    Remain   = find(~Del);
    if isempty(Remain)
        keyboard
    end
    Temp     = sort(Distance(Remain,Remain),2);
    [~,Rank] = sortrows(Temp);
    Del(Remain(Rank(1))) = true;
end
end