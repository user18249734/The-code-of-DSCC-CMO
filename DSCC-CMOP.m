classdef DSCC-CMOP < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% local auxiliary task-based constrained multiobjective optimization

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)

            Population1 = Problem.Initialization(); % Main task
            Fitness1    = CalFitness(Population1.objs,Population1.cons);
            Population2 = Problem.Initialization(); % Global auxiliary task
            Fitness2   = CalFitness(Population2.objs);
            Population4 = Problem.Initialization(); % Local auxiliary task
            Fitness4   = CalFitness(Population4.objs,Population4.cons);
            

            cons = Population1.cons;
            cons(cons<0) = 0;
            cons =sum(cons,2);
            index =find(cons>0);
            if isempty(index)
                VAR0 = 0;
            else
                VAR0 =  mean(cons(index));
            end

            cnt = 0; 
            flag = 0; 

           %% Generate the weight vectors  
            [RW,WVs] = UniformPoint(Problem.N,Problem.M); 
            T = ceil(WVs/10);  
            nr = ceil(WVs/100);  
            W = zeros(size(RW));
            for i = 1:WVs
                W(i, :) = RW(WVs - i + 1, :);
            end

            %% Detect the neighbours of each solution  
            B = pdist2(W,W); 
            [~,B] = sort(B,2);
            B = B(:,1:T); 
 
            %% Angle  
            angle=acos(1-pdist2(W,W,'cosine'));
            temp_angle=angle;
            temp_angle(logical(eye(size(temp_angle))))=inf;
            theta_min=min(temp_angle');
            theta_min=theta_min';
            theta=theta_min.*0.5;
            Population = [Population1,Population2]
            arch = ArchiveUpdate(Population,Problem.N);  
            Population = arch; 
            Fes = Problem.FE;           
            while Algorithm.NotTerminated(Population1)
                cnt =cnt +1; 
                if flag == 0
                    std_obj(cnt,:) = std(Population2.objs,[],1); 
                    if cnt>100  
                        if  sum(std(std_obj(cnt-100:cnt,:),[],1)<0.5) == Problem.M 
                            flag = 1;
                        end
                    end
                end


                for i = 1:Problem.N
                    ii = i;  
                    if i > WVs
                        ii = randi([1,WVs]); 
                    end                   
                    
                    if rand < 0.2  % Perform global selection with a probability of 0.2.
                        PP{i} = randperm(WVs);
                        P = PP{i};
                        E(i,:)=P(1:2); 
                    else     %  Perform neighborhood selection with a probability of 0.8.
                        PP{i} = B(ii,randperm(size(B,2)));
                        P = PP{i};
                        E(i,:)=P(1:2);
                    end
                end
                
                MatingPool = TournamentSelection(2,Problem.N,Fitness1);
                Offspring1 = OperatorGAhalf(Problem,[Population1(MatingPool)]); 


                if flag == 0  
                    Offspring2 = OperatorGAhalf(Problem,[Population2(E(:,1)),Population2(E(:,2))]); 
                else          
                    Offspring2 = OperatorDE(Problem,Population2,Population2(E(:,1)),Population2(E(:,2)));
                end

               Population3 = Offspring2;  
               % for each solution 
               Z = min(Population2.objs,[],1);
               for i = 1 : Problem.N
                   P = PP{i};
                   Z = min(Z,Population3(i).obj);
                   PopObj2=Population2(P).objs-repmat(Z,length(P),1);
                   Angle0 = acos(1 - pdist2(real(PopObj2),W(P,:),'cosine'));
                   Angle2 = diag(Angle0);
                   NewAngle2=acos(1-pdist2(real(Population3(i).objs-Z),W(P,:),'cosine'));
                   NewAngle2=NewAngle2';
                   
                   %
                   g_old2 = max(abs(Population2(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                   g_new2 = max(repmat(abs(Population3(i).obj-Z),length(P),1).*W(P,:),[],2);
                   
                   if (flag == 1)    %second stage         
                       CVO2 = sum(max(0,Population3(i).con));
                       CVP2 = sum(max(0,Population2(P).cons),2);
                       case1=NewAngle2<= theta(P) & Angle2<= theta(P) & (CVP2>CVO2 | (CVP2==CVO2 & g_old2>=g_new2));%both in niche
                       case2=NewAngle2<= theta(P) & Angle2> theta(P) & (CVP2>CVO2 | (CVP2==CVO2 & g_old2>=g_new2));
                       case3=NewAngle2<= theta(P) & Angle2> theta(P);
                       indices = mod(find([case1 ; case2 ; case3],nr),length(P));
                       if indices > 0
                           Population2(P(indices)) = Population3(i);
                       end
                   end
               end

                if length(Population4) <=1   
                    Offspring4 = [];
                else
                    MatingPool = TournamentSelection(2,min(length(Population4),Problem.N/2),Fitness4);
                    Offspring4 = OperatorGA(Problem,[Population4(MatingPool)]);
                end

                %% 

                if flag == 0  %first stage
                    [Population1,~] = EnvironmentalSelection([Population1,Offspring2,Offspring4],Problem.N,true);  
                    [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1],Problem.N,true);  
                    [Population2,~] = EnvironmentalSelection([Population2,Offspring1,Offspring2,Offspring4],Problem.N,false);
                else  %second stage
                    [Population1,~] = EnvironmentalSelection([Population1,Population2,Offspring4],Problem.N,true);  
                    [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1],Problem.N,true); 
                end

                [Population4,Fitness4] = EnvironmentalSelection_LAT([Population4,Offspring1,Offspring2,Offspring4],Problem.N,VAR0);

                
                cons = Offspring1.cons;
                cons(cons<0) = 0;
                cons = sum(cons,2);
                index = find(cons>0);
                if isempty(index)
                    VAR0 = 0;
                else
                    tanhFE = 8 * ((Problem.FE-Fes) / (Problem.maxFE-Fes)) - 4;
                    K = 0.5-tanh(tanhFE)/2;
                    VAR0 = K*mean(cons(index));
                end
            end
        end
    end
end
