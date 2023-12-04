clc;
clear;

for func_no=1:30    % No of Objective function from CEC2014
    Nf=func_no;
    
    for run_no=1:5
        tic
        Max_iter=1000; % Maximum numbef of iterations
        
        NFE=0;

        SearchAgents_no=30; % Number of search agents
        
        lb=-100;
        ub=100;
        dim=30;
        
        % initialize position vector and score for the leader
        Leader_pos=zeros(1,dim);
        Leader_score=inf; %change this to -inf for maximization problems
        
        
        %Initialize the positions of search agents
        Positions=initialization(SearchAgents_no,dim,ub,lb);
        
        Convergence_curve=zeros(1,Max_iter);
        fhd=str2func('cec14_func');%if you want to use form  Cost Functions CEC 2014
        
        func_num=Nf;
        iter=0;% Loop counter
        
        % Main loop
        while iter<Max_iter
            for i=1:size(Positions,1)
                
                % Return back the search agents that go beyond the boundaries of the search space
                Flag4ub=Positions(i,:)>ub;
                Flag4lb=Positions(i,:)<lb;
                Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                
                % Calculate objective function for each search agent
                fitness=feval(fhd,Positions(i,:)',func_num);%if you want to use form  Cost Functions CEC 2014
                NFE=NFE+1;
                % Update the leader
                if fitness<Leader_score % Change this to > for maximization problem
                    Leader_score=fitness; % Update alpha
                    Leader_pos=Positions(i,:);
                end
                
            end
            
            a=2-iter*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
            
            % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
            a2=-1+iter*((-1)/Max_iter);
            
            % Update the Position of search agents
            for i=1:size(Positions,1)
                r1=rand(); % r1 is a random number in [0,1]
                r2=rand(); % r2 is a random number in [0,1]
                
                A=2*a*r1-a;  % Eq. (2.3) in the paper
                C=2*r2;      % Eq. (2.4) in the paper
                
                
                b=1;               %  parameters in Eq. (2.5)
                l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
                
                p = rand();        % p in Eq. (2.6)
                
                for j=1:size(Positions,2)
                    
                    if p<0.5
                        if abs(A)>=1
                            rand_leader_index = floor(SearchAgents_no*rand()+1);
                            X_rand = Positions(rand_leader_index, :);
                            D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                            Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                            
                        elseif abs(A)<1
                            D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                            Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                        end
                        
                    elseif p>=0.5
                        
                        distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                        % Eq. (2.5)
                        Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                        
                    end
                    
                end
            end
            iter=iter+1;
            
            Convergence_curve(iter)=Leader_score;
        end
        
        run_times_WOA(run_no,func_no)=toc;

        % Results
        
        Cost_Rsult(1, run_no)=Leader_score;
        Rsult(run_no,:)= Convergence_curve;
    end
    
    Mean(func_no)=mean(Cost_Rsult);
    Best(func_no)=min(Cost_Rsult);
    Std(func_no)=std(Cost_Rsult);
    fprintf('Function #%2.0f,    Mean %7.4e,   Best %7.4e,   Std %7.4e\n',func_no,Mean(func_no),Best(func_no),Std(func_no))
end
