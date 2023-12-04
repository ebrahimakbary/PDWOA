clc;
clear;

for func_no=1:30    % No of Objective function from CEC2014
    Nf=func_no;

    for run_no=1:5
        tic
        Max_iter=1000; % Maximum number of iterations

        NFE=0;
        
        Crr=0.8;
        SearchAgents_no=30; % Number of search agents
        nPop=SearchAgents_no;
        
        lb=-100;
        ub=100;
        dim=30;
        
        % initialize position vector and score for the leader
        Leader_pos=zeros(1,dim);
        Leader_score=inf; %change this to -inf for maximization problems
        
        
        %Initialize the positions of search agents
        Pbest_pos=initialization(SearchAgents_no,dim,ub,lb);
        Positions1=initialization(SearchAgents_no,dim,ub,lb);
        Convergence_curve=zeros(1,Max_iter);
        fhd=str2func('cec14_func'); %if you want to use from Cost Functions CEC 2014
        
        func_num=Nf;
        Pbest_fit=zeros(1,nPop);
        fitness1=Pbest_fit;
        for i=1:size(Pbest_pos,1)
            Pbest_fit(i)=feval(fhd,Pbest_pos(i,:)',func_num);%if you want to use form  Cost Functions CEC 2014
            NFE=NFE+1;
            fitness1(i)=inf;
        end
        % Main loop
        iter=0;% Loop counter
        while iter<Max_iter
            for i=1:size(Pbest_pos,1)
                
                % Return back the search agents that go beyond the boundaries of the search space
                Flag4ub=Positions1(i,:)>ub;
                Flag4lb=Positions1(i,:)<lb;
                Positions1(i,:)=(Positions1(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                
                % Calculate objective function for each search agent
                %         fitness=fobj(Positions(i,:));
                fitness1(i)=feval(fhd,Positions1(i,:)',func_num);%if you want to use form  Cost Functions CEC 2014
                NFE=NFE+1;
                
                if fitness1(i)<Pbest_fit(i)
                    Pbest_pos(i,:)=Positions1(i,:);
                    Pbest_fit(i) =fitness1(i) ;
                end
                % Update the leader
                if Pbest_fit(i)<Leader_score % Change this to > for maximization problem
                    Leader_score=Pbest_fit(i); % Update alpha
                    Leader_pos=Pbest_pos(i,:);
                end
                
            end
            
            a=2-iter*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
            
            % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
            a2=-1+iter*((-1)/Max_iter);
            
            % Update the Position of search agents
            for i=1:size(Pbest_pos,1)
                r1=rand(); % r1 is a random number in [0,1]
                r2=rand(); % r2 is a random number in [0,1]
                
                A=2*a*r1-a;  % Eq. (2.3) in the paper
                C=2*r2;      % Eq. (2.4) in the paper

                b=1;               %  parameters in Eq. (2.5)
                l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)

                p = rand();        % p in Eq. (2.6)

                for j=1:size(Pbest_pos,2)
                    
                    if p<0.5
                        if abs(A)>=1
                            rand_leader_index = floor(SearchAgents_no*rand()+1);
                            X_rand = Pbest_pos(rand_leader_index, :);
                            D_X_rand=abs(C*X_rand(j)-Pbest_pos(i,j)); % Eq. (2.7)
                            Positions1(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)

                        elseif abs(A)<1
                            D_Leader=abs(C*Leader_pos(j)-Pbest_pos(i,j)); % Eq. (2.1)
                            Positions1(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                        end
                        
                    elseif p>=0.5
                        
                        distance2Leader=abs(Leader_pos(j)-Pbest_pos(i,j));
                        % Eq. (2.5)
                        Positions1(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                    end
                end
            end
            M1=randperm(nPop,3);
            M1(M1==i)=[];
            m1=M1(1);
            m2=M1(2);
            for jj=1:size(Pbest_pos,2)
                if rand>Crr
                    Positions1(i,jj)= Pbest_pos(i,jj)+rand*(Leader_pos(jj)-Pbest_pos(i,jj))+rand*(Pbest_pos(m1,jj)-Pbest_pos(m2,jj));
                end
            end
            iter=iter+1;
            
            Convergence_curve(iter)=Leader_score-(func_num*100);
        end
        
        run_times_PDWOA(run_no,func_no)=toc;
        
        % Results
        
        Cost_Rsult(1, run_no)=Leader_score-(func_num*100);
        Rsult(run_no,:)= Convergence_curve;
    end
    
    Mean(func_no)=mean(Cost_Rsult);
    Best(func_no)=min(Cost_Rsult);
    Std(func_no)=std(Cost_Rsult);
    fprintf('Function #%2.0f,    Mean %7.4e,   Best %7.4e,   Std %7.4e\n',func_no,Mean(func_no),Best(func_no),Std(func_no))
    
end
