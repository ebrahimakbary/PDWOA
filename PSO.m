clc;
clear;

for func_no=1:30    % No of Objective function from CEC2014
    Nf=func_no;
    
    for run_no=1:5
        tic
        Max_iter=1000; % Maximum number of iterations
        
        NFE=0;
        
        VarMin=-100;             % Decision Variables Lower Bound
        VarMax= -VarMin;             % Decision Variables Upper Bound
        
        npop=30;
        nvar=30;
        
        w=1;
        wdamp=0.99;
        
        c1=2;
        c2=2;
        
        xmin=VarMin;
        xmax=VarMax;
        dx=xmax-xmin;
        
        vmax=0.1*dx;
        
        fhd=str2func('cec14_func');%if you want to use from Cost Functions CEC 2014
        
        func_num=Nf;

        MM=zeros(Max_iter,1);
        iter=0;
        while iter<Max_iter
            iter=iter+1;
            if iter==1
                gbestcost=inf;
                for i=1:npop
                    velocity(i,:)=zeros(1,nvar);
                    position(i,:)=xmin+(xmax-xmin)*rand(1,nvar);
                    %          cost(i,:)=Cost(position(i,:), Ne);
                    cost(i,:) =feval(fhd,position(i,:)',func_num);%if you want to use form  Cost Functions CEC 2014
                    NFE=NFE+1;
                    
                    pbest(i,:)=position(i,:);
                    pbestcost(i,:)=cost(i,:);
                    
                    if pbestcost(i,:)<gbestcost
                        gbest=pbest(i,:);
                        gbestcost=pbestcost(i,:);
                    end
                end
                MM(iter)=gbestcost-(func_num*100);
            else
                for i=1:npop
                    velocity(i,:)=w*velocity(i,:)...
                        +c1*rand*(pbest(i,:)-position(i,:))...
                        +c2*rand*(gbest-position(i,:));
                    
                    velocity(i,:)=min(max(velocity(i,:),-vmax),vmax);
                    
                    position(i,:)=position(i,:)+velocity(i,:);
                    
                    position(i,:)=min(max(position(i,:),xmin),xmax);
                    
                    %            cost(i,:)=Cost(position(i,:), Ne);
                    cost(i,:) =feval(fhd,position(i,:)',func_num);%if you want to use form  Cost Functions CEC 2014
                    NFE=NFE+1;
                    if cost(i,:)<pbestcost(i,:)
                        pbest(i,:)=position(i,:);
                        pbestcost(i,:)=cost(i,:);
                        
                        if pbestcost(i,:)<gbestcost
                            gbest=pbest(i,:);
                            gbestcost=pbestcost(i,:);
                        end
                    end
                end
            end
            MM(iter)=gbestcost-(func_num*100);
            w=w*wdamp;
        end
        run_times_PSO(run_no,func_no)=toc;
        
        Rsult(run_no, :)=MM;
        Cost_Rsult(1, run_no)=gbestcost-(func_num*100);
    end
    
    Mean(func_no)=mean(Cost_Rsult);
    Best(func_no)=min(Cost_Rsult);
    Std(func_no)=std(Cost_Rsult);
    fprintf('Function #%2.0f,    Mean %7.4e,   Best %7.4e,   Std %7.4e\n',func_no,Mean(func_no),Best(func_no),Std(func_no))
    
end
