function [ps,pf]=CAMODE_SR(func_name,VRmin,VRmax,n_obj,NP,Max_Gen,N_rank)
%% 初始化参数
n_var=size(VRmin,2);    %Obtain the dimensions of decision space
Max_FES=Max_Gen*NP;     %Maximum fitness evaluations
CR=0.5;%交叉参数
count(200,3)=0;
%% Initialize population
VRmin=repmat(VRmin,NP,1);
VRmax=repmat(VRmax,NP,1);
pos=VRmin+(VRmax-VRmin).*rand(NP,n_var); %initialize the positions of the individuals
%% Evaluate the population
fitness=zeros(NP,n_obj);
for i=1:NP
    fitness(i,:)=feval(func_name,pos(i,:));
end
fitcount=NP;            % count the number of fitness evaluations
pop=[pos,fitness];
%% 初始化外部存档
Gobal_EXA=[];
Local_EXA=[];
[pop,gbest,~]=crowd_rank(pop,n_var,n_obj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pop2=pop;
gbest_index=gbest(:,n_var+n_obj+3);
pop2(gbest_index,:)=[];
Gobal_EXA=[Gobal_EXA;gbest];
EXA=gbest;
Local_EXA=[Local_EXA;pop2];
%%循环
for I=1:Max_Gen
    %% 种群聚类-- 动态聚类---选择
    [cluster,num]=AP1(pop,n_obj, n_var);
    %% 子种群进化
    tempop1=[];
    for ij=1:num
        subsize=size(cluster{ij,1},1);
        subpop=cluster{ij,1};
        Meta=[];
        trialfitness=[];
        if subsize<5 %%DE变异至少5个不同的个体
            %% 随机产生个体
            for subpop_i=1:subsize
                Meta(subpop_i,:)=VRmin(1,1:n_var)+(VRmax(1,1:n_var)-VRmin(1,1:n_var)).*rand(1,n_var);
            end
        else
            [subpop,Nbest,~]=crowd_rank(subpop(:,1:n_var+n_obj),n_var,n_obj);
            for subpop_i=1:subsize
                F=0.5;
                dx=randperm(subsize);
                dx(find(dx==subpop_i))=[];
                count(I,1)=count(I,1)+1;
                Meta(subpop_i,:)= subpop(dx(1),1:n_var)+F*(subpop(dx(2),1:n_var)-subpop(dx(3),1:n_var));
                sigma=tanh(I./Max_Gen);
                if rand<sigma
                    if subpop(subpop_i,n_var+n_obj+2)==max(subpop(:,n_var+n_obj+2))
                        count(I,2)=count(I,2)+1;
                        nbestpop=Gobal_EXA;
                        poptemp=repmat(subpop(subpop_i,1:n_var+n_obj),size(nbestpop,1),1);
                        distance_var=sqrt(sum((nbestpop(:,1:n_var)-poptemp(:,1:n_var)).*(nbestpop(:,1:n_var)-poptemp(:,1:n_var)),2));%%计算当前个体与种群中其他个体的欧式距离-distance from the j_t
                        [sortmiddlepop,index_sub]=sort(distance_var);
                        Guider=nbestpop(index_sub(end-4:end),:);
                        g_index=randperm(size(Guider,1),1);
                        guider=Guider(g_index,:);
                        Meta(subpop_i,:)=guider(1,1:n_var)+F*(subpop(dx(1),1:n_var)-subpop(dx(2),1:n_var));
                    end
                end
                Meta(subpop_i,:) = boundConstraint(Meta(subpop_i,:), subpop(subpop_i,:), [VRmin(1,:);VRmax(1,:)]);
            end
        end
        trial=subpop(:,1:n_var);
        %% 交叉操作
        for i=1:subsize
            r=randperm(n_var,1);
            for j=1:n_var
                if rand<=CR||j==r
                    trial(i,j)=Meta(i,j);
                end
            end
        end
        %% 评估种群
        for i=1:subsize
            trialfitness(i,:)=feval(func_name,trial(i,:));
            fitcount=fitcount+1;
        end
        trialpop=[trial,trialfitness];
        tempop1{ij,1}=[trialpop;subpop(:,1:n_var+n_obj)];
    end
    
    %% 非支配解排序---环境选择
    if I<0.75*Max_Gen
        GBA=[];
        tem_dele=[];
        for ij=1:num
            tempop=tempop1{ij,1};
            %% 在每个子种群中选择
            subpop1=tempop;
            [temp_pop,~,~]=Pareto_front_rank(subpop1(:,n_var+1:n_var+n_obj), n_obj);
            subpop1(:,n_var+n_obj+1)=temp_pop(:,end);      %%存储个体的排名值
            tempindex=find(subpop1(:,n_var+n_obj+1)==1);
            GBA{ij,1}=subpop1(tempindex,:);
            subpop1(tempindex,:)=[];
            tem_dele{ij,1}=subpop1(:,1:n_var+n_obj+1);
        end
        tempGBA=cell2mat(GBA);
        tempdele=cell2mat(tem_dele);
        %% 非支配解排序---环境选择
        if size(tempGBA,1)>NP
            tempGBA=non_domination_scd_sort(tempGBA(:,1:n_var+n_obj),n_obj,n_var); %%直接SCD不行
            tempGBA=sortrows(tempGBA,-(n_var+n_obj+3));  %%决策CD 直接SCD不行 全局和局部
            pop=tempGBA(1:NP,1:n_var+n_obj);             %%选择压力较大
            %             if N_rank>1
            %                 tempGBA=sortrows(tempGBA,-(n_var+n_obj+3));  %%决策CD 直接SCD不行 全局和局部
            %                 pop=tempGBA(1:NP,1:n_var+n_obj);             %%选择压力较大
            %
            %             else
            %                 [~,index2] = sortrows([tempGBA(:,n_var+n_obj+1) -tempGBA(:,n_var+n_obj+3)]); %%只有全局
            %                 tempGBA=tempGBA(index2,:);
            %                 pop=tempGBA(1:NP,1:n_var+n_obj);
            %             end
            %               [tempGBA2,Gbest,~]=crowd_rank(tempGBA,n_var,n_obj); %%直接SCD不行
            %               tempGBA3=sortrows(tempGBA2,(n_var+n_obj+2));%%决策CD   直接SCD不行
            % %               [~,index2] = sortrows([tempGBA2(:,n_var+n_obj+1) tempGBA2(:,n_var+n_obj+2)]);
            % %                tempGBA3=tempGBA2(index2,:);
            %                pop=tempGBA3(1:NP,1:n_var+n_obj);
        else
            remain=NP-size(tempGBA,1);
            [tempGBA1,Gbest,~]=crowd_rank(tempGBA,n_var,n_obj);          %%方案1还可以，但是比较复杂
            Gbest_index=Gbest(:,n_var+n_obj+3);
            tempGBA1(Gbest_index,:)=[];
            tempGBA=[Gbest;tempGBA1];
            
            data1=tempGBA1(:,1:n_var);
            RS1=0.2*prod(max(data1(:,1:n_var)-min(data1(:,1:n_var))).^(1./n_var));
            N1=size(data1,1);
            RS1=repmat(RS1,N1,1);
            
            data2=Gbest(:,1:n_var);
            RS2=0.2*prod(max(data2(:,1:n_var)-min(data2(:,1:n_var))).^(1./n_var));
            N2=size(data2,1);
            RS2=repmat(RS2,N2,1);
            Rs=[RS2; RS1];
            
            spop=[];
            n1=size(tempdele,1);
            number=1:1:n1;
            tempdele(:,n_var+n_obj+3)=number'; %%存储个体的索引
            for ii=1:size(tempGBA,1)
                popsort=tempdele;
                dist=zeros(size(popsort,1),1);
                dist(1:size(popsort,1),:)=sum((ones(size(popsort,1),1)*tempGBA(ii,1:n_var)-popsort(:,1:n_var)).^2,2)<Rs(ii);
                spop=[spop;popsort(dist==1,:)];
            end
            spop=unique(spop,'rows','stable');  %%过滤掉相同的个体
            %% 计算密度
            subsize=size(spop,1);
            if subsize>remain
                [spop,Gbest1,~]=crowd_rank(spop,n_var,n_obj); %%方案2还可以，与方案1差不多
                spop=sortrows(spop,(n_var+n_obj+2));          %%目标和决策
                pop=[tempGBA(:,1:n_var+n_obj);spop(1:remain,1:n_var+n_obj)];
            else
                indexspop=spop(:,n_var+n_obj+3);
                tempdele(indexspop,:)=[];
                tepop=[tempGBA(:,1:n_var+n_obj);spop(:,1:n_var+n_obj)];
                remaining=NP-size(tepop,1);
                tempdele=non_domination_scd_sort(tempdele(:,1:n_var+n_obj),n_obj,n_var);
                tempdele=sortrows(tempdele,-(n_var+n_obj+3)); %%决策SCD
                newpop2=tempdele(1:remaining,:);
                pop=[tepop;newpop2(:,1:n_var+n_obj)];
            end
        end
    else
        for ij=1:num
            tempop=tempop1{ij,1};
            parent=cluster{ij,1};
            subsize=size(parent,1);
            tempop=non_domination_scd_sort(tempop(:,1:n_var+n_obj),n_obj,n_var);
            cluster{ij,1}= tempop(1:subsize,1:n_var+n_obj);   %%SCD
        end
        pop=cell2mat(cluster);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    end
    
    [pop,best,~]=crowd_rank(pop,n_var,n_obj);
    pop2=pop;
    best_index1=best(:,n_var+n_obj+3);
    pop2(best_index1,:)=[];
    %% UPdate Gobal_EXA
    tempGobal_EXA=[Gobal_EXA(:,1:n_var+n_obj);best(:,1:n_var+n_obj)];  %%改了这个地方
    tempGobal_EXA=unique(tempGobal_EXA,'rows','stable');%%过滤掉相同的个体
    tempGobal_EXA=non_domination_scd_sort(tempGobal_EXA,n_obj,n_var);
    tempGobal_EXA1=tempGobal_EXA(tempGobal_EXA(:,n_var+n_obj+1)==1,1:n_var+n_obj);
    
    if size(tempGobal_EXA1,1)>NP
        Gobal_EXA=tempGobal_EXA1(1:NP,1:n_var+n_obj);
    else
        Gobal_EXA=tempGobal_EXA1;
    end
    
    
    %              I
    % %       figure(9)
    % %     %     %     % %     if n_var==2
    % %                              figure(I)
    %      plot(pop(:,1),pop(:,2),'ro');
    % %                         figure(I+50)
    %                          plot(pop(:,3),pop(:,4),'bo');
    %     %     else
    %     %         %         figure(I)
    %     %         plot3(pop(:,1),pop(:,2),pop(:,3),'ro');
    %     %         figure(I+50)
    %     %         plot(pop(:,4),pop(:,5),'ro');
    %     %     end
    %     pause(0.01)
    if fitcount>Max_FES
        break;
    end
end
%% Output
ps = pop(:,1:n_var);
pf = pop(:,1+n_var:n_var+n_obj);
