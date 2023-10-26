function [subpop,GBA,distance_var]=crowd_rank(subpop,n_var,n_obj)
%% 计算距离
Dis=[];
distance_var=[];
distance_obj=[];
subsize=size(subpop,1);
Crowd=[];
small_dis=[];
% N_ns=5-ceil(gen/Max_Gen.*6)+1; %
N_ns=5;
number=1:1:subsize;
for subpop_i=1:subsize
    poptemp=repmat(subpop(subpop_i,1:n_var+n_obj),subsize,1); %%当前个体与整个种群的距离
    distance_var{subpop_i,1}=sqrt(sum((subpop(:,1:n_var)-poptemp(:,1:n_var)).*(subpop(:,1:n_var)-poptemp(:,1:n_var)),2));%%计算当前个体与种群中其他个体的欧式距离-distance from the j_th individual to pop in ***decision space* [bb1,cc1]=sort(distance_var);
    distance_obj{subpop_i,1}=sqrt(sum((subpop(:,n_var+1:n_obj+n_var)-poptemp(:,n_var+1:n_obj+n_var)).*(subpop(:,n_var+1:n_obj+n_var)-poptemp(:,n_var+1:n_obj+n_var)),2));%%计算当前个体与种群中其他个体的欧式距离-distance from the j_th individual to pop in ***decision space* [bb1,cc1]=sort(distance_var);
    distance_var{subpop_i,1}=distance_var{subpop_i,1}./max(distance_var{subpop_i,1}); %%归一化决策空间距离
    distance_obj{subpop_i,1}=distance_obj{subpop_i,1}./max(distance_obj{subpop_i,1}); %%归一化目标空间距离
    %%排序个体，获取近邻个体的距离及个体
    [sorted_dis1,cc1]=sort(distance_var{subpop_i,1}); %%排序归一化的距离
    [sorted_dis2,cc2]=sort(distance_obj{subpop_i,1}); %%排序归一化的距离
    tem_small_dis1=sorted_dis1(1:N_ns); %%选取前N-ns个个体，构成个体的领域
    tem_small_dis2=sorted_dis2(1:N_ns); %%选取前N-ns个个体，构成个体的领域
    Crowd_var(subpop_i,1)=sum(tem_small_dis1); %%拥挤指标
    Crowd_obj(subpop_i,1)=sum(tem_small_dis2);
    %%解被选取的半径
    small_dis=[small_dis;tem_small_dis1];
end
Dis=mean(small_dis);%%选择过程中解被选择的半径--决策空间 记得是正交化之后的解的距离
CD=[];%%越小的CD值越好
Mvar=sum(Crowd_var)./subsize;
Mobj=sum(Crowd_obj)./subsize;
for subpop_i=1:subsize
    Crowd_v(subpop_i,1)=Crowd_var(subpop_i,1)./Mvar; %%拥挤指标
    Crowd_o(subpop_i,1)= Crowd_obj(subpop_i,1)./Mobj;
    CD(subpop_i,1)=1./(1+Crowd_v(subpop_i,1)+Crowd_o(subpop_i,1));
end
subpop(:,n_var+n_obj+2)=CD;  %%存储个体的拥挤度---可以进行改进
%%环境选择过程中解被选取的半径

%% 获取当前个体的排名在子种群中
[temp_subpop_rank,~,~]=Pareto_front_rank(subpop(:,n_var+1:n_var+n_obj),n_obj); %%对子种群进行非支配排名
subpop(:,n_var+n_obj+1)=temp_subpop_rank(:,n_obj+1); %%存储个体的排名值
subpop(:,n_var+n_obj+3)=number'; %%存储个体的索引
tempindex=find(subpop(:,n_var+n_obj+1)==1);
GBA=subpop(tempindex,:);
end
