 function [C,num]= AP1(tempop,n_obj, n_var)
    C=[];
    index=[];
    x=tempop(:,1:n_var);
    nsize=size(tempop,1);
%     %% 存储个体的索引
%     index=(1:1:nsize);
%     tempop(:,n_var+n_obj+1)=index';
    %% 种群聚类
    % generate similarity matri
    m     = size(x,1);
    s     = zeros(m);
    
    o     = nchoosek(1:nsize,2);      % set up all possible pairwise comparisons
    xx    = x(o(:,1),:)';           % point 1
    xy    = x(o(:,2),:)';           % point 2
    d2    = (xx - xy).^2;           % distance squared
    d     = -sqrt(sum(d2));         % distance
    
    k     = sub2ind([m m], o(:,1), o(:,2) );    % prepare to make square
    s(k)  = d;
    s     = s + s';
    di = 1:(m+1):m*m;         %index to diagonal elements
    
      s(di) = min(d);
%     s(di) = mean(d);
%     s(di) = median(d);
    %% clustering
    options.StallIter = 10;
    % options.OutputFcn = @(a,r) affprop_plot(a,r,x,'k.');
    %
    % figure
    ex = affprop(s, options );
    u=unique(ex );
    for k = 1:length(u)
        t=ex==u(k);
%         C{k,1}=tempop(t,1:n_var+n_obj);
        C{k,1}=tempop(t,:);
    end
    num=length(C);
 end