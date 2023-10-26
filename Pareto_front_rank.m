function [x,index_of_fronts,sorted_based_on_front,F]=Pareto_front_rank(x,n_obj)
% Input  x      num*n_obj
%        n_obj
% Output sorted_based_on_front

N=size(x,1);%
front = 1;
F(front).f = [];
individual = [];
for i = 1 : N
    % Number of individuals that dominate this individual
    individual(i).n = 0;
    % Individuals which this individual dominate
    individual(i).p = [];
    for j = 1 : N
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : n_obj
            if (x(i,k) < x(j,k))
                dom_less = dom_less + 1;
            elseif (x(i,k) == x(j,k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= n_obj
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= n_obj
            individual(i).p = [individual(i).p j];
        end
    end
    if individual(i).n == 0
        x(i,n_obj + 1) = 1;
        F(front).f = [F(front).f i];
    end
end
while ~isempty(F(front).f)
    Q = [];
    for i = 1 : length(F(front).f)
        if ~isempty(individual(F(front).f(i)).p)
            for j = 1 : length(individual(F(front).f(i)).p)
                individual(individual(F(front).f(i)).p(j)).n = ...
                    individual(individual(F(front).f(i)).p(j)).n - 1;
                if individual(individual(F(front).f(i)).p(j)).n == 0
                    x(individual(F(front).f(i)).p(j),n_obj  + 1) = ...
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
        end
    end
    front =  front + 1;
    F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,n_obj  + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
end

end







