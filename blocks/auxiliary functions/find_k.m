function [k,mid_val]=find_k(y,s)
% find the k defined in Nikolov's paper with respect to y which is a
% d-dimensional vector
    accuracy=1e-10;
    sumtemp=sum(y); % temporary variable storting the value sum(y((k+1):end))
    mid_val=sumtemp/s;
    
    if mid_val >= y(1)+accuracy
        k=0;
        return
    end
    
    for k=1:(s-1)
        sumtemp=sumtemp-y(k);
        mid_val=sumtemp/(s-k);
        if mid_val >= y(k+1)-accuracy && mid_val < y(k)+accuracy
            return
        end
    end
end
