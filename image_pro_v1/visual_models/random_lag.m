function R2=random_lag(R,lag,addto)
    lag_fit=0; c=0; MAX_Iterations=1000; L=length(R)+addto;
    while lag_fit==0,
        c=c+1;
        if c>MAX_Iterations,
            % Prevent endless repitition as problem can be impossible...
            error('Specified lag is too long!\n');
        end
        R2=randperm(L); 
        for ii=1:L,
            idx=find(R2(ii)==R);        % Find matching element in R1
            if ~isempty(idx)
                distance_diff(ii)=ii+L-idx; % Find distance to element
                % If distance is less then required lag, try again
                if distance_diff(ii)<=lag, break; end
            end
        end
        if min(distance_diff)>=lag, lag_fit=1; end % Distance always exceeds lag, R2 is found
    end
    fprintf('Randomization Facts\n')
    fprintf(['\tIterations: ' num2str(c) '\n']);
    fprintf(['\tMin Dstnce: ' num2str(min(distance_diff)) '\n']);
    fprintf(['\tAvg Dstnce: ' num2str(mean(distance_diff)) '\n']);
    fprintf(['\tMax Dstnce: ' num2str(max(distance_diff)) '\n']);
    
end