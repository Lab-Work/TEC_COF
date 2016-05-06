function [k_trans] = mapping(k_ori, original_interval, target_interval)
            % This function maps the values from original interval to a target interval
            % input:
            %       k_ori: the original values to be mapped
            %       original_interval: nx2 matrix. Each row is an interval to be mapped
            %           to the corresponding interval in target_interval.
            %       target_interval: nx2 matrix. Each row is the target interval.
            % output:
            %       k_trans: the mapped values which has the same dimension as k_ori
            % Example:
            % k_trans = mapping(k,[0 k_c; k_c k_m],[0 0.5*k_m; 0.5*k_m k_m]);
            % mapping (0~k_c)(k_c~k_,m) ==> (0~0.5*k_m)(0.5*k_m~k_m)
            
            k_trans = 0.0*k_ori;
            for i = 1:size(original_interval,1)
                isInDomain = (k_ori>=original_interval(i,1) & k_ori<=original_interval(i,2)+10e-6);
                k_trans(isInDomain) = target_interval(i,1) + ...
                    (k_ori(isInDomain)-original_interval(i,1))*(target_interval(i,2)-target_interval(i,1))...
                    /(original_interval(i,2)-original_interval(i,1));
            end
            
            
end