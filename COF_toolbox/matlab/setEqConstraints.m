

classdef setEqConstraints
    % This class sets the equality constraints (mass conservation) at junctions of a network
    % Yanning Li, Feb 17, 2016
    
    properties (Access = public)
        
       EqMatrix;  % The equality constraints matrix
        
    end
    
    
    methods (Access = public)
        
        function self = setEqConstraints(network, dv_index, dv_index_max)
            % This function saves the quality constraints.
            % input:
            %       - network, a initNetwork object with the network information
            %       - dv_index, struct, the decision variable index in CP
            %       - dv_index_max: the maximum number of decision variabels
            % output:
            %       - the equality constraints matrix is saved in self.EqMatrix
            
            self.EqMatrix = zeros(0,dv_index_max);
            
            for junc = network.junc_labels'
                
                juncStr = sprintf('junc_%d',junc);
                
                num_steps = length(network.network_junc.(juncStr).T);
                
                array = zeros(num_steps, dv_index_max );
                
                if strcmp(network.network_junc.(juncStr).type_junc,'diverge') ||...
                   strcmp(network.network_junc.(juncStr).type_junc,'offrampjunc')
                    
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).inlabel);
                    array(1:num_steps,...
                          dv_index.(linkStr).downstream(1,1):dv_index.(linkStr).downstream(2,1)) =...
                          diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).outlabel(1));
                    array(1:num_steps,...
                          dv_index.(linkStr).upstream(1,1):dv_index.(linkStr).upstream(2,1)) =...
                          -diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).outlabel(2));
                    array(1:num_steps,...
                          dv_index.(linkStr).upstream(1,1):dv_index.(linkStr).upstream(2,1)) =...
                          -diag(ones(num_steps,1));

                    self.EqMatrix = [self.EqMatrix; array];
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'merge') ||...
                       strcmp(network.network_junc.(juncStr).type_junc,'onrampjunc') 
                    
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).inlabel(1));
                    array(1:num_steps,...
                          dv_index.(linkStr).downstream(1,1):dv_index.(linkStr).downstream(2,1)) =...
                          diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).inlabel(2));
                    array(1:num_steps,...
                          dv_index.(linkStr).downstream(1,1):dv_index.(linkStr).downstream(2,1)) =...
                          diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).outlabel);
                    array(1:num_steps,...
                          dv_index.(linkStr).upstream(1,1):dv_index.(linkStr).upstream(2,1)) =...
                          -diag(ones(num_steps,1));
 
                    self.EqMatrix = [self.EqMatrix; array];
  
                elseif strcmp(network.network_junc.(juncStr).type_junc,'connection')
                    
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).inlabel);
                    array(1:num_steps,...
                          dv_index.(linkStr).downstream(1,1):dv_index.(linkStr).downstream(2,1)) =...
                          diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).outlabel);
                    array(1:num_steps,...
                          dv_index.(linkStr).upstream(1,1):dv_index.(linkStr).upstream(2,1)) =...
                          -diag(ones(num_steps,1));
                                          
                    self.EqMatrix = [self.EqMatrix; array];
                    
                end
                
            end
            
        end
        
        
        
    end
    
end
    
    
