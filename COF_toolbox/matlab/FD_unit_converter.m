% FD unit converter
% input:
%       para: the parameter struct, with at least vf, w, kc or q_max
%       from_unit: string, 'm' or 'mile', 'km'
%       to_unit: string, 'm' or 'mile', 'km'
%               'm': all in m, s; 'mile', all in mile, hour; 
%               'km' all in km, hour
% output:
%       converted_para: the converted para with all fields converted

function converted_para = FD_unit_converter(para, from_unit, to_unit)
            
    converted_para = para;
    
    if isfield(para, 'vf') && isfield(para, 'w') && isfield(para, 'kc_pl')
        
        if strcmp(from_unit, 'mile') && strcmp(to_unit, 'm')
        
            converted_para.vf = para.vf*1609/3600;  
            converted_para.w = para.w*1609/3600;    
            converted_para.kc_pl = para.kc_pl/1609; 
            converted_para.qmax_pl = converted_para.kc_pl*converted_para.vf; 
            converted_para.km_pl = converted_para.kc_pl*...
                (converted_para.w-converted_para.vf)/converted_para.w;
            
            % other optional fields;
            % if not exist, then add some default values
            if isfield(para, 'qon_max')
                converted_para.qon_max = para.qon_max/3600;
            else
                converted_para.qon_max = converted_para.qmax_pl;
            end
            
            if isfield(para, 'qoff_max')
                converted_para.qoff_max = para.qoff_max/3600;
            else
                converted_para.qoff_max = converted_para.qmax_pl;
            end
            
            if isfield(para, 'v_min')
                converted_para.v_min = para.v_min*1609/3600;
            else
                converted_para.v_min = 0;
            end
            
            if isfield(para, 'v_max')
                converted_para.v_max = para.v_max*1609/3600;
            else
                converted_para.v_max = 75*1609/3600;
            end
            
            return
        
        elseif strcmp(from_unit, 'mile') && strcmp(to_unit, 'km')
        
            converted_para.vf = para.vf*1.609;  
            converted_para.w = para.w*1.609;    
            converted_para.kc_pl = para.kc_pl/1.609; 
            converted_para.qmax_pl = converted_para.kc_pl*converted_para.vf; 
            converted_para.km_pl = converted_para.kc_pl*...
                (converted_para.w-converted_para.vf)/converted_para.w;
            
            % other optional fields;
            % if not exist, then add some default values
            if isfield(para, 'v_min')
                converted_para.v_min = para.v_min*1.609;
            else
                converted_para.v_min = 0;
            end
            
            if isfield(para, 'v_max')
                converted_para.v_max = para.v_max*1.609;
            else
                converted_para.v_max = 75*1.609;
            end
            
            return
        
        elseif strcmp(from_unit, 'm') && strcmp(to_unit, 'mile')
        
            converted_para.vf = para.vf*3600/1609;  
            converted_para.w = para.w*3600/1609;    
            converted_para.kc_pl = para.kc_pl*1609; 
            converted_para.qmax_pl = converted_para.kc_pl*converted_para.vf; 
            converted_para.km_pl = converted_para.kc_pl*...
                (converted_para.w-converted_para.vf)/converted_para.w;
            
            % other optional fields;
            % if not exist, then add some default values
            if isfield(para, 'qon_max')
                converted_para.qon_max = para.qon_max*3600;
            else
                converted_para.qon_max = converted_para.qmax_pl;
            end
            
            if isfield(para, 'qoff_max')
                converted_para.qoff_max = para.qoff_max*3600;
            else
                converted_para.qoff_max = converted_para.qmax_pl;
            end
            
            if isfield(para, 'v_min')
                converted_para.v_min = para.v_min*3600/1609;
            else
                converted_para.v_min = 0;
            end
            
            if isfield(para, 'v_max')
                converted_para.v_max = para.v_max*3600/1609;
            else
                converted_para.v_max = 75;
            end
            
            return
        
        elseif strcmp(from_unit, 'm') && strcmp(to_unit, 'km')
        
            converted_para.vf = para.vf*3600/1000;  
            converted_para.w = para.w*3600/1000;    
            converted_para.kc_pl = para.kc_pl*1000; 
            converted_para.qmax_pl = converted_para.kc_pl*converted_para.vf; 
            converted_para.km_pl = converted_para.kc_pl*...
                (converted_para.w-converted_para.vf)/converted_para.w;
            
            % other optional fields;
            % if not exist, then add some default values
            if isfield(para, 'qon_max')
                converted_para.qon_max = para.qon_max*3600;
            else
                converted_para.qon_max = converted_para.qmax_pl;
            end
            
            if isfield(para, 'qoff_max')
                converted_para.qoff_max = para.qoff_max*3600;
            else
                converted_para.qoff_max = converted_para.qmax_pl;
            end
            
            if isfield(para, 'v_min')
                converted_para.v_min = para.v_min*3600/1000;
            else
                converted_para.v_min = 0;
            end
            
            if isfield(para, 'v_max')
                converted_para.v_max = para.v_max*3600/1000;
            else
                converted_para.v_max = 75*1.609;
            end
            
            return
        
        elseif strcmp(from_unit, 'km') && strcmp(to_unit, 'm')
        
            converted_para.vf = para.vf*1000/3600;  
            converted_para.w = para.w*1000/3600;    
            converted_para.kc_pl = para.kc_pl/1000; 
            converted_para.qmax_pl = converted_para.kc_pl*converted_para.vf; 
            converted_para.km_pl = converted_para.kc_pl*...
                (converted_para.w-converted_para.vf)/converted_para.w;
            
            % other optional fields;
            % if not exist, then add some default values
            if isfield(para, 'qon_max')
                converted_para.qon_max = para.qon_max/3600;
            else
                converted_para.qon_max = converted_para.qmax_pl;
            end
            
            if isfield(para, 'qoff_max')
                converted_para.qoff_max = para.qoff_max/3600;
            else
                converted_para.qoff_max = converted_para.qmax_pl;
            end
            
            if isfield(para, 'v_min')
                converted_para.v_min = para.v_min*1000/3600;
            else
                converted_para.v_min = 0;
            end
            
            if isfield(para, 'v_max')
                converted_para.v_max = para.v_max*1000/3600;
            else
                converted_para.v_max = 75*1609/3600;
            end
            
            return   
            
        elseif strcmp(from_unit, 'km') && strcmp(to_unit, 'mile')
        
            converted_para.vf = para.vf/1.609;  
            converted_para.w = para.w/1.609;    
            converted_para.kc_pl = para.kc_pl*1.609; 
            converted_para.qmax_pl = converted_para.kc_pl*converted_para.vf; 
            converted_para.km_pl = converted_para.kc_pl*...
                (converted_para.w-converted_para.vf)/converted_para.w;
            
            if isfield(para, 'v_min')
                converted_para.v_min = para.v_min/1.609;
            else
                converted_para.v_min = 0;
            end
            
            if isfield(para, 'v_max')
                converted_para.v_max = para.v_max/1.609;
            else
                converted_para.v_max = 75;
            end
            
            return
        
        end
        
        
    end
    
    
    if isfield(para, 'vf') && isfield(para, 'w') && isfield(para, 'qmax_pl')
        
        para.kc_pl = para.qmax_pl/para.vf;
        
        converted_para = FD_unit_converter(para, from_unit, to_unit);
        return
        
    end
    
    
    if isfield(para, 'qmax_pl') && isfield(para, 'w') && isfield(para, 'kc_pl')
        
        para.vf = para.qmax_pl/para.kc_pl;
        
        converted_para = FD_unit_converter(para, from_unit, to_unit);
        return
        
    end
        
        

    
    

end