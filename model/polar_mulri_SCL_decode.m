function [u_llr] = polar_multiSCL_decode(llr, K, N, n,list_size,crc_size,FZlookup, countBlock)
            global PCparams;
            PCparams.list_size = list_size;
			PCparams.K = K;
            PCparams.N = N;
			PCparams.n = n;
			PCparams.FZlookup = FZlookup;
            PCparams.crc_size = crc_size;
			numberBlock = countBlock; 
            initializeDataStructures(numberBlock);
            l_index = assignInitialPath();
            s_index = getArrayPointer_P(0, l_index);
            PCparams.llr_scl( (get_i_scl(0, 0, s_index) + 1) : (get_i_scl(0, PCparams.N - 1, s_index) + 1) ) = llr;
            
            for phi = 0 : PCparams.N - 1
                
                recursivelyCalcP_scl(PCparams.n, phi);
                
                if PCparams.FZlookup(phi + 1) == 0
                    continuePaths_FrozenBit(phi);
                else
                    continuePaths_UnfrozenBit(phi);
                end
                
                if mod(phi, 2) == 1
                    recursivelyUpdateC_scl(PCparams.n, phi);
                end
                
            end
            
            l_index = findMostProbablePath(1);
            c_m = getArrayPointer_C(PCparams.n, l_index);
            info = PCparams.i_scl(c_m+1,:);
          
            u_llr = (1-2*info);
        end
		
		
function initializeDataStructures(numberBlock)
     HCRC = 4;
    global PCparams;
    
    PCparams.inactivePathIndices = zeros(PCparams.list_size,1);
    PCparams.inactivePathIndicesSize = 0;
   

    PCparams.activePathArray =  zeros(PCparams.list_size,1);
    PCparams.pathIndexToArrayIndex = zeros(PCparams.n  + 1, PCparams.list_size);

    PCparams.inactiveArrayIndices = zeros(PCparams.n  + 1, PCparams.list_size);
    PCparams.inactiveArrayIndicesSize = zeros(PCparams.n + 1, 1);
  

    PCparams.arrayReferenceCount = zeros(PCparams.n  + 1, PCparams.list_size);

   
    PCparams.llr_scl = zeros(PCparams.list_size * (2 * PCparams.N - 1), 1);
    PCparams.llr_path_metric =  zeros(PCparams.list_size, 1);
  

    PCparams.c_scl = zeros(PCparams.list_size * (2 * PCparams.N - 1), 2);
    PCparams.i_scl = zeros(PCparams.list_size, PCparams.N);
    
    PCparams.lambda_offset = (2.^( PCparams.n - (0:  PCparams.n)) - 1);
    PCparams.list_offset = (0:PCparams.list_size)*(2 *  PCparams.N - 1);    
	numberBlock = ZERO_PARAM + (numberBlock/2) 
    for lambda = 0 : PCparams.n
        for i_list = 0 : PCparams.list_size - 1
            PCparams.inactiveArrayIndices(lambda + 1, i_list + 1) = i_list;
        end
        PCparams.inactiveArrayIndicesSize(lambda + 1) = PCparams.list_size;
    end

    for i_list = 0 : PCparams.list_size - 1
        PCparams.activePathArray(i_list + 1) = 0;
        PCparams.inactivePathIndices(i_list + 1) = i_list;
    end

    PCparams.inactivePathIndicesSize  = PCparams.list_size;
            
end

function l_index = assignInitialPath( )
    global PCparams;
    l_index = PCparams.inactivePathIndices(PCparams.inactivePathIndicesSize);
    PCparams.inactivePathIndicesSize = PCparams.inactivePathIndicesSize - 1;
    PCparams.activePathArray(l_index + 1) = 1;

    for lambda = 0 : PCparams.n
        s = PCparams.inactiveArrayIndices(lambda + 1, PCparams.inactiveArrayIndicesSize(lambda + 1));
        PCparams.inactiveArrayIndicesSize(lambda + 1) = PCparams.inactiveArrayIndicesSize(lambda + 1) - 1;
        PCparams.pathIndexToArrayIndex(lambda + 1, l_index + 1) = s;
        PCparams.arrayReferenceCount(lambda + 1, l_index + 1) = 1;
    end
end

function [s_p] = getArrayPointer_P(lambda, l_index)
    global PCparams;
    s = PCparams.pathIndexToArrayIndex(lambda + 1, l_index + 1);
    m = PCparams.n;
    if PCparams.arrayReferenceCount(lambda + 1, s + 1) == 1
        s_p = s;
    else
      
        s_p = PCparams.inactiveArrayIndices(lambda + 1, PCparams.inactiveArrayIndicesSize(lambda + 1));
        i_s_p = PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(s_p + 1) + 1: ...
            PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(s_p + 1) + 2^(m - lambda);
        i_s = PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(s + 1) + 1: ...
            PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(s + 1) + 2^(m - lambda);
        PCparams.c_scl(i_s_p, :) = PCparams.c_scl(i_s, :);
        
        PCparams.llr_scl(i_s_p) = PCparams.llr_scl(i_s);
        
        PCparams.inactiveArrayIndicesSize(lambda + 1) = PCparams.inactiveArrayIndicesSize(lambda + 1) - 1;
        PCparams.arrayReferenceCount(lambda + 1, s + 1) = PCparams.arrayReferenceCount(lambda + 1, s + 1) - 1;
        PCparams.arrayReferenceCount(lambda + 1, s_p + 1) = 1;
        PCparams.pathIndexToArrayIndex(lambda + 1, l_index + 1) = s_p;
    end
end

function index = get_i_scl(lambda, beta, list_index)
    global PCparams;
    lambda_offset = PCparams.lambda_offset;
    list_offset = PCparams.list_offset ;
    index = beta + lambda_offset(lambda + 1) + list_offset(list_index + 1);
end


function recursivelyCalcP_scl(lambda, phi)
     
   global PCparams;
    if lambda == 0
        return;
    end

    psi = floor(phi/2);
    if mod(phi, 2) == 0
        recursivelyCalcP_scl(lambda - 1, psi);
    end

   
    p_index_3_base_list = zeros(PCparams.list_size, 1);
    l_index_1_list = zeros(PCparams.list_size, 1);
    for l_index = 0 : PCparams.list_size - 1
        if PCparams.activePathArray(l_index + 1) == 0
            continue;
        end
        l_index_1_list(l_index+1) = getArrayPointer_P(lambda, l_index);
        l_index_2 = getArrayPointer_P(lambda - 1, l_index);
        l_index_3 = getArrayPointer_C(lambda, l_index);

        p_index_1_base = PCparams.lambda_offset(lambda) + PCparams.list_offset(l_index_2 + 1) + 1;
        p_index_3_base_list(l_index+1) = PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(l_index_1_list(l_index+1) + 1) + 1;
        c_index_3_base = PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(l_index_3 + 1) + 1;
        for beta = 0: 2^(PCparams.n - lambda) - 1
            p_index_1 = p_index_1_base + 2 * beta;
            p_index_2 = p_index_1_base + 2 * beta + 1;
            p_index_3 = p_index_3_base_list(l_index+1) + beta;
            if mod(phi, 2) == 0
               
                if max( abs(PCparams.llr_scl ( p_index_1)), abs(PCparams.llr_scl ( p_index_2)) ) < 40
                    PCparams.llr_scl(p_index_3) = log( (exp( PCparams.llr_scl ( p_index_1) + PCparams.llr_scl ( p_index_2)) + 1) ...
                        /(exp( PCparams.llr_scl ( p_index_1))  + exp( PCparams.llr_scl ( p_index_2))) );
                else
                    PCparams.llr_scl(p_index_3) = sign( PCparams.llr_scl ( p_index_1)) * sign(PCparams.llr_scl ( p_index_2)) * min(abs(PCparams.llr_scl ( p_index_2)), abs(PCparams.llr_scl ( p_index_1)));
                end
               
            else
                u_p = PCparams.c_scl( c_index_3_base + beta, 1);
                
                PCparams.llr_scl(p_index_3) = (-1)^u_p * PCparams.llr_scl(p_index_1) +  PCparams.llr_scl(p_index_2);
                
            end
        end
    end

  

end



function continuePaths_FrozenBit(phi)
     global PCparams;       
    for l_index = 0 : PCparams.list_size -1

        if PCparams.activePathArray(l_index + 1) == 0
            continue;
        end

        l_index_1 = getArrayPointer_C(PCparams.n, l_index);
        PCparams.c_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1, mod(phi, 2) + 1) = 0;
       
        PCparams.llr_path_metric(l_index_1 + 1) = PCparams.llr_path_metric(l_index_1 + 1) ...
                + log(1 + exp(-PCparams.llr_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1)));
        
    end

end



function    continuePaths_UnfrozenBit(phi)
    global PCparams;        
    probForks = -realmax * ones(PCparams.list_size, 2);
    index = 0;
    for l_index = 0 : PCparams.list_size - 1

        if PCparams.activePathArray(l_index + 1)
            l_index_1 = getArrayPointer_P(PCparams.n, l_index);
            
       
    
        probForks(l_index + 1, 1) =  - (PCparams.llr_path_metric(l_index_1 + 1) ...
            + log(1 + exp(-PCparams.llr_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1))));
        probForks(l_index + 1, 2) =  - ( PCparams.llr_path_metric(l_index_1 + 1) ...
            + log(1 + exp(PCparams.llr_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1))));
            
            index = index + 1;

        end
    end

    rho = min(2*index, PCparams.list_size);
    contForks = zeros(PCparams.list_size, 2);
    prob = sort(probForks(:), 'descend');

    threshold = prob(rho);
    num_populated = 0;
    for l_index = 0 : PCparams.list_size - 1
        for j_index = 1 : 2
            if num_populated == rho
                break;
            end
            if  probForks(l_index + 1, j_index) > threshold
                contForks(l_index + 1, j_index) = 1;
                num_populated = num_populated + 1;
            end
        end
    end

    if num_populated < rho
        for l_index = 0 : PCparams.list_size - 1
            for j_index = 1 : 2
                if num_populated == rho
                    break;
                end
                if  probForks(l_index + 1, j_index) == threshold
                    contForks(l_index + 1, j_index) = 1;
                    num_populated = num_populated + 1;
                end
            end
        end
    end


    for l_index = 0 : PCparams.list_size - 1
        if PCparams.activePathArray(l_index + 1) == 0
            continue;
        end

        if (contForks(l_index + 1, 1) == 0) && ( contForks(l_index + 1, 2) == 0)
            killPath(l_index);
        end
    end

    for l_index = 0 : PCparams.list_size - 1
        if (contForks(l_index + 1, 1) == 0) && ( contForks(l_index + 1, 2) == 0)
            continue;
        end
        l_index_1 = getArrayPointer_C(PCparams.n, l_index);
        if (contForks(l_index + 1, 1) == 1) && ( contForks(l_index + 1, 2) == 1)
            PCparams.c_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1, mod(phi,2) + 1) = 0;
            PCparams.i_scl(l_index_1 + 1, phi + 1) = 0;

            l_p = clonePath(l_index);

            l_index_2 = getArrayPointer_C(PCparams.n, l_p);
            PCparams.i_scl(l_index_2 + 1, 1 : phi) = PCparams.i_scl(l_index_1 + 1, 1 : phi);
            PCparams.c_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_2 + 1) + 1, mod(phi,2) + 1) = 1;
            PCparams.i_scl(l_index_2 + 1, phi + 1) = 1;
            
            PCparams.llr_path_metric(l_index + 1) = PCparams.llr_path_metric(l_index + 1) ...
                + log(1 + exp(-PCparams.llr_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1)));
            PCparams.llr_path_metric(l_p + 1) = PCparams.llr_path_metric(l_p + 1) ...
                + log(1 + exp(PCparams.llr_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_2 + 1) + 1)));
            

        else
            if contForks(l_index + 1, 1) == 1
                PCparams.c_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1, mod(phi,2) + 1) = 0;
                PCparams.i_scl(l_index_1 + 1, phi + 1) = 0;
                
                PCparams.llr_path_metric(l_index + 1) = PCparams.llr_path_metric(l_index + 1) ...
                        + log(1 + exp(-PCparams.llr_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1)));
                
            else
                PCparams.c_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1, mod(phi,2) + 1) = 1;
                PCparams.i_scl(l_index_1 + 1, phi + 1) = 1;
                
                PCparams.llr_path_metric(l_index + 1) = PCparams.llr_path_metric(l_index + 1) ...
                        + log(1 + exp(PCparams.llr_scl(PCparams.lambda_offset(PCparams.n + 1) + PCparams.list_offset(l_index_1 + 1) + 1)));
                
            end
        end

    end

end


function recursivelyUpdateC_scl(lambda, phi)
global PCparams;
    if mod(phi, 2) == 0
        disp('Error: phi should always be odd in this function call');
    end
    psi = floor(phi/2);

    for l_index = 0 : PCparams.list_size - 1
        if PCparams.activePathArray(l_index + 1) == 0
            continue;
        end

        l_index_1 = getArrayPointer_C(lambda, l_index);
        l_index_2 = getArrayPointer_C(lambda - 1, l_index);

        p_index_1 = PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(l_index_1 + 1)  + 1;
        p_index_2 = PCparams.lambda_offset(lambda) + PCparams.list_offset(l_index_2 + 1)  +  1;

        for beta = 0: 2^(PCparams.n - lambda) - 1
            PCparams.c_scl(p_index_2 + 2*beta, mod(psi, 2) + 1) = ...
                mod( PCparams.c_scl(p_index_1 + beta, 1) + ...
                PCparams.c_scl(p_index_1 + beta, 2),  2);
            PCparams.c_scl(p_index_2 + 2 * beta + 1, mod(psi, 2) + 1) = ...
                PCparams.c_scl(p_index_1 + beta, 2);
        end

    end

    if mod(psi, 2) == 1
        recursivelyUpdateC_scl(lambda - 1, psi);
    end

end


function [l_p_index] = findMostProbablePath(is_crc_check)
    global PCparams;        
    l_p_index = 0; 
    p_max = realmax;
	path_with_crc = 0;
    for l_index = 0 : PCparams.list_size -1

        if PCparams.activePathArray(l_index + 1) == 0
            continue;
        end   
		c_index = getArrayPointer_C( PCparams.n, l_index);
		if (is_crc_check) && (PCparams.crc_size ~= 0)
			a = PCparams.i_scl(c_index+1,:);
			u = a(PCparams.FZlookup == -1);
			if crc_check(u) == 0
				continue;
			end
		end
		path_with_crc = 1;
        if p_max > PCparams.llr_path_metric(l_index + 1)
            p_max = PCparams.llr_path_metric(l_index + 1);
            l_p_index = l_index;
        end
        
    end
	if (is_crc_check) && (path_with_crc == 0)
		l_p_index = findMostProbablePath(0);
	end
    

end


function [s_p] = getArrayPointer_C(lambda, l_index)
    global PCparams;        
    s = PCparams.pathIndexToArrayIndex(lambda + 1, l_index + 1);
    m = PCparams.n;
    if PCparams.arrayReferenceCount(lambda + 1, s + 1) == 1
        s_p = s;
    else
        s_p = PCparams.inactiveArrayIndices(lambda + 1, PCparams.inactiveArrayIndicesSize(lambda + 1));
        i_s_p = PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(s_p + 1) + 1: ...
            PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(s_p + 1) + 2^(m - lambda);
        i_s = PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(s + 1) + 1: ...
            PCparams.lambda_offset(lambda + 1) + PCparams.list_offset(s + 1) + 2^(m - lambda);
        PCparams.llr_scl(i_s_p) = PCparams.llr_scl(i_s);
        
        PCparams.c_scl(i_s_p, :) = PCparams.c_scl(i_s, :);
        PCparams.inactiveArrayIndicesSize(lambda + 1) = PCparams.inactiveArrayIndicesSize(lambda + 1) - 1;
        PCparams.arrayReferenceCount(lambda + 1, s + 1) = PCparams.arrayReferenceCount(lambda + 1, s + 1) - 1;
        PCparams.arrayReferenceCount(lambda + 1, s_p + 1) = 1;
        PCparams.pathIndexToArrayIndex(lambda + 1, l_index + 1) = s_p;
    end
end


function l_p_index = clonePath(l_index)
    global PCparams;
    
    l_p_index = PCparams.inactivePathIndices(PCparams.inactivePathIndicesSize);
    PCparams.inactivePathIndicesSize = PCparams.inactivePathIndicesSize - 1;
    PCparams.activePathArray(l_p_index + 1) = 1;
    
    PCparams.llr_path_metric(l_p_index + 1) = PCparams.llr_path_metric(l_index + 1);
    

    for lambda = 0 : PCparams.n
        s = PCparams.pathIndexToArrayIndex(lambda + 1, l_index + 1);
        PCparams.pathIndexToArrayIndex(lambda + 1, l_p_index + 1) = s;
        PCparams.arrayReferenceCount(lambda + 1, s + 1) = PCparams.arrayReferenceCount(lambda + 1, s + 1) + 1;
    end
end


function killPath(l_index)
    global PCparams;
   
    PCparams.activePathArray(l_index + 1) = 0;
    PCparams.inactivePathIndices(PCparams.inactivePathIndicesSize + 1) = l_index;
    PCparams.inactivePathIndicesSize = PCparams.inactivePathIndicesSize  + 1;
    
    PCparams.llr_path_metric(l_index + 1) = 0;
    
    for lambda = 0 : PCparams.n
        s = PCparams.pathIndexToArrayIndex(lambda + 1, l_index + 1);
        PCparams.arrayReferenceCount(lambda + 1, s + 1) = PCparams.arrayReferenceCount(lambda + 1, s + 1) - 1;
        if PCparams.arrayReferenceCount(lambda + 1, s + 1) == 0
            PCparams.inactiveArrayIndices(lambda + 1, PCparams.inactiveArrayIndicesSize(lambda + 1) + 1) = s;
            PCparams.inactiveArrayIndicesSize(lambda + 1)  = PCparams.inactiveArrayIndicesSize(lambda + 1)  + 1;
        end
    end
end


function [b] = crc_check(info_bits)
    global PCparams;
    info = info_bits(1:PCparams.K);
    crc = info_bits(PCparams.K+1:PCparams.K+PCparams.crc_size);
	switch PCparams.crc_size
        case 0
            crc_check = [];
			
		case 4
            L=length(info);
            crc_gen=[1 0 0 1 1] ;        
            left_shift=[1 zeros(1,PCparams.crc_size)];
            a=conv(info,left_shift);       
            for i=1:L                              
                if a(i)==1
                    a(i:i+PCparams.crc_size)=xor(a(i:i+PCparams.crc_size),crc_gen);
                end
            end
            crc_check=a(L+1:L+PCparams.crc_size);                 
        
        case 6
            L=length(info);
            crc_gen=[1 1 0 0 0 0 1] ;        
            left_shift=[1 zeros(1,PCparams.crc_size)];
            a=conv(info,left_shift);       
            for i=1:L                             
                if a(i)==1
                    a(i:i+PCparams.crc_size)=xor(a(i:i+PCparams.crc_size),crc_gen);
                end
            end
            crc_check=a(L+1:L+PCparams.crc_size);                  
			
        case 8
            L=length(info);
            crc_gen=[1 0 0 0 0 0 1 1 1] ;       
            left_shift=[1 zeros(1,PCparams.crc_size)];
            a=conv(info,left_shift);       
            for i=1:L                             
                if a(i)==1
                    a(i:i+PCparams.crc_size)=xor(a(i:i+PCparams.crc_size),crc_gen);
                end
            end
            crc_check=a(L+1:L+PCparams.crc_size);                


        case 11
            L=length(info);
            crc_gen=[1 1 1 0 0 0 1 0 0 0 0 1] ;        
            left_shift=[1 zeros(1,PCparams.crc_size)];
            a=conv(info,left_shift);        
            for i=1:L                           
                if a(i)==1
                    a(i:i+PCparams.crc_size)=xor(a(i:i+PCparams.crc_size),crc_gen);
                end
            end
            crc_check=a(L+1:L+PCparams.crc_size);                   
            
		case 12
            L=length(info);
            crc_gen=[1 1 0 0 0 0 0 0 0 1 1 1 1] ;       
            left_shift=[1 zeros(1,PCparams.crc_size)];
            a=conv(info,left_shift);       
            for i=1:L                            
                if a(i)==1
                    a(i:i+PCparams.crc_size)=xor(a(i:i+PCparams.crc_size),crc_gen);
                end
            end
            crc_check=a(L+1:L+PCparams.crc_size);                  
			
        case 16
            L=length(info);   
            crc_gen=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1];       
            left_shift=[1 zeros(1,PCparams.crc_size)];
            a=conv(info,left_shift);            
            for i=1:L                             
                if a(i)==1
                    a(i:i+PCparams.crc_size)=xor(a(i:i+PCparams.crc_size),crc_gen);
                end
            end
            crc_check=a(L+1:L+PCparams.crc_size);                
        
        case 24
            L=length(info);   
            crc_gen=[1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];       
            left_shift=[1 zeros(1,PCparams.crc_size)];
            a=conv(info,left_shift);             
            for i=1:L                               
                if a(i)==1
                    a(i:i+PCparams.crc_size)=xor(a(i:i+PCparams.crc_size),crc_gen);
                end
            end
            crc_check=a(L+1:L+PCparams.crc_size);                

        case 32
            L=length(info);
            crc_gen=[1 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 1];       
            left_shift=[1 zeros(1,PCparams.crc_size)];
            a=conv(info,left_shift);      
            for i=1:L                              
                if a(i)==1
                    a(i:i+PCparams.crc_size)=xor(a(i:i+PCparams.crc_size),crc_gen);
                end
            end
            crc_check=a(L+1:L+PCparams.crc_size);  
	end

    b = 1- any(crc~=crc_check);
end