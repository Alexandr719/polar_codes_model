function y=pencode(u,myfrozenlookup,crc_size,bitreversedindices,F_kron_n) 


x = myfrozenlookup;
switch crc_size
        case 0
            crc_code = [];
		case 4
			L=length(u);
            crc_gen=[1 0 0 1 1] ;        
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);       
            for i=1:L                               
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                  
        case 6
			L=length(u);
            crc_gen=[1 1 0 0 0 0 1] ;       
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);       
            for i=1:L                             
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                  
			
        case 8
            L=length(u);
            crc_gen=[1 0 0 0 0 0 1 1 1] ;        
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);       
            for i=1:L                             
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                  
        case 11
            L=length(u);
            crc_gen=[1 1 1 0 0 0 1 0 0 0 0 1] ;        
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);       
            for i=1:L                             
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                 
			
		case 12
		L=length(u);
		crc_gen=[1 1 0 0 0 0 0 0 0 1 1 1 1] ;        
		left_shift=[1 zeros(1,crc_size)];
		a=conv(u,left_shift);      
		for i=1:L                            
			if a(i)==1
				a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
			end
		end
		crc_code=a(L+1:L+crc_size);                  


        case 16
            L=length(u);   
            crc_gen=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1];       
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);               
            for i=1:L                             
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                 
        case 24
            L=length(u);   
            crc_gen=[1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];       
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);              
            for i=1:L                             
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                   

        case 32
            L=length(u);
            crc_gen=[1 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 1];        
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);       
            for i=1:L                               
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);  
end
    
u = [u crc_code];

x (x == -1) = u;
x = x(bitreversedindices+1);
y = mod(x*F_kron_n,2);

end