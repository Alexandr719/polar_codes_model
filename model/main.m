
clear; tic;


polar_K = 80;
polar_N = 256;
polar_n = log2(polar_N);

SCAN_decoder_iter = 1;
BP_decoder_iter = 40;
SCL_list_size = 8;
crc_size = 0;

Rc = polar_K/polar_N;

construction_method = 0;
design_snr_dB = 0;
sigma = 0.9;

F = [1 0;1 1];
BB=1;
for ii=1:polar_n
    BB = kron(BB,F);
end
F_kron_n = BB;
bitreversedindices = zeros(1,polar_N);
for index = 1 : polar_N
    bitreversedindices(index) = bin2dec(wrev(dec2bin(index-1,polar_n)));
end
switch construction_method
    case 0
        constructed_code_file_name = sprintf('constructedPolarCode\\PolarCode_block_length_%d_designSNR_%.2fdB_method_BhattaBound.txt',polar_N,design_snr_dB);
    case 1
        constructed_code_file_name = sprintf('constructedPolarCode\\PolarCode_block_length_%d_designSNR_%.2fdB_method_MC.txt',polar_N,design_snr_dB);
    case 2
        constructed_code_file_name = sprintf('constructedPolarCode\\PolarCode_block_length_%d_sigma_%.2f_method_GA.txt',polar_N,sigma);
    otherwise
        error('The range of construction_method is from 0 to 2!');
end

indices = load(constructed_code_file_name);
FZlookup = zeros(1,polar_N);
if crc_size == 0
    FZlookup(indices(1:polar_K)) = -1;
else
    FZlookup(indices(1:polar_K+crc_size)) = -1;
end






Rm = 1;
ebn0 = 0:0.5:8;
ebn0_num = 10.^(ebn0/10);
SNR = ebn0 + 10*log10(Rc*Rm)+10*log10(2); 
noise_sigma = 1./(10.^(SNR/10));
min_simBits_errors = 100;
min_simBits_num = 1e5;
max_frame_num = 1e7;
max_simBits_num = 1e7;
FER = zeros(1,length(ebn0));
BER = zeros(1,length(ebn0));
bpsk_FER=zeros(1,length(ebn0));       
bpsk_BER=zeros(1,length(ebn0));  



for j = 1:length(ebn0)
	tt=tic();
	for l = 1:max_frame_num
        de_bpsk = zeros(1,polar_N);
        u=randi(2,1,polar_K)-1; 
		x=pencode(u,FZlookup,crc_size,bitreversedindices,F_kron_n);
        tx_waveform=2*x-1;
        noise=sqrt(noise_sigma(j))*randn(1,polar_N);
        rx_waveform = tx_waveform+noise;
        de_bpsk(rx_waveform>0)=1;
        
        nfails = sum(de_bpsk ~= x);
        bpsk_FER(j) = bpsk_FER(j) + (nfails>0);
        bpsk_BER(j) = bpsk_BER(j) + nfails;
        
	
        initial_llr = -2*rx_waveform/noise_sigma(j);
	

       [u_llr] = polar_SCL_decode(initial_llr, polar_K, polar_N, polar_n,SCL_list_size,crc_size,FZlookup);


        if crc_size
            uhat_crc_llr = u_llr(FZlookup == -1)';
            uhat_llr = uhat_crc_llr (1:polar_K);
        else
            uhat_llr = u_llr(FZlookup == -1)';
        end
		uhat = zeros(1,polar_K);
        uhat(uhat_llr<0) =1;

		nfails = sum(uhat ~= u);
        FER(j) = FER(j) + (nfails>0);
        BER(j) = BER(j) + nfails;

        if (l*polar_K>=min_simBits_num && BER(j)>=min_simBits_errors  || l*polar_K >=max_simBits_num)
            break;
        end
	end

        FER(j) = FER(j)/l;
        BER(j) = BER(j)/(polar_K*l);
        
        bpsk_BER(j) = bpsk_BER(j)/(polar_N*l);
        bpsk_FER(j) = bpsk_FER(j)/l;
        fprintf('EbN0 = %.2fdB, BER = %.7f \n',ebn0(j),BER(j));
        fprintf('Total time taken: %.2f sec (%d samples) \n',toc(tt),l);
		if(BER(j)<1e-5)
			break;
		end
end


figure;
semilogy(ebn0(1:j),BER(1:j),'bs-','LineWidth',1.5,'MarkerSize',6)
xlabel('Eb/No (dB)')
hold on;

semilogy(ebn0(1:j),FER(1:j),'gs:','LineWidth',1.5,'MarkerSize',6)
hold on;

semilogy(ebn0(1:j),bpsk_BER(1:j),'rv-','LineWidth',1.5,'MarkerSize',6)
hold on;

semilogy(ebn0(1:j),bpsk_FER(1:j),'rv:','LineWidth',1.5,'MarkerSize',6)


grid on;


BPSK_BER_legendname = 'BPSK BER';
BPSK_FER_legendname = 'BPSK FER';
polar_BER_legend_name = sprintf('BER,PolarCode(%d,%d),SCL(%d,%d)',polar_K,polar_N,SCL_list_size,crc_size);
polar_FER_legend_name = sprintf('FER,PolarCode(%d,%d),SCL(%d,%d)',polar_K,polar_N,SCL_list_size,crc_size);

legend(polar_BER_legend_name,polar_FER_legend_name,BPSK_BER_legendname,BPSK_FER_legendname);
toc

