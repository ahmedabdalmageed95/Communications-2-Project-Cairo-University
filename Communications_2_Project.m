SNR=-4:4;
FSK_Errors=[];
FSK_theoritical=[];
PSK_theoritical=[];
PSK_Errors=[];

for SNR_ind=1:9 %iterating over SNR values from -4 to 4
    total_PSK_error=0;
     total_FSK_error=0;
      total_th_PSK_error=0;
     total_th_FSK_error=0;
   
    for q=1:20 %20 realizations for each SNR value
Tb=40;
bit_no=100;
t=0:(Tb*bit_no)-1; %for 100 bits time
gen_seq=round(rand(1,bit_no)); %generating random sequence of bits

inputBits=[];
FSK_modulated_signal=[];
PSK_modulated_signal=[];
Nnode=2;
A=sqrt((10^(SNR(SNR_ind)/10))*(2*Nnode/Tb));
%calculating theoretical error
FSK_th=0.5*erfc(sqrt((10^(SNR(SNR_ind)/10))/2)); %calculating theoretical error for BFSK
                                                                                        %equivalent to 
                                                                                        %Eb=Tb*(A^2);
                                                                                        % FSK_th=0.5*erfc(sqrt(Eb/(4*Nnode)));

                                                                                        
PSK_th=0.5*erfc(sqrt((10^(SNR(SNR_ind)/10))));    %calculating theoretical error for BPSK 
                                                                                        %equivalent to 
                                                                                        %Eb=Tb*(A^2);
                                                                                        % PSK_th=0.5*erfc(sqrt(Eb/(2*Nnode)));
No=4;
N1=1;
w1=2*pi*(N1+No)/Tb;
w2=2*pi*(No-N1)/Tb;
wc=4*2*pi/Tb;
for i=1:length(gen_seq) % 1 is represented by 1 0 is represented by -1 one bit lasts for Tb , Tb=40
    for j=1:Tb
     if gen_seq(i)==1
        inputBits=[inputBits  1];
    
    elseif gen_seq(i)==0
         inputBits=[inputBits  -1];
    end
    end
end


for i=1:length(inputBits) % BFSK modulation
    if inputBits(i)==1
        FSK_modulated_signal=[FSK_modulated_signal  A*cos(w1*t(i))];
    
    elseif inputBits(i)==-1
         FSK_modulated_signal=[FSK_modulated_signal  A*cos(w2*t(i))];
    end
end


for i=1:length(inputBits)% BPSK modulation
    if inputBits(i)==1
        PSK_modulated_signal=[PSK_modulated_signal  A*cos(wc*t(i))];
    
    elseif inputBits(i)==-1
         PSK_modulated_signal=[PSK_modulated_signal  -A*cos(wc*t(i))];
    end
end





noise=wgn(1,length(FSK_modulated_signal),Nnode/2,'linear'); %generating noise
FSK_ch_signal=FSK_modulated_signal+noise; %adding noise to BFSK modulated signal
PSK_ch_signal=PSK_modulated_signal+noise; %adding BPSK to modulated signal

FSK_filter1=A*cos(w1*(-1*(1:Tb))); %2 filters for BFSK because 2 frequencies
FSK_filter2=A*cos(w2*(-1*(1:Tb)));

PSK_filter1=A*cos(wc*(-1*(1:Tb))); % a single filter for BPSK is enough

FSK_filter_FSK_ouput1=conv(FSK_ch_signal,FSK_filter1); %output of 1st BFSK filter
FSK_filter_FSK_ouput2=conv(FSK_ch_signal,FSK_filter2); %output of 2nd BFSK filter

PSK_filter_PSK_ouput1=conv(PSK_ch_signal,PSK_filter1); %output of BPSK filter


FSK_ouput=[];
PSK_output=[];
FSK_ctr=0;
PSK_ctr=0;
% decision making for BFSK
% if output of the first filter - output of 2nd filter > 0 then output is 1
% if output of the first filter - output of 2nd filter < 0 then output is 0
%if output of first filter > output of second filter output is 1
%else output is 0
for i=Tb:Tb:length(FSK_filter_FSK_ouput1) %taking sample every Tb
  
    if FSK_filter_FSK_ouput1(i)>FSK_filter_FSK_ouput2(i)
        FSK_ouput=[FSK_ouput 1];
    else 
         FSK_ouput=[FSK_ouput 0];
  
    end
end 
 % decision making for BPSK   
%assuming equiprobable ,therefore Threshold=0
for i=Tb:Tb:length(PSK_filter_PSK_ouput1) %taking sample every Tb
    if PSK_filter_PSK_ouput1(i)>0
        PSK_output=[PSK_output 1];
    else if PSK_filter_PSK_ouput1(i)< 0
         PSK_output=[PSK_output 0];
        end
    end
end 
%calculating number of error bits
for i=1:length(FSK_ouput)
    if FSK_ouput(i)~=gen_seq(i) %calculating number of error bits for BFSK output
       FSK_ctr=FSK_ctr+1;
    
    end
    if PSK_output(i)~=gen_seq(i)%calculating number of error bits for BPSK output
       PSK_ctr=PSK_ctr+1;
    
    end
end
FSK_error=FSK_ctr/length(gen_seq); %calculating error rate for BFSK 
PSK_error=PSK_ctr/length(gen_seq); %calculating error rate for BPSK 

total_PSK_error=total_PSK_error+PSK_error; %accumulating values of realizations of error rate for BPSK
total_FSK_error=total_FSK_error+FSK_error; %accumulating values of realizations of error rate for BFSK
total_th_PSK_error=total_th_PSK_error+PSK_th; %accumulating values of realizations of theoretical error rate for BPSK
total_th_FSK_error=total_th_FSK_error+FSK_th;  %accumulating values of realizations  of theoretical error rate for BFSK
    end
PSK_avg_error=total_PSK_error/20; % calculating average value for error rates for each SNR value for BPSK
FSK_avg_error=total_FSK_error/20; % calculating average value for error rates for each SNR value for BFSK
PSK_avg_th_error=total_th_PSK_error/20; %calculating average value for theoretical error rates for each SNR value for BPSK
FSK_avg_th_error=total_th_FSK_error/20; %calculating average value for theoretical error rates for each SNR value for BFSK
PSK_Errors=[PSK_Errors PSK_avg_error];
FSK_Errors=[FSK_Errors FSK_avg_error];
FSK_theoritical=[FSK_theoritical FSK_avg_th_error];
PSK_theoritical=[PSK_theoritical PSK_avg_th_error];
end
%plotting theoretical versus practical curves
 subplot(3,1,1); %plotting theoretical versus practical curves for BFSK
 plot (SNR,FSK_Errors,'r');
 hold on;
 plot (SNR,FSK_theoritical);
 title ('BFSK theoritical vs BFSK real')
 legend ('practical','theoritecal');
 subplot(3,1,2); %plotting theoretical versus practical curves for BPSK
plot(SNR,PSK_Errors,'r');
 hold on;
 plot (SNR,PSK_theoritical);
  title ('BPSK theoritical vs BPSK real')
  legend ('practical','theoritecal');
subplot(3,1,3); %plotting BFSK versus BPSK
 plot(SNR,FSK_Errors,'r')
  hold on;
  plot(SNR,PSK_Errors);
  title('BFSK vs BPSK')
  legend('BFSK','BPSK')