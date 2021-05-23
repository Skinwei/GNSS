%% function : aquire gps signal
%% date : 20210513
%% author: skingwei
%%%%%%%%%%%%%%55555555555555555555555555%%%%%%%%%
% PRN Doppler
% 21 3740
%%%%%%%%%5555555555555555555555555%%%%%%%%%%%%5555
clear
close all

filepath = 'OutCACodeResult/';
filename1 = 'GPSCACode.txt';
%% parameter
Fs = 16.368E6;  %采样率
Fc = 4.092E6;   %理论中频数据

Fcc = 1.023E6;    %码率
GPS_CACODE_NUM = 1023;
initDopp = -500;  %Hz
slew_count = Fs/Fcc; %%本地码率倍数
FFT_N = 4092;
FFT_Len = Fs/1E3;

%% read data
gpsIFdata = importdata('gps_data_20ms.mat');
gpsIFdata = gpsIFdata';

f1 = [filepath,filename1];
d1= importdata(f1); 
%%
L = length(gpsIFdata);
temp_Coh_PRN = zeros(FFT_Len*2-1,41); %ALLGPS 
FFTCode_Signal = zeros(FFT_Len,1);
FFTCode_CA = zeros(FFT_Len,1);
IFData = complex(gpsIFdata(:,1),0);

global reg_aquired regcount
reg_aquired = zeros(200,4);
regcount = 1;
%%
aqired_PRN_Rsult = cell(32,2);  %%corr dopp quiredresult
%% 并行码捕获
fcount = 1;
initdopp_a = (-40:40)*100;   %%多普勒范围

tic
for j = 1:32
    CA_Code = d1(:,j);
    codeValueIndex = ceil(Fcc*(1:1023*Fs/Fcc)/Fs);
    CA_Code_Resample = CA_Code(codeValueIndex); % resample C/A code
    FFTCode_CA = complex(CA_Code_Resample,0);
tic
for m = -40:40   %%-4000Hz~4000Hz

%init parameter
    FFT_count = 1;
    FFT_Cohcount = 1;
    initdopp = m*100;
    temp_Coh1 = zeros(FFT_Len*2-1,1);
    temp_Coh2 = zeros(FFT_Len*2-1,1);
    
    %%NCO
    fncos = initdopp + Fc;
    dt = 1/Fs;
    t= (0:L-1)*dt;
    sin_x = sin(2*pi*fncos*t);
    cos_x = cos(2*pi*fncos*t);
    x_real = IFData.*sin_x';
    x_imag = IFData.*cos_x';
    
    %%20ms
    for i = 1:20   
        temp_index = 1+(i-1)*16.368E3:16.368E3*i;
        FFTCode_Signal = complex(x_real(temp_index),x_imag(temp_index));  
       %归一化结果，并行码相关       
       [tempFFTCode_outVec,x] =  xcorr(FFTCode_Signal,FFTCode_CA,'coeff');       
       FFT_Cohcount = FFT_Cohcount + 1; 
       %相干积分
       if FFT_Cohcount >= 11
           temp_Coh2 = temp_Coh2 + tempFFTCode_outVec;
       else
           temp_Coh1 = temp_Coh1 + tempFFTCode_outVec;
       end
%        break;
    end
    
    %%捕获
    [aquiredresult,temp_Coh]= aquire_Prn(initdopp,temp_Coh1,temp_Coh2,x,j);  %% PRN Doppp codephase FFT-SNR
    
    %%存储
    temp_Coh_PRN(:,m+41) = temp_Coh;
  
end
if FFT_Cohcount >= 21
       disp('处理结束')
       disp(j)         
end
toc
% result
temp_Coh_PRNm = abs(temp_Coh_PRN);
%%%%%%%%%%%%%%%%%%%%%%%%% result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% figure
% figure(fcount)
% fcount= fcount + 1;
% plot(x/16,temp_Coh_PRNm,'DisplayName','temp_Coh_PRNm')
% xlabel('chip')
% ylabel('Amplitude')
% titlename = ['GPS',num2str(j),'-','Result'];
% title(titlename)
% grid on


%%save result
aqired_PRN_Rsult{j,1} = temp_Coh_PRNm;
aqired_PRN_Rsult{j,2} = initdopp_a;

%三维图
h = figure(fcount);
fcount= fcount + 1;
mesh(initdopp_a,x,aqired_PRN_Rsult{j,1})
titlename = ['GPS',num2str(j),'-','Result'];
title(titlename)
xlabel('dopp/Hz')
ylabel('codephase/码片')
zlabel('Amplitude')
savefig(h,titlename);
close(h)

end
toc
save aqired_PRN_Rsult.mat aqired_PRN_Rsult
save aquiredresult.mat aquiredresult
%%
% 由于太卡了就一个一个画吧 X= m Y = n  Z = n*m
% mesh(initdopp_a,x,aqired_PRN_Rsult{9})
% title('GPS-Result')
% xlabel('dopp/Hz')
% ylabel('codephase/码片')
% zlabel('Amplitude')
%%
%% 捕获结果
function [aquiredresult,temp_Coh]= aquire_Prn(dopp,temp_Coh1,temp_Coh2,x,j)

 global reg_aquired regcount
 thd_quired = 6.0; %%捕获阈值
 %%result 1
 temp = abs(temp_Coh1);
 [~,index_max1] = max(temp,[],1);  %%按列
 sum1 = sum(temp);
 sum2 = sum(temp(index_max1-15:index_max1+15));
 FFT_SNR1 = temp(index_max1)*(32735-31)/(sum1-sum2);
 
 %%result 2
 temp = abs(temp_Coh2);
 [~,index_max2] = max(temp,[],1);  %%按列
 sum1 = sum(temp);
 sum2 = sum(temp(index_max2-15:index_max2+15));
 FFT_SNR2 = temp(index_max2)*(32735-31)/(sum1-sum2);
 
 if FFT_SNR1 > FFT_SNR2
     FFT_SNR = FFT_SNR1;
     temp_Coh = temp_Coh1;
     c = index_max1;
 else
     FFT_SNR = FFT_SNR2;
     temp_Coh = temp_Coh2;
     c = index_max2;
 end
 
 %%捕获判断
 if FFT_SNR > thd_quired
    reg_aquired(regcount,:) = [j,dopp,x(c)/16,FFT_SNR];
    regcount = regcount+1;
 end

 aquiredresult = reg_aquired;

end






