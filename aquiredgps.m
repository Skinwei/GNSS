%% function : aquire gps signal
%% date : 20210507
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
temp_Coh_PRN1 = zeros(FFT_Len*2-1,32); %ALLGPS
temp_Coh_PRN2 = zeros(FFT_Len*2-1,32); %ALLGPS
FFTCode_Signal = zeros(FFT_Len,1);
FFTCode_CA = zeros(FFT_Len,1);
IFData = complex(gpsIFdata(:,1),0);
global reg_aquired regcount
reg_aquired = zeros(200,4);
regcount = 1;
%%

%% 并行码捕获
fcount = 1;
for m = 16   %%-4000Hz~4000Hz
for j = 14
% if j== 21 || j==15 || j==27 || j==9
%         continue
% end
% init parameter
tic
    FFT_count = 1;
    FFT_Cohcount = 1;
    CarrierLoop(Fs,Fc);
    initDopp = m*100;
    Set_CarrierFreq(initDopp,Fc,Fs);
    CodeLoop(j,d1);     %GPS 1-32
    
    temp_Coh1 = zeros(FFT_Len*2-1,1);
    temp_Coh2 = zeros(FFT_Len*2-1,1);
    for i = 1:L

    %载波剥离,暂时不考虑多普勒
    ComplexData = ComplexDDC(IFData(i));
    %%
    %码剥离,暂时不考虑多普勒
    [tempchar,Code_Signal, Code_CA] = Code_Correlate(ComplexData,slew_count);
    if tempchar == 1
        FFTCode_Signal(FFT_count,1) = Code_Signal;
        FFTCode_CA(FFT_count,1) = complex(Code_CA,0);
    
        if FFT_count < FFT_Len
            FFT_count= FFT_count +1;
        else 
           FFT_count = 1;
    %        tic
           [tempFFTCode_outVec,x] =  xcorr(FFTCode_Signal,FFTCode_CA,'coeff'); %归一化结果，并行码相关
    %        toc 
           %%相干积分
%            temp_Coh = temp_Coh + tempFFTCode_outVec;
%             temp_Coh = tempFFTCode_outVec;

           FFT_Cohcount = FFT_Cohcount + 1;
           if FFT_Cohcount >= 11
               temp_Coh2 = temp_Coh2 + tempFFTCode_outVec;
               temp_Coh1 = temp_Coh1 + tempFFTCode_outVec;
           else
               temp_Coh1 = temp_Coh1 + tempFFTCode_outVec;
           end
    %        break;
            if FFT_Cohcount >= 21
                disp('处理结束')
                disp(j)
            end
            
        end    
    end
    end
    temp_Coh_PRN1(:,j) = temp_Coh1;
    temp_Coh_PRN2(:,j) = temp_Coh2;
toc
end
% result
temp_Coh_PRNm1 = abs(temp_Coh_PRN1);
temp_Coh_PRNm2 = abs(temp_Coh_PRN2);
% search
aquiredresult = aquire_Prn(initDopp,temp_Coh_PRNm1,x);  %% PRN Doppp codephase FFT-SNR

% figure
figure(fcount)
fcount= fcount + 1;
plot(x/16,temp_Coh_PRNm1,'DisplayName','temp_Coh_PRNm')
xlabel('chip')
ylabel('Amplitude')
titlename = ['GPS',num2str(j),'-','Result'];
title(titlename)
disp("inidopp")
disp(initDopp)
grid on

aquiredresult = aquire_Prn(initDopp,temp_Coh_PRNm2,x);  %% PRN Doppp codephase FFT-SNR
% figure
figure(fcount)
fcount= fcount + 1;
plot(x/16,temp_Coh_PRNm2,'DisplayName','temp_Coh_PRNm')
xlabel('chip')
ylabel('Amplitude')
titlename = ['GPS',num2str(j),'-','Result'];
title(titlename)
disp("inidopp")
disp(initDopp)
grid on
end

%% func: 载波剥离
%% return : 正交剥离后的结果
%%//-----------------------中频数据下变频----------------------------//
function Complexdata = ComplexDDC(IFData)
    global costable sintable
    global  Carrier_phase Carrier_phaseStep
	%%Multiply 小心这里没处理好
	Carrier_phase = Carrier_phase + Carrier_phaseStep;			%相位步进
    if Carrier_phase >= pow2(32)
       Carrier_phase = mod(Carrier_phase,pow2(32));
    end
	phaseIndex = bitshift(Carrier_phase,-27,'uint32');
    phaseIndex = mod(phaseIndex,32)+1;
    %index is different to cpp
	TempComplexData = complex(costable(phaseIndex), sintable(phaseIndex));						%%NCO
    
	TempComplexData = IFData * TempComplexData;         %载波剥离
	TempComplexData = complex(bitshift(real(TempComplexData),-4,'int64'), bitshift(imag(TempComplexData),-4,'int64'));
   
    %%return
    Complexdata = TempComplexData;

end


%% //-----------------------载波频率设定----------------------------//
function Set_CarrierFreq(DopplerFreq,GPS_FREQ_COMPENSATE,SystemClockFre)
    global  Carrier_phaseStep
	double tempData;
	tempData = DopplerFreq + GPS_FREQ_COMPENSATE;
	tempData = pow2(32) * tempData / SystemClockFre;
	Carrier_phaseStep = floor(tempData);
end


%% //-----------------------初始 载波频率设定----------------------------//
function CarrierLoop(SystemClockFre,GPS_FREQ_COMPENSATE)
    global costable sintable
    global  Carrier_phase Carrier_phaseStep
    %% 正弦波表
    costable = [252, 247, 233, 210, 178, 140, 96, 49, 0, -49, -96, -140, -178, -210, -233, -247, -252, -247, -233, -210, -178, -140, -96, -49, 0, 49, 96, 140, 178, 210, 233, 247 ];
    sintable = [ 0, 49, 96, 140, 178, 210, 233, 247, 252, 247, 233, 210, 178, 140, 96, 49, 0, -49, -96, -140, -178, -210, -233, -247, -252, -247, -233, -210, -178, -140, -96, -49 ];

	double tempData;
	Carrier_phase = 0;
	%% 计算标称中频所对应的M
	tempData = GPS_FREQ_COMPENSATE;
	tempData = pow2(32) * tempData / SystemClockFre;
	Carrier_phaseStep = round(tempData);
    
end

%% func: 码剥离，码片重采样
%% return : 码剥离后的结果
%% 
%% result 重采样结果
function [temp_char,FFTCode_Signal, FFTCode_CA]= Code_Correlate(DDCData,slew_count)
    global Code_phaseIndex
    global CA_code count_cc 
    global div_count
    div_count = div_count + 1;
    if div_count>0  %重采样 1
        div_count =0;
        FFTCode_Signal = DDCData;
        %%重采样
        FFTCode_CA = CA_code(Code_phaseIndex);
        if count_cc < slew_count
             count_cc =  count_cc + 1;
        else
             count_cc = 1;
             Code_phaseIndex = Code_phaseIndex+1;        
             if Code_phaseIndex >1023
                 Code_phaseIndex = 1;
             end

        end
        temp_char = 1;
    else
        temp_char = 0;   
        FFTCode_Signal = 0;
        FFTCode_CA = 0;
    end

end


%% 初始设定
function CodeLoop(PRN,CA)
    global CurrentCodePhase   Code_phaseIndex 
    global CA_code count_cc
    global div_count
    div_count =0;
    CurrentCodePhase = 0;
    Code_phaseIndex = 1;
    count_cc = 1;

    CA_code = CA(:,PRN);
%     temp = CA(:,PRN);
%     temp2 = temp(122:end);
%     temp3 = temp(1:121);
%     CA_code = [temp2 ;temp3];
    
end

%% 捕获结果
function aquiredresult = aquire_Prn(dopp,recorr,x)

 global reg_aquired regcount
 b= max(recorr,[],1);  %%按列

 index = find(b>0.2);  %%捕获阈值
 if isempty(index)
 else
    L =length(index);
    PRN = index;
    initdopp = dopp;
    for i = 1:L
        temp = recorr(:,PRN(i));
        c = find(recorr(:,PRN(i))==b(PRN(i)));
        sum1 = sum(temp);
        sum2 = sum(temp(c-15:c+15));
        
        snr = temp(c)*(32735-31)/(sum1-sum2);
        
        reg_aquired(regcount,:) = [PRN(i),initdopp,x(c)/16,snr];
        regcount = regcount+1;
    end
    
 end
 aquiredresult = reg_aquired;

end






