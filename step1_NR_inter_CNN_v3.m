%step 1 产生NR、NR+LTE、NR+wifi、NR+蓝牙干扰信号，并数据预处理，准备送入CNN网络
%CIR=0.5,AWGN从-20db到20db,每隔5个db取一次。
clear;
clc
clear;
N_data=1000;%4000个点
% dlrefwaveform=zeros(N_data,307200);
c=2; %决定cir
nfft=128;
rxSampleRate=3.84e6;
N_t=rxSampleRate*0.01/nfft;%300
% Spec_Time=zeros(nfft,N_t,2,4*N_data);%128*300=38
% 
% \400

% data_label=zeros(4*N_data,1);
load('PSSCH_MCS1.mat')
load('PSSCH_MCS5.mat')
load('PSSCH_MCS13.mat')
load('PSSCH_MCS15.mat')

for SNRdB = -20:5:20%加噪
    tic;
    Spec_Time=zeros(nfft,N_t,2,4*N_data);%128*300=38400
%     data_label=zeros(4*N_data,1);
for i=1:N_data*5
 ran=unidrnd(4);
 switch ran
     case 1
         Waveform_sidelink=out_1.simout;
         sl_rxSampleRate=length(out_1.tout)-1;
     case 2
         Waveform_sidelink=out_5.simout;
         sl_rxSampleRate=length(out_5.tout)-1;
     case 3
         Waveform_sidelink=out_13.simout;
         sl_rxSampleRate=length(out_13.tout)-1;
     case 4
         Waveform_sidelink=out_15.simout;
         sl_rxSampleRate=length(out_15.tout)-1;
 end
 sl_waveform_0=[Waveform_sidelink; Waveform_sidelink; Waveform_sidelink];
 Point_temp=unidrnd(floor(length(sl_waveform_0)*0.06))-1;%产生的sidelink信号中随机取一段信号来用，这是这段信号起点。
 sl_waveform_1=sl_waveform_0(1+Point_temp:floor(sl_rxSampleRate*0.01)+Point_temp);%取10ms的长度信号
 sl_waveform_2=resample(sl_waveform_1,rxSampleRate,sl_rxSampleRate);%将信号下采到3.84Msample/s.那么信号长度是38400
% dlrefwaveinfo.Info.Nfft/ocr
%%
E=sum(sl_waveform_2.*conj(sl_waveform_2))/length(sl_waveform_2);
sl_waveform=sl_waveform_2/sqrt(E);

%%
if i<=N_data*1  %sidelink信号
    
    
%     rxSampleRate=7.68e6;
%     data_label(i)=1;
    rxWaveform=awgn(sl_waveform,SNRdB,-10*log10(double(nfft)));  %加噪声
    Data_label(i)=cellstr(['sidelink']);
    
elseif i<=N_data*2  %sidelink+NR信号   
    %  NR
    fr1testmodels = hNRReferenceWaveformGenerator.FR1TestModels;
    dlnrref=fr1testmodels(ran);% Model name and properties
    
    bws =["10MHz","20MHz","50MHz","40MHz"];
    ran=unidrnd(4);
    bw=bws(ran);% Channel bandwidth
    
    scss =["15kHz","30kHz","60kHz"];
    ran=unidrnd(3);
    scs=scss(ran);% Subcarrier spacing
    
    dm  = "FDD";  % Duplexing mode
    ncellid = unidrnd(1007);  % NCellID
    sv      = "15.2.0";  % TS 38.141-x version (NR-TM only)
    
    % Create generator object for the above NR-TM/PDSCH FRC reference model
    dlrefwavegen = hNRReferenceWaveformGenerator(dlnrref,bw,scs,dm,ncellid,sv);
    % Generate waveform
    [dlrefwaveform0,dlrefwaveinfo,dlresourceinfo] = generateWaveform(dlrefwavegen);
    dlrefwaveform1=[dlrefwaveform0;dlrefwaveform0;dlrefwaveform0];
    Point_temp=unidrnd(length(dlrefwaveform1)*0.6)-1;%产生的5G信号中随机取一段信号来用，这是这段信号起点。
    dlrefwaveform1=dlrefwaveform1(1+Point_temp:dlrefwaveinfo.Info.SampleRate*0.01+Point_temp);%取10ms的长度信号
    ocr=dlrefwaveinfo.Info.SampleRate/rxSampleRate;
    dlrefwaveform=resample(dlrefwaveform1,1,ocr);%将信号下采到3.84Msample/s.那么信号长度是38400
    % dlrefwaveinfo.Info.Nfft/ocr
    %%
%     E=sum(dlrefwaveform.*conj(dlrefwaveform))/length(dlrefwaveform);
%     dlrefwaveform=dlrefwaveform/sqrt(E);
%     rxWaveform0 = awgn(dlrefwaveform,SNRdB,-10*log10(double(nfft)));
%     rxWaveform=rxWaveform0;
    
     % sidelink+NR干扰信号相加
     sl_nr_waveform=sl_waveform+c*dlrefwaveform;
     E=sum(sl_nr_waveform.*conj(sl_nr_waveform))/length(sl_nr_waveform);
     sl_nr_waveform=sl_nr_waveform/sqrt(E);
     
     rxWaveform = awgn(sl_nr_waveform,SNRdB,-10*log10(double(nfft)));  %加噪声
%      data_label(i)=2;
     Data_label(i)=cellstr(['sidelink with NR inter']);
    
    
  
elseif i<=N_data*3  %sidelink+LTE信号
    %  LTE
     cfg = struct('RC', ['R.',num2str(unidrnd(8))], ...
        'DuplexMode', 'FDD', ...
        'NCellID', unidrnd(504-1), ...
        'TotSubframes', 30, ...
        'NumCodewords', 1, ...
        'Windowing', 0, ...
        'AntennaPort', 1);
     cfg.OCNGPDSCHEnable = 'On';
     cfg.OCNGPDCCHEnable = 'On';
     cfg.PDSCH.TxScheme = 'Port0';
     cfg.PDSCH.RNTI = 1;
     cfg.PDSCH.Rho = 0;
     cfg.PDSCH.RVSeq = [0 1 2 3];
     cfg.PDSCH.NHARQProcesses = 7;
     cfg.PDSCH.PMISet = 1;
     cfg = lteRMCDL(cfg);
     % input bit source:
     in = [1; 0; 0; 1];
     % waveform generation:
     [waveform0, grid, cfg] = lteRMCDLTool(cfg, in);
     Point_temp=unidrnd(length(waveform0)*0.6);%产生的信号中随机取一段信号来用，这是这段信号起点。
     waveform_lte0=waveform0(1+Point_temp:cfg.SamplingRate*0.01+Point_temp);%取10ms的长度信号
     ocr=cfg.SamplingRate/rxSampleRate;
     waveform_lte1=resample(waveform_lte0,rxSampleRate,cfg.SamplingRate);%将信号下采到3.84Msample/s.那么信号长度是38400
     
     E=sum(waveform_lte1.*conj(waveform_lte1))/length(waveform_lte1);
     waveform_lte1=waveform_lte1/sqrt(E); 
     % sidelink+LTE干扰信号相加
     sl_lte_waveform=sl_waveform+c*waveform_lte1;
     E=sum(sl_lte_waveform.*conj(sl_lte_waveform))/length(sl_lte_waveform);
     sl_lte_waveform=sl_lte_waveform/sqrt(E);
     
     rxWaveform = awgn(sl_lte_waveform,SNRdB,-10*log10(double(nfft)));  %加噪声
%      data_label(i)=2;
     Data_label(i)=cellstr(['sidelink with LTE inter']);
    
elseif i<=N_data*4  %sidelink+wifi信号
    % wifi
    temp=unidrnd(3);
    if(temp)==1
        cbw=20;
    elseif(temp)==2
        cbw=40;
    else
        cbw=80; 
    end
     vhtCfg = wlanVHTConfig('ChannelBandwidth',['CBW',num2str(cbw)], ...
    'NumUsers', 1, ...
    'NumTransmitAntennas', 1, ...
    'NumSpaceTimeStreams', [1], ...
    'SpatialMapping', 'Direct', ...
    'STBC', false, ...
    'MCS', unidrnd(7), ...
    'ChannelCoding', 'BCC', ...
    'APEPLength', 1024, ...
    'GuardInterval', 'Long', ...
    'GroupID', 63, ...
    'PartialAID', 275);
    numPackets = 100;
    % input bit source:
    in = randi([0, 1], 1000, 1);
    % waveform generation:子载波间隔为87.125khz
    waveform0 = wlanWaveformGenerator(in, vhtCfg, 'NumPackets', numPackets, 'IdleTime', 0, 'ScramblerInitialization', 93, 'WindowTransitionTime', 1e-07); 
    Fs_wlan = wlanSampleRate(vhtCfg); 	
    wave_wifi=resample(waveform0,rxSampleRate,Fs_wlan);%将信号下采到3.42Msample/s
    
    Point_temp=unidrnd(rxSampleRate*0.015)-1;%产生的信号中随机取一段信号来用，这是这段信号起点。
    waveform_wifi1=waveform0(1+Point_temp:rxSampleRate*0.01+Point_temp);%取10ms的长度信号
    
%     waveform_wifi1=resample(waveform_wifi0,rxSampleRate,Fs_wlan);%将信号下采到7.68Msample/s.那么信号长度是76800
    
    E=sum(waveform_wifi1.*conj(waveform_wifi1))/length(waveform_wifi1);
    waveform_wifi1=waveform_wifi1/sqrt(E);
    %sidelink+wifi
    sl_wifi_waveform=sl_waveform+c*waveform_wifi1;
    E=sum(sl_wifi_waveform.*conj(sl_wifi_waveform))/length(sl_wifi_waveform);
   sl_wifi_waveform=sl_wifi_waveform/sqrt(E);
%     SNRdB = randi(22);%加点噪声
    rxWaveform = awgn(sl_wifi_waveform,SNRdB,-10*log10(nfft)); 
%     rxSampleRate=wlanSampleRate(vhtCfg);
%     data_label(i)=3;
    Data_label(i)=cellstr(['sidelink with wifi inter']);  
else   %sidelink+蓝牙
    %蓝牙
    load('wave_blue.mat');
    waveform_blue=resample(waveform_blue,1,2);%信号下采，变成3.84m sample/s
    waveform_blue1=[waveform_blue;waveform_blue;waveform_blue];
    Point_temp=unidrnd(2*length(waveform_blue)-2);%产生的信号中随机取一段信号来用，这是这段信号起点。
    waveform_blue2=waveform_blue1(1+Point_temp:length(waveform_blue)+Point_temp);%取下采后的10ms
    
    E=sum(waveform_blue2.*conj(waveform_blue2))/length(waveform_blue2);
    waveform_blue2=waveform_blue2/sqrt(E);
    %sidelink+蓝牙
    sl_blue_waveform=sl_waveform+c*waveform_blue2;
    E=sum(sl_blue_waveform.*conj(sl_blue_waveform))/length(sl_blue_waveform);
    sl_blue_waveform=sl_blue_waveform/sqrt(E);
    
%     SNRdB = randi(22);%加点噪声
    rxWaveform = awgn(sl_blue_waveform,SNRdB,-10*log10(nfft)); 
%     data_label(i)=4;
    Data_label(i)=cellstr(['sidelink with bluetooth inter']);
end

spec_time0=spectrogram(rxWaveform,ones(nfft,1),0,nfft,'centered',rxSampleRate,'yaxis','MinThreshold',-130);
spec_time_real=real(spec_time0);
spec_time_imag=imag(spec_time0);
spec_time1(:,:,1)=spec_time_real;
spec_time1(:,:,2)=spec_time_imag;
    
Spec_Time(:,:,:,i)=spec_time1;
end
% Spec_Time2=Spec_Time;
% data_label2=data_label;
% Data_label2=Data_label;

% save('Data_nr_with_inter_12-13-300.mat','Spec_Time','Data_label');
save(['CIR_-6_Data_sidelink_with_inter_',num2str(SNRdB),'dB'],'Spec_Time','Data_label');
clear Spec_Time Data_label;
toc;
disp(['SNR',num2str(SNRdB),'dB ## data gen done!'])
end