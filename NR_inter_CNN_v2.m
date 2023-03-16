%step 1 产生NR、NR+LTE、NR+wifi、NR+蓝牙干扰信号，并数据预处理，准备送入CNN网络
%CIR=1,AWGN从-20db到20db,每隔5个db取一次。
clear;
N_data=4000;
% dlrefwaveform=zeros(N_data,307200);
nfft=128;
rxSampleRate=3.84e6;
N_t=rxSampleRate*0.01/nfft;%300
% Spec_Time=zeros(nfft,N_t,2,4*N_data);%128*300=38400
% data_label=zeros(4*N_data,1);

for SNRdB = -20:5:20%加点噪声
    tic;
    Spec_Time=zeros(nfft,N_t,2,4*N_data);%128*300=38400
%     data_label=zeros(4*N_data,1);
for i=1:N_data*4
 ran=unidrnd(8);
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
% dlrefwaveform=resample(dlrefwaveform1,1,ocr);%将信号下采到3.84Msample/s.那么信号长度是38400
dlrefwaveform=resample(dlrefwaveform1,rxSampleRate,dlrefwaveinfo.Info.SampleRate);%将信号下采到3.84Msample/s.那么信号长度是38400
% dlrefwaveinfo.Info.Nfft/ocr
%%
E=sum(dlrefwaveform.*conj(dlrefwaveform))/length(dlrefwaveform);
dlrefwaveform=dlrefwaveform/sqrt(E);

%%
if i<=N_data*1  %NR信号
    
    rxWaveform0 = awgn(dlrefwaveform,SNRdB,-10*log10(double(nfft)));
    rxWaveform=rxWaveform0;
%     rxSampleRate=7.68e6;
%     data_label(i)=1;
    Data_label(i)=cellstr(['NR']); %,num2str(data_label(i))
elseif i<=N_data*2  %NR+LTE信号
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
     % NR+LTE干扰信号相加
     nr_lte_waveform=dlrefwaveform+waveform_lte1;
     E=sum(nr_lte_waveform.*conj(nr_lte_waveform))/length(nr_lte_waveform);
     nr_lte_waveform=nr_lte_waveform/sqrt(E);
     
     rxWaveform = awgn(nr_lte_waveform,SNRdB,-10*log10(double(nfft)));  %加噪声
%      data_label(i)=2;
     Data_label(i)=cellstr(['NR with LTE inter']);
    
elseif i<=N_data*3  %NR+wifi信号
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
    % waveform generation:
    waveform0 = wlanWaveformGenerator(in, vhtCfg, 'NumPackets', numPackets, 'IdleTime', 0, 'ScramblerInitialization', 93, 'WindowTransitionTime', 1e-07); 
    Fs_wlan = wlanSampleRate(vhtCfg); 	
    wave_wifi=resample(waveform0,rxSampleRate,Fs_wlan);%将信号下采到3.42Msample/s
    
    Point_temp=unidrnd(rxSampleRate*0.015)-1;%产生的信号中随机取一段信号来用，这是这段信号起点。
    waveform_wifi1=waveform0(1+Point_temp:rxSampleRate*0.01+Point_temp);%取10ms的长度信号
%     waveform_wifi1=resample(waveform_wifi0,rxSampleRate,Fs_wlan);%将信号下采到7.68Msample/s.那么信号长度是76800
    
    E=sum(waveform_wifi1.*conj(waveform_wifi1))/length(waveform_wifi1);
    waveform_wifi1=waveform_wifi1/sqrt(E);
    %NR+wifi
    nr_wifi_waveform=dlrefwaveform+waveform_wifi1;
    E=sum(nr_wifi_waveform.*conj(nr_wifi_waveform))/length(nr_wifi_waveform);
    nr_wifi_waveform=nr_wifi_waveform/sqrt(E);
%     SNRdB = randi(22);%加点噪声
    rxWaveform = awgn(nr_wifi_waveform,SNRdB,-10*log10(nfft)); 
%     rxSampleRate=wlanSampleRate(vhtCfg);
%     data_label(i)=3;
    Data_label(i)=cellstr(['NR with wifi inter']);  
else   %NR+蓝牙
    %蓝牙
    load('wave_blue.mat');
    waveform_blue=resample(waveform_blue,1,2);%信号下采，变成3.84m sample/s
    waveform_blue1=[waveform_blue;waveform_blue;waveform_blue];
    Point_temp=unidrnd(2*length(waveform_blue)-2);%产生的信号中随机取一段信号来用，这是这段信号起点。
    waveform_blue2=waveform_blue1(1+Point_temp:length(waveform_blue)+Point_temp);%取下采后的10ms
    
    E=sum(waveform_blue2.*conj(waveform_blue2))/length(waveform_blue2);
    waveform_blue2=waveform_blue2/sqrt(E);
    %NR+蓝牙
    nr_blue_waveform=dlrefwaveform+waveform_blue2;
    E=sum(nr_blue_waveform.*conj(nr_blue_waveform))/length(nr_blue_waveform);
    nr_blue_waveform=nr_blue_waveform/sqrt(E);
    
%     SNRdB = randi(22);%加点噪声
    rxWaveform = awgn(nr_blue_waveform,SNRdB,-10*log10(nfft)); 
%     data_label(i)=4;
    Data_label(i)=cellstr(['NR with bluetooth inter']);
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
save(['Data_nr_with_inter_',num2str(SNRdB),'dB'],'Spec_Time','Data_label');
clear Spec_Time Data_label;
toc;
disp(['SNR',num2str(SNRdB),'dB ## data gen done!'])
end