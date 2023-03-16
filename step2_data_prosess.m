%step2 数据打乱，分出训练集和验证集
per_num=1000;
% ii=6000;
% load('data_nr_with_inter.mat')
% XTrain={};
% YTrain={};
% XTest={};
% YTest={};
Spec_Time;
% data_label=data_label';
nums=size(Spec_Time,4);
% for i=1:nums
%     if i<=per_num
%         Data_label(i)=cellstr(['NR ']); %,num2str(data_label(i))
%     elseif i<=2*per_num
%         Data_label(i)=cellstr(['NR with LTE inter ']); %,num2str(data_label(i))
%     elseif i<=3*per_num
%         Data_label(i)=cellstr(['NR with wifi inter ']);%,num2str(data_label(i))
%     else
%         Data_label(i)=cellstr(['NR with bluetooth inter ']);%,num2str(data_label(i))
%     end
% end


Ocr=0.7;%训练数据的比例
num_train=floor(nums*Ocr);
R=randperm(nums); %将数据打乱
XTrain=Spec_Time(:,:,:,R(1:num_train));
YTrain=Data_label(R(1:num_train));
XTest=Spec_Time(:,:,:,R(num_train+1:end));
YTest=Data_label(R(num_train+1:end));


% XTrain=cat(4,Spec_Time(:,:,:,1:0.8*ii),Spec_Time(:,:,:,ii+1:1.8*ii),Spec_Time(:,:,:,2*ii+1:2.8*ii),Spec_Time(:,:,:,3*ii+1:3.8*ii));
% YTrain=[data_label(1:0.8*ii);data_label(ii+1:1.8*ii);data_label(2*ii+1:2.8*ii);data_label(3*ii+1:3.8*ii)];
% XTest=cat(4,Spec_Time(:,:,:,0.8*ii+1:ii),Spec_Time(:,:,:,1.8*ii+1:2*ii),Spec_Time(:,:,:,2.8*ii+1:3*ii),Spec_Time(:,:,:,3.8*ii+1:4*ii));
% YTest=[data_label(0.8*ii+1:ii);data_label(1.8*ii+1:2*ii);data_label(2.8*ii+1:3*ii);data_label(3.8*ii+1:4*ii)];

YTrain= categorical(YTrain);
YTest= categorical(YTest);



