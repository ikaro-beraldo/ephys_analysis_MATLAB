function sleep_parameters = sleep_parameters(All_Sort)

%% Sleep parameters - Rafaela

e=size(All_Sort,1);

% - sleep efficiency (percentage of total sleep time during the recording time)

sleep_parameters.Sleep_E=((sum(All_Sort==1) + sum (All_Sort==2))/size(All_Sort,1))*100;


% - latency to sleep (time lag between the onset of recording and the first sleep period)

threshold_sleep=12; % 3 consecutive windows of 10s

aux_latency1=zeros(size(All_Sort));

aux_latency1(All_Sort==2)=1;


 % Light phase

S=find(diff([0;aux_latency1])==1);
E=find(diff([0;aux_latency1])==-1);
    
  
    if size(S,1)==size(E,1)
        
        W=[S,E,E-S];
        
    else 
        
        W=[S(1:end-1),E,E-S(1:end-1)];
        
    end
   
 sleep1=W((W(:,3)>threshold_sleep),1);
 sleep_parameters.latency_sleep1=(sleep1(1)*10); % in seconds
 
 % - Latency to the first REM episode (time lag between the first sleep period and the first episode of REM);
 
 aux_REM=find(All_Sort==1);
 
 aux2=(find(aux_REM>sleep1(1))); % in seconds
 
 sleep_parameters.latency_REM1=((aux_REM(aux2(1))-sleep1(1))*10); % in seconds

% sleep-wake state transitions

% NREM_WK

clear aux1
aux1=All_Sort;
aux1(All_Sort==1)=NaN;

t1=sum(diff(aux1)==1); %  NREM to WK-A

% REM_WK

aux2=All_Sort;
aux2(All_Sort==2)=NaN;

t2=sum(diff(aux2)==2); %  NREM_WK to WK-Q

sleep_parameters.NREM_wake=((t1/size(All_Sort,1))*6)*60; % Transitions/hour 
sleep_parameters.REM_wake=((t2/size(All_Sort,1))*6)*60; % Transitions/hour 
sleep_parameters.Sleep_wake_Total=(((t1+t2)/size(All_Sort,1))*6)*60; % Transitions/hour 


