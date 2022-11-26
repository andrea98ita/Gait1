clear all; close all; clc;

%% Load Data

[FileName,PathName] = uigetfile( 'D:\PhD - Backup\Analysis\Protocollo Pre FOG\*.csv');
filename = fullfile(PathName, FileName);
data = readtable(filename);
patient=data.patient(1,1);
task=data.task(1,1);
sensor = data.mac; sensorData = {[],[],[],[]};
sensorName = {'Right ankle', 'Left ankle', 'Lower back', 'Wrist'}; 
sensorList = {'CD:98:11:6E:DB:A3','F9:C2:5B:8A:BB:96','C9:CE:E1:82:70:C3','E8:43:CA:7C:8C:39'};
for i = 1:length(sensor)
   sId = find(contains(sensorList,sensor{i}));
   timeStamp = datetime(table2array(data(i,1)),'InputFormat','dd-MM-yyyy HH:mm:ss.S', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'); d = timeofday(timeStamp); timeStamp = seconds(d);
   sensorData{sId} = [sensorData{sId}; [table2array(data(i,2:10)) timeStamp]];
end
%  taskp{1}=[1 2 3 4;3 4 5 7];
%  taskp{2}=[ 2 2 2;3 3 3 ];
%  pz.pz1=taskp
%  paz1.task=swingtime_dx %QUESTO CHE HO FATTO E' LA SVOLTA!!
if patient==9004082173 || patient==9002208059 || patient==9001981881 ||patient==9001673389 ||patient==9000608670 ||patient==9000565268
    temp=sensorData{3}(:,1);
    temp2=sensorData{3}(:,4);
    sensorData{3}(:,1)=sensorData{3}(:,2);
    sensorData{3}(:,2)=temp;
    sensorData{3}(:,4)=sensorData{3}(:,5);
    sensorData{3}(:,5)=temp2;
end
FUSE=imufilter('SampleRate',60);
[orientation,angularVelocity] = FUSE(sensorData{1}(:,1:3),sensorData{1}(:,4:6))
figure()
plot(sensorData{3}(:,1:3))
title('Acc tronco');

trunkfilt=[];
[b,a] = butter(3,[0.5 2]/30,'bandpass');  %% rivedere questo filtraggio
trunkfilt(:,1)=filtfilt(b, a, sensorData{3}(:,3));
figure()
plot(trunkfilt(:,1))
[b,a] = butter(3,[0.5 2]/30,'bandpass');  %% rivedere questo filtraggio
trunkfilt(:,2)=filtfilt(b, a, sensorData{3}(:,1));
figure()
plot(trunkfilt(:,2))
 
y=skewness(sensorData{1}(:,6));
y2=skewness(sensorData{2}(:,6));
figure()
plot(sensorData{1}(:,6))
figure()
plot(sensorData{2}(:,6))

figure()
plot(sensorData{3}(:,4:6))
title('Velang tronco');
figure()
plot(sensorData{1}(:,1:3))
title('Acc caviglia');

shankfilt=[];
[b1,a1] = butter(5,10/(60/2),'low');
shankfilt(:,1:3)=filtfilt(b1, a1, sensorData{1}(:,4:6));
[b,a]=butter(2,2.3/30,'low');
shankfilt(:,3)=filtfilt(b,a,shankfilt(:,3));
figure()
plot(-shankfilt(:,3));


figure()
plot(-sensorData{1}(:,6))
hold on
plot(sensorData{2}(:,6))
title('Vel ang caviglia')

shankfilt=[];
[b1,a1] = butter(5,1/(60/2),'low');
shankfilt(:,1)=filtfilt(b1, a1, sensorData{2}(:,4));
[b,a] = butter(5,14/(60/2),'low');
shankfilt(:,2:3)=filtfilt(b, a, sensorData{2}(:,5:6));


shankfiltdx=[];
[b1,a1] = butter(5,1/(60/2),'low');
shankfiltdx(:,1)=filtfilt(b1, a1, sensorData{1}(:,4));
[b,a] = butter(5,14/(60/2),'low');
shankfiltdx(:,2:3)=filtfilt(b, a, sensorData{1}(:,5:6));



figure()
plot(shankfilt(:,3))
hold on
plot(-shankfiltdx(:,3));


figure()
plot(shankfilt(:,3))
hold on
plot(shankfilt(:,1))

figure()
plot(-shankfiltdx(:,3))
hold on
plot(shankfiltdx(:,1))

figure()
plot(sensorData{1}(:,4));

% [b,a] = butter(2,0.3/30,'low');
% trunkturn=filtfilt(b, a,sensorData{3}(:,4:6));
% figure()
% plot((abs(trunkturn(:,1))));
% hold on
% plot((abs(trunkturn(:,3))));
% figure()
% plot(((trunkturn(:,1))));
% hold on
% plot(((trunkturn(:,3))));
% 
% turnindin=[];
% turnindfi=[];
% [mturntr,indturntr]=findpeaks(abs(trunkturn(:,1)),'MinPeakHeight',std(abs(trunkturn(:,1))));
% for i=1:length(mturntr)
%     flag=1;
%     trovato=0;
% for k=1:length(trunkturn)
%     if abs(k-indturntr(i))<200
%     if abs(trunkturn(k,1))>(0.5*mturntr(i)) && flag ~=0
%         turnindin=[turnindin,k];
%         flag=0;
%         trovato=1;
%     end
%    
%     if abs(trunkturn(k,1))<0.5*mturntr(i) && trovato==1
%         turnindfi=[turnindfi,k];
%         break
%     end
%     end
% end
% end
%         
% figure()
% findpeaks(abs(trunkturn(:,1)),'MinPeakHeight',std(abs(trunkturn(:,1))));


%% Plot signals
figure,
for i = 1:length(sensorData)
    % Right ankle - angular velocity
    if i == 1
        subplot(4,1,i), plot(-sensorData{i}(:,4:6)), axis tight, title(sensorName{i})
    % Left ankle - angular velocity
    elseif i == 2
        subplot(4,1,i), plot(sensorData{i}(:,4:6)), axis tight, title(sensorName{i})
    % Lower back - acceleration
    elseif i == 3
        subplot(4,1,i), plot(sensorData{i}(:,1:3)-mean(sensorData{i}(:,1:3))), axis tight, title(sensorName{i})
    % Wrist - acceleration
    else
        subplot(4,1,i), plot(sensorData{i}(:,1:3)-mean(sensorData{i}(:,1:3))), axis tight, title(sensorName{i})
    end
end
%% VELOCITA ANGOLARI CAVIGLIA
figure()
plot(-sensorData{1}(:,6)); title('Angular velocity di ',sensorName{1});
figure()
plot(sensorData{2}(:,6));title('Angular velocity di ',sensorName{2});


figure()
plot(sensorData{3}(:,1)); title('Accelerazione x di ',sensorName{3});
figure()
plot(sensorData{3}(:,2));title('Accelerazione y di ',sensorName{3});
figure()
plot(-sensorData{3}(:,3)); title('Acclerazione z di ',sensorName{3});



%figure, plot(sensorData{1}(:,end)-sensorData{1}(1,end),-sensorData{1}(:,6),'linewidth',1.2), hold on, plot(sensorData{2}(:,end)-sensorData{2}(1,end),sensorData{2}(:,6),'linewidth',1.2), hold on, plot(sensorData{3}(:,end)-sensorData{3}(1,end),-100*normalize(sensorData{3}(:,1)),'linewidth',1.2),...
  %legend('R Ankle', 'L Ankle','Back'), axis tight, xlabel('time (s)'), set(gca,'fontsize',16)



%% Filtraggio velocità angolari caviglia SX
shankfilt=[];
[b1,a1] = butter(5,1/(60/2),'low');
shankfilt(:,1)=filtfilt(b1, a1, sensorData{2}(:,4));
[b,a] = butter(5,4/(60/2),'low');
shankfilt(:,2:3)=filtfilt(b, a, sensorData{2}(:,5:6));
% figure()
% findpeaks(-shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))))
% figure()
% findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))))


%% CAVIGLIAVELANG FINDPEAKS
% figure()
% findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))));
% title('MASSIMI di VelAngZ CAVIGLIA SX');
[m,ind]=findpeaks(-shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
[m2,ind2]=findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
med2=mean(m2);
med1=mean(m);
if med2>med1
 [m,ind]=findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);   
 figure()
 findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
 title('MASSIMI di VelAngZ CAVIGLIA SX');
 figure()
 findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
 title('Capire quando turn per Caviglia SX');
 hold on
 findpeaks(abs(shankfilt(:,1)),'MinPeakHeight',std(abs(shankfilt(:,1))));
 figure()
 findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
 hold on 
 plot(abs(shankfilt(:,1)));
 title('x2Capire quando turn per Caviglia SX ');
 [m1,ind1]=findpeaks(-shankfilt(:,3));
 figure()
 findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
 hold on 
 plot(shankfilt(:,1));
 title('Capire rotazioni')
else
  [m,ind]=findpeaks(-shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);  
  figure()
  findpeaks(-shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
  title('MASSIMI di VelAngZ CAVIGLIA SX');
  figure()
  findpeaks(-shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
  title('Capire quando turn per Caviglia SX');
  hold on
  findpeaks(abs(shankfilt(:,1)),'MinPeakHeight',std(abs(shankfilt(:,1))));
  figure()
  findpeaks(-shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
  hold on 
  plot(abs(shankfilt(:,1)));
  title('x2Capire quando turn per Caviglia SX ');
  [m1,ind1]=findpeaks(shankfilt(:,3));
  figure()
  findpeaks(-shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))),'MinPeakDistance',50);
  hold on 
  plot(shankfilt(:,1));
  title('Capire rotazioni')
end

%hold on
%findpeaks(-shankfilt(:,3));
sogliam=mean(m);
 
%% CAPIRE TURN
% figure()
% findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))));
% title('Capire quando turn per Caviglia SX');
% hold on
% findpeaks(abs(shankfilt(:,1)),'MinPeakHeight',std(abs(shankfilt(:,1))));
[mturn,indturn]=findpeaks(abs(shankfilt(:,1)),'MinPeakHeight',std(abs(shankfilt(:,1))));

% figure()
% findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))));
% hold on 
% plot(abs(shankfilt(:,1)));
% title('x2Capire quando turn per Caviglia SX ');
% [m1,ind1]=findpeaks(shankfilt(:,3));%picchi per Hs e toe off

% figure()
% findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))));
% hold on 
% plot(shankfilt(:,1));
% title('Capire rotazioni')


%% Rivedere !!!! luigi_p_9000815039_s_1_t_3.
to_sx=[];
hs_sx=[];
stepsx=[];
tmp=0;
tmpr=[];
stdturnsx=std(mturn);
diffcontactsx=diff(ind);
meancontactsx=mean(diffcontactsx);
for i=1:length(ind)
    flag=1;
%     if abs((shankfilt(ind(i),3))) < (0.70*sogliam)
%     %abs(3*shankfilt(ind(i),1))>m((i))  || abs((shankfilt(ind(i),3))) < (0.70*sogliam)
%       flag=0;
%    end
      for k=1:length(indturn)
        if ind(i)-30<=indturn(k) && ind(i)+30>=indturn(k)
            %if mturn(k)>=tmp
                if  (m(i)-mturn(k))/m(i)<0.50|| mturn(k)>m(i)  
                  %||mturn(k)>2*std(mturn) mturn(k)>2*mean(mturn)%|| mturn(k)/m(i)>0.60 %|| mturn(k)<0.9*tmpr abs(mturn(k)-m(i))<0.82*m(i)
                    flag=0;
                    tmpr=[tmpr,ind(i)];
                    break
                end
            %end
            %tmp=mturn(k);
        end
        
            
      end
    if(flag~=0)
    for j=1:length(ind1)
        if (ind1(j)>ind(i) && flag==1)
            to_sx=[to_sx,ind1(j-1)];
            hs_sx=[hs_sx,ind1(j)];
            stepsx=[stepsx,ind(i)];
            flag=0;
        end
     end
    end
end

% endstepsx=[];
% diffstepsx=diff(stepsx);
% %meantotr=mean(difftotr);
% for i=1:length(diffstepsx)-1
%     if diffstepsx(i+1)>2*diffstepsx(i)
%         endstepsx=[endstepsx,i+1];
%     end
% end
endstepsx=[];
for i=1:length(tmpr)
    for k=1:length(stepsx)
        if stepsx(k)>tmpr(i)
            endstepsx=[endstepsx,k-1];
            break
        elseif stepsx(k+1)-stepsx(k)>3*meancontactsx
            endstepsx=[endstepsx,k];
        end
    end
end
endstepsx=unique(endstepsx);

% 
% endstepsx=[];
% diffstepsx=diff(stepsx);
% stepmeansx=mean(diffstepsx);
% for i=1:length(stepsx)-1
%     if stepsx(i+1)-stepsx(i)>stepmeansx
%         endstepsx=[endstepsx,i];
%     end
% end
shsx=struct('stephsx',[]);
stsx=struct('steptosx',[]);
%s(1).step=hs_sx(1:10);
%diff(s(1).step)
%s(2).step=hs_sx(11:18);
%step1=[];
%step1=hs_sx(1:10);
for i=1:length(endstepsx)
    if i==1
        shsx(i).stephsx=hs_sx(1:endstepsx(i));
        stsx(i).steptosx=to_sx(1:endstepsx(i));
    else
    shsx(i).stephsx=hs_sx(endstepsx(i-1)+1:endstepsx(i));
    stsx(i).steptosx=to_sx(endstepsx(i-1)+1:endstepsx(i));
    end
end
% Aggiusto quando ne trovo uno e ce ne sono altre
if endstepsx(end)~=length(hs_sx)
    shsx(end+1).stephsx=hs_sx(endstepsx(i)+1:length(hs_sx));
    stsx(end+1).steptosx=to_sx(endstepsx(i)+1:length(to_sx));
end
for i=1:length(shsx)
    if length(shsx(i).stephsx)<4
        shsx(i)=[];
    end
        if length(stsx(i).steptosx)<4
            stsx(i)=[];
        end
end


figure()
findpeaks(shankfilt(:,3),'MinPeakHeight',100);
hold on
findpeaks(-shankfilt(:,3))
title('Controllare heel strike e toe off presi');

%% ANALISI CAVIGLIA DX
shankfiltdx=[];
[b1,a1] = butter(5,1/(60/2),'low');
shankfiltdx(:,1)=filtfilt(b1, a1, sensorData{1}(:,4));
[b,a] = butter(5,4/(60/2),'low');
shankfiltdx(:,2:3)=filtfilt(b, a, sensorData{1}(:,5:6));

% figure()
% findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))));
% title('MASSIMI di VelAngX caviglia dx');

[mdx,inddx]=findpeaks(shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
[mdx2,inddx2]=findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
media2=mean(mdx2);
media1=mean(mdx);
if media2>media1
    [mdx,inddx]=findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    figure()
    findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    title('MASSIMI di VelAngX caviglia dx');
    figure()
    findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    hold on
    findpeaks(abs(shankfiltdx(:,1)),'MinPeakHeight',std(abs(shankfiltdx(:,1))))
    title('Capire turn soggetto a caviglia dx');
    [m1dx,ind1dx]=findpeaks(shankfiltdx(:,3));%%PER HS E TOE OFF
    figure()
    findpeaks(shankfiltdx(:,3));
    figure()
    findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    hold on 
    plot(abs(shankfiltdx(:,1)));
    title('x2Capire turn soggetto a caviglia dx');
    figure()
    findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    hold on 
    plot(shankfiltdx(:,1))  
else
    [mdx,inddx]=findpeaks(shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    figure()
    findpeaks(shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    title('MASSIMI di VelAngX caviglia dx');
    figure()
    findpeaks(shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    hold on
    findpeaks(abs(shankfiltdx(:,1)),'MinPeakHeight',std(abs(shankfiltdx(:,1))))
    title('Capire turn soggetto a caviglia dx');
    [m1dx,ind1dx]=findpeaks(-shankfiltdx(:,3)); %%PER HS E TOE OFF
    figure()
    findpeaks(-shankfiltdx(:,3));
    figure()
    findpeaks(shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    hold on 
    plot(abs(shankfiltdx(:,1)));
    title('x2Capire turn soggetto a caviglia dx');
    figure()
    findpeaks(shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    hold on 
    plot(shankfiltdx(:,1))
    figure()
    findpeaks(shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))),'MinPeakDistance',50);
    hold on
    findpeaks(-shankfiltdx(:,3));
end
sogliadx=mean(mdx);
% figure()
% findpeaks(-shankfiltdx(:,3))


% figure()
% findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))));
% hold on
% findpeaks(abs(shankfiltdx(:,1)),'MinPeakHeight',std(abs(shankfiltdx(:,1))))
% title('Capire turn soggetto a caviglia dx');
[mturndx,indturndx]=findpeaks(abs(shankfiltdx(:,1)),'MinPeakHeight',std(abs(shankfiltdx(:,1))));

%[m1dx,ind1dx]=findpeaks(-shankfiltdx(:,3)); %%PER HS E TOE OFF
%plot(abs(shankfiltdx(:,1)))
%plot((shankfiltdx(:,1)))
% figure()
% findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))));
% hold on 
% plot(abs(shankfiltdx(:,1)));
% title('x2Capire turn soggetto a caviglia dx');
% 
% figure()
% findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))));
% hold on 
% plot(shankfiltdx(:,1))


% todx=[];
% hsdx=[];
% tmpdx=0;
% tmpr=0;
% cont=1;
% l=1;
% stdturndx=std(mturndx);
%  for i=1:length(inddx)
%     for j=1:length(ind1dx)
%         if (ind1dx(j)>inddx(i))
%            todx=[todx,ind1dx(j-1)];
%            hsdx=[hsdx,ind1dx(j)];
%            break
%         end
%     end
%  end
% hstotdx=[];
% cont=1;
% j=0;
% flag=0;
% for i=1:length(hsdx)
%     for k=1:length(indturndx)
%         if hsdx(i)-30<=indturndx(k) && hsdx(i)+30>=indturndx(k)
%             flag=0;
%             cont=0;
%            if mturndx(k)>mdx(i) || mturndx(k)>2*stdturndx
%             flag=flag+1;
%             cont=cont+1;
%             j=1;
%             break
%            else
%                cont=tempory;
%                j=j+1;
%                break
%            end
%         elseif(flag>0) 
%             cont=tempory+1;
%             flag=0;
%             
%       end
% %             if(flag>0)
% %                 cont=cont+1;
% %                 j=1;
% %                 break
%                 
% %             end
%             
%     end
%        if(flag==0)
%             j=j+1;
%             hstotdx(cont,j)=hsdx(i);
%             %j=j+1;
%             tempory=cont;
%             end
% end
%     
%         
%  
% %     flag=0;
% %     if abs((shankfiltdx(inddx(i),3))) < (0.80*sogliadx)
% %         %abs(3*shankfiltdx(inddx(i),1))>mdx((i)) || abs((shankfiltdx(inddx(i),3))) < (0.70*sogliadx)
% %         flag=0;
% %     end
%       for k=1:length(indturndx)
%           for i=1:length(inddx)
%               flag=1;
%         if inddx(i)-30<=indturndx(k) && inddx(i)+30>=indturndx(k)
%             %if mturn(k)>=tmp
%                 if mturndx(k)>mdx(i) || mturndx(k)>2*stdturndx
%                  %abs(mturndx(k)-mdx(i))<0.82*mdx(i) || mturndx(k)>mdx(i) 
%                     %||mturndx(k)>0.90*tmpr
%                     flag=0;
%                     tmpr=mturndx(k);
%                     cont=cont+1;
%                     %break
%                 end
%         end
%                   for j=1:length(ind1dx)
%                     if (ind1dx(j)>inddx(i) && flag==1)
%                         todx(cont,l)=ind1dx(j-1);
%                         hsdx(cont,l)=ind1dx(j);
%                         l=l+1;
%                         flag=0; %mi serve mentre scorre il for con k
%                     end
%                  end
%             %end
%             %tmp=mturn(k);
%         end
% end
% 
% A=[];
% A(1,1)=1;
% A(1,2)=2;
% 
%    if flag~=0
%     for j=1:length(ind1dx)
%         if (ind1dx(j)>inddx(i) && flag==1)
%             todx=[todx,ind1dx(j-1)];
%             hsdx=[hsdx,ind1dx(j)];
%             flag=0;
%          end
%     end
%    end
%    end

to_dx=[];
hs_dx=[];
stepdx=[];
tmpr=[];
stdturndx=std(mturndx);
diffcontact=diff(inddx);
meancontact=mean(diffcontact);
for i=1:length(inddx)-1
    flag=1;
%     if abs((shankfiltdx(inddx(i),3))) < (0.70*sogliadx)
%     %abs(3*shankfilt(ind(i),1))>m((i))  || abs((shankfilt(ind(i),3))) < (0.70*sogliam)
%       flag=0;
%    end
      for k=1:length(indturndx)
        if inddx(i)-30<=indturndx(k) && inddx(i)+30>=indturndx(k)
            %if mturn(k)>=tmp
                if (mdx(i)-mturndx(k))/mdx(i)<0.50|| mturndx(k)>mdx(i) 
                  %mturndx(k)>2*std(mturndx)  mturndx(k)>2*mean(mturndx %mturndx(k)/mdx(i)>0.60%|| mturn(k)<0.9*tmpr
                    %abs(mturn(k)-m(i))<0.82*m(i) || inddx(i+1)-inddx(i)>3*meancontact
                    flag=0;
                    tmpr=[tmpr,inddx(i)];
                    break
                end
            %end
            %tmp=mturn(k);
        end
      end
    if(flag~=0)
    for j=1:length(ind1dx)
        if (ind1dx(j)>inddx(i) && flag==1)
            to_dx=[to_dx,ind1dx(j-1)];
            hs_dx=[hs_dx,ind1dx(j)];
            stepdx=[stepdx,inddx(i)];
            flag=0;
        end
     end
    end
end
% for i=1:length(stepdx)-1
% if stepdx(i+1)-stepdx(i)>3*meancontact
%       endstepdx=[endstepdx,i-1];
% end
% end
endstepdx=[];
for i=1:length(tmpr)
    for k=1:length(stepdx)-1
        if stepdx(k)>=tmpr(i) 
            endstepdx=[endstepdx,k-1];
            break
        elseif stepdx(k+1)-stepdx(k)>3*meancontact
               endstepdx=[endstepdx,k]; 
               break
        end
    end
end
endstepdx=unique(endstepdx);
% endstepdx=[];
% diffstepdx=diff(stepdx);
% %meantotr=mean(difftotr);
% for i=1:length(diffstepdx)-1
%     if diffstepdx(i+1)>2*diffstepdx(i)
%         endstepdx=[endstepdx,i+1];
%     end
% end


figure()
plot(-shankfilt(:,3))
hold on
plot(shankfiltdx(:,3))
title('Blu sinistra,arancio dx');

% endstepdx=[];
% diffstepdx=diff(stepdx);
% stepmeandx=mean(diffstepdx);
% for i=1:length(stepdx)-1
%     if stepdx(i+1)-stepdx(i)>stepmeandx
%         endstepdx=[endstepdx,i];
%     end
% end
shdx=struct('stephsdx',[]);
stdx=struct('steptodx',[]);
%s(1).step=hs_sx(1:10);
%diff(s(1).step)
%s(2).step=hs_sx(11:18);
%step1=[];
%step1=hs_sx(1:10);
for i=1:length(endstepdx)
    if i==1
        shdx(i).stephsdx=hs_dx(1:endstepdx(i));
        stdx(i).steptodx=to_dx(1:endstepdx(i));
    else
    shdx(i).stephsdx=hs_dx(endstepdx(i-1)+1:endstepdx(i));
    stdx(i).steptodx=to_dx(endstepdx(i-1)+1:endstepdx(i));
    end
end
    if endstepdx(end)~=length(hs_dx)
    shdx(end+1).stephsdx=hs_dx(endstepdx(i)+1:length(hs_dx));
    stdx(end+1).steptodx=to_dx(endstepdx(i)+1:length(to_dx));
    end
for i=1:length(shdx)
    if length(shdx(i).stephsdx)<4
        shdx(i)=[];
    end
        if length(stdx(i).steptodx)<4
            stdx(i)=[];
        end
end


% figure()
% findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))));
% hold on 
% plot(shankfiltdx(:,1));

figure()
findpeaks(-shankfiltdx(:,3),'MinPeakHeight',std(abs(shankfiltdx(:,3))));
hold on
findpeaks(shankfiltdx(:,3))
title('Controllare heel strike e toe off presi');




%% Parameter spatio temporal caviglia 
% %Parametername=struct('paziente',{},'stridetimeleft',[], 'stridetimeright',[],'steptime',[],'swingtimeleft',[], 'swingtimeright',[],'stancetimeleft',[], 'stancetimeright',[],'ssright',[], 'ssleft',[],'dbleft',[],'dbright',[]);
% paziente={FileName,sId};
% name={FileName}
% Parameter.FileName=struct('stridetimeleft',[], 'stridetimeright',[],'steptime',[],'swingtimeleft',[], 'swingtimeright',[],'stancetimeleft',[], 'stancetimeright',[],'ssright',[], 'ssleft',[],'dbleft',[],'dbright',[]);
% Parametername.paziente{1:2}(:,1)
% stridetimeleft=[];
% stridetimeright=[];
% stancetimesx=[];
% stancetimedx=[];
% steptime=[];
% steptimedx=abs((hsdx-hs_sx))/60;
% steptimesx=abs((hs_sx-hsdx))/60;
% swingtimesx=abs(hs_sx-to_sx)/60; 
% swingtimedx=abs(hsdx-todx)/60;
% % stancetimesx=abs(to_sx-hs_sx)/60;
% % stancetimedx=abs(todx-hsdx)/60; 
% doublesupport=stancetimesx+stancetimedx; 
% singlesupportdx=swingtimesx;
% singlesupportsx=swingtimedx;
% for i=1:length(hs_sx)-1
%     stridetimeleft=[stridetimeleft,abs((hs_sx(i+1)-hs_sx(i)))/60]; 
%     stridetimeright=[stridetimeright,abs((hsdx(i+1)-hsdx(i)))/60];
%     stancetimesx=[stancetimesx,abs((to_sx(i+1)-hs_sx))/60];%dubbio passo successivo
%     stancetimedx=[stancetimedx,abs((todx(i+1)-hsdx))/60;];
% end
% cadence=60*length(hs_dx)/(length(shankfilt)/60);
figure()
plot(-shankfilt(:,3))
hold on
plot(shankfiltdx(:,3))

steptime=struct('step',[]);
doublesupport1piede=struct('piede',[]);
doublesupport2piede=struct('piede',[]);
swingtime_dx=struct('dx',[]);
swingtime_sx=struct('sx',[]);
stancetime_dx=struct('dx',[]);
stancetime_sx=struct('sx',[]);
stridetime_sx=struct('sx',[]);
stridetime_dx=struct('dx',[]);
for i=1:length(shdx)
    for j=1:length(shdx(i).stephsdx)
        swingtime_dx(i).dx(j)=(shdx(i).stephsdx(j)-stdx(i).steptodx(j))/60;
    end
end
for i=1:length(shdx)
    for j=1:length(shsx(i).stephsx)
        swingtime_sx(i).sx(j)=(shsx(i).stephsx(j)-stsx(i).steptosx(j))/60;
        %steptime(i).step(j)=abs(shdx(i).stephsdx(j)-shsx(i).stephsx(j))/60;
    end
end
% meanswingtimesx=[];
% meanswingtimedx=[];
% for i=1:length(shdx)
% meanswingtimedx(i)=mean(swingtime_dx(i).dx);
% meanswingtimesx(i)=mean(swingtime_sx(i).sx);
% end
% tohspiede=struct('piede',[]);
% tohspiede2=struct('piede',[]);
% for i=1:length(shdx)
%     tohspiede(i).piede=[shdx(i).stephsdx stsx(i).steptosx];
%     tohspiede(i).piede=sort(tohspiede(i).piede);
%     tohspiede2(i).piede=[shsx(i).stephsx stdx(i).steptodx];
%     tohspiede2(i).piede=sort(tohspiede2(i).piede);
% end 
% for i=1:length(tohspiede)
%     if tohspiede(i).piede(1)==shdx(i).stephsdx(1)
%         doublesupport1piede(i).piede=diff(
    




hstotpiede=struct('piede',[]);
for i=1:length(shdx)
    hstotpiede(i).piede=[shdx(i).stephsdx shsx(i).stephsx];
    hstotpiede(i).piede=sort(hstotpiede(i).piede);
end   
for i=1:length(hstotpiede)
    steptime(i).step=diff(hstotpiede(i).piede)/60;
end
% meansteptime=[];
% for i=1:length(hstotpiede)
%     meansteptime(i)=mean(steptime(i).step);
% end
singlesupport_dx=struct('dx',[]);
singlesupport_sx=struct('sx',[]);
for i=1:length(swingtime_sx)
    for j=1:length(swingtime_sx(i).sx)
        singlesupport_dx(i).dx(j)=swingtime_sx(i).sx(j);
        singlesupport_sx(i).sx(j)=swingtime_dx(i).dx(j);
     end
end
% meansinglesupportdx=[];
% meansinglesupportsx=[];
% for i=1:length(shdx)
% meansinglesupportdx(i)=mean(singlesupport_dx(i).dx);
% meansinglesupportsx(i)=mean(singlesupport_sx(i).sx);
% end
% for i=1:length(shdx)
%     for j=1:length(stdx(i).steptodx)-1     
%     doublesupport1piede(i).piede(j)=(stdx(i).steptodx(j+1)-shsx(i).stephsx(j))/60;
%     end
% end
% 
% for i=1:length(shdx)
%     for j=1:length(stsx(i).steptosx)-1     
%     doublesupport2piede(i).piede(j)=(stsx(i).steptosx(j+1)-shdx(i).stephsdx(j))/60;
%     end
% end
figure()
plot(shankfilt(:,3))
hold on
plot(-shankfiltdx(:,3))
stancetime_dx=struct('dx',[]);
stancetime_sx=struct('sx',[]);
stridetime_sx=struct('sx',[]);
stridetime_dx=struct('dx',[]);
for i=1:length(shdx)
    for j=1:length(shdx(i).stephsdx)-1
    stancetime_dx(i).dx(j)=(stdx(i).steptodx(j+1)-shdx(i).stephsdx(j))/60;
    end
end
for i=1:length(shdx)
    for j=1:length(shsx(i).stephsx)-1
    stancetime_sx(i).sx(j)=(stsx(i).steptosx(j+1)-shsx(i).stephsx(j))/60;
    %stridetime_dx(i).dx(j)=diff(shdx(i).stephsdx)/60;
    %stridetime_sx(i).sx(j)=diff(shsx(i).stephsx)/60;
    end
end
for i=1:length(shdx)
stridetime_dx(i).dx=diff(shdx(i).stephsdx)/60;
stridetime_sx(i).sx=diff(shsx(i).stephsx)/60;
end

% meanstancedx=[];
% meanstancesx=[];
% meanstridedx=[];
% meanstridesx=[];
% for i=1:length(shdx)
%     meanstancedx(i)=mean(stancetime_dx(i).dx);
%     meanstancesx(i)=mean(stancetime_sx(i).sx);
%     meanstridedx(i)=mean(stridetime_dx(i).dx);
%     meanstridesx(i)=mean(stridetime_sx(i).sx);
% end
% doublesuppmean=[];
% for i=length(meanstancedx)
% doublesuppmean(i)=meanstancedx(i)-meanswingtimesx(i);
% end
%[corrreg,lags]=xcorr(abs(shankfilt(:,3)));


%% Stride regularity
figure()
plot(shankfilt(:,3));
hold on
plot(shankfiltdx(:,3));

%% Tronco sensore
% figure()
% plot(sensorData{3}(:,1));
% figure()
% plot(sensorData{3}(:,2));
% figure()
% plot(sensorData{3}(:,3))
% figure()
% plot(sensorData{3}(:,4:6));

% trunkfilt=[];
% [b,a] = butter(3,0.1/30,'high');
% trunkfilt=filtfilt(b, a, sensorData{3}(:,1:2));
% [b,a] = butter(3,2/30,'low');
% trunkfilt=filtfilt(b, a, trunkfilt(:,1:2));
% figure()
% plot(trunkfilt(:,1))

% [mtrunk,indtrunk]=findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
% [mtrunk1,indtrunk1]=findpeaks(trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
% meantrunk=mean(mtrunk);
% meantrunk1=mean(mtrunk1);
figure()
plot(sensorData{3}(:,3))
trunkfilt=[];
[b,a] = butter(3,[0.5 2]/30,'bandpass');  %% rivedere questo filtraggio
trunkfilt(:,3)=filtfilt(b, a, sensorData{3}(:,3));
figure()
plot(trunkfilt(:,3))

figure()
findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));

figure()
findpeaks(7*trunkfilt(:,3),'MinPeakHeight',std(abs(7*trunkfilt(:,3))));
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
hold on
plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)
title('Vedere picchi Heel strike su blu e se usare questo da soolo è ok per individuare tutti picchi')




figure()
plot(sensorData{3}(:,3))
figure()
plot(trunkfilt(:,3))

%  if meantrunk>meantrunk1
%  [mtrunk,indtrunk]=findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
%  figure()
%  findpeaks(-3*trunkfilt(:,1),'MinPeakHeight',std(abs(3*trunkfilt(:,1))));
%  hold on
% findpeaks(3*trunkfilt(:,3),'MinPeakHeight',std(abs(3*trunkfilt(:,3))));
%  hold on
% plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
% hold on
% plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)
%  title('Vedere Picchi Heel Strike su arancione')
%  figure()
%  findpeaks(trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
%  hold on
%  findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
%  title('Blu individua toe off che deve essere dopo arancione')
%  [mg2,indg2]=findpeaks(trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
%  else
%  [mtrunk,indtrunk]=findpeaks(trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
%  figure()
%  findpeaks(3*trunkfilt(:,1),'MinPeakHeight',std(abs(3*trunkfilt(:,1))));
%  hold on
%  findpeaks(3*trunkfilt(:,3),'MinPeakHeight',std(abs(3*trunkfilt(:,3))));
%  hold on
% plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
% hold on
% plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)
%  title('Vedere picchi Heel strike su arancione')
%  figure()
%  findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
%  hold on 
%  findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
%  title('Blue individua Toe off che deve essere dopo arancione')
%  [mg2,indg2]=findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
%  end
[mg,indg]=findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
%sogliatr=mean(mtrunk);
soglia2=mean(mg);
%FORSE POTREBBE ANDARE
% figure()
% findpeaks((abs(trunkfilt(:,1))),'MinPeakHeight',std(((abs(trunkfilt(:,1))))));

% [b,a] = butter(3,[0.5 5]/30,'bandpass');  %% rivedere questo filtraggio
% trunkfilt(:,3)=filtfilt(b, a,trunkfilt(:,3));
% figure()
% plot(trunkfilt(:,3))

%trunkfilt=detrend(sensorData{3}(:,1:3));


%Rivedere questo filtraggio
 trunkvert=[];
 [b,a] = butter(4,5/30,'low');
 trunkvert=filtfilt(b, a, sensorData{3}(:,1));
 figure()
 plot(trunkvert(:,1))
%VECCHIO FILTRAGGIO PER VERTICALE!!
  trunkvert=[];
 [b,a] = butter(3,[0.5 20]/30,'bandpass');
 trunkvert=filtfilt(b, a, sensorData{3}(:,1));
  figure()
 plot(trunkvert(:,1))
figure()
plot(sensorData{3}(:,1));
% [b,a] = butter(3,[0.5 3]/30,'bandpass');  %% rivedere questo filtraggio
% trunkfilt(:,3)=filtfilt(b, a, sensorData{3}(:,3));
% 
 figure()
 plot(-trunkvert(:,1));
 hold on
 plot(trunkfilt(:,3));
 
 figure()
 findpeaks((trunkvert(:,1)),'MinPeakHeight',std(abs(trunkvert(:,1))))

%% Aggiustare scelta dei campioni vicini perchè non mi sta bene
[b,a] = butter(2,0.2/30,'low');
trunkturn=filtfilt(b, a,sensorData{3}(:,4:6));
figure()
plot((abs(trunkturn(:,1))));
figure()
findpeaks(abs(trunkturn(:,1)),'MinPeakHeight',std(abs(trunkturn(:,1))));
hold on
plot((abs(trunkturn(:,2))));
hold on
plot((abs(trunkturn(:,3))));
figure()
plot(((trunkturn(:,1))));
hold on
plot(((trunkturn(:,3))));
figure()
plot(sensorData{3}(:,4:6))

figure()
plot(sensorData{3}(:,1:3))

turnindin=[];
turnindfi=[];
figure()
tmp=0;
findpeaks(abs(trunkturn(:,1)),'MinPeakHeight',std(abs(trunkturn(:,1))));
[mturntr,indturntr]=findpeaks(abs(trunkturn(:,1)),'MinPeakHeight',std(abs(trunkturn(:,1))));
[mturntr1,indturntr1]=findpeaks(abs(trunkturn(:,3)),'MinPeakHeight',std(abs(trunkturn(:,3))));
[mturntr2,indturntr2]=findpeaks(abs(trunkturn(:,2)),'MinPeakHeight',std(abs(trunkturn(:,2))));

meantrturn=mean(mturntr);
meantrturn1=mean(mturntr1);
meantrturn2=mean(mturntr2);
if (meantrturn> meantrturn1 && meantrturn>meantrturn2)
    figure()
    findpeaks(abs(trunkturn(:,1)),'MinPeakHeight',std(abs(trunkturn(:,1))));
    [mturntr,indturntr]=findpeaks(abs(trunkturn(:,1)),'MinPeakHeight',std(abs(trunkturn(:,1))));
    o=1;
end
 if (meantrturn1> meantrturn && meantrturn1>meantrturn2)
    figure()
    findpeaks(abs(trunkturn(:,3)),'MinPeakHeight',std(abs(trunkturn(:,3))));
    [mturntr,indturntr]=findpeaks(abs(trunkturn(:,3)),'MinPeakHeight',std(abs(trunkturn(:,3))));
    o=3;
 end
if (meantrturn2> meantrturn && meantrturn2>meantrturn1)
    figure()
    findpeaks(abs(trunkturn(:,2)),'MinPeakHeight',std(abs(trunkturn(:,2))));
    [mturntr,indturntr]=findpeaks(abs(trunkturn(:,2)),'MinPeakHeight',std(abs(trunkturn(:,2))));
    o=2;
 end
for i=1:length(mturntr)
    flag=1;
    trovato=0;
%     if (indturntr(i)-tmp)<200 && i~=1
%         flag=0;
%         trovato=0;
%     end
for k=1:length(trunkturn)
  if abs(k-indturntr(i))<350 
    if abs(trunkturn(k,o))>(0.7*mturntr(i)) && flag ~=0
        turnindin=[turnindin,k];
        flag=0;
        trovato=1;
    end
   
    if abs(trunkturn(k,o))<0.7*mturntr(i) && trovato==1 && k<indturntr(i)
        turnindfi=[turnindfi,k];
        break       
    end
    if k==length(trunkturn) && mturntr(i)==mturntr(end) && trovato==1 && abs(trunkturn(k,o))>0.7*mturntr(i)
        turnindfi=[turnindfi,k];
    end
  end
%tmp=indturntr(i);
end
 tmp=indturntr(i);
end
Fine=vertcat(turnindin,turnindfi);
Fine=Fine';
% [b,a] = butter(2,1/30,'low');
% trunkfilt(:,2:3)=filtfilt(b, a, trunkfilt(:,2:3));
% [b,a] = butter(2,1/30,'low');
% trunkfilt(:,1)=filtfilt(b, a, trunkfilt(:,1));

% [b,a] = butter(3,[0.5 20]/30,'bandpass');
% trunkturn=filtfilt(b, a, sensorData{3}(:,4:6));
% [b,a] = butter(2,0.5/30,'low');
% trunkturn=filtfilt(b, a,sensorData{3}(:,4:6));
% figure()
% plot((abs(trunkturn(:,1))));
% hold on
% plot((abs(trunkturn(:,3))));
% hold on
% plot(trunkfilt(:,3))

% figure()
% plot(sensorData{3}(:,4:6));
% 
% 
% figure()
% plot(3*(-trunkfilt(:,1)))
% hold on
% plot(3*(trunkfilt(:,3)))
% hold on
% plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
% hold on
% plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)
% 
% 
% figure()
% plot(-sensorData{3}(:,1))
% hold on
% plot(2*sensorData{3}(:,3))
% hold on 
% plot(normalize(shankfilt(:,3),'range'))


% figure()
% plot(sensorData{3}(:,4:6))
% [b,a] = butter(3,1/30,'low');
% trunkturn=filtfilt(b, a, sensorData{3}(:,4:6));
% figure()
% plot(normalize(trunkturn(:,1:3)))


% figure()
% findpeaks((-trunkfilt(:,1)),'MinPeakHeight',std((abs(trunkfilt(:,1)))));
% hold on
% findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
% hold on
% plot((trunkturn(:,1:3)));

% figure()
% findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
% hold on
% plot(normalize(trunkturn(:,2)))
% title('Vedere da qui HS')


% figure()
% findpeaks(-trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
% title('Vedere da qui TO')
% [mtrunk,indtrunk]=findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));

[mg,indg]=findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
%[mto,indto]=findpeaks(-trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
[mto,indto]=findpeaks(-trunkfilt(:,3)); %%MEGLIO COSI LI PRENDO TUTTI SONO SICURO
%Credo che è nell'armonica verticale da pdf cosi pare!![mg2,indg2]=findpeaks(-trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
%sogliatr=mean(mtrunk);
figure()
findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
hold on
findpeaks(-trunkfilt(:,3))

figure()
findpeaks(-trunkfilt(:,3),'MinPeakHeight',(abs(trunkfilt(:,3))));

figure()
findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
hold on
findpeaks(-trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));

figure()
findpeaks(7*trunkfilt(:,3),'MinPeakHeight',std(abs(7*trunkfilt(:,3))));
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
hold on
plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)
title('Vedere picchi Heel strike su blu')

figure()
findpeaks(-5*trunkfilt(:,3),'MinPeakHeight',std(abs(5*trunkfilt(:,3))));
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
hold on
plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)
title('Vedere picchi TOE OFF su blu')


%Fine=vertcat(turnindin,turnindfi);
%% Calcolo heelstrike e toe off tronco
[mg,indg]=findpeaks(trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
[mg,indg]=findpeaks(trunkfilt(:,3),'MinPeakHeight',0);
figure()
findpeaks(trunkfilt(:,3));
%[mto,indto]=findpeaks(-trunkfilt(:,3),'MinPeakHeight',std(abs(trunkfilt(:,3))));
[mto,indto]=findpeaks(-trunkfilt(:,3)); %%MEGLIO COSI LI PRENDO TUTTI SONO SICURO
hstr=[];
flag=1;
for i=1:length(indg)
%     flag=1;
%     if  mg(i) < (0.70*soglia2)
%         flag=0;
%     end
%         if (flag~=0)
          hstr=[hstr,indg(i)];
         
end

% totr=[];
% conf_to=0;
% for i=1:length(indto)
%     for m=1:length(hstr)
%         if abs(hstr(m)-indto(i))<40 
%            totr=[totr,indto(i)];
%            break
%         end
%     end
% end
totr=[];
%Metodo2
for i=1:length(hstr)
    for m=2:length(indto)
        if indto(m)>hstr(i)
           totr=[totr,indto(m-1)];
           break
        end
    end
end
endsteptronco1=[];
dim=size(hstr,2);
      for j=1:size(Fine,1)
          for k=1:size(Fine,2)-1
              for i=1:dim
              if (hstr(i)>Fine(j,k)) && (hstr(i) < Fine(j,k+1))   
                   endsteptronco1=[endsteptronco1,i-1];
                   break
              end
          end
      end
  end
remove=[];
  for i=1:dim
      for j=1:size(Fine,1)
          for k=1:size(Fine,2)-1
              if (hstr(i)>Fine(j,k)) && (hstr(i) < Fine(j,k+1))   
                   remove=[remove,i];
                   %endsteptronco1=[endsteptronco1,i-1];
              end
          end
      end
  end
 hstr(remove(1:end))=[]; 
 
 
endsteptronco2=[];
dim=size(totr,2);
      for j=1:size(Fine,1)
          for k=1:size(Fine,2)-1
              for i=1:dim
              if (totr(i)>Fine(j,k)) && (totr(i) < Fine(j,k+1))   
                   endsteptronco2=[endsteptronco2,i-1];
                   break
              end
          end
      end
end
remove=[];
  for i=1:dim
      for j=1:size(Fine,1)
          for k=1:size(Fine,2)-1
              if (totr(i)>Fine(j,k)) && (totr(i) < Fine(j,k+1))   
                   remove=[remove,i];
              end
          end
      end
  end
  totr(remove(1:end))=[];
  endsteptronco1=unique(endsteptronco1);
  endsteptronco2=unique(endsteptronco2);
  

  

  
% endsteptr=[];
% diffhstr=diff(hstr);
% meanhstr=mean(diffhstr);
% for i=1:length(hstr)-1
%     if hstr(i+1)-hstr(i)>meanhstr
%         endsteptr=[endsteptr,i];
%     end
% end

% endsteptr=[];
% diffhstr=diff(hstr);
% %meantotr=mean(difftotr);
% for i=1:length(diffhstr)-1
%     if diffhstr(i+1)>3*diffhstr(i)
%         endsteptr=[endsteptr,i+1];
%     end
% end

hs_tronco=struct('hstrunk',[]);
for i=1:length(endsteptronco1)
    if i==1
        hs_tronco(i).hstrunk=hstr(1:endsteptronco1(i));
    else
    hs_tronco(i).hstrunk=hstr(endsteptronco1(i-1)+1:endsteptronco1(i));
    end
end
    if endsteptronco1(end)~=length(hstr)
    hs_tronco(end+1).hstrunk=hstr(endsteptronco1(i)+1:length(hstr));
    end
for i=1:length(hs_tronco)
    if length(hs_tronco(i).hstrunk)<6
        hs_tronco(i).hstrunk=[];
    end
end
  
% 
% endsteptr_to=[];
% difftotr=diff(totr);
% %meantotr=mean(difftotr);
% for i=1:length(difftotr)-1
%     if difftotr(i+1)>3*difftotr(i)
%         endsteptr_to=[endsteptr_to,i+1];
%     end
% end
to_tronco=struct('totrunk',[]);
for i=1:length(endsteptronco2)
    if i==1
        to_tronco(i).totrunk=totr(1:endsteptronco2(i));
    else
    to_tronco(i).totrunk=totr(endsteptronco2(i-1)+1:endsteptronco2(i));
    end
end
    if endsteptronco2(end)~=length(totr)
    to_tronco(end+1).totrunk=totr(endsteptronco2(i)+1:length(totr));
    end
for i=1:length(to_tronco)
    if length(to_tronco(i).totrunk)<6
        to_tronco(i).totrunk=[];
    end
end
 


figure()
findpeaks(3*trunkfilt(:,3),'MinPeakHeight',std(abs(3*trunkfilt(:,3))));
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
hold on
plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)
title('Vedere picchi Heel strike su blu')
steptime_trunk=struct('trunk',[]);
swingtime_tr=struct('trunk',[]);
stancetime_tr=struct('trunk',[]);
stridetime_tr=struct('trunk',[]);
singlesupporttrunk=struct('trunk',[]);
doublesupptrunk=struct('trunk',[]);
for i=1:length(hs_tronco)
    for j=1:length(hs_tronco(i).hstrunk)
        swingtime_tr(i).trunk(j)=(hs_tronco(i).hstrunk(j)-to_tronco(i).totrunk(j))/60;
        steptime_trunk(i).trunk=diff(hs_tronco(i).hstrunk)/60;
    end
end
meanswingtimetronco=[];
meansteptimetronco=[];
for i=1:length(hs_tronco)
    meanswingtimetronco(i)=mean(swingtime_tr(i).trunk);
    meansteptimetronco(i)=mean(steptime_trunk(i).trunk);
end
    
for i=1:length(swingtime_tr)
    for j=1:length(swingtime_tr(i).trunk)
        singlesupporttrunk(i).trunk(j)=swingtime_tr(i).trunk(j);
     end
end

for i=1:length(hs_tronco)
    for j=1:length(to_tronco(i).totrunk)-1
        doublesupptrunk(i).trunk(j)=(to_tronco(i).totrunk(j+1)-hs_tronco(i).hstrunk(j))/60;
    end
end
meandoublesupptronco=[];
for i=1:length(hs_tronco)
 meandoublesupptronco(i)=mean(doublesupptrunk(i).trunk);
end

%Controlla stance!!!!!!!!!!!!!
for i=1:length(hs_tronco)
    for j=1:length(hs_tronco(i).hstrunk)-2
            %if to_tronco(i).totrunk(1)<hs_tronco(i).hstrunk(1)
                stancetime_tr(i).trunk(j)=(to_tronco(i).totrunk(j+2)-hs_tronco.hstrunk(j))/60;
            %else
              %stancetime_tr(i).trunk(j)=(to_tronco(i).totrunk(j+1)-hs_tronco.hstrunk(j));
            %end
    end
end
meanstancetimetronco=[];
for i=1:length(hs_tronco)
    meanstancetimetronco(i)=mean(stancetime_tr(i).trunk);
end
             
for i=1:length(hs_tronco)
    for j=1:length(hs_tronco(i).hstrunk)-2
    stridetime_tr(i).trunk(j)=(hs_tronco(i).hstrunk(j+2)-hs_tronco(i).hstrunk(j))/60;
    end
end
meanstridetimetronco=[];
for i=1:length(hs_tronco)
meanstridetimetronco(i)=mean(stridetime_tr(i).trunk);
end

 %% Calcolo step length 
start=struct('vert',[]);
finish=struct('vert',[]);
 for i=1:length(hs_tronco) 
     start(i).vert= hs_tronco(i).hstrunk(1);
     finish(i).vert=to_tronco(i).totrunk(end);
 end
 vel=struct('tronco',[]);
 pos=struct('tronco',[]);
 
 for i=1:length(hs_tronco)
  vel(i).tronco=cumtrapz(abs(trunkvert(start(i).vert(1):finish(i).vert(1))));
  pos(i).tronco=cumtrapz(abs(vel(i).tronco));
  [b,a]=butter(4,0.1/30,'high');
  pos(i).tronco=filtfilt(b,a,pos(i).tronco);
 end
 figure()
 plot(abs(pos(2).tronco));
 figure()
 plot(trunkvert(1:5))
  
  
  
  
  
  
  
  
  
  
 






% hstr=[];
% steptrunk=[];
% difft=diff(indtrunk);
% for i=1:length(indtrunk)
%     flag=1;
%   for k=1:length(indg)
%     if abs(trunkfilt(indtrunk(i),1)) < (0.70*sogliatr)
%         flag=0;
%         prova=[prova,i];
%     end
%     for m=1:length(difft)
%       if (flag~=0) & ((abs(indtrunk(i)-indg(k))<20) | abs(indg(k)-indtrunk(i)) < 0.5*(difft(m)))
%           %& (indtrunk(i)-indg(k)>0) Quando è in corrispondenza o poco dopo
%           %non me lo fa prendere per quello forse la rimuovo qui
%           %for j=1:length(Fine)-1
%               %if (indg(k)<Fine(j,:)) & (indg(k) < Fine(j+1,:))
%                hstr=[hstr,indg(k)];
%                steptrunk=[steptrunk,k];
%       end
%    end
%   end
% end
% hstr1=[];
% steptrunk=[];
% prova1=[];
% for i=1:length(indtrunk)
%     flag=1;
%    for k=1:length(indg)
%     if abs(trunkfilt(indtrunk(i),1)) < (0.65*sogliatr)
%         flag=0;
%         prova1=[prova1,i];
%     end
%       if (flag~=0) & (abs(indtrunk(i)-indg(k))<15) & (indtrunk(i)-indg(k)>0)
%           %for j=1:length(Fine)-1
%               %if (indg(k)<Fine(j,:)) & (indg(k) < Fine(j+1,:))
%                hstr1=[hstr1,indg(k)];
%                steptrunk=[steptrunk,k];
%               end
%           end
% end
% dim=size(hstr,2);
% remove=[];
%   for i=1:dim
%       for j=1:size(Fine,1)
%           for k=1:size(Fine,2)-1
%               if (hstr(i)>Fine(j,k)) && (hstr(i) < Fine(j,k+1))   
%                    remove=[remove,i];
%               end
%           end
%       end
%   end
%   hstr(remove(1:end))=[]; 
% totr=[];
% difftrunk=diff(hstr);
%   for j=1:length(indg2)
%       trovato=0;
%       for i=1:length(hstr)
%           if indg2(j)> hstr(i) && trovato ~=1 && abs(indg2(j)-hstr(i))<30
%               totr=[totr,indg2(j)];
% %               trovato=1;
% %           end
% %       end
% %   end
% %     for j=1:length(ind1dx)
% %         if (ind1dx(j)>inddx(i) && flag==1)
% %             todx=[todx,ind1dx(j-1)];
% %             hsdx=[hsdx,ind1dx(j)];
% %             flag=0;
% %         end
% 
% 
% 
% %mtrunk,indtrunk=findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));      
%  figure()
%  plot(-trunkfilt(:,1))
%  hold on
%  plot(trunkfilt(:,3))
%  %hold on 
% % plot(normalize(trunkturn))
%  
%  figure()
% [coeff,f]=cwt(-trunkfilt(:,1),'bump',60);
% 
% passi=-trunkfilt(:,1)'.*abs(coeff(38,:));
% figure()
% plot(passi)
% %signal=idwt(coeff,zeros(size(coeff)),'sym4');
% figure()
% plot(normalize(-trunkfilt(:,1)))
% hold on
% plot(normalize(passi));
% 
% 
% 
% 
% 
% 
% 
% 
% 
% plot(coeff)
% int=coeff(5,:);
% abs(int);
% figure()
% plot(-trunkfilt(:,1));
% hold on
% plot(abs(int))
% 
% 
% 
% 
% trunkfilt(:,3)=detrend(trunkfilt(:,3));
% trunkfilt(:,3)=cumtrapz(trunkfilt(:,3));
% [b,a] = butter(3,0.5/30,'high');
% trunkfilt(:,3)=filtfilt(b, a, trunkfilt(:,3));
% figure()
% cwt(trunkfilt(:,3),'bump',60)
% [coeff,f]=cwt(trunkfilt(:,3),'bump',60);
% figure()
% plot(trunkfilt(:,3));
% hold on 
% plot(abs(coeff(40,:)))
% 
% [psi,x]=gauswavf(25,30,,6000,1);
% figure()
% plot(x,psi);

%% STRIDE REGULARITY AND STEP REGULARITY
[b1,a1]=butter(4,3/30,'low');
trunkfilt(:,1)=filtfilt(b1, a1, trunkfilt(:,1)); %filtraggio per correlazione

[r,lags]=xcorr(trunkfilt(:,1),'unbiased');
figure()
plot(lags,r);
figure()
findpeaks(r,'MinPeakHeight',std(abs(r)))


%% Step length
vel1=cumtrapz(trunkfilt(:,1));
posiz=cumtrapz(vel1);

[b,a]=butter(4,0.1/30,'high');
posiz1filt=filtfilt(b,a,posiz);

figure()
plot(posiz1filt)



%[h,f]=freqz(b,a);
%figure()
%plot(f,20*log(abs(h)));
% d = designfilt('bandpassiir','FilterOrder',8, ...
%     'HalfPowerFrequency1',0.5,'HalfPowerFrequency2',20, ...
%     'SampleRate',60);
% sos = ss2sos(A,B,C,D);
% fvt = fvtool(sos,d,'Fs',60);
% legend(fvt,'butter','designfilt')
%[b,a] = butter(3,1/30,'low');
%trunkfilt=filtfilt(b, a, trunkfilt(:,1:3));

figure()
plot(trunkfilt(:,3));
title('Acc AP')

figure()
plot(-trunkfilt(:,1));
title('Acc Vert')


figure()
findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
title('MASSIMI di AccX TRONCO');
[mtr,indtr]=findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
%sogliadx=mean(mdx);

figure()
findpeaks(-trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
hold on
findpeaks(abs(trunkfiltgiro(:,1)),'MinPeakHeight',std(abs(trunkfiltgiro(:,1))))
title('Capire turn soggetto tronco');
[mturtr,indturntr]=findpeaks(abs(trunkfiltgiro(:,1)),'MinPeakHeight',std(abs(trunkfiltgiro(:,1))));

[m1tr,ind1tr]=findpeaks(trunkfilt(:,1)); %%PER HS E TOE OFF
%plot(abs(shankfiltdx(:,1)))
%plot((shankfiltdx(:,1)))
figure()
findpeaks(trunkfilt(:,1),'MinPeakHeight',std(abs(trunkfilt(:,1))));
hold on 
plot((trunkfiltgiro(:,1)));
title('z2Capire turn soggetto tronco');













% trunkfilt(:,3)=cumtrapz(trunkfilt(:,3));
% 
% figure()
% plot(trunkfilt(:,3));
% 
% [b,a] = butter(3,0.5/30,'high');
% trunkfilt=filtfilt(b, a, trunkfilt(:,1:3));
% figure()
% plot(trunkfilt(:,3));
% 
% % Gaussian wavelet
% [psi,x] = gauswavf(0,length(trunkfilt),1);
% figure()
% plot(x,psi)
% 
% trunkcwt=[];
% [w,f]=cwt(trunkfilt(:,3),'gauss',60);
% figure()
% plot(trunkfilt(:,3))
% hold on
% plot(abs(15*w(20,:)))
% 
% 
% 
% 
% trunkcwt=gauswavf(trunkfilt(:,3),'gauss',60);







%% Polso

figure()
plot(sensorData{4}(:,1:3))

figure()
plot(sensorData{4}(:,4:6))


polsot=[];
[b,a]=butter(10,0.5/30,'low');
polsot=filtfilt(b,a,sensorData{4}(:,4));
figure()
plot(abs(polsot(:,1)))
figure()
findpeaks(abs(polsot(:,1)),'MinPeakHeight',std(abs(polsot(:,1))));


polsofilt=[];
[b,a] = butter(3,0.2/30,'high');  %% rivedere questo filtraggio
polsofilt(:,1:3)=filtfilt(b, a, sensorData{4}(:,4:6));
[b,a] = butter(3,2/30,'low');  %% rivedere questo filtraggio
polsofilt(:,1:3)=filtfilt(b, a, polsofilt(:,1:3));
figure()
plot(sensorData{4}(:,4:6))
figure()
plot(polsofilt(:,1:3))

figure()
plot(polsofilt(:,1))
hold on
plot(shankfilt(:,3),'linewidth',1.5)
hold on
plot(-shankfiltdx(:,3),'linewidth',1.5)

figure()
plot(polsofilt(:,2))
hold on
plot(shankfilt(:,3),'linewidth',1.5)
hold on
plot(-shankfiltdx(:,3),'linewidth',1.5)



figure()
plot(polsofilt(:,3))
hold on
plot(shankfilt(:,3),'linewidth',1.5)
hold on
plot(-shankfiltdx(:,3),'linewidth',1.5)

figure()
plot(sensorData{4}(:,1:3))
polsofiltacc=[];

[b,a] = butter(3,0.2/30,'high');  %% rivedere questo filtraggio
polsofiltacc(:,1:3)=filtfilt(b, a, sensorData{4}(:,1:3));
[b,a] = butter(3,2/30,'low');  %% rivedere questo filtraggio
polsofiltacc(:,1:3)=filtfilt(b, a, polsofiltacc(:,1:3));
figure()
plot(sensorData{4}(:,1:3))
figure()
plot(polsofiltacc(:,1:3))

figure()
plot(polsofiltacc(:,1))
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
hold on
plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)

figure()
plot(polsofiltacc(:,2))
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
hold on
plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)



figure()
plot(polsofiltacc(:,3))
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
hold on
plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)




























figure()
plot(sensorData{4}(:,1:3))
hold on 


% polsot2=[];
% [b,a]=butter(10,5/30,'high');
% polsot2=filtfilt(b,a,sensorData{4}(:,1:3));
% [b,a]=butter(10,1/30,'low');
% polsot2=filtfilt(b,a,polsot2(:,1:3));
% figure()
% plot(abs(polsot2(:,1:3)))
% polsopicchi=[];
 figure()
 plot(sensorData{4}(:,1)-mean(sensorData{4}(:,1)))
% polsopicchi=sensorData{4}(:,1)-mean(sensorData{4}(:,1));
 figure()
 plot(sensorData{4}(:,1))
 
 figure()
 plot(sensorData{4}(:,4))
hold on
plot(shankfilt(:,3),'linewidth',1.5)
hold on
plot(-shankfiltdx(:,3),'linewidth',1.5)
polso=[];
[b,a]=butter(10,3/30,'low');
polso=filtfilt(b,a,sensorData{4}(:,5:6));
[b,a]=butter(10,1/30,'high');
polso=filtfilt(b,a,polso(:,1:2))
figure()
plot(polso(:,1:2))

figure()
findpeaks(3*polso(:,1),'MinPeakHeight',std(abs(3*polso(:,1))))
hold on
plot(shankfilt(:,3))

figure()
plot(3*polso(:,1))
hold on
plot(shankfilt(:,3))

figure()
findpeaks(7*polso(:,2),'MinPeakHeight',std(abs(7*polso(:,2))))
hold on
plot(shankfilt(:,3),'linewidth',1.5)
hold on
plot(-shankfiltdx(:,3),'linewidth',1.5)

figure()
findpeaks(7*polso(:,3),'MinPeakHeight',std(abs(7*polso(:,3))))
hold on
plot(shankfilt(:,3))
 
 figure()
 plot(sensorData{4}(:,1:3))
 
figure()
plot((polso(:,1)))
hold on
plot(normalize(polsoturn,'range'));
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)
hold on
plot(normalize(-shankfiltdx(:,3),'range'),'linewidth',1.5)


polsopicchi=[];
[b,a]=butter(10,[0.5 20]/30,'bandpass');
polsopicchi=filtfilt(b,a,sensorData{4}(:,1));
[b,a]=butter(5,5/30,'low');
polsopicchi=filtfilt(b,a,polsopicchi);
[b,a]=butter(4,1/30,'high');
[b,a]=butter(5,1/30,'low');
polsoturn=filtfilt(b,a,sensorData{4}(:,6));

figure()
plot(4*(-polsopicchi(:,1)))
hold on
plot(normalize(polsoturn,'range'));
hold on
plot(normalize(shankfilt(:,3),'range'),'linewidth',1.5)


figure()
plot((sensorData{4}(:,6))-mean(sensorData{4}(:,6)));
figure()
plot(sensorData{4}(:,6));

[b,a]=butter(5,1/30,'low');
polsoturn=filtfilt(b,a,sensorData{4}(:,6));
figure()
plot(polsoturn);


figure()
findpeaks(polsopicchi(:,1),'MinPeakHeight',std(abs(polsopicchi(:,1))));
title('Capire quando turn per Caviglia SX');
hold on
findpeaks(abs(shankfilt(:,1)),'MinPeakHeight',std(abs(shankfilt(:,1))));
[mturn,indturn]=findpeaks(abs(shankfilt(:,1)),'MinPeakHeight',std(abs(shankfilt(:,1))));

figure()
findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))));
hold on 
plot(abs(shankfilt(:,1)));
title('x2Capire quando turn per Caviglia SX ');
[m1,ind1]=findpeaks(-shankfilt(:,3)); %picchi per Hs e toe off

figure()
findpeaks(shankfilt(:,3),'MinPeakHeight',std(abs(shankfilt(:,3))));
hold on 
plot(shankfilt(:,1));
title('Capire rotazioni')


































