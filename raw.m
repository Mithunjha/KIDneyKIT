%% hard coded - peak detection
clc;
clear all;
close all;
%% file opening %%
fid=fopen('C:/Users/Mithunjha/Desktop/hackx_semi/ecgraw5.txt');
data= textscan(fid,'%s');
fclose(fid);
data=data{1};
N=cellfun(@(x)str2double(x),data);
N(isnan(N))=[];
t=[1:length(N)];
%% smoothing %%
%%F1=smooth(N,10);
%%F2=smooth(N,'sgolay');
%%F3=sgolayfilt(N,7,21);
F=smooth(t,N,20,'lowess');%7
%% plot unfiltered %%
figure;
plot(t,N);
title('Ecg unfiltered');
ylabel('readings');
xlabel('samples)');
xlim([0 2000]);
ylim([250 700]);
grid on; 
%% plot filtered %%
figure;
plot(t,F);
title('Ecg matlab filtered');
ylabel('readings');
xlabel('samples');
xlim([0 2000]);
ylim([250 700]);
grid on;
%% Peakdetection R,P,T%%
[R_sample,R_points]=findpeaks(F,'MinPeakHeight',500,'MinPeakDistance',300);%550
[P_sample,low_P_points]=findpeaks(F,'MinPeakHeight',370);%360
[T_sample,low_T_points]=findpeaks(F,'MinPeakHeight',380); %360
l=1;
for i=1:length(low_P_points)
    if (F(low_P_points(i))>350 & F(low_P_points(i))<400)   %% Excludes R and other peaks
        P_points_temp(l)=low_P_points(i);
        l=l+1;
end
end
j=1;
for i=1:length(low_T_points)
    if (F(low_T_points(i))>350 & F(low_T_points(i))<500) %350,400  %% Exlcludes R and other peaks
        T_points_temp(j)=low_T_points(i);
        j=j+1;
    end
end
q=1;w=1;ind1=1;ind2=1;
for i=1:length(P_points_temp)
    if P_points_temp(i)<R_points(1)
        ind1=ind1+1;
    end
end
for i=1:length(T_points_temp)
    if T_points_temp(i)<R_points(1)        %rejects all the points before 1st peak
        ind2=ind2+1;
    end
end
ppoint=P_points_temp(ind1-1:end);%ind1-1
tpoint=T_points_temp(ind2:end);         % p and t points without first few points
for u=1:length(ppoint)
    if (q<=length(R_points))
        if (R_points(q)> ppoint(u)) & (R_points(q)-ppoint(u))<150
            P_points(q)=ppoint(u);
            q=q+1;
        end
    end
end
for i=1:length(tpoint)
    if (w<=length(R_points))
        if (R_points(w)<tpoint(i)) & (tpoint(i)-R_points(w)<150)
            T_points(w)=tpoint(i);
            w=w+1;
        end
    end
end
%% plot inverted %%
figure;
plot(t,-F);
title('Ecg inverted');
ylabel('readings');
xlabel('samples');
xlim([0 2000]);
ylim([-700 -250]);
grid on;
%% peak detection Q S %%
[Q_sample,low_Q_points]=findpeaks(-F,'MinPeakHeight',-360);%320
[S_sample,low_S_points]=findpeaks(-F,'MinPeakHeight',-360);%320/330
c=1;
d=1;
for u=1:length(low_Q_points)
    if (c<=length(R_points))
        if (R_points(c)> low_Q_points(u)) & (R_points(c)-low_Q_points(u))<30
            Q_points(c)=low_Q_points(u);
            c=c+1;
        end
    end
end
for i=1:length(low_S_points)
    if (d<=length(R_points))
        if (R_points(d)<low_S_points(i)) & (low_S_points(i)-R_points(d)<40)
            S_points(d)=low_S_points(i);
            d=d+1;
        end
    end
end
%% plotting detected points  %%
figure;
hold on;
plot(t,F);
plot(R_points,F(R_points),'rv','MarkerFaceColor','r');
plot(P_points,F(P_points),'rv','MarkerFaceColor','b');
plot(T_points,F(T_points),'rv','MarkerFaceColor','g');
plot(Q_points,F(Q_points),'r.','MarkerEdgeColor','m','MarkerSize',10);
plot(S_points,F(S_points),'r.','MarkerEdgeColor','k','MarkerSize',10);
xlim([0 2000]);
ylim([250 700]);
grid on;
legend('ECG Signal','R-waves','P-waves','T-waves','Q-waves','S-waves');
xlabel('Samples');
title('PQRST detection in filtered ECG Signal');
hold off;
%% Slope plotting %%
dx=mean(diff(t));
dy=gradient(F,dx);
figure;
%plot(t,F);
plot(t,dy,'-r');
hold on;
xlim([0 2000]);
ylim([-50 50]);
plot(R_points,0,'r.','MarkerEdgeColor','k','MarkerSize',10);
plot(P_points,0,'r.','MarkerEdgeColor','b','MarkerSize',10);
plot(T_points,0,'r.','MarkerEdgeColor','g','MarkerSize',10);
grid on;
legend('Slope of ECG Signal','R-waves','P-waves','T-waves');
xlabel('Samples');
title('slope in filtered ECG Signal');
hold off;
%% T right slope value %%
slope_val_mat=zeros(length(T_points),0);
K=1;
for i=1:length(T_points)
    slope_value=0;
    v=(length(t)-T_points(i));
    h=1;
    for j=1:20 %% add dy values until it starts to increase
        if dy(T_points(i)+j-1)>=dy(T_points(i)+j)
            slope_value= slope_value + dy(T_points(i)+j-1);
            h=h+1;
        end
    end
    slope_val_mat(K)=slope_value/h;
    K=K+1; 
end
max(slope_val_mat)-min(slope_val_mat)
