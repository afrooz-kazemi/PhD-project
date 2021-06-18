%Rainflow Counting
%Spn1My is flap bending moments for root region calculated for all winds in 30s.
tt = Spn1My(:,1);
t = seconds(tt{:, 1});
rf = [];
%Rainflow Counting for winds 3.5-24.5m/s.
for i = 2:23
x = Spn1My(:,i);
xx = x{:, 1};
r = rainflow (timetable(t , xx));
%Rayleigh distribution to find wind speed probapility.
speedRef = 42.5; %wind turbine class is ll.
speedAve = 0.2*speedRef; %IEC 61400
RayleighScaleP2 = 2*speedAve^2/pi;
speedHub = i+1.5;
probability = speedHub/RayleighScaleP2*exp(-(speedHub^2/(2*RayleighScaleP2)));
%change number of cycles from 30s to 1 year.
r(:,1) = r(:, 1)*probability*2*60*24*365;
rf = [rf; r];
end
%create Markov matrix
rf = MarkovMatrix(rf);
%unique and sorte Rainflow Counting 
rfCounting = unique(rf(:,[2 3]), 'rows');
rfCounting(:,2:end+1) = rfCounting(:,1:end);
rfCounting(:,1) = 0;
for i = 1:size(rf,1)
    for j = 1:size(rfCounting,1)
    if rf(i,2) == rfCounting(j,2) && rf(i,3)==rfCounting(j,3)
        rfCounting(j,1) = rfCounting(j,1)+rf(i,1);
    end
    end
end
%Rainflow Counting finishes.

for i = 1:size(rfCounting,1)
    rfCounting(i,4) = (2*rfCounting(i,3)-rfCounting(i,2))./(2*rfCounting(i,3)+rfCounting(i,2));
    rfCounting(i,5) = rfCounting(i,2)./(1-rfCounting(i,4));
end
%convert flatwise bending moment(kN.m) to stress11(Mpa) using transfer function(Mpa/kN.m) calculated by Abaqus.
%tension side transfer functions for all plys.
ff = TF1Flatwise(:,1);
tF = ff{:, 1};
for k = 1:length(tF)
    rfCountingply{k}(:,3) = rfCounting(:,3)*abs(tF(k));
    rfCountingply{k}(:,2) = rfCounting(:,2)*abs(tF(k));
    rfCountingply{k}(:,1) = rfCounting(:,1);
    rfCountingply{k} = MarkovMatrix(rfCountingply{k});
end
%convert finishes.

%static strength of the UD laminate on-axis.
UTS = 780.47;
UCS = -539.11;
%input data for three stress ratio.
%T-T
R01 = 0.1;
Test1 = sortrows(TestData1,2);
aa = Test1(:,1);
sigmaA1 = aa{:, 1};
nn = Test1(:,2);
N1 = nn{:, 1};
[c1, s1] = sendeckyj_parameter(UTS, sigmaA1, N1);
%T-C loading
R02 = -1;
Test2 = sortrows(TestData2,2);
aa = Test2(:,1);
sigmaA2 = aa{:, 1};
nn = Test2(:,2);
N2 = nn{:, 1};
[c2, s2] = sendeckyj_parameter(UTS, sigmaA2, N2);
%C-C loading
R03 = 10;
Test3 = sortrows(TestData3,2);
aa = Test3(:,1);
sigmaA3 = aa{:, 1};
nn = Test3(:,2);
N3 = nn{:, 1};
[c3, s3] = sendeckyj_parameter(abs(UCS), sigmaA3, N3);
%data finishes.

%estimate number of cycles to failure for rfCountingply loads.
for k = 1:size(rfCountingply,2)
    for i = 1:size(rfCountingply{k},1)
        rfCountingply{k}(i,4) = (2*rfCountingply{k}(i,3)-rfCountingply{k}(i,2))./(2*rfCountingply{k}(i,3)+rfCountingply{k}(i,2));%find stress ratio for each row.
        if R01<=rfCountingply{k}(i,4) && rfCountingply{k}(i,4)<1
            x = [UTS rfCountingply{k}(i,3)];
            y = [0 rfCountingply{k}(i,2)/2];
            p1 = polyfit(x,y,1);
            x = [0 1+R01];
            y = [0 1-R01];
            p2 = polyfit(x,y,1);
            x_intersect = fzero(@(x) polyval(p2-p1,x),0);
            rfCountingply{k}(i,5) = (1/c1)*(c1-1+(UTS/(2*x_intersect/(1+R01)))^(1/s1));
        end
        if R02<=rfCountingply{k}(i,4) && rfCountingply{k}(i,4)<R01
            for j = 1:floor(UTS)
                amp2(j) = j;
                mean2(j) = 0;
                NoCycle2(j) = (1/c2)*(c2-1+(UTS/(j*2/(1-R02)))^(1/s2));
                amp1(j) = UTS*(1-R01)/2/((NoCycle2(j)*c1-c1+1)^s1);%amp1 which have the same nember of cycles to failure.
                mean1(j) = (1+R01)/(1-R01)*amp1(j);
            end
            for j = 1:floor(UTS)-1
                polygonx = [mean2(j),mean2(j+1),mean1(j+1),mean1(j)];
                polygony = [amp2(j),amp2(j+1),amp1(j+1),amp1(j)];
                in(j) = inpolygon(rfCountingply{k}(i,3),rfCountingply{k}(i,2)/2,polygonx,polygony);
                if in(j) == 1
                    x = [mean1(j) mean2(j)];
                    y = [amp1(j) amp2(j)];
                    p1 = polyfit(x,y,1);
                    x = [0 1+rfCountingply{k}(i,4)];
                    y = [0 1-rfCountingply{k}(i,4)];
                    p2 = polyfit(x,y,1);
                    x_intersect1 = fzero(@(x) polyval(p2-p1,x),0);
                    y_intersect1 = polyval(p2,x_intersect1);
                    x = [mean1(j+1) mean2(j+1)];
                    y = [amp1(j+1) amp2(j+1)];
                    p1 = polyfit(x,y,1);
                    x = [0 1+rfCountingply{k}(i,4)];
                    y = [0 1-rfCountingply{k}(i,4)];
                    p2 = polyfit(x,y,1);
                    x_intersect2 = fzero(@(x) polyval(p2-p1,x),0);
                    y_intersect2 = polyval(p2,x_intersect2);
                    amplitude = (amp2(j+1)-amp2(j))*((rfCountingply{k}(i,2)/2-y_intersect1)/(y_intersect2-y_intersect1))+amp2(j);
                    rfCountingply{k}(i,5) = (1/c2)*(c2-1+(UTS/(2*amplitude/(1-R02)))^(1/s2));
                    continue
                end
            end
            if all(in(:)==0)
                x = [mean1(1) mean2(1)];
                y = [amp1(1) amp2(1)];
                p1 = polyfit(x,y,1);
                x = [0 1+rfCountingply{k}(i,4)];
                y = [0 1-rfCountingply{k}(i,4)];
                p2 = polyfit(x,y,1);
                x_intersect = fzero(@(x) polyval(p2-p1,x),0);
                y_intersect = polyval(p2,x_intersect);
                amplitude = amp2(1)*(((rfCountingply{k}(i,2)/2)-y_intersect)/y_intersect)+amp2(1);
                rfCountingply{k}(i,5) = (1/c2)*(c2-1+(UTS/(2*amplitude/(1-R02)))^(1/s2));
            end
        end
        if R02>rfCountingply{k}(i,4) || rfCountingply{k}(i,4)>R03
            for j = 1:floor(UTS)
                amp2(j) = j;
                mean2(j) = 0;
                NoCycle2(j) = (1/c2)*(c2-1+(UTS/(j*2/(1-R02)))^(1/s2));
                amp3(j) = abs(UCS*(1-R03))/2/((NoCycle2(j)*c3-c3+1)^s3);%amp3 which have the same nember of cycles to failure.
                mean3(j) = abs((1+R03)/(1-R03))*amp3(j);
            end
            for j = 1:floor(UTS)-1
                polygonx = [mean2(j),mean2(j+1),mean3(j+1),mean3(j)];
                polygony = [amp2(j),amp2(j+1),amp3(j+1),amp3(j)];
                in(j) = inpolygon(rfCountingply{k}(i,3),rfCountingply{k}(i,2)/2,polygonx,polygony);
                if in(j) == 1
                    x = [mean3(j) mean2(j)];
                    y = [amp3(j) amp2(j)];
                    p1 = polyfit(x,y,1);
                    x = [0 1+rfCountingply{k}(i,4)];
                    y = [0 1-rfCountingply{k}(i,4)];
                    p2 = polyfit(x,y,1);
                    x_intersect1 = fzero(@(x) polyval(p2-p1,x),0);
                    y_intersect1 = polyval(p2,x_intersect1);
                    x = [mean3(j+1) mean2(j+1)];
                    y = [amp3(j+1) amp2(j+1)];
                    p1 = polyfit(x,y,1);
                    x = [0 1+rfCountingply{k}(i,4)];
                    y = [0 1-rfCountingply{k}(i,4)];
                    p2 = polyfit(x,y,1);
                    x_intersect2 = fzero(@(x) polyval(p2-p1,x),0);
                    y_intersect2 = polyval(p2,x_intersect2);
                    amplitude = (amp2(j+1)-amp2(j))*((rfCountingply{k}(i,2)/2-y_intersect1)/(y_intersect2-y_intersect1))+amp2(j);
                    rfCountingply{k}(i,5) = (1/c2)*(c2-1+(UTS/(2*amplitude/(1-R02)))^(1/s2));
                    continue
                end
            end
            if all(in(:)==0)
                x = [mean3(1) mean2(1)];
                y = [amp3(1) amp2(1)];
                p1 = polyfit(x,y,1);
                x = [0 1+rfCountingply{k}(i,4)];
                y = [0 1-rfCountingply{k}(i,4)];
                p2 = polyfit(x,y,1);
                x_intersect = fzero(@(x) polyval(p2-p1,x),0);
                y_intersect = polyval(p2,x_intersect);
                amplitude = amp2(1)*(((rfCountingply{k}(i,2)/2)-y_intersect)/y_intersect)+amp2(1);
                rfCountingply{k}(i,5) = (1/c2)*(c2-1+(UTS/(2*amplitude/(1-R02)))^(1/s2));
            end
        end
        if 1<=rfCountingply{k}(i,4) && rfCountingply{k}(i,4)<R03
            x = [UCS rfCountingply{k}(i,3)];
            y = [0 rfCountingply{k}(i,2)/2];
            p1 = polyfit(x,y,1);
            x = [0 1+R03];
            y = [0 1-R03];
            p2 = polyfit(x,y,1);
            x_intersect = fzero(@(x) polyval(p2-p1,x),0);
            rfCountingply{k}(i,5) = (1/c3)*(c3-1+(UCS/(2*x_intersect/(1+R03)))^(1/s3));
        end        
    end
end

for year = 1:200
    for k = 1:size(rfCountingply,2)
        rfCountingply{k}(:,1) = rfCountingply{k}(:,1)*year;
        for i = 1:size(rfCountingply{k},1)
            rfCountingply{k}(i,6) = rfCountingply{k}(i,1)./rfCountingply{k}(i,5);
            if i == 1
                rfCountingply{k}(1,7) = rfCountingply{k}(1,6);
            else
                rfCountingply{k}(i,7) = rfCountingply{k}(i,6)+rfCountingply{k}(i-1,7);
            end
            if rfCountingply{k}(i,7) >= 1
                fprintf('FAIL\n');
                return
            end
        end
    end
end

function f = MarkovMatrix(rainflow)
minRange = floor(min(rainflow(:,2)));
maxRange = ceil(max(rainflow(:,2)));
minMean = floor(min(rainflow(:,3)));
maxMean = ceil(max(rainflow(:,3)));
d = 1; %unit distance for Markov matrix
aveRange(1) = minRange+ d/2;
aveMean(1) = minMean + d/2;
i = 1;
while aveRange(1) + d*i <= maxRange
    aveRange(i+1) = aveRange(1) + d*i;
    i = i+1;
end
i = 1;
while aveMean(1) + d*i <= maxMean
    aveMean(i+1) = aveMean(1) + d*i;
    i = i+1;
end
for i = 1:size(rainflow,1)
    [~,idx] = min(abs(rainflow(i,2)-aveRange));
    rainflow(i,2) = aveRange(idx);
    [~,idx] = min(abs(rainflow(i,3)-aveMean));
    rainflow(i,3) = aveMean(idx);
end
f = rainflow;
end

function [c, s] = sendeckyj_parameter(SStrength, sigmaA, N)
index = find(N<10000,1,'last');
LowCP = polyfit(log10(N(1:index)),sigmaA(1:index),1);
HighCP = polyfit(log10(N(index+1:length(N))),log10(sigmaA(index+1:length(N))),1);
c0=linspace(0.01,1,101);
s0=linspace(0.01,1,101);

for i = 1:length(c0)
    for j = 1:length(s0)
        for k = 1:length(N)
            if N(k) < 10000
                Diff(k) = abs((LowCP(1).*log10(N(k))+LowCP(2))-(SStrength./((1-c0(i)+c0(i).*N(k)).^s0(j))))^2;
            else
                Diff(k) = abs((10.^(HighCP(1).*log10(N(k))+HighCP(2)))-(SStrength./((1-c0(i)+c0(i).*N(k)).^s0(j))))^2;
            end
        end
        SSE(i,j) = sqrt(mean(Diff));
    end
end
MinSSE = min(SSE(:));
[row,col] = find(SSE==MinSSE);
c = c0(row);
s = s0(col);

for l = 1:2
    m = c - (1/(10^(l+1)));
    n = c + (1/(10^(l+1)));
    o = s - (1/(10^(l+1)));
    p = s + (1/(10^(l+1)));
    c0 = linspace(m,n,21); %New smaller range for C
    s0 = linspace(o,p,21); %New smaller range for S
    for i = 1:length(c0)
        for j = 1:length(s0)
            for k = 1:length(N)
                if N(k) < 10000
                    Diff(k) = abs((LowCP(1).*log10(N(k))+LowCP(2))-(SStrength./((1-c0(i)+c0(i).*N(k)).^s0(j))))^2;
                else
                    Diff(k) = abs((10.^(HighCP(1).*log10(N(k))+HighCP(2)))-(SStrength./((1-c0(i)+c0(i).*N(k)).^s0(j))))^2;
                end
            end
            SSE(i,j) = sqrt(mean(Diff));
        end
    end
    MinSSE = min(SSE(:));
    [row,col] = find(SSE==MinSSE);
    c = c0(row);
    s = s0(col);
end
end