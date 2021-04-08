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
for i = 1:length(rf)
    for j = 1:length(rfCounting)
    if rf(i,2) == rfCounting(j,2) && rf(i,3)==rfCounting(j,3)
        rfCounting(j,1) = rfCounting(j,1)+rf(i,1);
    end
    end
end
%Rainflow Counting finishes.

%convert flatwise bending moment(kN.m) to stress11(Mpa) using transfer function(Mpa/kN.m) calculated by Abaqus.
%tension side transfer functions for all plys.
ff = TF1Flatwise(:,1);
tF = ff{:, 1};
for k = 1:length(tF)
    rfCountingply{k}(:,3) = rfCounting(:,3)*tF(k);
    rfCountingply{k}(:,2) = rfCounting(:,2)*tF(k);
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
aa = TestData1(:,1);
sigmaA1 = aa{:, 1};
nn = TestData1(:,2);
N1 = nn{:, 1};
rr = TestData1(:,3);
sigmaR1 = rr{:, 1};
[c1, s1] = sendeckyj_parameter(sigmaA1, N1, sigmaR1);
%T-C loading
R02 = -1;
aa = TestData2(:,1);
sigmaA2 = aa{:, 1};
nn = TestData2(:,2);
N2 = nn{:, 1};
rr = TestData2(:,3);
sigmaR2 = rr{:, 1};
[c2, s2] = sendeckyj_parameter(sigmaA2, N2, sigmaR2);
%C-C loading
R03 = 10;
aa = TestData3(:,1);
sigmaA3 = aa{:, 1};
nn = TestData3(:,2);
N3 = nn{:, 1};
rr = TestData3(:,3);
sigmaR3 = rr{:, 1};
[c3, s3] = sendeckyj_parameter(sigmaA3, N3, sigmaR3);
%data finishes.

%estimate number of cycles to failure for rfCountingply loads.
for k = 1:length(rfCountingply)
    for i = 1:length(rfCountingply{k})
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
                amp2(1)*(((rfCountingply{k}(i,2)/2)-y_intersect)/y_intersect)+amp2(1)
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

%estimate nou and residual strength for rfCountingply loads.
for year = 1:200
    for k = 1:length(rfCountingply)
        rfCountingply{k}(:,1) = rfCountingply{k}(:,1)*year;
        for i = 1:length(rfCountingply{k})%check for the first failure criteria.
            if rfCountingply{k}(i,1) >= rfCountingply{k}(i,5)
                fprintf('FAIL\n');
                return
            end
        end
        for i = 1:length(rfCountingply{k})
            amp3 = 0; amp2 = amp3; amp1 = amp2;
            mean3 = 0; mean2 = mean3; mean1 = mean2;
            nou3 = 0; nou2 = nou3; nou1 = nou2;
            in = 0;
            if R02<=rfCountingply{k}(i,4) && rfCountingply{k}(i,4)<=1
                for j = 1:floor(UTS/((rfCountingply{k}(i,1)*c2-c2+1)^s2))%j=maxstress is considered.
                    amp2(j) = j;
                    mean2(j) = 0;
                    nou2(j) = round(schaff_nou(UTS,c2,s2,rfCountingply{k}(i,1),amp2(j),R02)*10)/10; 
                end
                for j = 1:floor(UTS/((rfCountingply{k}(i,1)*c2-c2+1)^s2))-1
                    if nou2(j+1) < nou2(j)
                        nou2(j+1) = nou2(j);
                    end
                end
                nou = 1:0.1:max(nou2);
                indic2 = zeros(1, numel(nou)+1);
                for n = 1:numel(nou)
                    if ismember(nou(n), nou2) == 1
                        indic2(n) = min(find(abs(nou2 - nou(n))<1e-6));
                    else
                        if n == 1
                            indic2(n) = 1;
                        else
                            indic2(n) = indic2(n-1);
                        end
                    end
                end
                indic2(n+1) = length(amp2);
                for j = 1:floor(min((UTS/((rfCountingply{k}(i,1)*c1-c1+1)^s1)), 2*UTS/(1+R01)))
                    amp1(j) = j*(1-R01)/2;
                    mean1(j) = j*(1+R01)/2;
                    nou1(j) = round(schaff_nou(UTS,c1,s1,rfCountingply{k}(i,1),amp1(j),R01)*10)/10;
                end
                if range(nou1) == 0
                    rfCountingply{k}(i,6) = nou1(1);
                else
                    for j = 1:floor(min((UTS/((rfCountingply{k}(i,1)*c1-c1+1)^s1)), 2*UTS/(1+R01)))-1
                        if nou1(j+1) < nou1(j)
                            nou1(j+1) = nou1(j);
                        end
                    end
                    indic1 = zeros(1, numel(nou)+1);
                    for n = 1:numel(nou)
                        if ismember(nou(n), nou1) == 1
                            indic1(n) = min(find(abs(nou1 - nou(n))<1e-6));
                        else
                            if n == 1
                                indic1(n) = 1;
                            else
                                indic1(n) = indic1(n-1);
                            end
                        end
                    end
                    indic1(n+1) = length(amp1);
                    for n = 1:length(indic2)
                        x = [mean2(indic2(n)), mean1(indic1(n))];
                        y = [amp2(indic2(n)), amp1(indic1(n))];
                        p(n,:) = polyfit(x,y,1);
                    end
                    for n = 1:length(indic2)-1
                        if p(n,1)<0 && p(n+1,1)
                            polygonx = [mean2(indic2(n)),-p(n,2)/p(n,1),-p(n+1,2)/p(n+1,1),mean2(indic2(n+1))];
                            polygony = [amp2(indic2(n)),0,0,amp2(indic2(n+1))];
                            in(n) = inpolygon(rfCountingply{k}(i,3),rfCountingply{k}(i,2)/2,polygonx,polygony);
                            if in(n) == 1
                                rfCountingply{k}(i,6) = nou(n);
                                continue
                            end
                        else
                            in(n) = 0;
                        end
                    end
                    if all(in(:)==0)
                        rfCountingply{k}(i,6) = nou(1);
                    end
                end
            else
                for j = 1:floor(UTS/((rfCountingply{k}(i,1)*c2-c2+1)^s2))%j=maxstress is considered.
                    amp2(j) = j;
                    mean2(j) = 0;
                    nou2(j) = round(schaff_nou(UTS,c2,s2,rfCountingply{k}(i,1),amp2(j),R02)*10)/10; 
                end
                for j = 1:floor(UTS/((rfCountingply{k}(i,1)*c2-c2+1)^s2))-1
                    if nou2(j+1) < nou2(j)
                        nou2(j+1) = nou2(j);
                    end
                end
                nou = 1:0.1:max(nou2);
                indic2 = zeros(1, numel(nou)+1);
                for n = 1:numel(nou)
                    if ismember(nou(n), nou2) == 1
                        indic2(n) = min(find(abs(nou2 - nou(n))<1e-6));
                    else
                        if n == 1
                            indic2(n) = 1;
                        else
                        indic2(n) = indic2(n-1);
                        end
                    end
                end
                indic2(n+1) = length(amp2);
                for j = 1:floor(min((abs(UCS)/((rfCountingply{k}(i,1)*c3-c3+1)^s3)), 2*abs(UCS)/abs(1+R03)))
                    amp3(j) = abs(j*(1-R03)/2);
                    mean3(j) = abs(j*(1+R03)/2);
                    nou3(j) = round(schaff_nou(UCS,c3,s3,rfCountingply{k}(i,1),amp3(j),R03)*10)/10;
                end
                if range(nou3) == 0
                    rfCountingply{k}(i,6) = nou3(1);
                else
                    for j = 1:floor(min((abs(UCS)/((rfCountingply{k}(i,1)*c3-c3+1)^s3)), 2*abs(UCS)/abs(1+R03)))-1
                        if nou3(j+1) < nou3(j)
                            nou3(j+1) = nou3(j);
                        end
                    end
                    indic3 = zeros(1, numel(nou)+1);
                    for n = 1:numel(nou)
                        if ismember(nou(n), nou3) == 1
                            indic3(n) = min(find(abs(nou3 - nou(n))<1e-6));
                        else
                            if n == 1
                                indic3(n) = 1;
                            else
                                indic3(n) = indic3(n-1);
                            end
                        end
                    end
                    indic3(n+1) = length(amp3);
                    for n = 1:length(indic2)
                        x = [mean2(indic2(n)), mean3(indic3(n))];
                        y = [amp2(indic2(n)), amp3(indic3(n))];
                        p(n,:) = polyfit(x,y,1);
                    end
                    for n = 1:length(indic2)-1
                        if p(n,1)<0 && p(n+1,1)
                            polygonx = [mean2(indic2(n)),-p(n,2)/p(n,1),-p(n+1,2)/p(n+1,1),mean2(indic2(n+1))];
                            polygony = [amp2(indic2(n)),0,0,amp2(indic2(n+1))];
                            in(n) = inpolygon(rfCountingply{k}(i,3),rfCountingply{k}(i,2)/2,polygonx,polygony);
                            if in(n) == 1
                                rfCountingply{k}(i,6) = nou(n);
                                continue
                            end
                        else
                            in(n) = 0;
                        end       
                    end
                    if all(in(:)==0)
                        rfCountingply{k}(i,6) = nou(1);
                    end
                end
            end %nou estimation finishes, residual strength estimation starts.
            if R02<=rfCountingply{k}(i,4) && rfCountingply{k}(i,4)<=1
                if i == 1 || abs(UTS - rfCountingply{k}(i-1,8))<1e-6
                    rfCountingply{k}(i,7) = 0;
                else
                    rfCountingply{k}(i,7) = (((UTS-rfCountingply{k}(i-1,8))/(UTS-rfCountingply{k}(i,2)))^(1/rfCountingply{k}(i,6)))*rfCountingply{k}(i,5);%neff
                end
                rfCountingply{k}(i,8) = UTS-(UTS-rfCountingply{k}(i,2)/(1-rfCountingply{k}(i,4)))*((rfCountingply{k}(i,7)+rfCountingply{k}(i,1))/rfCountingply{k}(i,5))^rfCountingply{k}(i,6);
                rfCountingply{k}(i,9) = rfCountingply{k}(i,8)*abs(UCS)/UTS;%strength convrsion
            else
                if i == 1 || abs(abs(UCS)-rfCountingply{k}(i-1,9))<1e-6
                    rfCountingply{k}(i,7) = 0;
                else
                    rfCountingply{k}(i,7) = (((abs(UCS)-rfCountingply{k}(i-1,9))/(abs(UCS)-abs(rfCountingply{k}(i,2))))^(1/rfCountingply{k}(i,6)))*rfCountingply{k}(i,5);
                end
                rfCountingply{k}(i,9) = abs(UCS)-(abs(UCS)-abs(rfCountingply{k}(i,2)/(1-rfCountingply{k}(i,4))))*((rfCountingply{k}(i,7)+rfCountingply{k}(i,1))/rfCountingply{k}(i,5))^rfCountingply{k}(i,6);
                rfCountingply{k}(i,8) = rfCountingply{k}(i,9)*UTS/abs(UCS);
            end
            if abs(rfCountingply{k}(i,2)/(1-rfCountingply{k}(i,4))) >= rfCountingply{k}(i,8) || abs(rfCountingply{k}(i,2)/(1-rfCountingply{k}(i,4))) >= rfCountingply{k}(i,9)%check for second failure criteria.
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
for i = 1:length(rainflow)
    [~,idx] = min(abs(rainflow(i,2)-aveRange));
    rainflow(i,2) = aveRange(idx);
    [~,idx] = min(abs(rainflow(i,3)-aveMean));
    rainflow(i,3) = aveMean(idx);
end
f = rainflow;
end

function [c, s] = sendeckyj_parameter(sigmaA, N, sigmaR)
c0=linspace(0,1,101);
s0=linspace(0,1,101);
for i=1:length(c0)
    for j=1:length(s0)
        for k=1:length(sigmaA)
            sigmaE(k) = sigmaA(k)*(((sigmaR(k)/sigmaA(k)).^(1/s0(j)) + (N(k)-1)*c0(i)).^s0(j));%Equivalent static strength.
        end
        wbp = wblfit(sigmaE);%find weibull parameters.
        alpha(i,j) = wbp(2);
    end
end
max_alpha = max(alpha(:));
[row,col] = find(alpha==max_alpha);
c = c0(row);
s = s0(col);
for l=1:2
    m = c - (1/(10^(l+1)));
    n = c + (1/(10^(l+1)));
    o = s - (1/(10^(l+1)));
    p = s + (1/(10^(l+1)));
    c0 = linspace(m,n,21);
    s0 = linspace(o,p,21);
    for i=1:length(c0)
        for j=1:length(s0)
            for k=1:length(sigmaA)
                sigmaE(k) = sigmaA(k)*(((sigmaR(k)/sigmaA(k)).^(1/s0(j)) + (N(k)-1)*c0(i)).^s0(j));
            end
            wbp = wblfit(sigmaE);
            alpha(i,j) = wbp(2);
        end
    end
max_alpha = max(alpha(:));
[row, col] = find(alpha==max_alpha);
c = min(c0(row),1);
s = min(s0(col),1);
end
end

function f = schaff_nou(U,c,s,n,amp,R)
U = abs(U);
sLoad = abs(2*amp/(1-R));
sNumber = 1/c*(c-1+(U/sLoad)^(1/s));
sResidual= sLoad*((U/sLoad)^(1/s)-c*(n-1))^s;
if U <= sResidual
    f = 1;
else
f = log10((U-sResidual)/(U-sLoad))/log10(n/sNumber);
end
if f < 1
    f = 1;
end
end
