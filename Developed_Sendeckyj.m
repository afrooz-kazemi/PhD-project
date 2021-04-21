TestData = sortrows(Test,2);
aa = TestData(:,1);
SigmaA = aa{:, 1}; %Maximum applied cyclic stress
bb = TestData(:,2);
NofC = bb{:, 1}; %Number of cycles

SStrength = 292.2;

index = find(NofC<10000,1,'last');
LowCP = polyfit(log10(NofC(1:index)),SigmaA(1:index),1);
HighCP = polyfit(log10(NofC(index+1:length(NofC))),log10(SigmaA(index+1:length(NofC))),1);

c = linspace(0.01,1,100);
s = linspace(0.01,1,100);

for i = 1:length(c)
    for j = 1:length(s)
        for k = 1:length(NofC)
            if NofC(k) < 10000
                Diff(k) = abs((LowCP(1).*log10(NofC(k))+LowCP(2))-(SStrength./((1-c(i)+c(i).*NofC(k)).^s(j))))^2;
            else
                Diff(k) = abs((10.^(HighCP(1).*log10(NofC(k))+HighCP(2)))-(SStrength./((1-c(i)+c(i).*NofC(k)).^s(j))))^2;
            end
        end
        SSE(i,j) = sqrt(mean(Diff));
    end
end

MinSSE = min(SSE(:));
[row,col] = find(SSE==MinSSE);
ParameterC = c(row);
ParameterS = s(col);

for l = 1:2
    m = ParameterC - (1/(10^(l+1)));
    n = ParameterC + (1/(10^(l+1)));
    o = ParameterS - (1/(10^(l+1)));
    p = ParameterS + (1/(10^(l+1)));
    c = linspace(m,n,21); %New smaller range for C
    s = linspace(o,p,21); %New smaller range for S
    for i = 1:length(c)
        for j = 1:length(s)
            for k = 1:length(NofC)
                if NofC(k) < 10000
                    Diff(k) = abs((LowCP(1).*log10(NofC(k))+LowCP(2))-(SStrength./((1-c(i)+c(i).*NofC(k)).^s(j))))^2;
                else
                    Diff(k) = abs((10.^(HighCP(1).*log10(NofC(k))+HighCP(2)))-(SStrength./((1-c(i)+c(i).*NofC(k)).^s(j))))^2;
                end
            end
            SSE(i,j) = sqrt(mean(Diff));
        end
    end
    MinSSE = min(SSE(:));
    [row,col] = find(SSE==MinSSE);
    ParameterC = c(row);
    ParameterS = s(col);
end


