aa = TestData(:,1);
Smax = aa{:, 1}; %Maximum applied cyclic stress
bb = TestData(:,2);
NofC = bb{:, 1}; %Number of cycles
cc = TestData(:,3);
Sr = cc{:, 1}; %Residual strength

m=0;
n=1;
o=0;
p=1;

c=linspace(m,n,101); %Constraint for C
s=linspace(o,p,101); %Constraint for S

for i=1:length(c)
    for j=1:length(s)
        for k=1:length(Smax)
            sigmaE(k) = Smax(k)*(((Sr(k)/Smax(k)).^(1/s(j)) + (NofC(k)-1)*c(i)).^s(j)); %Equivalent static strength Equation (1)
        end
        parmHat=wblfit(sigmaE); %Wibul from equivalent static strength data
        alpha(i,j)=parmHat(2); %Find Weibul shape parameter
    end
end

max_alpha=max(alpha(:));
[row,col]=find(alpha==max_alpha);
new_c=c(row);
new_s=s(col);

for l=1:2
    m=new_c - (1/(10^(l+1)));
    n=new_c + (1/(10^(l+1)));
    o=new_s - (1/(10^(l+1)));
    p=new_s + (1/(10^(l+1)));
    c=linspace(m,n,21); %New smaller range for C
    s=linspace(o,p,21); %New smaller range for S
    for i=1:length(c)
        for j=1:length(s)
            for k=1:length(Smax)
                sigmaE(k) = Smax(k)*(((Sr(k)/Smax(k)).^(1/s(j)) + (NofC(k)-1)*c(i)).^s(j));
            end
            parmHat=wblfit(sigmaE);
            alpha(i,j)=parmHat(2);
        end
    end
max_alpha=max(alpha(:));
[row,col]=find(alpha==max_alpha);
new_c=min(c(row),1);
new_s=min(s(col),1);
end

for k=1:length(Smax)
    sigmaE(k) = Smax(k)*(((Sr(k)/Smax(k)).^(1/s(col)) + (NofC(k)-1)*c(row)).^s(col));
end
parmHat=wblfit(sigmaE);            
beta=parmHat(1); %Wibul scale parameter