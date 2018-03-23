%parametric CGCI

alpha = 0.05;
num_windows=5;

% load the prewhitened time series
% y0-y4 αντιστοιχο�?ν στις prewhitened χ�?ονοσει�?ές για την 1η εως την 5η 6ετία
y0 = load('C:\Users\xx\Desktop\ergasiaD2\prewhite0.dat');
y1 = load('C:\Users\xx\Desktop\ergasiaD2\prewhite1.dat');
y2 = load('C:\Users\xx\Desktop\ergasiaD2\prewhite2.dat');
y3 = load('C:\Users\xx\Desktop\ergasiaD2\prewhite3.dat');
y4 = load('C:\Users\xx\Desktop\ergasiaD2\prewhite4.dat');
[n,m]=size(y0);

% τα ονόματα για όλα τα αγαθά
nameM = textread('C:\Users\xx\Desktop\ergasiaD2\WBCommodityUSDnames.dat','%s');

% διάνυσμα με τα επιλεγμένα αγαθά
iV=[6 7 8 9 10 11 12 13 15 17 18 19 24 25 26 27 28 31 32 33];
names = nameM(iV,:);

A=[y0 y1 y2 y3 y4]
for k=1:5
for i=0,20,40,60,80
        y{k,:}=A(:,1+i:20+i)
end
end
% plot τιχ prewhitened χ�?ονοσει�?ές για την 1η 6ετία
figure(1)
clf
plot(y0,'.-')
xlabel('month t')
ylabel('y(t)')
legend(names)
title("prewhitened time series 1st window of 6 years")

% CGCI conditional Granger Causality Index
% β�?ίσκει άμεσες συσχετίσεις
%% [CGCIM,pCGCIM] = CGCinall(xM,m,maketest)
% CGCinall computes the conditional Granger Causality index (GCI) for all 
% time series pairs (Xi,Xj) in the presence of the rest time series in the
% vector time series given in X, for both directions (Xi->Xj|(X-{Xi,Xj} and 
% Xj->Xi|(X-{Xi,Xj}). 
% INPUTS
% - xM          : the vector time series  
% - m           : Order of the restricted and unrestricted AR model 
% - maketest    : If 1 make parametric test and give out the p-values
% OUTPUTS
% - CGCIM       : The matrix KxK of the conditional Granger Causality
%                 indexes, (i,j) for CGCI (Xi->Xj) 
% - pCGCIM      : The p-values of the parametric significance test for the
%                 values in CGCIM 

% θα δο�?με το δείκτη α�?χικά για τάξη 1
s=2   % ? ???? ??? ?? ??????? 1 ? 2

figno=1
for s=1:2
    for k=1:5
        [cgci,pcgci]=CGCinall(y{k,:},s,1);
tit2txt = sprintf('parametric p-value, CGCI_{X->Y}(%d)',s);
tmp01 = pcgci < alpha;  %???????? logical ??????, ???? 1=true ???? ? p_value<a 
                          %?????? ???????????? H0: ??? ??????? ??????? ???? ????????? ??? ?????????
h2 = plotnetworktitle(tmp01,[],names,tit2txt,figno+1);
figno=figno+1;
h3 = plotcolormap(pcgci,names,tit2txt,figno,1)
end
end

        
[cgci01,pcgci01]=CGCinall(y0,s,1);
[cgci02,pcgci02]=CGCinall(y1,s,1);
[cgci03,pcgci03]=CGCinall(y2,s,1);
[cgci04,pcgci04]=CGCinall(y3,s,1);
[cgci05,pcgci05]=CGCinall(y4,s,1);

%colormaps
figno = 2;
tit2txt = sprintf('parametric p-value, CGCI_{X->Y}(%d)',s);
tmp01 = pcgci01 < alpha;  %???????? logical ??????, ???? 1=true ???? ? p_value<a 
                          %?????? ???????????? H0: ??? ??????? ??????? ???? ????????? ??? ?????????
h2 = plotnetworktitle(tmp01,[],names,tit2txt,figno);
figno=figno+1;
h3 = plotcolormap(pcgci01,names,tit2txt,figno,1)


########################################
# in another way
A=[y0 y1 y2 y3 y4]
for k=1:5
for i=0,20,40,60,80
        y{k,:}=A(:,1+i:20+i)
end
end

figno=1
for s=1:2
    for k=1:5
        [cgci,pcgci]=CGCinall(y{k,:},s,1);
tit2txt = sprintf('parametric p-value, CGCI_{X->Y}(%d)',s);
tmp01 = pcgci < alpha;  %???????? logical ??????, ???? 1=true ???? ? p_value<a 
                          %?????? ???????????? H0: ??? ??????? ??????? ???? ????????? ??? ?????????
h2 = plotnetworktitle(tmp01,[],names,tit2txt,figno+1);
figno=figno+1;
h3 = plotcolormap(pcgci,names,tit2txt,figno,1)
end
end
##apla den mporw na kalesw pisw ta parathyra tou cgci, pcgci gia kathe 6etia