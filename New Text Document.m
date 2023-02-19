T = readtable('Debudi_xl.csv');
B = [];
A = [];
for i = 2:82
    B = [B , T(:,i)];
end
y = T(:,2);

z = table2array(y);
c = table2array(B);
start = 0;
en = 0;
for i = 1:502
    if ((c(i,1) <= 1650) && (c(i,1) >= 950))
        A = [A ; c(i,(2:end))];
        if (start == 0)
           start = i;
        end
        if (en < i)
            en = i;
        end
    end
end
X = z(start:en);
%%% its the comment to check git%%%
len = size(A);
l_1 = len(1);
l_2 = len(2);

figure(26)
plot(X,A);
title('RAW SPCTRA PLOT')


 GA = detrend(A);
 
figure(1)
plot(X, GA);
title('DETRENDED DATA PLOT')
hold on;

% 
% % % % %SNV DATA PLOT(GB)
m=mean(A);
s=std(A);
for i=1:l_1
    for j=1:l_2
        GB(i,j)=(A(i,j)-m(j))/s(j);
       
    end
end

figure(2)
plot(X, GB);
title('SNV DATA PLOT')
hold on;

for i=1:l_1
    for j=1:l_2
        GC(i,j)=(A(i,j)-m(j));
       
    end
end

figure(3)
plot(X, GC);
title('MEAN-CENTERING DATA PLOT');

% % %%%%% MSC DATA PLOT(GD)
GD = [];
m = mean(A');
n = size(m);
n_1 = n(1,2);

p = ones(1,n_1);
s = [p' m'];
psudoinverse = inv(s'*s)*s';

q = size(A);
for i = 1:q(2)
   w = psudoinverse * A(:,i); 
   GD = [GD  ((A(:,i) - p' * w(1))./w(2))];
end

figure(4)
plot(X, GD);
title('MSC DATA PLOT')

% % % %%%% MIN-MAX (GE)
for i=1:l_2
    r(i)=min(A(:,i));
    h(i)=max(A(:,i));
end

   for j=1:l_1
       for i=1:l_2
        GE(j,i)=(A(j,i)-r(i))/(h(i)-r(i));     
       end
   end

figure(5)
plot(X, GE);
title('MIN-MAX DATA PLOT')
hold on;
    
%  %  % % MEAN-CENTERING +DETREND DATA PLOT(GF)
   
% % DETRENDED DATA PLOT
GF=detrend(GC);

figure(6)
plot(X, GF);
title('MEAN-CENTERING +DETREND DATA PLOT')
hold on;

% % % MEAN-CENTERING + SNV DATA PLOT(GG)

% %%% SNV DATA PLOT
m=mean(GC);
s=std(GC);
for i=1:l_1
    for j=1:l_2
        GG(i,j)=(GC(i,j)-m(j))/s(j);
       
    end
end

figure(7)
plot(X, GG);
title('MEAN-CENTERING + SNV DATA PLOT')
hold on;
 
% % % MEAN-CENTERING + MSC DATA PLOT(GH)

% %%%% MSC DATA PLOT

GH = [];
m = mean(GC');
n = size(m);
n_1 = n(1,2);

p = ones(1,n_1);
s = [p' m'];
psudoinverse = inv(s'*s)*s';

q = size(GC);
for i = 1:q(2)
   w = psudoinverse * GC(:,i); 
   GH = [GH  ((GC(:,i) - p' * w(1))./w(2))];
end

figure(8)
plot(X, GH);
title('MEAN-CENTERING + MSC DATA PLOT')
hold on;

 
% % % % MEAN-CENTERING + MIN_MAX DATA PLOT(GI)

% %%%% MIN-MAX
for i=1:l_2
    r(i)=min(GC(:,i));
    h(i)=max(GC(:,i));
end

   for j=1:l_1
       for i=1:l_2
        GI(j,i)=(GC(j,i)-r(i))/(h(i)-r(i));     
   end
   end

figure(9)
plot(X, GI);
title('MEAN-CENTERING + MIN_MAX DATA PLOT')
hold on;
% % % % MSC +DETREND DATA PLOT(GJ)

% % DETRENDED DATA PLOT
GJ=detrend(GD);

figure(10)
plot(X, GJ);
title('MSC + DETREND DATA PLOT')
hold on;
 
% % % % MSC +SNV DATA PLOT(GK)
% 
% %SNV DATA PLOT
m=mean(GD);
s=std(GD);
for i=1:l_1
    for j=1:l_2
        GK(i,j)=(GD(i,j)-m(j))/s(j);
       
       end
end

figure(11)
plot(X, GK);
title('MSC +SNV DATA PLOT')
hold on;

% % % % MSC + MEAN-CENTERING DATA PLOT(GL)

% % MEAN-CENTERING DATA PLOT
m=mean(GD);

for i=1:l_1
    for j=1:l_2
        GL(i,j)=(GD(i,j)-m(j));
       
    end
end

figure(12)
plot(X, GL);
title('MSC + MEAN-CENTERING DATA PLOT')
hold on;

% % % MSC + MIN-MAX DATA PLOT(GM)

% %MIN-MAX 
for i=1:l_2
    r(i)=min(GD(:,i));
    h(i)=max(GD(:,i));
end

   for j=1:l_1
       for i=1:l_2
        GM(j,i)=(GD(j,i)-r(i))/(h(i)-r(i));     
       end
   end
  
figure(13)
plot(X, GM);
title('MSC + MIN-MAX DATA PLOT')
hold on;
   
   % % % SNV + DETREND DATA PLOT(GN)

   % DETRENDED DATA PLOT
GN=detrend(GB);

figure(14)
plot(X, GN);
title('SNV + DETREND DATA PLOT')
hold on;

% % % SNV + MSC DATA PLOT(GO)

%MSC DATA PLOT

GO = [];
m = mean(GB');
n = size(m);
n_1 = n(1,2);

p = ones(1,n_1);
s = [p' m'];
psudoinverse = inv(s'*s)*s';

q = size(GB);
for i = 1:q(2)
   w = psudoinverse * GB(:,i); 
   GO = [GO  ((GB(:,i) - p' * w(1))./w(2))];
end

figure(15)
plot(X, GO);
title('SNV + MSC DATA PLOT')
hold on;

% % % % SNV + MEAN-CENTERING DATA PLOT(GP)
 
% % MEAN-CENTERING DATA PLOT
m=mean(GB);

for i=1:l_1
    for j=1:l_2
        GP(i,j)=(GB(i,j)-m(j));
       
    end
end

figure(16)
plot(X, GP);
title('SNV + MEAN-CENTERING DATA PLOT')
hold on;

% % % % SNV + MIN-MAX  DATA PLOT(GQ)

% % % %MIN-MAX (GD)
for i=1:l_2
    r(i)=min(GB(:,i));
    h(i)=max(GB(:,i));
end

   for j=1:l_1
       for i=1:l_2
            GQ(j,i)=(GB(j,i)-r(i))/(h(i)-r(i));     
       end
   end
  
figure(17)
plot(X, GQ);
title('SNV + MIN-MAX  DATA PLOT')
hold on;
   
% % % MIN-MAX+ MEAN-CENTERING DATA PLOT(GR)
  
% %MEAN-CENTERING DATA PLOT(I)
m=mean(GE);

for i=1:l_1
    for j=1:l_2
        GR(i,j)=(GE(i,j)-m(j));
       
    end
end

figure(18)
plot(X, GR);
title('MIN-MAX+ MEAN-CENTERING DATA PLOT')
hold on;


% % % % MIN-MAX+ SNV DATA PLOT(GS)

% % %SNV DATA PLOT
m=mean(GE);
s=std(GE);
for i=1:l_1
    for j=1:l_2
        GS(i,j)=(GE(i,j)-m(j))/s(j);
       
    end
end
 
figure(19)
plot(X, GS);
title('MIN-MAX+ SNV DATA PLOT')
hold on;
% % % % MIN-MAX+ MSC DATA PLOT(GT)
% 
% % %MSC DATA PLOT

GT = [];
m = mean(GE');
n = size(m);
n_1 = n(1,2);

p = ones(1,n_1);
s = [p' m'];
psudoinverse = inv(s'*s)*s';

q = size(GE);
for i = 1:q(2)
   w = psudoinverse * GE(:,i); 
   GT = [GT  ((GE(:,i) - p' * w(1))./w(2))];
end

figure(20)
plot(X, GT);
title('MIN-MAX+ MSC DATA PLOT')
hold on;

% % % % MIN-MAX+ DETREND DATA PLOT(GU)

% % DETRENDED DATA PLOT
GU=detrend(GE);

figure(21)
plot(X, GU);
title('MIN-MAX+ DETREND DATA PLOT')

% % % %  DETREND + MSC DATA PLOT(GV)
 
% % %MSC DATA PLOT

GV = [];
m = mean(GA');
n = size(m);
n_1 = n(1,2);

p = ones(1,n_1);
s = [p' m'];
psudoinverse = inv(s'*s)*s';

q = size(GA);
for i = 1:q(2)
   w = psudoinverse * GA(:,i); 
   GV = [GV  ((GA(:,i) - p' * w(1))./w(2))];
end

figure(22)
plot(X, GV);
title('DETREND + MSC DATA PLOT')

% % % %  DETREND + SNV DATA PLOT(GW)

%SNV DATA PLOT(H)
m=mean(GA);
s=std(GA);
for i=1:l_1
    for j=1:l_2
        GW(i,j)=(GA(i,j)-m(j))/s(j);
       
    end
end

figure(23)
plot(X, GW);
title('DETREND + SNV DATA PLOT')
hold on;

% % % % %  DETREND + MC DATA PLOT(GX)

% % MEAN-CENTERING DATA PLOT
m=mean(GA);

for i=1:l_1
    for j=1:l_2
        GX(i,j)=(GA(i,j)-m(j));
       
    end
end

figure(24)
plot(X, GX);
title('DETREND + MC DATA PLOT')

% % % %  DETREND + MIN_MAX DATA PLOT(GY)
% 
% % %MIN-MAX (K)
for i=1:l_2
    r(i)=min(GA(:,i));
    h(i)=max(GA(:,i));
end

   for j=1:l_1
       for i=1:l_2
            GY(j,i)=(GA(j,i)-r(i))/(h(i)-r(i));     
       end
   end
   
figure(25)
plot(X, GY);
title('DETREND + MIN-MAX DATA PLOT')

% %standardizing the dataset inputs for Principal component calculation
% %Z-Score standardization 
% m = mean(A);
% S = std(A);
% b = [];
% A = A';
% for i = 1:l_2
%     b = [b ; (1/S(i))*(A(i,:) - m(i)*ones(1,l_1))];
% end

%PCA on raw data
no_of_samples = 8;
num_samp_perclass = [10 10 10 10 10 10 10 10];

da = A';
[coeff,score,latent,tsquared,explained,mu] = pca(da);

lab = [];
test = 1;
for i = 1:80
   temp = strcat('sample',num2str(test));
   lab = [lab ; temp];
   if (mod(i,10) == 0)
      test = test + 1;
   end
end

pc1 = score(:,1);
pc2 = score(:,2);
figure(27)
gscatter(pc1,pc2,lab);
title('PCA ON RAW DATA');
xlabel(explained(1));
ylabel(explained(2));
Separability = separability(score(:,1:2)',no_of_samples,num_samp_perclass);
disp(Separability);

%PCA on DETRENDED data
da_1 = GA';
[coeff,score,latent,tsquared,explained,mu] = pca(da_1);

pc1 = score(:,1);
pc2 = score(:,2);
figure(28)
gscatter(pc1,pc2,lab);
xlabel(explained(1));
ylabel(explained(2));
title('PCA ON DETRENDED DATA');
Separability = separability(score(:,1:2)',no_of_samples,num_samp_perclass);
disp(Separability);

%PCA on SNV data
da_2 = GB';
[coeff,score,latent,tsquared,explained,mu] = pca(da_2);

pc1 = score(:,1);
pc2 = score(:,2);
figure(29)
gscatter(pc1,pc2,lab);
xlabel(explained(1));
ylabel(explained(2));
title('PCA ON SNV DATA');
Separability = separability(score(:,1:2)',no_of_samples,num_samp_perclass);
disp(Separability);

%PCA on MEAN-CENTERING DATA
da_3 = GC';
[coeff,score,latent,tsquared,explained,mu] = pca(da_3);

pc1 = score(:,1);
pc2 = score(:,2);
figure(30)
gscatter(pc1,pc2,lab);
xlabel(explained(1));
ylabel(explained(2));
title('PCA ON MEAN CENTERED DATA');
Separability = separability(score(:,1:2)',no_of_samples,num_samp_perclass);
disp(Separability);

%PCA on MSC DATA
da_4 = GD';
[coeff,score,latent,tsquared,explained,mu] = pca(da_4);

pc1 = score(:,1);
pc2 = score(:,2);
figure(31)
gscatter(pc1,pc2,lab);
xlabel(explained(1));
ylabel(explained(2));
title('PCA ON MSC DATA');
Separability = separability(score(:,1:2)',no_of_samples,num_samp_perclass);
disp(Separability);

%PCA on MSC _ SNV DATA
da_5 = GK';
[coeff,score,latent,tsquared,explained,mu] = pca(da_5);
  
pc1 = score(:,1);
pc2 = score(:,2);
figure(32)
gscatter(pc1,pc2,lab);
xlabel(explained(1));
ylabel(explained(2));
title('PCA ON MSC _ SNV DATA');
Separability = separability(score(:,1:2)',no_of_samples,num_samp_perclass);
disp(Separability);

%PCA ON min_max data
da_6 = GE';
[coeff,score,latent,tsquared,explained,mu] = pca(da_1);

pc1 = score(:,1);
pc2 = score(:,2);
figure(33)
gscatter(pc1,pc2,lab);
xlabel(explained(1));
ylabel(explained(2));
title('PCA ON MIN MAX DATA');
Separability = separability(score(:,1:2)',no_of_samples,num_samp_perclass);
disp(Separability);

%PCA ON MEAN _CENTERING + MSC
%PCA ON min_max data
da_7 = GH';
[coeff,score,latent,tsquared,explained,mu] = pca(da_1);

pc1 = score(:,1);
pc2 = score(:,2);
figure(34)
gscatter(pc1,pc2,lab);
xlabel(explained(1));
ylabel(explained(2));
title('PCA ON MEAN_CENTERING + MSC');
Separability = separability(score(:,1:2)',no_of_samples,num_samp_perclass);
disp(Separability);