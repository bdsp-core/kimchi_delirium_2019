clear all; clc; format compact; 

load rawData % T -- needs to be cleaned up before use
v = T.Properties.VariableNames; 
Y01 = T.cam3d; 
char(v); 

%% predict delirium using eeg
v(20:34) % - eeg findings
getThese = {'PDR'    'ThetaSlowing'    'ThetaFocal'  ...
    'ThetaGen' 'Delta'    'DeltaFocal'    'DeltaGen' ...
    'LPD'    'GPD'  'TPW'    'LRDA'    'GRDA'...    
    'LowVolt' 'Sporadic'}; 
ind = find(ismember(v,getThese)); 
Te = T(:,ind); 
L = v([20:32 34]);
y = T(:,44); 

%% ordinal regression
X = table2array(Te); 
Y = table2array(y); 
[~,jj]=sort(Y); 
Y = Y(jj); 
X = X(jj,:); 
% fill in nans
for i = 1:size(X,2); 
    x = X(:,i); 
    ind = find(isnan(x)); 
    m = nanmedian(x); 
    x(ind)=m; 
%     x = (x-mean(x))/std(x); 
    X(:,i)=x; 
end

%% put in order by least to most predictive in isolation
for i = 1:size(X,2); 
   x = X(:,i); 
   [r(i),p] = corr(x,Y,'Type','Spearman');
end
[ii,jj] = sort(r); 
X = X(:,jj); 
L = L(jj); 

%% change features
getThese = {'PDR'    ...
    'ThetaSlowing'  'ThetaGen' ...
    'Delta'    'DeltaGen' 'GRDA'...
    'LPD'    'GPD'  'TPW'    'LRDA' ...        
    'LowVolt' 'Sporadic'}; 

% getThese = {'PDR'    'ThetaGen'  'DeltaGen' ...
%     'LPD'    'GPD'  'TPW'    'LRDA'    'GRDA'...    
%     'LowVolt' 'Sporadic'}; 

ind = find(ismember(getThese,L))
L(ind)
X = X(:,ind); 
L = getThese; 
% temp1 = 
temp = X(:,9)|X(:,10); 
XX = [X temp]; 
LL = [L 'genSlow']; 
X = XX; 
L = LL; 
%%
[B,dev,stats] = mnrfit(X,Y+1,'model','ordinal');
[pihat,dlow,hi] = mnrval(B,X,stats,'model','ordinal');
pihat

%% put cases in order according to Y
[~,jj]=sort(Y); 
ph = pihat(jj,:);
yt = Y(jj); 
pt = zeros(size(ph)); 
for i = 1:length(yt); 
   ii = yt(i)+1; 
   pt(i,ii)=1; 
end

%% sort within severity scores
idx = [];
for i = 0:7
   i0 = find(Y==i);  
   temp = pihat(i0,8); 
   [ii,jj] = sort(temp); 
   idx = [idx i0(jj)']; 
end

X = X(idx,:); 
Y = Y(idx); 
pt = pt(idx,:); 
pihat = pihat(idx,:); 
ph=ph(idx,:); 


%%

figure(1); clf; 
subplot(211); 
bar(ph,'stacked')
axis tight
subplot(212); 
bar(pt,'stacked'); 
axis tight

%% try predicting actual value
yh = ph*([0:7]')
figure(2); clf; 
scatter(yt,yh); axis equal
axis([0 7 0 7])


%% 
figure(3); clf; 
subplot(7, 1,[2:7]); 
imagesc(1-X'); colormap gray; 
set(gca,'ytick',1:14); 
set(gca,'yticklabel',L); 
xlabel('Patient number'); 

% draw lines
for i = 0:size(X,1); 
    x = [i (i+1) (i+1) i i]; 
    y = [0 0 12 12 0]; 
    hold on; 
    plot(x,y,'b'); 
end

hold on
for i = 0:size(X,2); 
   x = [0 200 200 0 0]; 
   y = [i i (i+1) (i+1) i]+0.5; 
   plot(x,y,'k'); 
end

subplot(7,1,1); 
imagesc(1-Y'); 
set(gcf,'color','w')
set(gca,'ytick','','xtick',''); 

%% logistic regression
ind = find(Y01==1);
length(ind)/length(Y01); 
% find(
