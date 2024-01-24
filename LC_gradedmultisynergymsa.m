close all
clear all


tic

load relativelesions

BA_label = {'BA1','BA2','BA3','BA4','BA5','BA6','BA7','BA8','BA9','BA10',...
    'BA11','BA17','BA18','BA19','BA20','BA21','BA22','BA23','BA24','BA25',...
    'BA26','BA27','BA28','BA29','BA30','BA32','BA34','BA35','BA36','BA37',...
    'BA38','BA39','BA40','BA41','BA42','BA43','BA44','BA45','BA46','BA47','BA48'};

n_pat = size(relative_lesions,1);
n_reg = size(relative_lesions,2);

figure
subplot(131)
imagesc(relative_lesions)
title('graded lesion data, level of lesion')
colorbar
colormap(jet)
subplot(132)
imagesc(lesion_map)
title('graded lesion data, level of intactness')
colorbar
colormap(jet)
subplot(133)
imagesc(1-relative_lesions)
title('graded lesion data, level of intactness double check')
colorbar
colormap(jet)

RHO1 = corr(relative_lesions,'type','Pearson');
figure
imagesc(RHO1)
ax = gca;
ax.YTick = [1:1:size(BA_label,2)];
ax.YTickLabel = BA_label;
colorbar
ax = gca;
ax.XTick = [1:1:size(BA_label,2)];
ax.XTickLabel = BA_label;
ax.XTickLabelRotation = 90;
colorbar
colormap(jet)
title('correlation between lesion size')

%% 

new_RHO1 = tril(RHO1,-1);
figure
imagesc(new_RHO1)
colorbar
colormap(jet)
ax = gca;
ax.YTick = [1:1:size(BA_label,2)];
ax.YTickLabel = BA_label;
colorbar
ax = gca;
ax.XTick = [1:1:size(BA_label,2)];
ax.XTickLabel = BA_label;
ax.XTickLabelRotation = 90;
colorbar
colormap(jet)
title('correlation between lesion size')


dd_lc = double(new_RHO1<0.3 & new_RHO1>0);
max_pairs_lc = size(find(dd_lc),1);
pairs_lc = zeros(2,max_pairs_lc);

[r,c,v] = find(dd_lc);
pairs_lc(1,:) = r';
pairs_lc(2,:) = c';

max_pairs = 50;
pairs = pairs_lc(:,[3:6:max_pairs_lc]);



%sorting confusion matrix
[MATreordered,MATindices,MATcost]  = reorderMAT(RHO1,1000000,'line');
label_column=BA_label(MATindices);
label_row=BA_label(MATindices);
figure
imagesc(MATreordered)
ax = gca;
ax.XTick = [1:1:size(BA_label,2)];
ax.YTick = [1:1:size(BA_label,2)];
ax.XTickLabel = label_column;
ax.YTickLabel = label_row;
ax.XTickLabelRotation = 90;
set(gca,'FontSize',12,'FontWeight','bold')

colorbar
colormap(jet)
set(colorbar,'FontSize',12,'FontWeight','bold')
title('Correlation of lesion size','FontSize',16,'FontWeight','bold' )


t = sum(relative_lesions,2);
[Y,I] = sort(t,'descend');
relative_lesions_sorted  = relative_lesions(I,:);

figure
imagesc(relative_lesions_sorted)
title('graded lesion data sorted descend, level of lesion')
colormap(jet)

t = sum(lesion_map,2);
[Y,I] = sort(t,'ascend');
lesion_map_sorted  = lesion_map(I,:);
figure
imagesc(lesion_map_sorted)
title('graded lesion data sorted ascend, level of intactness')
colorbar
colormap(jet)

%thresholding-binarization
for i = 1:n_reg
    [r,c,v] = find(relative_lesions(:,i));
    med(i) = median(v); %median value of lesion
end

data_binary = zeros(n_pat,n_reg);
for i = 1:n_reg
    data_binary(:,i) = double((relative_lesions(:,i)<med(i))); %intactness
end

figure
subplot(121)
imagesc(1-relative_lesions)
colorbar
subplot(122)
imagesc(data_binary)
colorbar
colormap(jet)



n_groundT = max_pairs; %I have 50 GT models in this case
score = zeros(n_pat,max_pairs); %I have 50 GT models in this case, default is all zeros, which means bad performance and 1 means good

scorebin= zeros(n_pat,max_pairs);
for k = 1:n_groundT
    first = pairs(1,k);
    second = pairs(2,k);
    
    for h = 1:n_pat
       
            score(h,k) = lesion_map(h,first)*lesion_map(h,second);
        if data_binary(h,first)== 1 && data_binary(h,second) == 1  %try with double lesions
            scorebin(h,k) = 1;
        end
    end
end


%% ESTIMATED MSA - ACCURACY AND YOUDEN
accuracy_L1o = zeros(1,n_groundT); %I have one accuracy for every groundtruth, i.e. for every region
youden = zeros(1,n_groundT);
accuracy_new = zeros(1,n_groundT);
rmse_total = zeros(1,n_groundT);
p = length(data_binary(:,1));
tr = p-1;
combos = combntns(1:p,tr);

maxerror = 0.01; 
for k = 1:n_groundT
    %
    count = 0;
        realhealthy = 0;
        realneglect = 0;
        unrealhealthy = 0;
        unrealneglect = 0;
    
    for in = 1:p
        
        training_label_vector=score(combos(in,:),k);
        training_instance_matrix=lesion_map(combos(in,:),:); %i need to use intactness
        
        testing_label_vector=zeros(p-tr);
        testing_instance_matrix=lesion_map(p-in+1,:);
        
       t =  fitrtree(training_instance_matrix,training_label_vector);
        predicted_label = predict(t,testing_instance_matrix);
        
        e_bis(in)=predicted_label-score(p-in+1,k); %predicted-real
        
        if abs(e_bis(in)) <=maxerror
            count = count +1;
        else count = count;
        end
        
               if score(p-in+1,k)== 0
                    if predicted_label == 0
                        realneglect=realneglect+1;
                    else unrealhealthy=unrealhealthy+1;
                    end
                end
                if score(p-in+1,k)==1
                    if predicted_label == 1
                        realhealthy=realhealthy+1;
                    else unrealneglect=unrealneglect+1;
                    end
                end
        
    end
    rmse_total(k) = sqrt((sum(e_bis.^2))/p);
    accuracy_L1o(k) = (count/p)*100;
    
        accuracy_new(k) = ((realhealthy+realneglect)/p)*100;
        specificity=realhealthy/(realhealthy+unrealneglect);
        sensitivity=realneglect/(realneglect+unrealhealthy);
        youden(k)=sensitivity+specificity-1;
    
end



clear t

%% estimated MSA with bootstrap
% %
nb_bootest=100; %number of bootstrap samples
num_perms = 6500;
resultEpbootest = msaE_get_samples_cfg(num_perms,n_reg);%conf wise format
permsEpbootest = resultEpbootest.perms;
configurationsbootest = resultEpbootest.configs;
permsEpbootest.matrix;
predicted_label_svmbootest = zeros(size(configurationsbootest,1),n_groundT);


for k = 1:n_groundT %1 groundtruth model at a time
       clear t bigmatboot bigmat
    
   [sortdata, sortid]=sort(scorebin(:,k))

    clear lesion_mapre datatmp scoretmp scorere lesion_maptmp
     lesion_maptmp=lesion_map(sortid,:);
    scoretmp=int64(score(sortid,k)*10);
    scoretmp=double(scoretmp);
    lesion_mapre=lesion_maptmp;%(1:300,:);
    scorere(:,k)=scoretmp;%(1:300);
    
   
   
    imbal(k)=sum(scorere(:,k))/length(lesion_mapre(:,1));
    afindbal=find(sortdata==0);
    bfindbal=find(sortdata==1);
    maxminicase=max(afindbal);
   
      bigmat=[lesion_mapre, scorere(:,k)]; %matrice avec les pateints en ligne et en colone, les regions puis le score clinique
   
    p = length(lesion_mapre(:,1));
    for samplbootest=1:nb_bootest %pour tous les echantillons de bootstrap
        resultEpbootest = msaE_get_samples_cfg(num_perms,n_reg);%conf wise format
        permsEpbootest = resultEpbootest.perms;
        configurationsbootest = resultEpbootest.configs;
        clear nouvvect
         for pat=1:200 
            nouvvect(pat)=randperm(maxminicase,1);
        end
        for pat=201:581
            nouvvect(pat)=randperm(p-maxminicase,1)+maxminicase;
        end
        
        n_reg = size(lesion_mapre,2);
    


        bigmatboot(:,:)=bigmat([nouvvect],:); %matrice avec en ligne les patients avec remise du tirage bootstrap
     
        training_label_vectorbootest=bigmatboot(:,end);
        training_instance_matrixbootest=bigmatboot(:,[1:end-1]);
        testing_label_vectorbootest = zeros(size(configurationsbootest,1),1);
        testing_instance_matrixbootest = configurationsbootest;
        
        tbootest = fitrtree(training_instance_matrixbootest,training_label_vectorbootest);
        predicted_label_svmbootest = predict(tbootest,testing_instance_matrixbootest);
        
       msaDatabootest=struct('num_elements',n_reg,'configs',configurationsbootest,'perfs',predicted_label_svmbootest,'perms',permsEpbootest); %conf wise format
        
        cperfbootest = msaE_cfg2prm_perfs(msaDatabootest); %perm wise format
        msaData_Ebootest=struct('num_elements',n_reg,'max_perturbed',n_reg,'perms',permsEpbootest,'perfs',cperfbootest); %perm wise -%estimated MSA
        cvEpbootest = msaE_compute_CVs(msaData_Ebootest);
        
        sommabootest(samplbootest)=sum(abs(cvEpbootest.sh));
        m_cvEpbootest(samplbootest,:,k)=cvEpbootest.sh/sommabootest(samplbootest);
        s_cvEpbootest(samplbootest,:,k)=cvEpbootest.se/sommabootest(samplbootest);
        
        
        
        k
        
        
        
    end
    
    
    
end


essaiboot=nanmean(m_cvEpbootest);
essaibootsq=squeeze(essaiboot);
essaisdtboot=nanmean(s_cvEpbootest);
essaisdtbootsq=squeeze(essaisdtboot);

m_cvEp_normboot = essaibootsq/sum(sum(abs(essaibootsq)));
s_cvEp_normboot = essaisdtbootsq/sum(sum(abs(essaibootsq)));


TPboot = 0;
FPboot = 0;
TNboot = 0;
FNboot = 0;
matrixboot = m_cvEp_normboot;



errboot = 0;
goodboot = 0;
sommaboot = 0;


for k = 1:n_groundT
    for m = 1:n_reg
        if  m == pairs(1,k) | m == pairs(2,k)
            goodboot = goodboot + abs(matrixboot(m,k));
            if  matrixboot(m,k) ~= 0
                TPboot = TPboot+1; FNboot = FNboot;
            else
                FNboot = FNboot +1; TPboot = TPboot;
            end
        elseif  m ~= pairs(1,k) & m ~= pairs(2,k)
            
            errboot = errboot + abs(matrixboot(m,k));
            if  matrixboot(m,k) ==0% 0.001
                TNboot = TNboot+1; FPboot = FPboot;
            else
                FPboot = FPboot +1; TNboot = TNboot;
            end
            
        end
    end
end


TPboot
TNboot
FNboot
FPboot

accuboot = ((TPboot+TNboot)/(TPboot+ TNboot+ FPboot +FNboot))*100
sommaboot= TPboot+ TNboot+ FPboot +FNboot
sensitivityboot = TPboot/(TPboot+FNboot)*100
specificityboot = TNboot/(TNboot+FPboot)*100

sommaboot = 0;
for i =1:n_groundT
    
    if pairs(1,i) == pairs(2,i)
    sommaboot = sommaboot + abs(matrixboot(pairs(1,i),i));
    else 
        sommaboot = sommaboot + abs(matrixboot(pairs(1,i),i)) + abs(matrixboot(pairs(2,i),i));
    end
end

sommaboot
goodboot

err2boot = (sum(sum(abs(matrixboot))) - sommaboot)*100
errboot*100

matrix_nboot = zeros(2,n_groundT);

for k= 1:n_groundT
    matrix_nboot(:,k) = matrixboot(pairs(:,k),k);
end
% matrix_n=matrix_n(:,[1:4,6:67,69:77,79:98,100]);
% figure
% imagesc(matrix_n)
% colorbar
% colormap(jet)

   
matrix_nboot = matrix_nboot(:);
accu_nboot =((TPboot*(1-var(matrix_nboot))+TNboot*(1-errboot/100))/(TPboot+ TNboot+ FPboot +FNboot))*100

figure
imagesc(essaibootsq)
set(gca,'FontName','Arial','FontSize',16)
title('MSA bootstrap','FontSize',24)
ax = gca;
ax.XTick = [1:1:n_groundT];
ax.YTick = [1:1:n_reg];
ax.XTickLabel = BA_label;
ax.YTickLabel = BA_label;
ax.XTickLabelRotation = 90;
xlabel('lesioned','FontSize',16)
ylabel('inferred','FontSize',16)
colorbar
set(colorbar,'FontSize',12,'FontWeight','bold')
colormap(jet)


save('LC_doublesourcegradedmultisynergyMSA.mat')


