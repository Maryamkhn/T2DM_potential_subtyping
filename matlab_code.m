%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code related to the analyzes performed in paper entitled "Unraveling the Molecular Heterogeneity in type 2 Diabetes: A Potential     %
% Subtype Discovery Study Followed by Metabolic Modeling                                                                                           %
% Maryam Khoshnejat1, Kaveh Kavousi1, Ali Mohammad Banaei- Moghaddam, Ali Akbar Moosavi-Movahedi                                                   %
% Laboratory of Complex Biological Systems and Bioinformatics (CBB), Department of Bioinformatics, Institute of Biochemistry and Biophysics (IBB), %
% University of Tehran, Tehran, Iran                                                                                                               %
% Corresponding author: Kaveh Kavousi (kkavousi@ut.ac.ir)                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Matlab codes
initCobraToolbox

[Deseq_FPKM, genes]=xlsread('normalized_expression_data.xlsx');
Deseq_fpkm=genes;
for i=2:57821
for j=2:155
Deseq_fpkm{i,j}=Deseq_FPKM(i-1,j-1);
end
end

NGT_Deseq_FPKM_log2=log2(cell2mat(Deseq_fpkm(2:end,2:92))+1);
T2D_Deseq_FPKM_log2=log2(cell2mat(Deseq_fpkm(2:end,93:155))+1);


NGT_expression.gene=genes(2:end,1);
T2D_expression.gene=genes(2:end,1);

NGT_expression.rawValue=NGT_Deseq_FPKM_log2;
T2D_expression.rawValue=T2D_Deseq_FPKM_log2;
load model

 expressionRxns_NGT=zeros(5740,91);
 for i=1:91
 NGT_expression.value=NGT_expression.rawValue(:,i);
[expressionRxns_NGT(1:5740,i) ~] = mapExpressionToReactions(model, NGT_expression,'minsum');
 end
 
 expressionRxns_T2D=zeros(5740,91);;
 for i=1:63
 T2D_expression.value=T2D_expression.rawValue(:,i);
 [expressionRxns_T2D(1:5740,i) ~] = mapExpressionToReactions(model, T2D_expression,'minsum');
 end

 
%%%%%%%%%personalized metabolic modeling
%%%%%do these codes below for NGT too
Eflux_patient=model;
Eflux_patient=changeObjective (Eflux_patient, 'HMR_6916');
j=1;
errorbar_T2D=[];
for i=1:12
try
Eflux_patient = E_Flux(Eflux_patient, expressionRxns_T2D(:,i),max_expression);
Eflux_patient.lb(5740)=0.8*0.0014;
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9024',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9025',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9026',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9027',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9028',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9029',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9030',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9031',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9810',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9811',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9812',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9814',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9084',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9088',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9089',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9101',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9103',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9108',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9126',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9127',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9128',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9130',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9155',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9156',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9163',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9171',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9172',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9201',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9213',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9256',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9314',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9357',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9364',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9365',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9367',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9368',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9386',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9389',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9393',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9429',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9443',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9444',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9445',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9446',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9447',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9448',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9449',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9450',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9451',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9456',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9457',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9458',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9460',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9461',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9681',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9682',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9683',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9684',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9685',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9687',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9689',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9690',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9692',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9693',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9694',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9696',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9697',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9698',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9699',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9700',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9701',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9702',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9707',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9709',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9710',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9711',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9713',0,'u');
Eflux_patient = changeRxnBounds(Eflux_patient,'HMR_9715',0,'u');

[mintest,maxtest]= fluxVariability(Eflux_patient);
T2D_max_flux(1,j)=i;
T2D_max_flux(2:5741,j)=maxtest;
T2D_min_flux(1,j)=i;
T2D_min_flux(2:5741,j)=mintest;
j=j+1; 
catch
errorbar_T2D=[errorbar_T2D,i];
end
end
%%%%%%%%%% t-test 
for i=1:5740
[~,p_c1(i),~,stats_c1(i)]= ttest2(NGT_max_flux(i+1,:),c1_max_flux(i+1,:),'Vartype','unequal');
end
[Q_c1] = mafdr(p_c1,'BHFDR',true);





%%%%%%%%%%%%%%%%% Genetic algorithm/SVM
MZ2=(1:247)';
Y2=cell2mat(deseq2_norm_5_38(2:end,2:end));
grp2=deseq2_norm_5_38(1,2:end);
id2 = grp2idx(grp2);  


options2 = gaoptimset('PopulationSize',500,...
'PopulationType', 'bitstring',...
'Generations',100,...
'FitnessLimit',5,...
'UseParallel',true,...
'Display', 'iter');
            
                 
nVars2 = 247;                          % set the number of features
FitnessFcn2 = {@mysvm,Y2,id2};       % set the fitness function
[x,fval,exitflag,output,population,scores] = ga(FitnessFcn2,nVars2,options2);



%%%%%%%%%%%%fitness function
function classPerformance = mysvm(thePopulation,Y,id)
%svm as fitness score in ag
p1=sum(thePopulation);
FPKM=Y(find(thePopulation),:)';
k=10;
id2 = ismember(id,1);
groups=id2;
cvFolds = crossvalind('Kfold', groups, k);   %% get indices of 10-fold CV
cp = classperf(groups);                      %% init performance tracker
for i = 1:k                                  %% for each fold
testIdx = (cvFolds == i);                %% get indices of test instances
trainIdx = ~testIdx;                     %% get indices training instances
%% train an SVM model over training instances
svmModel = svmtrain(FPKM(trainIdx,:), groups(trainIdx),'Kernel_Function', @(u,v) mysigmoid(u,v,p1));
%% test using test instances
pred = svmclassify(svmModel, FPKM(testIdx,:), 'Showplot',false);
%% evaluate and update performance object
cp = classperf(cp, pred, testIdx);
end
classPerformance=(1-cp.CorrectRate)*100;
end



