data=cell(9,1);
data{1} = 'Ackley';
data{2} = 'EllipticE3';
data{3} = 'GriewankG3';
data{4} = 'RastriginR3';
data{5} = 'RastriginS3';
data{6} = 'RosenbrockR3';
data{7} = 'SchwefelS3';
data{8} = 'SphereS3';
data{9} = 'WeierstrassW3';

xmin=[-32,-100,-600,-5.12,-5,-30,-100,-100,-0.5];
xmax=[32,100,600,5.12,5,30,100,100,0.5];
true_min=[0,82,87,5,98,11,83,111,118];
errgoal=[0.001,500,1,1,1,1,0.01,0.001,1];
cutoff_length=400000;

funcs = numel(data);
numRuns=10;
seed=123;

fevDef = zeros(numRuns*funcs,1);
fevILS2 = zeros(numRuns*funcs,1);

k=1;

for i = 1 : numRuns
	for j = 1 : funcs

		% default
		c1=2.05;
		c2=2.05;
		rad=1;
		lbd=1.0;
		tmax=5;

		[~, ~, fevDef(k), ~] = MPSO(char(data(j)),xmin(j),xmax(j),true_min(j),errgoal(j),cutoff_length,seed,c1,c2,rad,lbd,tmax);
		
		if fevDef(k) > 400000
			fevDef(k) = 400000;
		end

		% ParamILS Run 2
		c1=2.1;
		c2=2.1;
		rad=2;
		lbd=0.8;
		tmax=12;
	
		[~, ~, fevILS2(k), ~] = MPSO(char(data(j)),xmin(j),xmax(j),true_min(j),errgoal(j),cutoff_length,seed,c1,c2,rad,lbd,tmax);
		
		if fevILS2(k) > 400000
			fevILS2(k) = 400000;
		end
		
		seed = seed + 1;
		k = k + 1;
		
	end
end

fevDef = log10(fevDef);
fevILS2 = log10(fevILS2);

figure
grid on
set(gca,'fontsize',12)
c1 = cdfplot(fevDef);
hold on
c2 = cdfplot(fevILS2);
set(c1,'color','r')
set(c1,'linewidth',2)
set(c2,'color','b')
set(c2,'linewidth',2)
title('CDF of FEvals for 10 ind. runs each on a set of 9 functions')
xlabel('log10 of FEvals')
ylabel('Proportion of trials')
l = legend('Default Config.','Run 2 Config.');
set(l,'Location','southeast');
