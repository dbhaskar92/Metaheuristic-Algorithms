
func_name='EllipticE2';
xmin=-100;
xmax=100;
true_min=-51;
errgoal=500;
cutoff_length=400000;

numRuns=10;
seed=123;

fevDef = zeros(numRuns,1);
fevILS0 = zeros(numRuns,1);
fevILS1 = zeros(numRuns,1);
fevILS2 = zeros(numRuns,1);

for i = 1 : numRuns

	% default
	c1=2.05;
	c2=2.05;
	rad=1;
	lbd=1.0;
	tmax=5;

	[optimum, gbest, fevDef(i), tElapsed] = MPSO(func_name,xmin,xmax,true_min,errgoal,cutoff_length,seed,c1,c2,rad,lbd,tmax);
	
	if fevDef(i) > 400000
		fevDef(i) = 400000;
	end
	
	% ParamILS Run 0
	c1=2.05;
	c2=2.15;
	rad=3;
	lbd=0.8;
	tmax=5;
	
	[optimum, gbest, fevILS0(i), tElapsed] = MPSO(func_name,xmin,xmax,true_min,errgoal,cutoff_length,seed,c1,c2,rad,lbd,tmax);
	
	if fevILS0(i) > 400000
		fevILS0(i) = 400000;
	end
	
	% ParamILS Run 1
	c1=2.075;
	c2=2.05;
	rad=5;
	lbd=0.5;
	tmax=5;
	
	[optimum, gbest, fevILS1(i), tElapsed] = MPSO(func_name,xmin,xmax,true_min,errgoal,cutoff_length,seed,c1,c2,rad,lbd,tmax);
	
	if fevILS1(i) > 400000
		fevILS1(i) = 400000;
	end

	% ParamILS Run 2
	c1=2.1;
	c2=2.1;
	rad=2;
	lbd=0.8;
	tmax=12;
	
	[optimum, gbest, fevILS2(i), tElapsed] = MPSO(func_name,xmin,xmax,true_min,errgoal,cutoff_length,seed,c1,c2,rad,lbd,tmax);
	
	if fevILS2(i) > 400000
		fevILS2(i) = 400000;
	end
    
    seed = seed + 1;

end

x = 1:1:400000;
y = 400000.*ones(400000,1);

figure
grid on
set(gca,'fontsize',11)
scatter(fevDef,fevILS0,25,'r','filled')
hold on
scatter(fevDef,fevILS1,25,'g','filled')
hold on
scatter(fevDef,fevILS2,25,'b','filled')
hold on
plot(x,x,'color','k')
hold on
plot(x,y,'color','k')
hold on
plot(y,x,'color','k')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([1,500000,1,500000])
title('Runlengths for High Conditioned Elliptic Function (10 ind. runs)')
xlabel('log10 of FEvals for Default Configuration')
ylabel('log10 of FEvals for ParamILS Configuration')
l = legend('Run 0','Run 1','Run 2');
set(l,'Location','southwest');
