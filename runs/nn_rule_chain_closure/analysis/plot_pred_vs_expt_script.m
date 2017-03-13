% tag,  chain_closure_vals (M), nupack_val, num_strands, is_symmetric
% Note that NUPACK give 1.62 kcal/mol for singlet G-C --> may assume 55 M as reference concentration? (intead of 1 M, with 3.4 kcal/mol init in Serra-Turner)
data = {...
{'stacked_pair/gc_gc',[ 5419.220950, 9569.737107, 2186.908872, 4224.555449],-1.72, 2, 'w','o'},...   % check  1.62 - 3.4 = -1.78?
{'bulge/gc_gac',[  3.736660, 38.194024, 21.961424], 2.02, 2,'k','s'},...                            % check  1.62 - 3.4 + 3.9 = 2.12
{'triloop/guuuc',[  0.294882, 0.370346, 0.799784, 0.408375], 4.12, 1, 'g','o'},...                   % check  4.1 (Table III, Serra-Turner)
{'tetraloop/guuuuc',[  0.488806, 0.148445, 1.450323, 1.280446], 4.20, 1, [0,0.8,0],'o'},...                 % check  4.9 - 0.7 = 4.2 (Table III, Serra-Turner, UUmismatch bonus)
{'pentaloop/guuuuuc',[  0.244949, 0.398444, 0.294847, 1.548277], 3.70, 1, [0,0.6,0],'o'},...               % check  4.4 - 0.7 = 3.7 (Table III, Serra-Turner)
{'decaloop/guuuuuuuuuuc',[  0.068007, 0.123614, 0.228074, 0.887496], 4.60, 1,[0,0.5,0],'o'},...           % check 5.2 + 1.75 * 0.606 * log( 10/9) - 0.7 = 4.6
{'two_way/guuuc_guuuc',[  3.897240, 1.272395], 2.32, 2,'c' ,'s'},...
{'two_way/guuc_guuuuc',[0.283563, 0.606822], 2.92, 2,[0 0.7 1],'s'},...
{'two_way/guc_guuuuuc',[0.091259, 0.380388], 5.52, 2,[0 0.4 1],'s'},...
{'two_way/gc_guuuuuuc',[0.188561, 0.103352], 6.62, 2,[0 0.2 1],'s'},...                                     % check  1.62 + 5.0 = 6.62 (Table III, Serra_Turner) -- no stacking
{'two_way/guc_guc',[ 12.845838 ,  61.138562 ], 2.42, 2,[0.6 0 0],'s'},...              % check: 
{'three_way/guc_guc_guc',[  0.147206, 0.141576, 0.314864,0.403610], 9.34, 3,[0.6 0 0],'v'},...              % check: 1.62*2 + 4.6 + 0.4*3 +0.1*3 = 9.34
{'four_way/guc_guc_guc_guc', [ 3.438419, 4.680015, 6.192666, 5.127130], 11.46, 4, [0.7 0 0],'d'},...
{'five_way/guc_guc_guc_guc_guc', [ 0.265101, 0.115717, 0.818074, 0.815670 ], 13.57, 5, [0.8 0 0],'p'},...
{'six_way/guc_guc_guc_guc_guc_guc',[  0.726702, 1.454485, 3.076583, 1.035727], 15.69, 6,'r','h'},...  % check: 1.62*5 + 4.6 + 0.4*6 +0.1*6 = 15.7
}


kT = 0.606; % kcal/mol
for i = 1:length( data )
    legends{i} = data{i}{1};
    delG_NUPACK = data{i}{3};
    num_strands = data{i}{4};
    C_eff_NUPACK = exp( -( delG_NUPACK - 1.62*( num_strands - 1 ) - 4.04 ) / kT )
    C_eff_FARFAR = data{i}{2};
    plot( C_eff_NUPACK * ones(1,length(C_eff_FARFAR)), C_eff_FARFAR,...
    	  'ko','markerfacecolor',data{i}{5},...
	  'marker',data{i}{6});
    hold on
end
plot( 10.^[-4:0.1:6], 10.^[-4:0.1:6],'k-' );
xlabel( 'C_{eff}(NUPACK,Serra-Turner,1995) (M)' ); ylabel( 'C_{eff}(FARFAR) (M)' );
hold off
legend( legends ,'interpreter', 'none', 'Location','NorthWest' );
set(gca,'xscale','log','yscale','log','xgrid','on','ygrid','on');
