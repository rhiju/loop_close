ADD_TRANSLATIONS = 1; % Do 6D translations & rotations (SE(3)) vs. 3D rotations (SO(3))
%[xyz, Ixyz] = get_aligned_xyz( 'xyz_G.txt' ); % use G nucleobase
[xyz, Ixyz] = get_aligned_xyz( 'xyz_O3p_P_O5p.txt' );

% this is goofy -- use P-centered coordinate system, 'stub' -- should not
% match exact expressions, but might be close in SE(3) case.
%xyz = get_Pcentered_xyz( 'xyz_O3p_P_O5p.txt' ); 

N = size(xyz,1);
Nval = 10000000;
M = get_uniform_random_matrices( Nval );
xyz_rot = tprod( xyz, [1 -1], M, [-1,2,3] ); % tensor product for speed.
xyz_transform = xyz_rot;
if ADD_TRANSLATIONS % from SO(3) to SE(3)
    TRANS_MAG = 4.0; % in Angstroms
    translations = repmat( TRANS_MAG * randn( 1, 3, Nval ), [N,1,1] ) ; 
    xyz_transform = xyz_transform + translations;
end
rmsd = squeeze( sqrt( sum( sum( ( xyz_transform - repmat( xyz, [1, 1, Nval] ) ).^2 / N ) ) ) );

% copied from analyze_chainbreak
% /Users/rhiju/Dropbox/projects/RNA/nn_rule_chain_closure/analysis/analyze_chainbreak.m
delx = 0.05;
x = [0:delx:100];
h = hist( rmsd, x ) / Nval / delx; % dp/dr.
loglog( x, h, 'o' ); hold on; 

% compare to expected. 
% Note, I derived this again from scratch, but this exactly matches my
% previous computation for SO(3) rotations:
% ~/src/rosetta/demos/pilot/rhiju/rb_entropy/rb_entropy_analysis_script.m.
pred_curve = N^(3/2) * x.^2 / (2*pi) / sqrt(Ixyz(1)+Ixyz(2)) / sqrt(Ixyz(1)+Ixyz(3)) / sqrt(Ixyz(2)+ Ixyz(3));
if ADD_TRANSLATIONS
    pred_curve = N^(3/2) * (1/sqrt(2*pi)/TRANS_MAG)^3 * x.^5 * pi/8 / sqrt(Ixyz(1)+Ixyz(2)) / sqrt(Ixyz(1)+Ixyz(3)) / sqrt(Ixyz(2)+ Ixyz(3));
end
plot( x, pred_curve, 'k' ); 
% new term to check
pred_curve_corrected = pred_curve .* ( 1 + (1/48) * x.^2 * N *( 1/(Ixyz(1)+Ixyz(2)) + 1/(Ixyz(1)+Ixyz(3)) + 1/(Ixyz(2)+Ixyz(3))));
if ADD_TRANSLATIONS
    pred_curve_corrected = pred_curve .* ( 1 + (1/96) * x.^2 * N *( 1/(Ixyz(1)+Ixyz(2)) + 1/(Ixyz(1)+Ixyz(3)) + 1/(Ixyz(2)+Ixyz(3))));
end
plot( x, pred_curve_corrected, 'b' );
hold off
legend( 'Measured', 'Predicted','Predicted with next order correction');
xlabel( 'RMSD (Angstroms)');
ylabel( 'Probability distribution (dP/dRMSD)');
title( 'RMSD distribution in SO(3)' )
if ADD_TRANSLATIONS; title( 'RMSD distribution in SE(3)' ); end;

