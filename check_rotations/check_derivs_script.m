% 'rotation-vector' space. Will look at points in sphere with radius pi.
%gv = [-pi:0.01:pi]; too fine.
gv = [-pi:0.05:pi];
% X,Y,Z = v_x, v_y, v_z (in angle notation)
[X,Y,Z] = meshgrid(gv);
E = ( X - pi/3 ).^4 + ( Y ).^2  + ( Z ).^4;

% draw it
figure(4);
clf;
FV = isosurface(X,Y,Z,E,0.1);
p = patch( FV );
p.FaceColor = 'red';
p.EdgeColor = 'none';
axis([-pi pi -pi pi -pi pi ]);
axis vis3d
camlight; lighting phong

% Now cycle through random rotations in this space
N_rot = 20;
M_random = get_uniform_random_matrices( N_rot );
V_random = SpinCalc( 'DCMtoEV', M_random, 1.0e-6, 0 );
V_random(:,4) = V_random(:,4) - 360*(V_random(:,4)>180);
dEdomega_analytical = [];
dEdomega_numerical  = [];
count = 0;
figure(5);  clf;
interpolation_scheme = 'cubic';
for j = 1:N_rot
    %Xq = 0.5; Yq = 0.0; Zq = 0.0;
    Xq = V_random(j,1)*V_random(j,4)*pi/180.0;
    Yq = V_random(j,2)*V_random(j,4)*pi/180.0;
    Zq = V_random(j,3)*V_random(j,4)*pi/180.0;

    Eq = interp3(X,Y,Z,E,Xq,Yq,Zq,interpolation_scheme);
    delq = 1.0e-8;
    dEdv(1) = (interp3(X,Y,Z,E,Xq+delq,Yq     ,Zq     ,interpolation_scheme)-Eq)/delq;
    dEdv(2) = (interp3(X,Y,Z,E,Xq     ,Yq+delq,Zq     ,interpolation_scheme)-Eq)/delq;
    dEdv(3) = (interp3(X,Y,Z,E,Xq     ,Yq     ,Zq+delq,interpolation_scheme)-Eq)/delq;
    [Xq, Yq, Zq, dEdv]

    % Used in analytical solution below
    V_norm = [Xq,Yq,Zq];
    V = norm( V_norm ); % angle of the rotation, in radians
    V_norm = V_norm/ V; % unit vector for axis

     % get Matrix associated with this V, too
     Vq = [Xq,Yq,Zq]*(180.0/pi);
     Vq(4) = norm(Vq);
     Vq(1:3) = Vq(1:3)/Vq(4);
     M = SpinCalc( 'EVtoDCM', Vq, 1.0e-8, 0 );

    % cycle through random axes of rotation for perturbation
    %omega = [0, 0, 1];
    N = 10;
    for i = 1:N
        count = count+1;
        omega = randn(1,3); % random direction
        omega = omega/norm(omega);

        domega = 1.0e-8;
        m_omega = vrrotvec2mat( [-omega, domega]); 

        M_perturb = M*m_omega;

        V_perturb = SpinCalc( 'DCMtoEV', M_perturb, 1.0e-8, 0 );
        Xp = V_perturb(:,1)*V_perturb(:,4)*(pi/180.0);
        Yp = V_perturb(:,2)*V_perturb(:,4)*(pi/180.0);
        Zp = V_perturb(:,3)*V_perturb(:,4)*(pi/180.0);
        Eq_perturb = interp3(X,Y,Z,E,Xp,Yp,Zp,interpolation_scheme);
        dEdomega_numerical(count) = ( Eq_perturb - Eq ) / domega;

        % compare to my analytical solution.
        F1 = dot(dEdv, V_norm ) * V_norm + ...
             (V/2)*cot( V/2 )*( dEdv - dot(dEdv,V_norm)*V_norm ) ...
             - (V/2)*cross(dEdv,V_norm);
        dEdomega_analytical(count) = dot( omega, F1 );
    end
    %Eq_perturb
    %V_norm
    %plot( dEdomega_numerical, dEdomega_analytical ,'o' ); 
    %xlabel( 'numerical');ylabel('analytical');pause
end
plot( dEdomega_numerical, dEdomega_analytical ,'o' );
hold on
plot( [-100 100],[-100 100], 'k' );
hold off
axis( [-100 100 -100 100] );
    xlabel( 'numerical');ylabel('analytical');

