setup_test_functions;

boundary = {'periodic','flat','linear','cubic'};

% check 1D
x = 0.1211;
[val, deriv_analytic, deriv_numeric] = get_numerical_deriv( F_1D, x, minval, binwidth, 'periodic' );
x = -1.23;
for i = 1:length( boundary )
[val, deriv_analytic, deriv_numeric] = get_numerical_deriv( F_1D, x, minval, binwidth, boundary{i} );
end
fprintf( '\n' );

% check 2D
x = [0.1211,0.5555];
[val, deriv_analytic, deriv_numeric] = get_numerical_deriv( F, x, minval, binwidth, 'periodic' );
fprintf( '\n' );

% check 3D
x = [0.1211,0.5555,0.4461];
[val, deriv_analytic, deriv_numeric] = get_numerical_deriv( F_3D, x, minval, binwidth, 'periodic' );

x = [-1.23,1.5,-0.233333];
for i = 1:length( boundary )
[val, deriv_analytic, deriv_numeric] = get_numerical_deriv( F_3D, x, minval, binwidth, boundary{i} );
end
fprintf( '\n' );

% check 4D
x = [0.1211,0.5555,0.4461,0.0022];
[val, deriv_analytic, deriv_numeric] = get_numerical_deriv( F_4D, x, minval, binwidth, 'periodic' );
