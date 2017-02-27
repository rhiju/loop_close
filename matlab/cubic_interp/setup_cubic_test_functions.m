%test_fn = @cubic_test_fn
test_fn = @periodic_test_fn

clf
minval = 0; binwidth = 0.1;
ctrs = [minval : binwidth : 0.9]; 
[X,Y] = ndgrid( ctrs, ctrs );
ctrs_fine = [-0.2:0.01:1.2];
[Xf,Yf] = ndgrid( ctrs_fine, ctrs_fine );
F  = test_fn( X, Y );
Ff = test_fn( Xf, Yf );

% for 1D scans
F_1D  = F(find(ctrs==0.0),:)';
Ff_1D = Ff(find(ctrs_fine==0.0),:)';
