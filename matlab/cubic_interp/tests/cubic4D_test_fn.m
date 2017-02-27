function F = cubic4D_test_fn( X , Y, Z, W )
F = 3*X + Y + 3*X.*X - X.*X.*X - 3 * Y.*Y.*Y - 2 * X.*Y + Y .* Z + W.*W.*W - W.*X + 0.5 * W;
