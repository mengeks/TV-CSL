# Load necessary library
library(splines)

# Generate a sequence of x values (1-dimensional variable)
x <- seq(0, 10, length.out = 100)

# Generate the natural cubic spline basis matrix with 4 degrees of freedom
ns_basis <- ns(x, df = 4)

# Plot each basis function
plot(x, ns_basis[,1], type = 'l', col = 'red', lwd = 2, ylim = c(-0.5, 1.5),
     ylab = 'Basis function value', xlab = 'x', main = 'Natural Cubic Spline Basis Functions')
lines(x, ns_basis[,2], col = 'blue', lwd = 2)
lines(x, ns_basis[,3], col = 'green', lwd = 2)
lines(x, ns_basis[,4], col = 'purple', lwd = 2)
legend('topright', legend = paste("B", 1:4, sep=""), col = c('red', 'blue', 'green', 'purple'), lwd = 2)
