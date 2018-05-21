

load( "Molecular Subgroup Adaptation Matrix.RData")

XX <- log2( RNA.TPM[,rownames( COEF.MATRIX )[-1]] + 0.25) - log2( 0.25 ) # rows = samples, columns = genes (Ensembl IDs) where RNA.TPM is unscaled TPM data

YY <- apply( XX , 2 , function(x) (x-min(x))/(max(x)-min(x))   ) # range scale data

RES <- cbind( 1 , YY ) %*% COEF.MATRIX # apply regularized regression coefficients

PRED <- colnames( COEF.MATRIX )[apply( RES , 1 , function(x) which.max(x))] 

table( PRED )









