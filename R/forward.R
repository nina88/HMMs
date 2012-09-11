
forward <- function(A, B, C, D, myArray, myDims){
	.Call( "forward", A, B, C, D, myArray, myDims, PACKAGE = "HMMs" )
}

