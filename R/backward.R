
backward <- function(A, B, myArray, myDims){
	.Call( "backward", A, B, myArray, myDims, PACKAGE = "HMMs" )
}


