/*
	NOT EFFICIENT! 
*/
let recursive = (arr, n) => {
	if(n === 0) return 1;
	return recursive_helper(arr, n, arr.length - 1);
}
let recursive_helper = (arr, n, i) => {
	if(n == 0) return 1;	// empty subset
	if(n < 0) return 0;		// nothing adds up to n
	if(i < 0) return 0; 	// index out of bounds

	if(n < arr[i])
		return recursive_helper(arr, n, i-1);
	else		
		return recursive_helper(arr, n - arr[i], i-1) + recursive_helper(arr, n, i - 1);

}




console.log(recursive([1,3,5,7],8))