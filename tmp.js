let solution = (arr, n) => {
	if(n == 0) return 1;
	return rec(arr, n, arr.length - 1);
}

let rec = (arr, n, i) => {
	if(n == 0) return 1;
	if(n < 0) return 0;
	if(i < 0) return 0;
	let count = 1;
	let tmp = n;
	while(Math.floor(n/arr[i])){
		count ++;
		
	} 
}