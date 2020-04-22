function longestIncreasingSubsequence(seq){
	let memo = new Array(seq.length);
	// The least longest increasing subsequence is 1 for each element
	memo.fill(1);
	for(let i = 1; i < seq.length; i++){
		for(let j = 0; j < i; j ++){
			if(memo[i] > memo[j]) continue;
			if(seq[j] < seq[i]) memo[i] = memo[j] + 1;
		}
	}
	return Math.max.apply(null, memo)
}
