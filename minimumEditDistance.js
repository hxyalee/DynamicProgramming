function minimumEditDistance(s1, s2){
	let memo = new Array(s1.length + 1);
	for(let i = 0; i < memo.length; i++){
		memo[i] = new Array(s2.length + 1);
		memo[i][0] = i;
	}	
	for(let i = 0; i < memo[0].length; i++) memo[0][i] = i;
	for(let i = 1; i < memo.length; i++){
		for(j = 1; j < memo[i].length; j++){
			if(s1.charAt(i - 1) === s2.charAt(j - 1)) 
				memo[i][j] = memo[i-1][j-1]
			else
				memo[i][j] = Math.min(memo[i-1][j], memo[i-1][j-1], memo[i][j-1]) + 1;
		}
	}	
	return memo[s1.length][s2.length]
}
