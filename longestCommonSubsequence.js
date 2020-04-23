// Given two strings, determine the longest commmon subsequence
// More in README.md
// Example input: "abcdefghij", "cdgi"
// Example output: "cdgi" => longest common subsequence since they are 'in order'


// Exponential time algorithm
function longestCommonSubsequenceRecursive(s1,s2, i, j){
	if(i == s1.length || j == s2.length) return 0;
	if(s1.charAt(i) === s2.charAt(j)) 
		return 1 + longestCommonSubsequenceRecursive(s1,s2, i+1,j+1)
	else 
		return Math.max(longestCommonSubsequenceRecursive(s1, s2, i+1, j), longestCommonSubsequenceRecursive(s1, s2, i, j+1))
}

// Dynamic programming
/*
 * s1 starts from character 1 ... n
 * s2 starts from character 1 ... m
 * if they match, the substring length is of length 1 + the previous diagonal
 * Otherwise, get the max of the previous i or j
 * the reason for +1 in the array and -1 in the char
 * is to allocate initial values for memoization
 * Time complexity: O(mn)
*/
function longestCommonSubsequence(s1,s2){
	let memo = new Array(s1.length + 1);
	for(let i = 0; i < memo.length; i++){
		memo[i]= new Array(s2.length + 1).fill(0);
	}
	for(let i = 1; i < memo.length; i++){
		for(let j = 1; j < memo[i].length; j++){
			// Look above for why i - 1?
			if(s1.charAt(i-1) === s2.charAt(j-1)) 
				memo[i][j] = 1  + memo[i-1][j-1]
			else
				memo[i][j] = Math.max(memo[i-1][j], memo[i][j-1])
		}
	}
	return memo[s1.length][s2.length]
}
console.log(longestCommonSubsequence('bd','abcd'))