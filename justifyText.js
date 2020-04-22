
/*
	['Joshua', 'Kim', 'Loves', 'To', 'Code']
	Greedy ==> put as many words as possible in each line
				if the word doesnt fit, put it in next line
			   Lines   Empty Spaces 	Empty^2
			=> Joshua Kim	0				0
			=> Loves To		2				4
			=> Code			6				36
							8				40

	DP 	==> Joshua			4				16
		==> Kim Loves		1				1
		==> To Code			3				9
							8				26 <- Better! since less empty space 
*/

let justify = (arr, maxlength) => {
	let nWords = arr.length;
	// initialise an empty 2D array to store the costs of empty spaces per words
	let justifyArray = new Array(nWords)
	for(let i = 0; i < nWords; i++){
		justifyArray[i] = new Array(nWords)
	}
	// get the cost of empty spaces per word if the words fit in one line
	// i is row
	for(let i = 0; i < nWords; i++){
		// j is column
		let prev = 0;
		for(let j = i; j < nWords; j ++){
			if(prev === 0) prev += arr[j].length;
			else prev += arr[j].length + 1;
			if(prev > maxlength) break;
			// exponent is essential to look if efficient
			justifyArray[i][j] = Math.pow(maxlength - prev, 2);
		}
	}

	let minCost = new Array(nWords);
	let result = new Array(nWords);
	// start from right to left
	for(let i = nWords - 1; i >= 0; i--){
		minCost[i] = justifyArray[i][nWords-1];
		result[i] = nWords
		for(let j = nWords - 1; j > i; j--){
			if(!justifyArray[i][j-1]) continue;
			if(!minCost[i] || minCost[i] > minCost[j] + justifyArray[i][j - 1]){
				minCost[i] = minCost[j] + justifyArray[i][j-1]
				result[i] = j
			}
		}
	}

	let ans = new Map()
	for(let i = 0; i < nWords; i++){
		if(ans.has(result[i])){
			ans.set(result[i], ans.get(result[i]) + ' ' + arr[i])
		} else{
			ans.set(result[i], arr[i])
		}
	}
	let seen = new Set()
	for(let i = 0; i < nWords; i++){
		if(seen.has(result[i])) continue;
		else{
			console.log(ans.get(result[i]))
			seen.add(result[i])
		}
	}
}

justify(['i','love','you', 'sike', 'youre','berkenshire?', 'a', 'fat', 'old', ',man'], 32)