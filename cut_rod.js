

// return the best obtainable price for a 
// rod of length n and prices of different pieces
const cut_rod = (price,n) => {
	let val = new Array(n+1).fill(0);
	for(let i = 1; i < n + 1; i++){
		let max = -Infinity;
		for(let j = 0; j < i; j++){
			max = Math.max(max, price[j] + val[i - j - 1])
		}
		val[i] = max;
	}
	return val[n]
}