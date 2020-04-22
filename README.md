# DYNAMIC PROGRAMMING
* Build an optimal solution to the problem from optimal soultions for subproblems
	
* Subproblems are chosen in a way which allows recursive construction of optimal solutions to problems from optimal solutions to smaller size problems

* Efficiency of DP comes from the fact that the set of subproblems needed to solve larger problems heavily overlap: each subproblem is solved only once and its solution is stored in a table for multiple use for solving many larger problems

## Activity Selection
* **Instance:** A list of activities *a<sub>i</sub>*, *1  <= i <= n* with starting times *s<sub>i</sub>* and finishing times *f<sub>i</sub>*. No two activities can take place simultaneously.
* **Task:** Find a subset of compatible activities of **maximal total duration**
> *Note: Greedy method was used to solve a similar problem in the past: finding a subset with the **largest possible number** of compatible activities, but Greedy method **does not** work with this problem.*
* Solution:
	* Start by sorting the activities by their finishing time into a non-decreasing sequencel; and hence will assume that *f<sub>1</sub> <= f<sub>2</sub> <=  ... <= f<sub>n</sub>*
	* For every *i<=n* solve the following subpproblems: <br />
	**Subproblem** P(*i*): find a subsequence *x<sub>i</sub>* of the sequence of activities *S<sub>i</sub> = < a<sub>1</sub>, a<sub>2</sub>, ... a<sub>i</sub> >* such that:
		1. *x<sub>i</sub>* consists of non-overlapping activities;
		2. *x<sub>i</sub>* ends with activity a*<sub>i</sub>*;
		3. *x<sub>i</sub>* is of maximal total duration among all subsequences of S<sub>i</sub> which satisfy a and b.
	* Let T(*i*) be the total duration of the optimal solution S(*i*) of the subproblem P(*i*).
	* For S(1), choose a<sub>1</sub>; thus, T(1) = f<sub>1</sub> - s<sub>1</sub>.
	* **Recursion** assuming that we have solved subproblems for all *j < i* and stored them in a table, let <br />
###### <div align="center"> T(*i*) = max{T(j) + f<sub>i</sub> - s<sub>i</sub> : *i < j & f<sub>j</sub> < s<sub>i</sub>*} </div>

## Longest incresing subsequence
* Given a sequence of *n* real numbers A[1...*n*], determine a subsequence (not essentitally contiguous) of maximum length in which the values in the subsequence are strictly increasing
* Solution:
	* For each 1 <= *i* <= n **Subproblem** P(*i*):<br /> *Find a subsequence of the sequence A[1...i] of maximum length in which the values are strictly increasing and which ends with A[i].*
	* Recursion: Assume all subproblems for *j < i* has been solved. <br /> Then, look for all A[*m*] such that *m < i* and such that *A[m] < A[i]*; <br /> Among those, pick *m* which produced the longest increasing subsequecne ending with *A[m]* and extend it with *A[i]* to obtain the longest increasing subsequence which ends with *A[i]*.
	> *Why does this prodice optimal solutions to subproblems?<br /> Truncating optimal solution for P(i) push produce optimal solution of the subproblem P(m), otherwise we could find a better solution for P(i) as well, by extending such better solution of P(m)*
	* Store in the i<sup>th</sup> slot of the table the length of the longest increasing subsequnce ending with A[*i*] and *m* such that the optimal solution for P(*i*) extends the optimal solution for P(*m*). <br /> So, for every *i <= n* the longest increasing subsequence of the sequence A[1...*i*] is found which ends with A[*i*]. <br /> Then, pick the longest one from such subsequence. <br /> **Time complexity: *O(n<sup>2</sup>)***
###### Code
```javascript 
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
```
## Making Change
* **Instance:** You are given *n* types of coin denominations of values *v<sub>1</sub> < v<sub>2</sub> < ... v<sub>n</sub>*(all integers). Assume *v<sub>1</sub> = 1*, so that you can always make change for any integer amount. Assume that you ahve an unlimited supply of coins of each denomination.
* **Task:** Give an algorithm which makes change for any given integer amount *C* with as few coins as possible.
* Solution:
	* Consider an optimal solution *S<sub>i</sub>* for amount *i<=C*
	* if *i>0*, then *S<sub>i</sub>* includes at least one coint, say of denomination v<sub>k</sub>.
	* Removing this coin must produce an optimal solution for the amount *i-v<sub>k</sub>,  S<sub>i-v<sub>k</sub></sub>* again by the *cut-and-paist* argument.
	* We do not know which coins *S<sub>i</sub>* includes,so we try all the available coins and then pick *k* for which *S<sub>i-v<sub>k</sub></sub>* uses the fewest number of coins.
	* For 0 <= *i* <= *C*, subproblem P(*i*) is to make change for amount *i* with as few coins as possible. Let *m(i)* be the number of coins required.
	* If *i* = 0, the solution is trivial: you dont need any coin *m(0)=0*.
	* Assume optimal solution for amounts *j < i* and find an optimal solution for the amount *i*. That is, <br />
	###### <div align="center"> m(*i*) = min{*m(i-v<sub>k</sub>) + 1* | *1<= k <= n, i - v<sub>k</sub> >= 0*}</div>
	> Dont forget the condition *i-v<sub>k</sub> >= 0* or else define m(*i*) = *infinity* for i < 0.
	* **Time complexity: *O(nC)*** *(pseudo-ploynomial)*
###### Code
``` javascript
function makeChange(coins, total){
	let memo = new Array(coins.length);
	for(let i = 0; i < memo.length; i++) memo[i] = new Array(total + 1).fill(0);
	for(let i = 0; i < memo.length; i++){
		for(let j = 0; j < memo[i].length; j++){
			if(i == 0) memo[i][j] = Math.floor(j/coins[i]);
			else{
				if(coins[i] > j) memo[i][j] = memo[i-1][j]
				else memo[i][j] = Math.min(memo[i][j-coins[i]] + 1, memo[i-1][j])
			}
		}
	}
	return memo[memo.length-1][memo[0].length-1]
}
```

## Balanced Partition
* Given a set *S* of *n* integers, partition *S* into two subsets *S<sub>1</sub>S<sub>2</sub>* such that |*s<sub>1</sub> - s<sub>2</sub>*|, where *s<sub>1</sub>* and *s<sub>2</sub>* denote the sum of the elements in each of the two subsets, is minimised. 
* Solution:
	* Let *s* be the total sum of all integers in the set; consider the Knapsack problem (without duplicates) with the knapsack of size *s*/2 and with each integer *x<sub>i</sub>* of both size and value equal to *x<sub>i</sub>*.
	* Claim:
		* The best packing of such knapsack produces optimally balaned partition, with *S<sub>1</sub>* being all the integers in the knapsack and *S<sub>2</sub>* all the integers left out of the knapsack.
	* Why? Since *s = s<sub>1</sub> + s<sub>2</sub>*, we obtain 2(*s*/2 - *s<sub>1</sub>*) = *s - 2*s<sub>1</sub>* = *s<sub>2</sub> - s<sub>1</sub>*.
	* Thus, minimising *s*/2-*s<sub>1</sub>* will minimise *s<sub>2</sub> - s<sub>1</sub>*.
	* So, all we have to do is find a subset of these numbers with the largest possible total sum which fits inside a knapsack of size *s*/2.
	