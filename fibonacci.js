/**
 * @param {number} n
 * @return {number}
 */
var climbStairs = function(n) {
    return fibonacci(n)
    
};


let fibonacci = (n) => {
    let arr = new Array(n);
    arr[0] = 1;
    arr[1] = 1;
    
    for(let i = 2; i <= n; i ++){
        arr[i] = arr[i-2] + arr[i -1];
    }
    return arr[n]
}