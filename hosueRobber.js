/**
 * @param {number[]} nums
 * @return {number}
 */
var rob = function(nums) {
    if (nums.length === 0) return 0;
    let length = nums.length;
    let memo = new Array(length + 1);
    // minimum value
    memo[0] = 0
    memo[1] = nums[0]
    for(let i = 1; i < length; i++){
        memo[i+1] = Math.max(memo[i], memo[i - 1] + nums[i]);
    }
    return memo[length]
};