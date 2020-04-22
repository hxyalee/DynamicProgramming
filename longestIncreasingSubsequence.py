class Solution:
    def lengthOfLIS(self, nums: List[int]) -> int:
        if not nums:
            return 0
        memo = [1] * len(nums)
        for i in range(1, len(nums)):
            for j in range(i):
                if memo[i] > memo[j]:
                    continue
                if nums[i] > nums[j]:
                    memo[i] = memo[j] + 1
        return max(memo)