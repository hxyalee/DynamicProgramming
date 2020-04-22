
# items in format
# [ { price : int,
#	 weight : int } ]
# limit is the maximum weight possible
def knapsack_main(items, limit):
	noItem = len(items)
	# + 1 to include zeros
	dp = [ [0 for x in range(limit + 1)] for y in range(noItem + 1) ]
	items.insert(0, {'price' : 0, 'weight' : 0})
	for i in range(1, noItem + 1):
		for j in range(1, limit + 1):
			curWeight = items[i]['weight']
			curPrice = items[i]['price']
			if(curWeight > j):
				dp[i][j] = dp[i-1][j]
			else:
				dp[i][j] = max(dp[i-1][j], dp[i-1][j-curWeight] + curPrice)

	for x in dp:
		print(x)

	print()
	print()
	return getPath(items, dp)

def getPath(items, dp):
	col = len(dp) - 1
	row = len(dp[0]) - 1
	path = list()
	while(col > 0):
		while(row > 0):
			if(dp[col - 1][row] < dp[col][row] and dp[col][row - 1] < dp[col][row]):
				path.append(items[col])
				row = row - items[col]['price']
				break
			elif(dp[col - 1][row] >= dp[col][row]):
				col -= 1
			elif(dp[col][row - 1] >= dp[col][row]):
				row -= 1
		print(path)
		col -= 1


if __name__ == '__main__':
	items = [
		{'price':1, 'weight':2},
		{'price':2, 'weight':3},
		{'price':5, 'weight':4},
		{'price':6, 'weight':5}
			]
	limit = 8
	print(knapsack_main(items, limit))

