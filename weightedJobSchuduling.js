class Job{
    constructor(s,e,p){
        this.start = s;
        this.end = e;
        this.val = p;
    }
}
// Find the latest job (in sorted array) that doesnt conflix with the job[i]
let latestNonConflict = (arr, i) => {
    for(let j = i - 1; j >= 0; j --){
        // does not overlap
        if(arr[j].end <= arr[i].start) 
            return j;
    }
    return -1;
}

var jobScheduling = function(startTime, endTime, profit) {
    let n = endTime.length;
    let job = new Array(n);
    // combine every input into a single Job class 
    for(let i = 0; i < n; i++)
        job[i] = new Job(startTime[i], endTime[i], profit[i]);
    job.sort((a,b) => a.end - b.end);
    
    let memo = new Array(n);
    memo[0] = job[0].val;
    
    for(let i = 1; i < n; i++){        
        let currPrice = job[i].val;
        let l = latestNonConflict(job, i);
        if(l != -1)
            currPrice += memo[l];
        // store maximum of including and excluding
        memo[i] = Math.max(currPrice, memo[i-1]);
    }
    return memo[n-1];
};
