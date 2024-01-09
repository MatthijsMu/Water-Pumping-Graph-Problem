# Water Pumping Graph Problem

## Introduction
 My solution for the second Project of the course ADS. The problem is a Dynamic graph problem, where one has to find an optimal route through a graph, turning on water pumps to maximize the total water pumped in a finite time window. I proposed a solution that employed:
 
 - A graph reduction phase: use Dijkstra's Algorithm to find shortest paths between all waterpumps, in order to reduce the state space for the dynamic programming phase. The graph reduction phase needs to calculate the shortest paths between all W waterpumps, and therefore it is inefficient to employ Floyd-Warshall's Algorithm, which would run in O(V^3) where V ~ 5000, whereas W <= 12 so running Dijkstra W times is O(W(W+E)log(V)).
 
 - Then use dynamic programming (a detailed explanation of the problem statement and how dynamic programming applies, can be found in the report).
 
 - During the dynamic programming, employ a branch-and-bound mechanism to prune partial solutions early if they are guaranteed to only lead to suboptimal solutions. This does not affect the worst-case complexity of the algorithm, but in practice it lead to lower wallclock times on most test cases, in particular for the large problems (V ~ 5000, t ~ 10000)
 
## What to find here

- The problem statement, providing a description of the problem and of the way input is given (also helpful if you want to know how the test cases are passed to the program)
 
- The solution is implemented in C++, see `WaterProblem.cpp`. This code is generously commented and this should make it completely self-explanatory. I encourage you to have a look at the code first.

- Then, a report that explains the algorithm and shows its correctness and complexity.

- A folder containing 27 test cases, which you can try on the compiled code. Non-sample cases are numbered 01 to 23. Each test case consists of two files, namely `<nr>-<size>.in` containing input in plain text, and `<nr>-<size>.out` containing the required output (maximum water pumped) in a single line of plain text. The test cases are formatted according to the problem statement. 

- A makefile to compile the code, because typing `g++ -o WaterProblem.cpp` is obviously too tedious.





## Interaction

You can compile the program using the accompanying makefile:
```bash
../Water-Pumping-Graph-Problem$ make
g++ -o WaterProblem WaterProblem.cpp
../Water-Pumping-Graph-Problem$
```

Then run it:

```
../Water-Pumping-Graph-Problem$ ./WaterProblem

```

The program waits for input, which you can provide by copying the entire contents of a `.in` file (hint: `Ctrl`+`a`) into the terminal (`Ctrl`+`Shift`+`v`) and pressing return once the buffer has been copied completely into stdin.
It should give the answer within 0.1 second for most test cases. You can check this against the `.ans` file.


