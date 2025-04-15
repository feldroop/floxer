This project was created for my master's thesis. It is experimental and I will not develop it further. If anyone in the future wishes to keep working on this project, here are a few notes on the state of the implementation:

The next steps to improve the performance of the program would be to integrate a faster and more memory efficient alignment algorithm into the prototype. The current one from the SeqAn3 library uses too much memory,
which also hurts the running time performance.

Initially, the program was parallelized using a single OpenMP for loop in the main function. I exchanged this apporach for a complex task-based solution to achieve more finely grained load balancing.
This ended up not significantly helping the performance and added a lot of complexity. I would return to the simple OpenMP solution as a first step when restarting the work on this project.

The `src/main` folder contains the source code for a number of small programs, apart form the main of the whole floxer program. Some of these programs could be useful useful in future for analyzing the reults of `floxer`,
others could probably be deleted.

The code I used to evaluate this prototype can be found [here](https://github.com/feldroop/msc-thesis-benchmark).
