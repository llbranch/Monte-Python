# Topics to remember 
[resource](https://superfastpython.com/multiprocessing-in-python/#What_Are_Processes)
1. Concurrent 
- Code that can run out of order 
2. Simultaneously 
- interchangeable with Parallel
3. Parallel 
- code capable to execute simultaneously 
4. Process 
- is a computer program or an instance of the python interpreter
- a process will have at least one thread (main thread)
- will terminate once all (non-background) threads are terminated
5. Thread
- always exists within a process and represents the manner in which the instructions or code is executed 
6. Cycle of a Process 
- new process 
    - creating an instance of the .Process class 
- running process
    - start() method 
    - also creates and starts the main thread to execute the code 
- terminated process 
    - cannot terminate until 
        - all non-deamon threads have terminated including main thread 
        - all non-deamon child processes have terminated, including main thread 

7. Blocked running process 
- can happen if the main thread is reading or writing to a file 

8. Child Process 
- to create: fork, spawn, fork server 
- depending on the creation it may or may not inherit properties of the parent process 
- forked process: may inherit a copy of the global variables from the parent process