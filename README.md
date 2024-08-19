# bidirectional_QCSP

File "bidirectional_QCSP.cpp": the C++ code of bidirectional_QCSP model, implemented through CPLEX.

File "bidirectional_QCSP_Benders.cpp": the C++ code of Benders approach, including the lower bound formulation and the unidirectional formulation, implemented through CPLEX.

File "test.zip": the compressed file of all test data, most of which are obtained from the benchmark instances in Kim and Park (2004) and Meisel and Bierwirth (2011), while the file "test - real instances" includes the data from real practice. All of the above instances are named by "n-b-q", where n is the number of tasks, b is the number of bays, and q is the mumber of QCs. In the txt file corresponding to each instance, the first "[]" records the number of tasks, the number of bays, the number of precedence pairs, the number of QCs, the travel time, safety margins. The next two "[]"s record the processing times and bay locations of each task, respectively. The following two "[]"s record the ready times, initial bay locations of each QC, respectively. The rest "[]"s record the detailed task precedence pairs. 


File "QCSP-Bid-Tables": summarized computational results, recording all the details.
