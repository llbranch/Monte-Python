

import multiprocessing
import multiprocessing_import_worker_function

if __name__ == '__main__': 
    jobs = []
    for i in range(5): 
        process_i = multiprocessing.Process(
            target=multiprocessing_import_worker_function.worker_function,
        )
        jobs.append(process_i)
        process_i.start()
        # print(type(process_i))

    