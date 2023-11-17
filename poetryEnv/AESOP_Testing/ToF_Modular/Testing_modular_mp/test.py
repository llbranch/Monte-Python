import mod
import multiprocessing
"""
Imports are only needed on the side of the function not on the main file 
"""

def other_run(): 
    print('worker')


if __name__ == '__main__':
    print(mod.pi)

    def run():
        jobs = []
        for i in range(5): 
            process_i = multiprocessing.Process(
                target=mod.worker
                )
            jobs.append(process_i)
            process_i.start()
            process_i.join()
        print('ran')
        return
    
    run()

    