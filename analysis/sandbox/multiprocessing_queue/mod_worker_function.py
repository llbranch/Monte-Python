# Module used for example file: multiprocessing_calling_module.py
import multiprocessing

def worker_function():
    """worker function"""
    name = multiprocessing.current_process().name
    print(name, ": Worker function")
    return # in a module we have the return 