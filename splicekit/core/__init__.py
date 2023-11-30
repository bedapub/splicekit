"""
core module of splicekit
"""

import os
import splicekit.core as core
import splicekit
from queue import *
from threading import *
import os
import multiprocessing

num_worker_threads = max(1, multiprocessing.cpu_count()-1)
q = Queue()

def worker():
    while True:
        task = q.get()
        os.system(task)
        q.task_done()

def mprocess(fname):
    tasks = []
    f = open(fname, "rt")
    r = f.readline()
    while r:
        r = r.replace("\r", "").replace("\n", "")
        if r!="" and len(r)>0:
            tasks.append(r)
        r = f.readline()
    f.close()
    print(f"splicekit | started multicore ({num_worker_threads}) processing: {fname}")
    for i in range(num_worker_threads):
        t = Thread(target=worker)
        t.daemon = True
        t.start()

    for task in tasks:
        q.put(task)

    q.join()
    print("splicekit | done multicore ({num_worker_threads}) processing: {fname}")

