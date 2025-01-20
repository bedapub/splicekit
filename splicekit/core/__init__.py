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
import re

# split by semicolons not inside quotes
def split_ignore_quoted(input_str):
    pattern = r'(?:[^";]|"(?:\\.|[^"])*")+'
    parts = re.findall(pattern, input_str)
    return [part.strip() for part in parts]

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

def smart_number_format(number):
    if abs(number) < 0.0001 or abs(number) >= 1e4:
        return f"{number:.4e}"
    else:
        return f"{number:.4f}"