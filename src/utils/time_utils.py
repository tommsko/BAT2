import time

def current_milli_time() -> int:
    return round(time.time() * 1000)

def current_sec_time() -> float:
    return time.time()