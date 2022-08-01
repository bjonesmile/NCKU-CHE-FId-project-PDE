import time

class MyTimer():
    timer_num = 0
    def __init__(self, name = None):
        if name is None:
            self.name = "timer"+str(MyTimer.timer_num)
            MyTimer.timer_num += 1
        else:
            self.name = name
            MyTimer.timer_num += 1
        self.basictime = time.process_time()
        self.basictime_real = time.time()
        self.loop = {"process": [], "real": []}

    def getTime(self, kind="process"):
        if kind == "process":
            ptime = time.process_time() - self.basictime
            print("process time: ",ptime)
            self.loop["process"].append(ptime)
            return ptime
        else:
            ptime_real = time.time() - self.basictime_real
            print("process time(real): ",ptime_real)
            self.loop["real"].append(ptime_real)
            return ptime_real

    def restart(self):
        print(f"{self.name} reset!")
        self.basictime = time.process_time()
        self.basictime_real = time.time()
        pass

if __name__ == "__main__":
    timer1 = MyTimer()
    timer2 = MyTimer()
    time.sleep(3)
    timer1.getTime(kind="process")
    timer2.getTime(kind="real")
    for i in range(10000):
        print("i:",i)
    timer1.getTime(kind="process")
    timer2.getTime(kind="real")
    timer1.restart()
    time.sleep(2)
    timer1.getTime(kind="real")
    timer2.getTime(kind="real")
    exit()
