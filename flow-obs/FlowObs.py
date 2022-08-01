from gams import *
import sys
import os
import threading
import numpy as np
import itertools
import matplotlib.pyplot as plt
import csv

class FlowObs():
    def __init__(self,dx,dy):
        self.ws = GamsWorkspace(working_directory = os.getcwd())
        self.dx = dx
        self.dy = dy
        self.drawing = 1
        
        self.jobGAMSfile = "flowobs.gms"

    def solve_result(self):
        db = self.ws.add_database()
        job = self.ws.add_job_from_file(self.jobGAMSfile)
        opt = self.ws.add_options()
        opt.defines["gdxincname"] = db.name
        opt.nlp = "conopt"
        opt.idcgdxoutput = "test.gdx"
        dx = db.add_parameter("dx", 0)
        dx.add_record().value = self.dx
        dy = db.add_parameter("dy", 0)
        dy.add_record().value = self.dy
        job.run(opt,databases = db)
        FId = job.out_db.get_variable("obj").first_record().level

        return FId

if __name__ == '__main__' :
    test = FlowObs(1,1)
    result = test.solve_result()

    print(result)

    exit(0)