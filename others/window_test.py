import tkinter as tk
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import dill as pickle

window = tk.Tk()

window.title('graph loader')
window.geometry('600x400')
window.configure(background='white')

header_label = tk.Label(window, text='3D interactive graph loader')
header_label.pack()

input_frame = tk.Frame(window)
input_frame.pack(side=tk.TOP)
input_label = tk.Label(input_frame, text='file name:')
input_label.pack(side=tk.LEFT)
input_entry = tk.Entry(input_frame)
input_entry.pack(side=tk.LEFT)

def click_show_btn_action():
    filename = input_entry.get()
    load_fig(filename)

def load_fig(picklefile):
    plt.figure().add_subplot(111, projection='3d')
    plotdata = pickle.load(open(picklefile, 'rb'))
    plotdata['plotfun'](plotdata)
    plt.show()

show_btn = tk.Button(window, text='load', command=click_show_btn_action)
show_btn.pack()

window.mainloop()