import tkinter as tk
from tkinter import ttk

import matplotlib

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import gexpmath as gm
import random
from datetime import datetime
from decimal import *

import numpy as np
import scipy.stats
import clipboard
from scipy import interpolate
from scipy.interpolate import interp1d
from collections import Counter
matplotlib.use("TkAgg")
experiments = []
bins = 0
probabilities = []
global_bin_max = 0
global_y_max = 0


def draw_histogram(spline=False):

    p = f.add_subplot(111)
    bins = int(bins_in.get())
    hist = np.histogram(experiments, bins, range=(0, bins))
    x = hist[1][:-1]
    y = hist[0]
    p.bar(x, y, width=0.1)
    p.plot(x, y, 'o', markersize=10)

    # p.plot(hist_dist)

    y_max = max(Counter(experiments).values())
    y_max = int(1.1 * y_max)
    # print(Counter(experiments))
    # print(hist)
    print(y_max)
    global global_y_max
    global global_bin_max

    if y_max > global_y_max:
        global_y_max = y_max

    if bins > global_bin_max:
        global_bin_max = bins

    p.set_xlim([-1, global_bin_max])
    p.set_ylim([0, global_y_max])
    p.set_xlabel('Value', fontsize=10)
    p.set_ylabel('Frequency', fontsize=10)

    canvas.draw()
    draw_histogram_btn.config(state='disabled')
    draw_spline_btn.config(state='normal')
    clear_figure_btn.config(state='normal')
    global subp
    subp = p


def draw_spline():
    p = f.add_subplot(111)
    bins = int(bins_in.get())
    hist = np.histogram(experiments, bins, range=(0, bins))
    x = hist[1][:-1]
    y = hist[0]

    max_bin = int(x[-1])
    cubic_interp = interp1d(x, y, kind='cubic')
    x_partition = np.linspace(0, max_bin, max_bin * 10)
    p.plot(x_partition, cubic_interp(x_partition))
    canvas.draw()
    draw_spline_btn.config(state='disabled')


def copy_experiments_to_clipboard():
    copy_string = ''
    for val in experiments:
        copy_string += str(val) + '\n'
    clipboard.copy(copy_string)


def generate_sample():
    global experiments
    experiments = []
    alpha = Decimal(alpha_in.get())
    u = Decimal(u_in.get())
    p = Decimal(p_in.get())

    N = int(N_in.get())
    bins = int(bins_in.get())

    global probabilities
    probabilities = []

    for i in range(bins):
        probabilities.append(gm.GNB_Probability(i, u, p, alpha))

    # for prob in probabilities:
    #   print(prob)

    # print(sum(probabilities))
    # print()
    for k, v in enumerate(probabilities):
        print(k, '-', v)

    random.seed(datetime.now())
    clean_expirements = []
    for i in range(N):
        pr = random.random()
        j = 0
        cumprob = 0
        # out_of_range = False
        while (j < bins):
            cumprob += probabilities[j]
            if (cumprob > pr):
                break
            j += 1
        if cumprob > pr:
            clean_expirements.append(j)
            experiments.append(j)

    # print(experiments)
    copy_to_clipboard_btn.config(state='normal')
    draw_histogram_btn.config(state='normal')

    filename = 'experiments_log_' + str(datetime.now()) + '.txt'
    log_file = open(filename, 'w')

    log_file.write('Params: u = ' + str(u) + '; p = ' + str(p) + '; alpha = ' + str(alpha) + '; N = ' + str(N) + '; bins = ' + str(bins) + '\n')
    for val in experiments:
        log_file.write(str(val) + '\n')

    info_str = ''
    for i, pr in enumerate(probabilities):
        info_str += 'P(k = ' + str(i) + ') = ' + str(pr)[0:5] + '\n'

    print(sum(probabilities))
    info_str += 'P(k > ' + str(bins - 1) + ') = ' + \
        str(1 - sum(probabilities))[0:5] + '\n'
    info_str += 'Mean = ' + str(np.mean(clean_expirements))[0:5] + '\n'
    info_str += 'Std = ' + str(np.std(clean_expirements))[0:5] + '\n'
    desc_text.delete(1.0, tk.END)
    desc_text.insert(1.0, info_str)


def clear_figure():
    global global_y_max
    global global_bin_max
    global_bin_max = 0
    global_y_max = 0
    f.clf()
    canvas.draw()
    clear_figure_btn.config(state='disabled')


window = tk.Tk()
window.title("GNB")
window.geometry("1400x900+0+0")

mainmenu = tk.Menu(window)
window.config(menu=mainmenu)

mainmenu.add_command(label='Generation')
mainmenu.add_command(label='Test and parameteres approximation')
# mainmenu.pack()

header = tk.Label(window, text='GNB Generation')
header.pack()


params = tk.Frame(window, height=400, width=200, pady=20, padx=20)
params.place(anchor=tk.NW)

alpha_in = tk.StringVar()
N_in = tk.StringVar()
u_in = tk.StringVar()
p_in = tk.StringVar()
bins_in = tk.StringVar()

alpha_label = tk.Label(params, text="\u03B1", width=3)
alpha_input = tk.Entry(params, width=5, textvariable=alpha_in)
u_label = tk.Label(params, text="u", width=3)
u_input = tk.Entry(params, width=5, textvariable=u_in)
p_label = tk.Label(params, text="p", width=3)
p_input = tk.Entry(params, width=5, textvariable=p_in)
N_label = tk.Label(params, text="N", width=3)
N_input = tk.Entry(params, width=5, textvariable=N_in)
bins_label = tk.Label(params, text="bins", width=3)
bins_input = tk.Entry(params, width=5, textvariable=bins_in)


# pois(E\lambda)
# GG(поверх ООБР с теми же параметрами)


params_widgets = [
    (alpha_label, alpha_input),
    (u_label, u_input),
    (p_label, p_input),
    (N_label, N_input),
    (bins_label, bins_input)
]

i = 0
for label, input_ in params_widgets:
    label.place(x=20, y=30 * i + 2)
    input_.place(x=50, y=30 * i)
    i += 1
i += 1

generate_GBN_btn = tk.Button(
    params, text=" Generate sample ", command=generate_sample)
generate_GBN_btn.place(x=20, y=30 * i)

draw_histogram_btn = tk.Button(
    params, text="  Draw histogram  ", command=draw_histogram)
draw_histogram_btn.config(state='disabled')
draw_histogram_btn.place(x=20, y=30 * (i + 1))

draw_spline_btn = tk.Button(
    params, text='     Draw spline     ', command=draw_spline)
draw_spline_btn.config(state='disabled')
draw_spline_btn.place(x=20, y=30 * (i + 2))

copy_to_clipboard_btn = tk.Button(
    params, text='Copy to clipboard', command=copy_experiments_to_clipboard)
copy_to_clipboard_btn.config(state='disabled')
copy_to_clipboard_btn.place(x=20, y=30 * (i + 3))

clear_figure_btn = tk.Button(
    params, text='          Clear          ', command=clear_figure)
clear_figure_btn.config(state='disabled')
clear_figure_btn.place(x=20, y=30 * (i + 4))


histogram_frm = tk.Frame(window, height=800, width=800)
histogram_frm.place(x=420, y=30)

f = Figure(figsize=(8, 8), dpi=100)
canvas = FigureCanvasTkAgg(f, master=histogram_frm)
canvas.get_tk_widget().grid(row=0, column=0)

dist_description = tk.Frame(window, height=400, width=200, pady=20, padx=20)
dist_description.place(x=0, y=350)

desc_text = tk.Text(dist_description, width=200)
desc_text.place(x=0, y=0)

window.mainloop()
