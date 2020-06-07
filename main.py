import tkinter as tk
from tkinter import ttk

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import gexpmath as gm
import random
from datetime import datetime
from decimal import *
D = Decimal

import numpy as np
import pandas as pd
import scipy.stats
import clipboard
from scipy import interpolate
from scipy.interpolate import interp1d
from collections import Counter
from tkinter.filedialog import askopenfilename

class App(tk.Tk):
    """docstring for App"""

    def __init__(self):
        tk.Tk.__init__(self)
        self.title("GNB")
        self.geometry("1400x900+0+0")
        container = tk.Frame(self)

        container.pack(side="top", fill="both", expand=True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        for F in [VisualFrame, ApproximationFrame]:
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        menubar = tk.Menu(container)

        pagemenu = tk.Menu(menubar, tearoff=0)
        pagemenu.add_command(label="Sampling & visualisation", command=lambda: self.show_frame(VisualFrame))
        pagemenu.add_command(label="Parameters approximation", command=lambda: self.show_frame(ApproximationFrame))
        menubar.add_cascade(label="Go to...", menu=pagemenu)

        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Desciption", command=lambda: self.show_description())
        helpmenu.add_command(label="Input files", command=lambda: self.show_inputs())
        menubar.add_cascade(label="Help", menu=helpmenu)

        tk.Tk.config(self, menu=menubar)

        self.show_frame(VisualFrame)

    def show_frame(self, frm):
        frame = self.frames[frm]
        frame.tkraise()

    def show_description(self):
        info_string = '''This app is for analysis of samples of discrete variables \u2265 0. 
        It's based on represantation of Generalized Negative Binomial distribution (GNB) through gamma-exponential function by Kudryavtsev A.A.
        More info on topic: http://www.mathnet.ru/links/4020ec82ae3281e43e80576e2c756072/ia632.pdf

        App by Egor Tsvetkov
        egor.tsvetkov98@gmail.com

        '''
        
        desc_window = tk.Toplevel(self)
        label = tk.Label(desc_window, text=info_string)
        label.pack()

    def show_inputs(self):
        info_string = '''For input files use .csv with one column and no header. It should like this:
        0
        1
        2
        3
        4
        5

        '''

        desc_window = tk.Toplevel(self)
        label = tk.Label(desc_window, text=info_string)
        label.pack()

class VisualFrame(tk.Frame):
    """
    Main frame with visualisation
    """

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        histogram_frm = tk.Frame(self, height=800, width=800)
        histogram_frm.place(x=240, y=-50)

        self.fig = Figure(figsize=(12, 7), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=histogram_frm)
        self.canvas.get_tk_widget().pack(side=tk.TOP)
        self.p = self.fig.add_subplot(111)

        header = tk.Label(self, text='GNB Generation')
        header.pack()

        params = tk.Frame(self, height=800, width=200, pady=20, padx=20)
        params.place(x=40, y=30)

        self.alpha_in = tk.StringVar()
        self.N_in = tk.StringVar()
        self.u_in = tk.StringVar()
        self.p_in = tk.StringVar()
        self.bins_in = tk.StringVar()

        alpha_label = tk.Label(params, text="\u03B8", width=3)
        alpha_input = tk.Entry(params, width=5, textvariable=self.alpha_in)

        u_label = tk.Label(params, text="v", width=3)
        u_input = tk.Entry(params, width=5, textvariable=self.u_in)

        p_label = tk.Label(params, text="q", width=3)
        p_input = tk.Entry(params, width=5, textvariable=self.p_in)

        N_label = tk.Label(params, text="N", width=3)
        N_input = tk.Entry(params, width=5, textvariable=self.N_in)

        bins_label = tk.Label(params, text="bins", width=3)
        bins_input = tk.Entry(params, width=5, textvariable=self.bins_in)

        params_widgets = [
            (u_label, u_input),
            (p_label, p_input),
            (alpha_label, alpha_input),
            (N_label, N_input),
            (bins_label, bins_input)
        ]

        row = 0
        for label, input_ in params_widgets:
            label.place(x=20, y=30 * row + 2)
            input_.place(x=50, y=30 * row)
            row += 1
        row += 1

        generate_GBN_btn      = tk.Button(params, text=" Generate sample ", command=self.generate_sample)
        draw_histogram_btn    = tk.Button(params, text="  Draw histogram  ", command=self.draw_histogram)
        draw_spline_btn       = tk.Button(params, text='     Draw spline     ', command=self.draw_spline)
        clear_figure_btn      = tk.Button(params, text='          Clear          ', command=self.clear_figure)
        copy_to_clipboard_btn = tk.Button(params, text='Copy to clipboard', command=self.copy_experiments_to_clipboard)
        export_to_csv         = tk.Button(params, text='   Export to csv    ', command=self.copy_experiments_to_clipboard)
        export_to_excel       = tk.Button(params, text='  Export to excel  ', command=self.copy_experiments_to_clipboard)
        draw_histogram_btn.config(state='disabled')
        draw_spline_btn.config(state='disabled')
        copy_to_clipboard_btn.config(state='disabled')
        clear_figure_btn.config(state='disabled')
        export_to_csv.config(state='disabled')
        export_to_excel.config(state='disabled')

        self.buttons = {
            "gen": generate_GBN_btn,
            "hist": draw_histogram_btn,
            "spline": draw_spline_btn,
            "clear": clear_figure_btn,
            "copy": copy_to_clipboard_btn,
            "csv": export_to_csv,
            "excel": export_to_excel
        }

        for btn in self.buttons.values():
            btn.place(x=20, y=30 * row)
            row += 1

        self.info_label = tk.Label(params)
        self.info_label.place(x=20, y=30 * row)
        self.y_max = 0
        self.x_max = 0

    def generate_sample(self):
        alpha = D(self.alpha_in.get())
        p = D(self.p_in.get())
        u = D(self.u_in.get())
        bins = int(self.bins_in.get())
        size = int(self.N_in.get())
        if alpha <= 0 or p <= 0 or u <= 0 or size <= 0 or bins <= 0:
            return

        self.dist = gm.GNBDistribution(u, p, alpha, bins)
        self.dist.generate_sample(size)
        self.enable_button("copy")
        self.enable_button("hist")
        self.enable_button("csv")
        self.enable_button("excel")

        info_str = ''
        for i, pr in enumerate(self.dist.probabilities):
            info_str += 'P(X = ' + str(i) + ') = ' + str(pr)[0:5] + '\n'

        info_str += 'P(X > ' + str(self.dist.max_n - 1) + ') = ' + str(abs(1 - sum(self.dist.probabilities)))[0:5] + '\n'
        self.info_label.config(text=info_str)

    def draw_histogram(self, density=False):
        self.disable_button("hist")

        hist = self.dist.make_hist_dict()
        x = hist["x"]
        y = hist["y"]

        if density:
            for i in range(len(y)):
                y[i] /= self.dist.sample_size
        self.p.bar(x, y, width=0.1, label=str(self.dist))
        self.p.plot(x, y, 'o', markersize=10)

        if self.y_max < max(y):
            self.y_max = max(y)
            self.p.set_ylim([0, self.y_max * 1.1])

        if self.x_max < max(x):
            self.x_max = max(x)
            self.p.set_xlim([-1, self.x_max + 1])

        self.p.set_xlabel('Value', fontsize=10)
        self.p.set_ylabel('Frequency', fontsize=10)

        self.p.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), shadow=True, ncol=4)

        self.canvas.draw()
        self.enable_button("spline")
        self.enable_button("clear")

    def draw_spline(self):
        hist = self.dist.make_hist_dict()
        x = hist["x"]
        y = hist["y"]
        cubic_interp = interp1d(x, y, kind='cubic')
        x_partition = np.linspace(0, self.x_max, self.x_max * 10)
        self.p.plot(x_partition, cubic_interp(x_partition))
        self.canvas.draw()
        self.disable_button("spline")

    def copy_experiments_to_clipboard(self):
        copy_string = ''
        for val in self.dist.sample:
            copy_string += str(val) + '\n'
        clipboard.copy(copy_string)

    def clear_figure(self):
        self.x_max = 0
        self.y_max = 0
        self.fig.clf()
        self.p = self.fig.add_subplot(111)
        self.canvas.draw()
        self.disable_button("spline")
        self.disable_button("clear")

    def enable_button(self, btn_name):
        self.buttons[btn_name].config(state='normal')

    def disable_button(self, btn_name):
        self.buttons[btn_name].config(state='disabled')


class ApproximationFrame(tk.Frame):
    """docstring for ApproximationFrame"""

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        header = tk.Label(self, text='GNB Parameters Approximation')
        header.pack()

        self.controller = controller
        input_frame = tk.Frame(self, width=340, height=800, pady=20, padx=20)
        input_frame.place(x=40, y=30)

        self.alpha_in = tk.StringVar()
        self.u_in = tk.StringVar()
        self.p_in = tk.StringVar()
        self.a_in = tk.StringVar()
        alpha_label = tk.Label(input_frame, text="\u03B8", width=3)
        alpha_input = tk.Entry(input_frame, width=5, textvariable=self.alpha_in)

        u_label = tk.Label(input_frame, text="v", width=3)
        u_input = tk.Entry(input_frame, width=5, textvariable=self.u_in)

        p_label = tk.Label(input_frame, text="q", width=3)
        p_input = tk.Entry(input_frame, width=5, textvariable=self.p_in)

        a_label = tk.Label(input_frame, text="\u03B1", width=3)
        a_input = tk.Entry(input_frame, width=5, textvariable=self.a_in)

        params_widgets = [
            (u_label, u_input),
            (p_label, p_input),
            (alpha_label, alpha_input),
            (a_label, a_input)
        ]

        row = 0
        for label, input_ in params_widgets:
            label.place(x=20, y=30 * row + 2)
            input_.place(x=50, y=30 * row)
            row += 1

        row += 1

        self.param_method = tk.StringVar()
        r1 = tk.Radiobutton(input_frame, text="Nelder-Mead simplex", value='l2opt', variable=self.param_method)
        r2 = tk.Radiobutton(input_frame, text="Maximum likelihood", value='ml', variable=self.param_method)
        r1.place(x=20, y=30 * row)
        r2.place(x=20, y=30 * (row + 1))
        r1.select()

        read_file_button = tk.Button(input_frame, text="            Open file            ", command=self.read_file)
        # read_file_button.place(x=20, y=30 * (row + 2))

        chi2_test_button = tk.Button(input_frame, text="       Pearson \u03C7\u00B2 test      ", command=self.pearson_test)
        chi2_test_button.place(x=20, y=30 * (row + 3))

        param_estimation_button = tk.Button(input_frame, text=" Parameters estimation ", command=self.param_estimation)
        param_estimation_button.place(x=20, y=30 * (row + 4))
        self.df = None

        output_frame = tk.Frame(self, height=800, width=800)
        output_frame.place(x=240, y=-50)

        self.fig = Figure(figsize=(12, 7), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=output_frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP)
        self.p = self.fig.add_subplot(111)

        self.output_label = tk.Label(input_frame)
        self.output_label.place(x=-10, y=30 * (row + 6))

    def read_file(self):

        v = tk.StringVar()
        self.controller.update()
        csv_file_path = askopenfilename()
        # self.controllers.destroy()
        if csv_file_path == '':
            return -1
        self.df = pd.read_csv(csv_file_path, header=None, names=['val'])
        v.set(csv_file_path)
        self.controller.update()
        print(self.df.head())
        return 1

    def pearson_test(self):
        fl = self.read_file()
        if fl == -1:
            return
        self.controller.update()
        sample = list(self.df['val'])
        v = D(self.u_in.get())
        q = D(self.p_in.get())
        theta = D(self.alpha_in.get())

        dist = gm.GNBDistribution(v, q, theta, max(sample) + 1)

        alpha = self.a_in.get() or '0.05'
        alpha = float(alpha)
        hist = make_hist(sample)
        result, chi2_stat, crit_val, pvalue = dist.chi_square(alpha, hist)
        suff = ' not' if result else ''
        result_string = 'Hypothesis is' + suff + ' rejected\n' + f'\u03C7\u00B2 = {chi2_stat:0.2f}\npvalue = {pvalue:0.3f}.'
        self.output_label.config(text=result_string)

    def param_estimation(self):
        fl = self.read_file()
        if fl == -1:
            return
        sample = list(self.df['val'])
        if self.param_method.get() == 'l2opt':
            hist = make_hist(sample, density=True)['y']
            res = gm.l2optimization(hist)
        else:
            hist = make_hist(sample)['y']
            res = gm.ml_gnb(sample)
        v = res[0]
        q = res[1]
        theta = res[2]
        dist = gm.GNBDistribution(v, q, theta, max(sample) + 1)
        # test = dist.chi_square(sample)
        result_string = f'v = {res[0]:0.5f}\nq = {res[1]:0.5f}\n\u03B8 = {res[2]:0.5f}'
        self.output_label.config(text=result_string)
        self.fig.clf()
        self.p = self.fig.add_subplot(111)
        est_probs = dist.probabilities
        x_arr = list(range(len(est_probs)))
        self.p.bar(np.array(x_arr) - 0.1, est_probs, width=0.2, label='Estimated')

        h2 = make_hist(sample, density=True)
        self.p.bar(np.array(h2['x']) + 0.1, h2['y'], width=0.2, label='Real')
        self.p.legend()
        self.canvas.draw()
        print('-' * 30)
        print(est_probs)
        print(h2['y'])
        print('-' * 30)


def make_hist(sample, density=False, no_outliers=False):
    if sample is None:
        raise(Exception("No sample in object. Call object.generate_sample(sample_size) before calling this function"))

    max_n = max(sample)
    hist = np.histogram(sample, max_n, range=(0, max_n))

    x = hist[1][:-1]
    y = hist[0]
    if no_outliers and max_n in x:
        x = x[:-1]
        y = y[:-1]

    if density:
        y = np.array(y) / sum(y)

    return {'x': x, 'y': y}

if __name__ == "__main__":
    app = App()
    print("running app")
    app.mainloop()
