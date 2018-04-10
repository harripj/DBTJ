#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 14:47:27 2017

@author: pjh523
"""

from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg,
                                                NavigationToolbar2QT)
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QSizePolicy, QMainWindow, QWidget,
                             QGridLayout, QDockWidget)

import matplotlib.pyplot as plt
from ..helper_functions import remove_line_from_axes

class MPLCanvas(FigureCanvasQTAgg):
    '''generic MPL figure canvas for use in GUI'''

    def __init__(self, fig_size=(10, 10)):
        # initialise figure
        self.fig, self.ax = plt.subplots()
        #self.fig.set_size_inches(fig_size)
        self.fig.tight_layout()

        # mouse in axes flag
        self.mia = False
        # button click location
        self.xclick_loc = 0
        self.yclick_loc = 0
        
        # holder for image
        self.image = None
        # holder for line
        self.line = None

        # initialise canvas
        FigureCanvasQTAgg.__init__(self, self.fig)

        # set canvas geometry rules
        FigureCanvasQTAgg.setSizePolicy(self,
                                        QSizePolicy.Expanding,
                                        QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

        # set axes GUI events
        FigureCanvasQTAgg.mpl_connect(self, 'axes_enter_event', self.ax_enter)
        FigureCanvasQTAgg.mpl_connect(self, 'axes_leave_event', self.ax_leave)
        FigureCanvasQTAgg.mpl_connect(self,
                                      'button_press_event', self.button_press)

    def ax_enter(self, event):
        '''flag when enters edge axes'''
        self.mia = True

    def ax_leave(self, event):
        '''flag when leaving edge axes'''
        self.mia = False

    def button_press(self, event):
        '''activated during click in axis'''
        if self.mia:
            self.xclick_loc = event.xdata
            self.yclick_loc = event.ydata

    def plot(self, x, y, cla=False, remove_line=None,
             xlim=False, ylim=False, **kwargs):
        '''plot x vs y on axes and update canvas
        kwargs are normal plt.plot kwargs, eg. label'''
        if cla:
            self.ax.cla()
        if ylim:
            self.ax.set_ylim(min(y), max(y))  # may not be numpy array
        if xlim:
            self.ax.set_xlim(min(x), max(x))
        if remove_line is not None:
            remove_line_from_axes(self.ax, remove_line)
        self.ax.plot(x, y, **kwargs)
        self.draw()
        
    def plot_image(self, image, cla=False, aspect=None):
        ''''''
        if cla:
            self.ax.cla()
        self.image = self.ax.matshow(image, aspect=aspect)
        self.draw()
        
    def update_image(self, image):
        ''''''
        self.image.set_array(image)
        #print(self.image.get_array())
        self.ax.draw_artist(self.image)
        self.update()

class QtWindowCanvas(QMainWindow):
    '''generic app window to display a figure with a toolbar'''

    def __init__(self, bottom_dock=True, left_dock=True):
        # initialsie window
        QMainWindow.__init__(self)
        self.main_widget = QWidget(self)

        # add dock widget
        if bottom_dock:
            self.bottom_dock = QDockWidget()
            self.addDockWidget(Qt.BottomDockWidgetArea, self.bottom_dock)
            self.bottom_dock.setFeatures(QDockWidget.DockWidgetFloatable)

        if left_dock:
            self.left_dock = QDockWidget()
            self.addDockWidget(Qt.LeftDockWidgetArea, self.left_dock)
            self.left_dock.setFeatures(QDockWidget.DockWidgetFloatable)

        # define layout
        self.grid = QGridLayout(self.main_widget)

        # create figure on window
        self.canvas = MPLCanvas()
        self.grid.addWidget(self.canvas, 0, 0, 1, 1)  # all rows

        # figure toolbar
        toolbar = NavigationToolbar2QT(self.canvas, self)
        self.grid.addWidget(toolbar, 1, 0, 1, 1)

        self.setCentralWidget(self.main_widget)
