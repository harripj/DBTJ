#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 14:45:49 2017

@author: pjh523
"""

import sys

from PyQt5.QtCore import Qt, QSize
from PyQt5.QtWidgets import (QApplication, QSlider, QGridLayout,
                             QWidget, QLabel, QTextEdit, QPushButton,
                             QFileDialog, QSplitter, QMessageBox)

from .Generic_QtWindow import QtWindowCanvas
from ..models import dbtj
from ..STS.spectra import Bias_Spectroscopy
from ..helper_functions import sort_ylim
from numpy import linspace


class QFloatSlider(QSlider):
    def __init__(self, float_range, offset):
        '''Basically a QSlider but with a multiplication factor, mfactor,
        which can be used to return float values'''
        QSlider.__init__(self, Qt.Horizontal)
        self.float_range = float_range
        self.offset = offset
        
    def getFloat(self):
        '''returns calibrated value'''
        int_range = self.maximum() - self.minimum()
        return self.offset + self.float_range * self.value()/int_range
        

class DBTJSlider:
    '''generic GUI with figure and sliders to show how fns change'''
    def __init__(self):
        self.app = QApplication(sys.argv)
        self.aw = QtWindowCanvas(bottom_dock=True, left_dock=False)
        self.aw.setWindowTitle('DBTJ Slider, by PJH')
        
        self.fn_bias = 3
        self.bias = linspace(-self.fn_bias, self.fn_bias, 256)

        self.setup_bottom_dock()
        #self.setup_left_dock()
        
        self.aw.canvas.ax.set_xlim(-self.fn_bias, self.fn_bias)
        self.aw.canvas.ax.set_xlabel('Bias/ V')
        self.aw.canvas.ax.set_ylabel('dIdV/ (A/V)')
        
        self.aw.canvas.fig.tight_layout()
        
        about_menu = self.aw.menuBar().addMenu('Help')
        about_menu.addAction('About').triggered.connect(self.show_about_window)

        self.aw.show()
        self.app.exec_()
        
    def show_about_window(self):
        QMessageBox.about(self.aw, 'About', 
                          'DBTJ Slider GUI created by Patrick Harrison.\n\n' +
                          'DBTJ model:\nHanna and Tinkham, Phys. Rev. B, 44 (11), 5919-5921, 1991.')        

    def setup_bottom_dock(self):
        ''''''
        splitter = QSplitter(Qt.Horizontal)
        #grid = QGridLayout(splitter)
        
        widget_left = QWidget()
        widget_right = QWidget()
        
        grid_left_dock = QGridLayout(widget_left)
        self.add_sliders_to_widget(widget_right)
        
        grid_left_dock.addWidget(QLabel('Bias'), 0, 0, 1, 1)
        
        self.bias_from = QTextEdit('-3')
        self.bias_from.setMaximumSize(QSize(50, 20))
        self.bias_from.textChanged.connect(self.bias_range_changed)
        grid_left_dock.addWidget(self.bias_from, 0, 1, 1, 1)
        
        grid_left_dock.addWidget(QLabel(':'), 0, 2, 1, 1)
        
        self.bias_to = QTextEdit('3')
        self.bias_to.setMaximumSize(QSize(50, 20))
        self.bias_to.textChanged.connect(self.bias_range_changed)
        grid_left_dock.addWidget(self.bias_to, 0, 3, 1, 1)
        
        grid_left_dock.addWidget(QLabel('V'), 0, 4, 1, 1)
        
        self.button_adj_axes = QPushButton('Adjust Axes')
        self.button_adj_axes.clicked.connect(self.adjust_axes)
        grid_left_dock.addWidget(self.button_adj_axes, 1, 0, 1, 0)
        
        self.button_open_file = QPushButton('Open .dat File')
        self.button_open_file.clicked.connect(self.open_file)
        grid_left_dock.addWidget(self.button_open_file, 2, 0, 1, -1)
        
        #widget_left.setLayout(self.grid_left_dock)
        splitter.addWidget(widget_left)#, 0, 0, 1, 1)
        splitter.addWidget(widget_right)#, 0, 1, 1, 1)
        splitter.setSizes([1, 10]) #ratio 1, 10
        self.aw.bottom_dock.setWidget(splitter)
    
    def adjust_axes(self):
        sort_ylim(self.aw.canvas.ax, padding=True)
        self.aw.canvas.draw()
    
    def open_file(self):
        file, filter = QFileDialog.getOpenFileName(
                directory='/Users/pjh523/Data/Machine_Data_Backups/STM/Patrick/HOPG_ZYB_Pt_923/sample_june', filter='*.dat')
        bs = Bias_Spectroscopy(file)
        self.aw.canvas.plot(bs.bias_data, bs.dIdV_data,
                            remove_line='_data',
                            color='b', label='_data')
    
    def bias_range_changed(self):
        try:
            bias_from = float(self.bias_from.toPlainText())
            bias_to = float(self.bias_to.toPlainText())
            self.bias = linspace(bias_from, bias_to, 256)
            self.aw.canvas.ax.set_xlim(bias_from, bias_to)
        except ValueError:
            print('Bias range must be a number.')
    
    def add_sliders_to_widget(self, widget):
        # add slider to dock
        grid_right_dock = QGridLayout(widget)

        self.r1, self.r1_label = self.create_slider_label_group('R1', 0, 100,
                                                                0,
                                                                1e9 - 1e5,
                                                                1e5,
                                                                grid_right_dock)

        self.r2, self.r2_label = self.create_slider_label_group('R2', 0, 100,
                                                                1,
                                                                1e9 - 1e5,
                                                                1e5,
                                                                grid_right_dock)

        self.c1, self.c1_label = self.create_slider_label_group('C1', 0, 100, 
                                                                2,
                                                                1e-18 - 1e-20,
                                                                1e-20,
                                                                grid_right_dock)

        self.c2, self.c2_label = self.create_slider_label_group('C2', 0, 100, 
                                                                3, 
                                                                1e-18 - 1e-20,
                                                                1e-20,
                                                                grid_right_dock)

        self.q0, self.q0_label = self.create_slider_label_group('Q0', 0, 100,
                                                                4,
                                                                0.5 - (-0.5),
                                                                -0.5,
                                                                grid_right_dock)

        self.t, self.t_label = self.create_slider_label_group('T', 0, 100,
                                                                5,
                                                                293 - 4,
                                                                4,
                                                                grid_right_dock)
        # set default slider Values, in values of index
        self.r1.setValue(60)
        self.r2.setValue(20)
        self.c1.setValue(25)
        self.c2.setValue(10)
        self.q0.setValue((self.q0.maximum()-self.q0.minimum())/2)
        self.t.setValue(0)
        #widget.setLayout(self.grid_bottom_dock)
        #self.aw.bottom_dock.setWidget(widget)

    def create_slider(self, float_range, offset):
        '''create a slider widget'''
        slider = QFloatSlider(float_range, offset)
        slider.valueChanged.connect(self.slider_changed)
        return slider
    
    def create_slider_label_group(self, label, smin, smax, row,
                                  float_range, offset, grid):
        '''create a label slider group'''
        slider = self.create_slider(float_range, offset)
        slider.setMinimum(smin)
        slider.setMaximum(smax)
        
        label = QLabel(label + ': {:.2g}'.format(slider.getFloat()))
        
        grid.addWidget(slider, row, 1, 1, 1)
        grid.addWidget(label, row, 0, 1, 1)
        return slider, label

    def slider_changed(self, val):
        '''catches slider changed signal'''
        self.r1_label.setText('R1 (Ohm): {:.2g}'.format(self.r1.getFloat()))
        self.r2_label.setText('R2 (Ohm): {:.2g}'.format(self.r2.getFloat()))
        self.c1_label.setText('C1 (Farad): {:.2g}'.format(self.c1.getFloat()))
        self.c2_label.setText('C2 (Farad): {:.2g}'.format(self.c2.getFloat()))
        self.q0_label.setText('Q0 (#e): {:.2g}'.format(self.q0.getFloat()))
        self.t_label.setText('T (K): {:.2g}'.format(self.t.getFloat()))
        
        I, dIdV = dbtj(self.bias, self.r1.getFloat(), self.r2.getFloat(),
                       self.c1.getFloat(), self.c2.getFloat(),
                       self.q0.getFloat(), T=self.t.getFloat())
        
        self.aw.canvas.plot(self.bias, dIdV, remove_line='_model', ylim=False,
                            color='r', label='_model')
