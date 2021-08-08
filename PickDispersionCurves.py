# A simple GUI 
# for extracting dispersion curves from dispersion spectrum
# TODO: outlier detection

import  tkinter
from matplotlib.backends.backend_tkagg import(
    FigureCanvasTkAgg, NavigationToolbar2Tk
)
from matplotlib.backend_bases import(
    button_press_handler, key_press_handler
)
from matplotlib import pyplot as plt
from matplotlib.path import Path

import numpy as np
import h5py
import os
import yaml

def find_closel(a,x):
    return int((a-x[0])/(x[1]-x[0]))

def find_closer(a,x):
    return int((a-x[0])/(x[1]-x[0]))+1

def inpolygon(xq, yq, xv, yv):
    """
    reimplement inpolygon in matlab
    :type xq: np.ndarray
    :type yq: np.ndarray
    :type xv: np.ndarray
    :type yv: np.ndarray
    """
    # 合并xv和yv为顶点数组
    vertices = np.vstack((xv, yv)).T
    # 定义Path对象
    path = Path(vertices)
    # 把xq和yq合并为test_points
    test_points = np.hstack([xq.reshape(xq.size, -1), yq.reshape(yq.size, -1)])
    # 得到一个test_points是否严格在path内的mask，是bool值数组
    _in = path.contains_points(test_points)
    # 得到一个test_points是否在path内部或者在路径上的mask
    _in_on = path.contains_points(test_points, radius=-1e-10)
    # 得到一个test_points是否在path路径上的mask
    _on = _in ^ _in_on

    return _in_on, _on

def get_yaml_data(filename):
    file = open(filename,'r',encoding='utf-8')
    file_data = file.read()
    file.close()
    data = yaml.load(file_data,Loader=yaml.FullLoader)
    return data

class DispersionCurve():
    # Class for saving dispersion curves of each dispersion spectrum
    def __init__(self):
        self.norders = 0
        self.dc = dict()
     
    
    def addorder(self,order):
        self.norders += 1
        self.dc[order] = dict()

    def add(self,order,x,y):
        # add dispersion curves for one order
        self.addorder(order)
        self.dc[order]['f'] = x
        self.dc[order]['c'] = y
    
class DispersionSpectrum():
    # simple GUI for extracting dispersion curves
    def __init__(self,f,c,data,outfile,master):
        self.f = f
        self.c = c
        self.nf = len(f)
        self.nc = len(c)
        self.data = data
        self.NumP = len(data)
        self.outfile = outfile
        if os.path.exists(outfile):
            self.dispersionCurves = get_yaml_data(outfile)
        else:
            self.dispersionCurves = dict()
        self.root = master
        self.indx = 0
        self.x = []
        self.y = []
        self.index = range(3,len(f),5)
        self.fp = f[self.index]
        self.F,self.C = np.meshgrid(self.f,self.c)
        
        self.fig, self.ax = plt.subplots(figsize=(8,6),dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig,master=self.root)
        self.canvas.get_tk_widget().pack(side=tkinter.TOP,fill=tkinter.BOTH,expand=1)
        toolbar = NavigationToolbar2Tk(self.canvas,self.root)
        toolbar.update()
        self.drawDS()
        #### set buttons
        self.buttonQuit = tkinter.Button(
            master=self.root, text = "Quit",command = self._quit
        )
        self.buttonNextPoint = tkinter.Button(
            master=self.root, text = "Next Point",command = self.nextPoint
        )
        self.buttonLastPoint = tkinter.Button(
            master=self.root, text = "Last Point",command = self.lastPoint
        )
        self.buttonPick = tkinter.Button(
            master=self.root, text = "Pick",command = self.Pick
        )
        self.buttonSearch = tkinter.Button(
            master=self.root, text = "Search",command = self.search
        )
        self.buttonRedraw = tkinter.Button(
            master=self.root, text = "Redraw",command = self.drawDS
        )
       
        #### buttons distribution
        self.buttonQuit.pack(side=tkinter.RIGHT,fill=tkinter.Y)
        self.buttonNextPoint.pack(side=tkinter.RIGHT,fill=tkinter.Y)
        self.buttonLastPoint.pack(side=tkinter.RIGHT,fill=tkinter.Y)
        self.buttonPick.pack(side=tkinter.RIGHT,fill=tkinter.Y)
        self.buttonSearch.pack(side=tkinter.RIGHT,fill=tkinter.Y)
        self.buttonRedraw.pack(side=tkinter.RIGHT,fill=tkinter.Y)
    
    

    def nextPoint(self):
        if self.indx < (self.NumP-1):
            self.indx = self.indx + 1
        else:
            self.indx = 0
        self.drawDS()

    def lastPoint(self):
        if self.indx > 0 :
            self.indx = self.indx - 1
        else:
            self.indx = self.NumP-1
        self.drawDS()

    def Pick(self):
        self.connect()

    def _quit(self):
        self.writein()
        self.root.quit()
        self.root.destroy()

    def drawDS(self):
        self.order = 0
        self.ax.clear()
        self.ax.pcolormesh(self.f,self.c,self.data[self.indx,:,:],vmin=0.0,vmax=1.0)
        self.ax.set_xlabel('Frequency')
        self.ax.set_ylabel('Phase Velocity')
        self.ax.set_title('order:'+str(self.order))
        if self.indx in self.dispersionCurves:
            for tmpdict in self.dispersionCurves[self.indx]:
                tmpdict = self.dispersionCurves[self.indx][tmpdict]
                self.ax.plot(tmpdict['f'],tmpdict['c'],'r.')
        self.canvas.draw()
        #self.dc = DispersionCurve()
        self.x = []
        self.y = [] 
        title =  str(self.indx)+'/'+str(self.NumP)
        self.root.title(title)
    
    def connect(self):
        self.canvas.mpl_connect("button_press_event",self.on_button_click)
        self.canvas.mpl_connect("key_press_event",self.on_key_press)
    
    def disconnect(self):
        self.canvas.mpl_disconnect("button_press_event")
        self.canvas.mpl_disconnect("key_press_event")

    def on_button_click(self,event):
        x = round(event.xdata,3)
        y = round(event.ydata,3)
        self.ax.plot(x,y,'k.')
        self.x.append(x)
        self.y.append(y)
        self.canvas.draw()

    def on_key_press(self,event):
        valid_range = ["{:d}".format(i) for i in range(10)]
        if event.key in valid_range:
            self.order = int(event.key)
            self.ax.set_title('order:'+str(self.order))
            self.x = []
            self.y = []
        else:
            info = "Invalid key"
            self.ax.set_title(info)
        self.canvas.draw()

    def search(self):
        self.disconnect()
        i1 = find_closer(min(self.x),self.fp)
        i2 = find_closel(max(self.x),self.fp)
        inon, on = inpolygon(self.F,self.C,self.x,self.y)
        inon = np.array(inon)
        inon = inon.reshape(self.nf,self.nc)
        tmp = inon.astype(float)*np.abs(np.squeeze(self.data[self.indx,:,:]))
        self.ax.pcolormesh(self.f,self.c,tmp,vmin=0,vmax=1)
        x = self.fp[i1:i2]
        y = []
        for i  in range(i1,i2):
            y.append(self.c[np.argmax(tmp[:,self.index[i]])])
        
        x = list(x)
        x = [float(t) for t in x]
        y = [float(t) for t in y]
        self.addDC(x,y)
        self.ax.plot(x,y,'k.')
        self.x.append(self.x[0])
        self.y.append(self.y[0])
        self.ax.plot(self.x,self.y,'r')
        self.x = []
        self.y = []
        self.canvas.draw()
        #self.dispersionCurves[str(self.loc[self.indx,:])] = self.dc.dc
        #print(self.dc.dc)

    def writein(self):
        with open(self.outfile,'w',encoding='utf-8') as f:
            yaml.dump(self.dispersionCurves,f)
    
    def addPoint(self):
        if self.indx in self.dispersionCurves:
            pass
        else:
            tmpdict = dict()
            tmpdict[0] = dict()
            tmpdict[0]['f'] = []
            tmpdict[0]['c'] = []
            self.dispersionCurves[self.indx] = tmpdict
    
    def addOrder(self):
        if self.order in self.dispersionCurves[self.indx]:
            pass
        else:
            self.dispersionCurves[self.indx][self.order] = dict()

    def addDC(self,x,y):
        self.addPoint()
        self.addOrder()
        self.dispersionCurves[self.indx][self.order]['f'] = x
        self.dispersionCurves[self.indx][self.order]['c'] = y


#Dir =  'G:\\FromOneDrive20100119\\GetdataPython\\GatheredData\\Circle35'
Dir =  '.'
filename = 'ds.h5'
outfile  = 'out.yml'
filepath = os.path.join(Dir, filename)
#outfile = os.path.join(Dir,outfile)
h5file = h5py.File(filepath, 'r')
data = h5file["ds"][:]
f = h5file["f"][:]
c = h5file["c"][:]
h5file.close()


root = tkinter.Tk()
DispersionSpectrum(f, c, data,  outfile, root)

root.mainloop()