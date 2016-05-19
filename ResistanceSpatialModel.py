
# coding: utf-8

# In[11]:

# Biologia de sistemas
# Proyecto final
# Resistencia a antibióticos

# Pablo Cárdenas Ramírez
# Juan Sebastián Andrade
# Laura Milena Forero Junco
# Andrea Velásquez

# imports
import numpy as np
import math
import time
import pylab
import copy
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation


# In[12]:

### Setup

# Constants
fps = 40; # framerate in frames/second
size = 2 + 100; # square board size. Remember that edges are not part of simulation (add two rows and columns)
numTurns = 1000; # number of iterations

repProb = [0, 0.3, 0.2]; # probability of reproduction for empty tiles, defectors, and cooperators
deathProb = [0, 0.05, 0.05]; # probability of death for empty tiles, defectors, and cooperators
abDedProb = [0.25, 0]; # probability of dying due to antibiotic if antibiotic is present for cells in tiles [without, with] resistance
abTime = 50; # time at which antibiotic is added
mutProb = 0.01; # probability that a reproduction event results in a cell of the opposite type

seedProbCell = 0.01; # probability that a tile is initialized with a cell
seedProbCoop = 0.5; # probability that a initialized cell is a cooperator

# Arrays and matrices
cells = np.array([[0 if np.random.random()>seedProbCell else 2 if np.random.random()<seedProbCoop else 1 for x in range(size)] for y in range(size)]); # initializaes cell matrix 
cells[0,:] = 0; # set borders to empty. borders are not used in sim
cells[cells.shape[0]-1,:] = 0;
cells[:,0] = 0;
cells[:,cells.shape[0]-1] = 0;

ab = np.array([[0 for x in range(size)] for y in range(size)]); # initializaes cell matrix
res = np.array([[0 for x in range(size)] for y in range(size)]); # initializaes resistance matrix

nc = []; # stores number of cooperators
nd = []; # stores number of defectors
mad= []; # matrix vector, stores system state at each turn

ncoop = 0; # init cell counters used each turn
ndef = 0;


# In[13]:

### Main loop

for t in range(numTurns): # every turn
    
    # We register the system state and update the resistance matrix first, then simulate death and reproduction.
    
    mad.append(copy.deepcopy(cells)); # add copy of cells array to store state at beginning of turn. Uses deepcopy to store values, not references
    nc.append(ncoop); # total number of cooperators at end of turn
    nd.append(ndef); # total number of defectors at end of turn
    ncoop = 0; # cooperator cell counter
    ndef = 0; # defector cell counter
    
    if t == abTime: # add antibiotic at given turn 
        ab[:,:] = 1;
        
    # Resistance matrix has to be updated before death checks
    res[:,:] = 0; # wipe resistance matrix
    for i in range(1,res.shape[0]-1): # every row except edges
        for j in range(1,res.shape[1]-1):# every column except edges
            if not(cells[i,j]-2): # if tile contains cooperator
                res[i-1:i+2,j-1:j+2] = 1; # update resistance matrix (this is why borders are excluded)
                ncoop = ncoop + 1; # update coop counter
            elif not(cells[i,j]-1): # else if defector
                ndef = ndef + 1; # update def counter
    
    for i in range(1,cells.shape[0]-1): # every row except edges
        for j in range(1,cells.shape[1]-1):# every column except edges
            if np.random.random() < deathProb[cells[i,j]]+abDedProb[res[i,j]]*ab[i,j]: # check if death
                cells[i,j] = 0; # die
                
            if np.random.random() < repProb[cells[i,j]]: # check if reproduction
                r = 0; # stores row modifier. Can be -1, 0 or 1
                c = 0; # stores col modifier. Can be -1, 0 or 1
                while (not(r) and not(c))==True: # while both modifiers are 0
                    r = int(math.floor(np.random.random()*3)-1); # recalculate modifiers until at least one of them isn't 0
                    c = int(math.floor(np.random.random()*3)-1);
                
                isMutant = np.random.random()<mutProb; # boolean, is or not mutant
                mutant = cells[i,j]-1 or 2; # stores opposite type
                cells[i+r,j+c] = cells[i,j]*(not isMutant) + mutant*isMutant; # new offspring (be it parental or mutant) occupies random neighboring tile
    
    if (t*100./numTurns)%10 == 0: # tracker
        print str(int(t*100./numTurns))+ '% de turnos simulados';
    
print 'Done simulating';


# In[ ]:

### Rendering

# Board animation
plt.switch_backend('tkagg'); # esoteric stuff. switches plotting backend or something.

fig, ax = plt.subplots(1); # make figure
fig.set_size_inches(fig.get_size_inches()[0]+1, fig.get_size_inches()[1], forward=True); # Makes figure wider to accomodate legend

# make a color map of fixed colors
colMap = matplotlib.colors.ListedColormap(['beige', 'y', 'darkolivegreen']); # list colors for color map
bounds = [0,1,2,3]; # assign colors to ranges (less than or equal to value)
norm = matplotlib.colors.BoundaryNorm(bounds, colMap.N); # produces normalized range for color map 

im = plt.imshow(mad[0][1:cells.shape[0]-1, 1:cells.shape[1]-1], interpolation="None", origin='upper', cmap=colMap, norm=norm); # make axisimage object with colormap out of cells matrix. Exclude borders.

ax.xaxis.tick_top(); # sets x axis on top position

cBar = plt.colorbar(im, cmap=colMap, norm=norm, boundaries=bounds, ticks=[]); # make a color bar with empty ticks

for i, lab in enumerate(['Empty','Defectors (R-)','Cooperators (R+)']): # iterate and add labels to color bar
    cBar.ax.text(1.5, (2*i+1)/6., lab, ha='left', va='center');
    
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5);
turnTxt = ax.text(0.83, 0.07, 't = 0', transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props); # place a text box in upper left in axes coords

# function to update figure
def updatefig(j):
    im.set_data(mad[j][1:cells.shape[0]-1, 1:cells.shape[1]-1]); # set data in the axesimage object
    turnTxt.set_text('t = '+str(j));
    return im, turnTxt # return artists set

# kick off animation
ani = animation.FuncAnimation(fig, updatefig, frames=range(numTurns), interval=1000./fps, blit=True);
plt.show();

# Static population graph
tSeries = range(0,numTurns);
c = []; # total cells
for i in range(0,len(nc)):
    c.append(nc[i]+nd[i]);
    
plt.axvline(x=abTime,ymin=0,ymax=max(c)+2,color="blue",label='Antibiotic',linewidth=2); # tracks antibiotic
        
plt.plot(tSeries,nc,label="Cooperators",color='darkolivegreen');
plt.plot(tSeries,nd,label="Defectors",color='y');
plt.plot(tSeries,c,label="Total",color='green');
plt.title("Population Dynamics");
plt.xlabel("Time (turns)");
plt.ylabel("Number of cells");
plt.legend(loc=1,prop={'size':10});
plt.show();

#Dynamic population graph
ant=[50] #Turn  (or turns) in which antibiotic is introduced

#Time vector
tm=np.linspace(0,len(c)-1,len(c))

#Basic graph parameters
fig = plt.figure()
ax1 = plt.axes(xlim=(0, max(tm)+2), ylim=(0, max(c)+2))
line, = ax1.plot([], [], lw=2)

#Vertical lines (mark antibiotic introduction)
for i in range(0,len(ant)):
    if i==0:
        plt.axvline(x=ant[i],ymin=0,ymax=max(c)+2,color="blue",label=u"Int. Antibi�tico")
    else:
        plt.axvline(x=ant[i],ymin=0,ymax=max(c)+2,color="blue")
#Axis labels  and title
plt.xlabel("Tiempo  (turnos)")
plt.ylabel(u"N�mero")
plt.title(u"Din�mica Poblacional")

#Creation of data lines (one for cooperators, one for defectors, one for both)
plotlays, plotcols,labels = [2], ["black","red","green"],["Defectoras","Cooperadoras","Total"]
lines = []
for index in range(3):
    lobj = ax1.plot([],[],lw=2,color=plotcols[index],label=labels[index])[0]
    lines.append(lobj)

#Define initial state 
def init():
    for line in lines:
        line.set_data([],[])
    return lines

#Define function which updates the data lines as time progresses

x1,y1 = 0,0
x2,y2 = 0,0
x3,y3 = 0,0

def animate(i):
    x1=tm[:i]
    y1=nd[:i]
    lines[0].set_data(x1,y1)

    x2=tm[:i]
    y2=nc[:i]
    lines[1].set_data(x2,y2)
    
    x3=tm[:i]
    y3=c[:i]
    lines[2].set_data(x3,y3)
    return lines,

#Kick off animation
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=range(len(tm)-360), interval=1)
plt.legend(loc=1,fontsize=15)
plt.show();


# In[267]:

print 2-1 or 2; # 


# In[292]:

print range(0,2);


# In[ ]:



