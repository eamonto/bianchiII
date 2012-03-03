#!/usr/bin/env python
import Tkinter
import os
from Tkconstants import *

'''
    Copyright (C) 2012  Edison Montoya, eamonto@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Up to date: 3 Mar 2012					
'''

def calculate(*args):
    try:
        f = open("param.txt","w")

        f.write(input_box[0].get()+"\n")
        f.write(input_box[1].get()+"\n")
        f.write(input_box[2].get()+"\n")
        f.write(input_box[3].get()+"\n")
        f.write(input_box[4].get()+"\n")
        f.write(input_box[5].get()+"\n")
        f.write(input_box[6].get()+"\n")
        f.write(input_box[7].get()+"\n")
        f.write(input_box[8].get()+"\n")
        f.write(input_box[9].get()+"\n")
        f.write(input_box[10].get()+"\n")
        f.write(input_box[11].get()+"\n")
        f.write(input_box[12].get()+"\n")

        f.write(input_box2[0].get()+"\n")
        f.write(input_box2[1].get()+"\n")
        f.write(input_box2[2].get()+"\n")
        f.write(input_box2[3].get()+"\n")
        f.write(input_box2[4].get()+"\n")
        f.write(input_box2[5].get()+"\n")
        f.write(input_box2[6].get()+"\n")

        f.write(input_box3[0].get()+"\n")
        f.write(input_box3[1].get()+"\n")
        f.write(input_box3[2].get()+"\n")
        f.write(input_box3[3].get()+"\n")
        f.write(input_box3[4].get()+"\n")

        f.close()

        print os.system("./exec.out param.txt "+input_box[13].get())

    except ValueError:
        pass


tk = Tkinter.Tk()

#Window Tittle
tk.title("Input Parameters")

#Principal Frame
frame = Tkinter.Frame(tk, relief=RIDGE, borderwidth=2)
frame.grid(column=0, row=0, sticky=(N, W, E, S))

frame2 = Tkinter.Frame(tk, relief=RIDGE, borderwidth=2)
frame2.grid(column=0, row=1, sticky=(N, W, E, S))

frame3 = Tkinter.Frame(tk, relief=RIDGE, borderwidth=2)
frame3.grid(column=0, row=2, sticky=(N, W, E, S))


#io variables 
l = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
input_box = [Tkinter.StringVar() for n in l]

l2 = [0,1,2,3,4,5,6]
input_box2 = [Tkinter.StringVar() for n in l2]

l3 = [0,1,2,3,4]
input_box3 = [Tkinter.StringVar() for n in l3]


#Input box
box_entry=[Tkinter.Entry(frame, width=7, textvariable=input_box[n]) for n in l]

box3_entry=[Tkinter.Entry(frame3, width=7, textvariable=input_box3[n]) for n in l3]



Tkinter.Radiobutton(frame2, text="RK4",   variable=input_box2[0], value="0").grid(column=2, row=0, sticky=W)
Tkinter.Radiobutton(frame2, text="RK-F",  variable=input_box2[0], value="1").grid(column=3, row=0, sticky=W)
Tkinter.Radiobutton(frame2, text="RK-CK", variable=input_box2[0], value="2").grid(column=4, row=0, sticky=W)

Tkinter.Radiobutton(frame2, text="I",   variable=input_box2[1], value="0.0").grid(column=2, row=1, sticky=W)
Tkinter.Radiobutton(frame2, text="II",  variable=input_box2[1], value="1.0").grid(column=3, row=1, sticky=W)

Tkinter.Radiobutton(frame2, text="Internal (N=V)",variable=input_box2[2], value="0").grid(column=2, row=2, sticky=W)
Tkinter.Radiobutton(frame2, text="Cosmic (N=1)",  variable=input_box2[2], value="1").grid(column=3, row=2, sticky=W)

Tkinter.Radiobutton(frame2, text="Classical",variable=input_box2[3], value="0").grid(column=2, row=3, sticky=W)
Tkinter.Radiobutton(frame2, text="Effective",variable=input_box2[3], value="1").grid(column=3, row=3, sticky=W)

Tkinter.Radiobutton(frame2, text="Off",  variable=input_box2[4], value="0").grid(column=2, row=4, sticky=W)
Tkinter.Radiobutton(frame2, text="On",   variable=input_box2[4], value="1").grid(column=3, row=4, sticky=W)

Tkinter.Radiobutton(frame2, text="Off", variable=input_box2[5], value="0").grid(column=2, row=5, sticky=W)
Tkinter.Radiobutton(frame2, text="On",  variable=input_box2[5], value="1").grid(column=3, row=5, sticky=W)

Tkinter.Radiobutton(frame2, text="Inflation", variable=input_box2[6], value="0").grid(column=2, row=6, sticky=W)
Tkinter.Radiobutton(frame2, text="Moduli (Cyclic)",    variable=input_box2[6], value="1").grid(column=3, row=6, sticky=W)


#Text Labels
Tkinter.Label(frame, text="mu1c1 =").grid( column=1, row=0, sticky=E)
Tkinter.Label(frame, text="mu2c2 =").grid( column=1, row=1, sticky=E)
Tkinter.Label(frame, text="mu3c3 =").grid( column=1, row=2, sticky=E)
Tkinter.Label(frame, text="p1 ="   ).grid( column=1, row=3, sticky=E)
Tkinter.Label(frame, text="p2 ="   ).grid( column=1, row=4, sticky=E)
Tkinter.Label(frame, text="p3 ="   ).grid( column=1, row=5, sticky=E)
Tkinter.Label(frame, text="phi ="  ).grid( column=1, row=6, sticky=E)
Tkinter.Label(frame, text="Tolerance error =").grid( column=1, row=7, sticky=E)
Tkinter.Label(frame, text="Max iterations =" ).grid( column=1, row=8, sticky=E)
Tkinter.Label(frame, text="Initial time ="   ).grid( column=1, row=9, sticky=E)
Tkinter.Label(frame, text="Final time ="     ).grid( column=1, row=10, sticky=E)
Tkinter.Label(frame, text="dt ="             ).grid( column=1, row=11, sticky=E)
Tkinter.Label(frame, text="Output every "    ).grid( column=1, row=12, sticky=E)
Tkinter.Label(frame, text="Output file "     ).grid( column=1, row=13, sticky=E)
Tkinter.Label(frame, text="                      ").grid(column=2, row=14, sticky=E)


Tkinter.Label(frame2, text="Integrator").grid(column=1, row=0, sticky=E)
Tkinter.Label(frame2, text="Bianchi").grid(column=1, row=1, sticky=E)
Tkinter.Label(frame2, text="Time").grid(column=1, row=2, sticky=E)
Tkinter.Label(frame2, text="Equations").grid(column=1, row=3, sticky=E)
Tkinter.Label(frame2, text=" Standard output").grid(column=1, row=4, sticky=E)
Tkinter.Label(frame2, text="Potential switch").grid(column=1, row=5, sticky=E)
Tkinter.Label(frame2, text="Select potential").grid(column=1, row=6, sticky=E)

Tkinter.Label(frame3, text="Field mass =").grid(column=1, row=0, sticky=E)
Tkinter.Label(frame3, text=" Field interaction =").grid(column=1, row=1, sticky=E)
Tkinter.Label(frame3, text="V0 =").grid(    column=1, row=2, sticky=E)
Tkinter.Label(frame3, text="sigma1 =").grid(column=1, row=3, sticky=E)
Tkinter.Label(frame3, text="sigma2 =").grid(column=1, row=4, sticky=E)
Tkinter.Label(frame3, text="").grid(column=1, row=5, sticky=E)
Tkinter.Label(frame3, text="").grid(column=2, row=5, sticky=E)


#Button
Tkinter.Button(frame3, text="Calculate", command=calculate).grid(column=2, row=6, sticky=E)


n=0
while (n<14):
    box_entry[n].grid(column=2, row=n, sticky=(W, E))
    n=n+1

n=0
while (n<5):
    box3_entry[n].grid(column=2, row=n, sticky=(W, E))
    n=n+1

#Initial Values
input_box[0].set(1.57079633) # mu1c1
input_box[1].set(1.57079633) # mu2c2
input_box[2].set(1.57079633) # mu3c3
input_box[3].set(100.0)    # p1
input_box[4].set(100.0)      # p2 
input_box[5].set(100.0)      # p3
input_box[6].set(0.0)        # phi
input_box[7].set(1.0e-7)     # error
input_box[8].set(100)        # iteration
input_box[9].set(0.0)        # init time
input_box[10].set(-30)       # final time
input_box[11].set(-5.0e-4)   # dt
input_box[12].set(20)        # output time
input_box[13].set("output")  # name of output file

input_box2[0].set(0)   # integrator
input_box2[1].set(1.0) # bianchi I or II
input_box2[2].set(1)   # time selection
input_box2[3].set(1)   # equations of motion
input_box2[4].set(0)   # stdout switch
input_box2[5].set(0)   # potential switch
input_box2[6].set(0)   # election of potential

input_box3[0].set(0.0) # mass         (inflationary potential)
input_box3[1].set(0.0) # interaction  (inflationary potential)
input_box3[2].set(0.0) # V0     (cyclic potential)
input_box3[3].set(0.0) # sigma1 (cyclic potential)
input_box3[4].set(0.0) # sigma2 (cyclic potential)


#Put the cursor on the input box
box_entry[0].focus()

#If the user press Enter or the button then call the calculate routine
tk.bind('<Return>', calculate)

#Run the window
tk.mainloop()
