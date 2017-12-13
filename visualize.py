#! /usr/bin/python

import sys
from math import sqrt
from tkinter import *

master = None
w = None
size, scale = 0, 0

def init_graphics(sz=700):
	global master, w, size, scale
	master = Tk()
	w = Canvas(master, width=sz, height=sz)
	size, scale = sz, sz / 100
	w.pack()

def grid():
	for y in range(1, 100, 1):
		w.create_line(y*scale, 0, y*scale, size, fill="#cccccc")

	for x in range(1, 100, 1):
		w.create_line(0, x*scale, size, x*scale, fill="#cccccc")

	for y in range(10, 100, 10):
		w.create_line(y*scale, 0, y*scale, size)

	for x in range(10, 100, 10):
		w.create_line(0, x*scale, size, x*scale)

def point(x, y, t='c', s=1.0, c='red'):
	x, y = float(x), float(y)

	if t == 'c':
		w.create_oval((y-s/2) * scale, (x-s/2) * scale, (y+s/2) * scale, (x+s/2)*scale, fill=c)
	elif t == 's':
		w.create_rectangle((y-s/2) * scale, (x-s/2) * scale, (y+s/2) * scale, (x+s/2)*scale, fill=c)
	elif t == 'd':
		pass
	else:
		raise ValueError('Unknown point type!')

def circle(x, y, s, o='blue'):
	x, y = float(x), float(y)
	w.create_oval((y-s) * scale, (x-s) * scale, (y+s) * scale, (x+s)*scale, outline=o, fill='')

def show():
	mainloop()

def visualize(instance_file):
	init_graphics()
	grid()

	with open(instance_file, 'r') as f:
		header = f.readline()
		max_walk = float(header.split(',')[2].strip().split(' ')[0])
		f.readline()	# empty

		stops = []
		studs = []
		school = f.readline().split()[1:]

		for l in f:
			if l == '\n':
				break
			stops.append(l.split()[1:])

		f.readline() 	# empty

		for l in f:
			if l == '\n':
				break
			studs.append(l.split()[1:])

	master.title(header)

	for stud in studs:
		point(*stud, c='blue')

	for stop in stops:
		point(*stop, s=1.5, t='s')
		circle(*stop, max_walk)

	point(*school, s=2, c='black')

	show()

if __name__ == "__main__":
	visualize(sys.argv[1])
