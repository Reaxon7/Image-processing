#This is a black-box system with unknown properties
#the syntax is the function y=bboxn(x), where x is the vector
#of the input signal and the y the vector of the output signal
#both vectors are assumed to start at time zero.
import numpy as np
# from scipy import fft


def bbox1(x):

    num_cols = len(x)

    if (num_cols < 1):
        print("Error: Input variable x invalid")


    y = np.linspace(0,0,num_cols)
    y[num_cols-1] = 1.4142*x[num_cols-1]

    for i in range(0,num_cols-1):
        y[i] = 0.74*y[i+1]+1.4142*x[i]

    return y

def bbox2(x):
    num_cols = len(x)

    if (num_cols < 1):
        print("Error: Input variable x invalid")

    y = np.real(np.fft.fft(x))

    return y

def bbox3(x):

    num_cols = len(x)
    filter_length = 10
    if (num_cols < 1):
        print("Error: Input variable x invalid")

    #pad input with zeroes before time zero (shift up)
    x2 = np.linspace(0,0,num_cols+filter_length-1)
    x2[filter_length-1:num_cols+filter_length-1] = x

    ##compute output
    y = np.linspace(0,0,num_cols)

    for i in range(0,num_cols):
        y[i] = np.median(x2[i:i+filter_length-1])

    return y
