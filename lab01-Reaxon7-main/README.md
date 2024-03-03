# Lab 1: Discrete and Continuous-Time Signals

The purpose of this lab is to illustrate the properties of continuous and discrete-time signals using digital computers and the [Python 3](https://www.python.org/) environment. A continuous-time signal takes on a value at every point in time, whereas a discrete-time signal is only defined at integer values of the “time” variable. However, while discrete-time signals can be easily stored and processed on a computer, it is impossible to store the values of a continuous-time signal for all points along a segment of the real line. In later labs, we will see that digital computers are actually restricted to the storage of quantized discrete-time signals. Such signals are appropriately known as digital signals.

How then do we process continuous-time signals? In this lab, we will show that continuous-time signals may be processed by first approximating them by discrete-time signals using a process known as sampling. We will see that proper selection of the spacing between samples is crucial for an efficient and accurate approximation of a continuous-time signal. Excessively close spacing will lead to too much data, whereas excessively distant spacing will lead to a poor approximation of the continuous-time signal. Sampling will be an important topic in future labs, but for now we will use sampling to approximately compute some simple attributes of both real and synthetic signals.

## Getting Started

1. Git clone the repository.

2. Lab instructions can be found in [`lab01_instructions.pdf`](lab01_instructions.pdf).

3. Open a new Terminal, run the command `jupyter notebook` and open `lab01_individual_report.ipynb`.

4. Complete the individual lab report `lab01_individual_report.ipynb`.

5. If you are in a group of more than 1, work with your group member(s) to complete the final lab report `lab01_final_report.ipynb`.

## Submission

```bash
git add lab01_report.ipynb
git commit -m "update lab report"
git push
```
