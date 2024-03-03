# Lab 5b: Digital Filter Design (Week 2)

## Getting Started

1. Clone the repository by running the command

    ```bash
    git clone https://github.com/Purdue-ECE438-Labs/lab05b-<github_username>.git  # using web URL
    # or
    git clone git@github.com:Purdue-ECE438-Labs/lab05b-<github_username>.git  # using SSH
    ```

2. Lab instructions can be found in [`lab05b_instructions.pdf`](lab05b_instructions.pdf).

3. Complete the individual lab report `lab05b_individual_report.ipynb`.

4. If you are in a group of more than 1, work with your group member(s) to complete the final lab report `lab05b_final_report.ipynb`.

## Special Instructions

### Inserting Diagrams

To insert an image (e.g., `image.png`) in your report
  
  1. Insert a new cell of type `Markdown`.
  2. Enter the following command `<img src="image.png">`.
  3. To scale the inserted image to 60% for example, enter `<img src="image.png" style="width:60%;height:60%;">`.

### Playing Audio

To play the audio, we need to `import IPython.display as ipd`, and then

```python
ipd.Audio(audio, rate=fs)
```

where `audio` is the audio signal, and `fs` is the sampling rate.

## Submission

### ⚠️ Attention ⚠️

In this lab, interactive plots will be generated. If you choose to __Restart & Run All__ after you finish everything, there is a chance that these interactive plots will not be displayed. Make sure you rerun these cells manually to make sure the plots are displayed correctly. 

### Individual Code Submission

In your terminal, change current directory to `lab05b-<github_username>`.

```bash
git add -A 
git commit -m "update lab report"
git push
```

### Group Lab Report Submission

1. Open `lab05b_final_report.ipynb`.

2. Export it as a PDF file.
    * Ctrl + P (Windows) or Cmd + P (MacOS)
    * Save it as a PDF file. You may want to try setting the scale lower if the plots are not fully visible.

3. ⚠️ Make sure all plots are correctly displayed.

4. Submit it to [Gradescope](https://www.gradescope.com/).
