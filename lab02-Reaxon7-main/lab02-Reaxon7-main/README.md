# Lab 2: Discrete-Time Systems

## Getting Started

1. Clone the repository by running the command

    ```bash
    git clone https://github.com/Purdue-ECE438-Labs/lab02-<github_username>.git  # using web URL
    # or
    git clone git@github.com:Purdue-ECE438-Labs/lab02-<github_username>.git  # using SSH
    ```

2. Lab instructions can be found in [`lab02_instructions.pdf`](lab02_instructions.pdf).

3. Complete the lab report `lab02_report.ipynb`.

## Special Instructions

### Inserting Diagrams

1. To draw a diagram, you have two options:
    1. Use a drawing tool such as [draw.io](https://app.diagrams.net/) and export it as an image file to `lab02-<github_username>` folder.
    2. Draw the diagram on a blank sheet of paper, scan it and upload the image file to `lab02-<github_username>` folder.

2. To insert it (e.g., `image.png`) in your report
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

In your terminal, change directory to `lab02-<github_username>`.

```bash
git add -A 
git commit -m "update lab report"
git push
```
