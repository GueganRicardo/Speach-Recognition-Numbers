# Digit Recognition Project

## Introduction
This project explores digit recognition through speech analysis. We aim to identify each digit from zero to ten using a dataset of spoken numbers by analyzing signals in both time and frequency domains and constructing a navigable decision tree.

## Time Domain Analysis
We began by examining the `.wav` files in the time domain to understand the characteristics of each digit. We normalized the signals for amplitude and timing alignment, making it easier to identify common patterns within each digit.

## Frequency Domain Analysis
We then transformed the signals into the frequency domain using Fourier transform, applying a highpass filter and a Blackman window to isolate useful information. We used spectral roll-off and peak localization metrics to distinguish digits.

## Short-Term Fourier Transform (STFT)
Finally, we used STFT to get a hybrid view of both time and frequency domains. We generated spectrograms for each digit and identified specific zones corresponding to different letters and syllables.

## Conclusion
The project allowed us to apply signal analysis techniques for digit detection. However, we identified potential overfitting issues and suggested more advanced classification systems like recurrent neural networks for future improvements.

## How to Use
1. Clone the repository.
2. Run the `ProjetoATD.m` file to analyze the signals.
3. Explore the results in both time and frequency domains.

## Contributors
- Raul Sofia
- Ricardo Guegan
